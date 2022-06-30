/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <cassert>

#include <fstream>
#include <numeric>
#include <regex>
#include <set>
#include <shared_mutex>
#include <stack>
#include <tuple>
#include <unordered_map>
#include <mutex>

#include <filesystem>

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/logic/tribool.hpp>

#include "cif++/Cif++.hpp"
#include "cif++/CifParser.hpp"
#include "cif++/CifUtils.hpp"
#include "cif++/CifValidator.hpp"

namespace ba = boost::algorithm;
namespace io = boost::iostreams;
namespace fs = std::filesystem;

namespace cif
{

CIFPP_EXPORT int VERBOSE = 0;

static const char *kEmptyResult = "";

// --------------------------------------------------------------------
// most internal data structures are stored as linked lists
// Item values are stored in a simple struct. They should be const anyway

struct ItemValue
{
	ItemValue *mNext;
	uint32_t mColumnIndex;
	char mText[4];

	ItemValue(const char *v, size_t columnIndex);
	~ItemValue();

	bool empty() const { return mText[0] == 0 or ((mText[0] == '.' or mText[0] == '?') and mText[1] == 0); }
	bool null() const { return mText[0] == '.' and mText[1] == 0; }
	bool unknown() const { return mText[0] == '?' and mText[1] == 0; }

	void *operator new(size_t size, size_t dataSize);
	void operator delete(void *p);
	void operator delete(void *p, size_t dataSize);
};

// --------------------------------------------------------------------

ItemValue::ItemValue(const char *value, size_t columnIndex)
	: mNext(nullptr)
	, mColumnIndex(uint32_t(columnIndex))
{
	assert(columnIndex < std::numeric_limits<uint32_t>::max());
	strcpy(mText, value);
}

ItemValue::~ItemValue()
{
	// remove recursion (and be paranoid)
	while (mNext != nullptr and mNext != this)
	{
		auto n = mNext;
		mNext = n->mNext;
		n->mNext = nullptr;
		delete n;
	}
}

void *ItemValue::operator new(size_t size, size_t dataSize)
{
	return malloc(size - 4 + dataSize + 1);
}

void ItemValue::operator delete(void *p)
{
	free(p);
}

void ItemValue::operator delete(void *p, size_t dataSize)
{
	free(p);
}

// --------------------------------------------------------------------
// itemColumn contains info about a column or field in a Category

struct ItemColumn
{
	std::string mName; // store lower-case, for optimization
	const ValidateItem *mValidator;
};

// itemRow contains the actual values for a Row in a Category

struct ItemRow
{
	~ItemRow();

	void drop(size_t columnIx);
	const char *c_str(size_t columnIx) const;

	std::string str() const
	{
		std::stringstream s;

		s << '{';
		for (auto v = mValues; v != nullptr; v = v->mNext)
		{
			s << mCategory->getColumnName(v->mColumnIndex)
			  << ':'
			  << v->mText;
			if (v->mNext != nullptr)
				s << ", ";
		}
		s << '}';

		return s.str();
	}

	ItemRow *mNext;
	Category *mCategory;
	ItemValue *mValues;
	uint32_t mLineNr = 0;
};

std::ostream &operator<<(std::ostream &os, const ItemRow &r)
{
	os << r.mCategory->name() << '[';
	for (auto iv = r.mValues; iv != nullptr; iv = iv->mNext)
	{
		os << iv->mText;
		if (iv->mNext)
			os << ',';
	}
	os << ']';

	return os;
}

// --------------------------------------------------------------------

ItemRow::~ItemRow()
{
	// remove recursive
	while (mNext != nullptr and mNext != this)
	{
		auto n = mNext;
		mNext = n->mNext;
		n->mNext = nullptr;
		delete n;
	}

	delete mValues;
}

void ItemRow::drop(size_t columnIx)
{
	if (mValues != nullptr and mValues->mColumnIndex == columnIx)
	{
		auto v = mValues;
		mValues = mValues->mNext;
		v->mNext = nullptr;
		delete v;
	}
	else
	{
		for (auto v = mValues; v->mNext != nullptr; v = v->mNext)
		{
			if (v->mNext->mColumnIndex == columnIx)
			{
				auto vn = v->mNext;
				v->mNext = vn->mNext;
				vn->mNext = nullptr;
				delete vn;

				break;
			}
		}
	}

#if DEBUG
	for (auto iv = mValues; iv != nullptr; iv = iv->mNext)
		assert(iv != iv->mNext and (iv->mNext == nullptr or iv != iv->mNext->mNext));
#endif
}

const char *ItemRow::c_str(size_t columnIx) const
{
	const char *result = kEmptyResult;

	for (auto v = mValues; v != nullptr; v = v->mNext)
	{
		if (v->mColumnIndex == columnIx)
		{
			result = v->mText;
			break;
		}
	}

	return result;
}

// --------------------------------------------------------------------

namespace detail
{

	ItemReference &ItemReference::operator=(const std::string &value)
	{
		if (mConst)
			throw std::logic_error("Attempt to write to a constant row");

		if (mRow.mData == nullptr)
			throw std::logic_error("Attempt to write to an uninitialized row");

		mRow.assign(mName, value, false);
		return *this;
	}

	const char *ItemReference::c_str() const
	{
		const char *result = kEmptyResult;

		if (mRow.mData != nullptr /* and mRow.mData->mCategory != nullptr*/)
		{
			//		assert(mRow.mData->mCategory);

			for (auto iv = mRow.mData->mValues; iv != nullptr; iv = iv->mNext)
			{
				if (iv->mColumnIndex == mColumn)
				{
					if (iv->mText[0] != '.' or iv->mText[1] != 0)
						result = iv->mText;

					break;
				}
			}
		}

		return result;
	}

	const char *ItemReference::c_str(const char *defaultValue) const
	{
		const char *result = defaultValue;

		if (mRow.mData != nullptr and mRow.mData->mCategory != nullptr)
		{
			for (auto iv = mRow.mData->mValues; iv != nullptr; iv = iv->mNext)
			{
				if (iv->mColumnIndex == mColumn)
				{
					// only really non-empty values
					if (iv->mText[0] != 0 and ((iv->mText[0] != '.' and iv->mText[0] != '?') or iv->mText[1] != 0))
						result = iv->mText;

					break;
				}
			}

			if (result == defaultValue and mColumn < mRow.mData->mCategory->mColumns.size()) // not found, perhaps the category has a default defined?
			{
				auto iv = mRow.mData->mCategory->mColumns[mColumn].mValidator;
				if (iv != nullptr and not iv->mDefault.empty())
					result = iv->mDefault.c_str();
			}
		}

		return result;
	}

	bool ItemReference::empty() const
	{
		return c_str() == kEmptyResult;
	}

	bool ItemReference::is_null() const
	{
		boost::tribool result;

		if (mRow.mData != nullptr and mRow.mData->mCategory != nullptr)
		{
			for (auto iv = mRow.mData->mValues; iv != nullptr; iv = iv->mNext)
			{
				if (iv->mColumnIndex == mColumn)
				{
					result = iv->mText[0] == '.' and iv->mText[1] == 0;
					break;
				}
			}

			if (result == boost::indeterminate and mColumn < mRow.mData->mCategory->mColumns.size()) // not found, perhaps the category has a default defined?
			{
				auto iv = mRow.mData->mCategory->mColumns[mColumn].mValidator;
				if (iv != nullptr)
					result = iv->mDefaultIsNull;
			}
		}

		return result ? true : false;
	}

	void ItemReference::swap(ItemReference &b)
	{
		Row::swap(mColumn, mRow.mData, b.mRow.mData);
	}

	std::ostream &operator<<(std::ostream &os, ItemReference &item)
	{
		os << item.c_str();
		return os;
	}

} // namespace detail

// --------------------------------------------------------------------
// Datablock implementation

Datablock::Datablock(const std::string_view name)
	: mName(name)
	, mValidator(nullptr)
	, mNext(nullptr)
{
}

Datablock::~Datablock()
{
	delete mNext;
}

auto Datablock::emplace(std::string_view name) -> std::tuple<iterator, bool>
{
	// LRU code

	std::shared_lock lock(mLock);

	bool isNew = true;

	auto i = begin();
	while (i != end())
	{
		if (iequals(name, i->name()))
		{
			isNew = false;

			if (i != begin())
			{
				auto n = std::next(i);
				mCategories.splice(begin(), mCategories, i, n);
			}

			break;
		}

		++i;
	}

	if (isNew)
	{
		mCategories.emplace(begin(), *this, std::string(name), mValidator);

		for (auto &cat : mCategories)
			cat.updateLinks();
	}

	return std::make_tuple(begin(), isNew);
}

Category &Datablock::operator[](std::string_view name)
{
	iterator i;
	std::tie(i, std::ignore) = emplace(name);
	return *i;
}

const Category &Datablock::operator[](std::string_view name) const
{
	using namespace std::literals;

	auto result = get(name);
	if (result == nullptr)
		// throw std::out_of_range("The category with name " + std::string(name) + " does not exist");
	{
		std::unique_lock lock(mLock);
		if (not mNullCategory)
			mNullCategory.reset(new Category(const_cast<Datablock&>(*this), "<null>", nullptr));
		result = mNullCategory.get();
	}
	return *result;
}

Category *Datablock::get(std::string_view name)
{
	std::shared_lock lock(mLock);

	for (auto &cat : mCategories)
	{
		if (iequals(cat.name(), name))
			return &cat;
	}

	return nullptr;
}

const Category *Datablock::get(std::string_view name) const
{
	std::shared_lock lock(mLock);

	for (auto &cat : mCategories)
	{
		if (iequals(cat.name(), name))
			return &cat;
	}

	return nullptr;
}

bool Datablock::isValid()
{
	std::shared_lock lock(mLock);

	if (mValidator == nullptr)
		throw std::runtime_error("Validator not specified");

	bool result = true;
	for (auto &cat : *this)
		result = cat.isValid() and result;
	return result;
}

void Datablock::validateLinks() const
{
	std::shared_lock lock(mLock);

	for (auto &cat : *this)
		cat.validateLinks();
}

void Datablock::setValidator(const Validator *v)
{
	std::shared_lock lock(mLock);

	mValidator = v;

	for (auto &cat : *this)
		cat.setValidator(v);
}

void Datablock::add_software(const std::string_view name, const std::string &classification, const std::string &versionNr, const std::string &versionDate)
{
	std::shared_lock lock(mLock);

	Category &cat = operator[]("software");
	auto ordNr = cat.size() + 1;
	// TODO: should we check this ordinal number???

	cat.emplace({{"pdbx_ordinal", ordNr},
		{"name", name},
		{"version", versionNr},
		{"date", versionDate},
		{"classification", classification}});
}

void Datablock::getTagOrder(std::vector<std::string> &tags) const
{
	std::shared_lock lock(mLock);

	for (auto &cat : *this)
		cat.getTagOrder(tags);
}

void Datablock::write(std::ostream &os)
{
	std::shared_lock lock(mLock);

	os << "data_" << mName << std::endl
	   << "# " << std::endl;

	// mmcif support, sort of. First write the 'entry' Category
	// and if it exists, _AND_ we have a Validator, write out the
	// audit_conform record.

	for (auto &cat : mCategories)
	{
		if (cat.name() == "entry")
		{
			cat.write(os);

			if (mValidator != nullptr)
			{
				Category auditConform(*this, "audit_conform", nullptr);
				auditConform.emplace({{"dict_name", mValidator->dictName()},
					{"dict_version", mValidator->dictVersion()}});
				auditConform.write(os);
			}

			break;
		}
	}

	for (auto &cat : mCategories)
	{
		if (cat.name() != "entry" and cat.name() != "audit_conform")
			cat.write(os);
	}
}

void Datablock::write(std::ostream &os, const std::vector<std::string> &order)
{
	std::shared_lock lock(mLock);

	os << "data_" << mName << std::endl
	   << "# " << std::endl;

	std::vector<std::string> catOrder;
	for (auto &o : order)
	{
		std::string cat, Item;
		std::tie(cat, Item) = splitTagName(o);
		if (find_if(catOrder.rbegin(), catOrder.rend(), [cat](const std::string &s) -> bool
				{ return iequals(cat, s); }) == catOrder.rend())
			catOrder.push_back(cat);
	}

	for (auto &c : catOrder)
	{
		auto cat = get(c);
		if (cat == nullptr)
			continue;

		std::vector<std::string> items;
		for (auto &o : order)
		{
			std::string catName, Item;
			std::tie(catName, Item) = splitTagName(o);

			if (catName == c)
				items.push_back(Item);
		}

		cat->write(os, items);
	}

	// for any Category we missed in the catOrder
	for (auto &cat : mCategories)
	{
		if (find_if(catOrder.begin(), catOrder.end(), [&](const std::string &s) -> bool
				{ return iequals(cat.name(), s); }) != catOrder.end())
			continue;

		cat.write(os);
	}
}

bool operator==(const cif::Datablock &dbA, const cif::Datablock &dbB)
{
	bool result = true;

	std::shared_lock lockA(dbA.mLock);
	std::shared_lock lockB(dbB.mLock);

	std::vector<std::string> catA, catB;

	for (auto &cat : dbA)
	{
		if (not cat.empty())
			catA.push_back(cat.name());
	}
	sort(catA.begin(), catA.end());

	for (auto &cat : dbB)
	{
		if (not cat.empty())
			catB.push_back(cat.name());
	}
	sort(catB.begin(), catB.end());

	// loop over categories twice, to group output
	// First iteration is to list missing categories.

	std::vector<std::string> missingA, missingB;

	auto catA_i = catA.begin(), catB_i = catB.begin();

	while (catA_i != catA.end() and catB_i != catB.end())
	{
		std::string nA = *catA_i;
		ba::to_lower(nA);

		std::string nB = *catB_i;
		ba::to_lower(nB);

		int d = nA.compare(nB);
		if (d > 0)
		{
			missingA.push_back(*catB_i);
			++catB_i;
		}
		else if (d < 0)
		{
			missingB.push_back(*catA_i);
			++catA_i;
		}
		else
			++catA_i, ++catB_i;
	}

	while (catA_i != catA.end())
		missingB.push_back(*catA_i++);

	while (catB_i != catB.end())
		missingA.push_back(*catB_i++);

	if (not(missingA.empty() and missingB.empty()))
	{
		if (cif::VERBOSE > 1)
		{
			std::cerr << "compare of datablocks failed" << std::endl;
			if (not missingA.empty())
				std::cerr << "Categories missing in A: " << ba::join(missingA, ", ") << std::endl
						  << std::endl;

			if (not missingB.empty())
				std::cerr << "Categories missing in B: " << ba::join(missingB, ", ") << std::endl
						  << std::endl;

			result = false;
		}
		else
			return false;
	}

	// Second loop, now compare category values
	catA_i = catA.begin(), catB_i = catB.begin();

	while (catA_i != catA.end() and catB_i != catB.end())
	{
		std::string nA = *catA_i;
		ba::to_lower(nA);

		std::string nB = *catB_i;
		ba::to_lower(nB);

		int d = nA.compare(nB);
		if (d > 0)
			++catB_i;
		else if (d < 0)
			++catA_i;
		else
		{
			if (not(*dbA.get(*catA_i) == *dbB.get(*catB_i)))
			{
				if (cif::VERBOSE > 1)
				{
					std::cerr << "Compare of datablocks failed due to unequal values in category " << *catA_i << std::endl;
					result = false;
				}
				else
					return false;
			}
			++catA_i;
			++catB_i;
		}
	}

	return result;
}

std::ostream &operator<<(std::ostream &os, const Datablock &data)
{
	// whoohoo... this sucks!
	const_cast<Datablock &>(data).write(os);
	return os;
}

// --------------------------------------------------------------------
//

namespace detail
{

	void KeyCompareConditionImpl::prepare(const Category &c)
	{
		mItemIx = c.getColumnIndex(mItemTag);

		auto cv = c.getCatValidator();
		if (cv)
		{
			auto iv = cv->getValidatorForItem(mItemTag);
			if (iv != nullptr and iv->mType != nullptr)
			{
				auto type = iv->mType;
				mCaseInsensitive = type->mPrimitiveType == DDL_PrimitiveType::UChar;
			}
		}
	}

	void KeyIsEmptyConditionImpl::prepare(const Category &c)
	{
		mItemIx = c.getColumnIndex(mItemTag);
	}

	void KeyMatchesConditionImpl::prepare(const Category &c)
	{
		mItemIx = c.getColumnIndex(mItemTag);
	}

} // namespace detail

// --------------------------------------------------------------------
//
//	class to compare two rows based on their keys.

class RowComparator
{
  public:
	RowComparator(Category *cat)
		: RowComparator(cat, cat->getCatValidator()->mKeys.begin(), cat->getCatValidator()->mKeys.end())
	{
	}

	template <typename KeyIter>
	RowComparator(Category *cat, KeyIter b, KeyIter e);

	int operator()(const ItemRow *a, const ItemRow *b) const;

	int operator()(const Row &a, const Row &b) const
	{
		return operator()(a.mData, b.mData);
	}

  private:
	typedef std::function<int(const char *, const char *)> compareFunc;

	typedef std::tuple<size_t, compareFunc> keyComp;

	std::vector<keyComp> mComp;
};

template <typename KeyIter>
RowComparator::RowComparator(Category *cat, KeyIter b, KeyIter e)
{
	auto cv = cat->getCatValidator();

	for (auto ki = b; ki != e; ++ki)
	{
		std::string k = *ki;

		size_t ix = cat->getColumnIndex(k);

		auto iv = cv->getValidatorForItem(k);
		if (iv == nullptr)
			throw std::runtime_error("Incomplete dictionary, no Item Validator for Key " + k);

		auto tv = iv->mType;
		if (tv == nullptr)
			throw std::runtime_error("Incomplete dictionary, no type Validator for Item " + k);

		using namespace std::placeholders;

		mComp.emplace_back(ix, std::bind(&ValidateType::compare, tv, _1, _2));
	}
}

int RowComparator::operator()(const ItemRow *a, const ItemRow *b) const
{
	assert(a);
	assert(b);

	int d = 0;
	for (auto &c : mComp)
	{
		size_t k;
		compareFunc f;

		std::tie(k, f) = c;

		const char *ka = a->c_str(k);
		const char *kb = b->c_str(k);

		d = f(ka, kb);

		if (d != 0)
			break;
	}

	return d;
}

// --------------------------------------------------------------------
//
//	class to keep an index on the keys of a Category. This is a red/black
//	tree implementation.

class CatIndex
{
  public:
	CatIndex(Category *cat);
	~CatIndex();

	ItemRow *find(ItemRow *k) const;

	void insert(ItemRow *r);
	void erase(ItemRow *r);

	// batch create
	void reconstruct();

	// reorder the ItemRow's and returns new head and tail
	std::tuple<ItemRow *, ItemRow *> reorder()
	{
		std::tuple<ItemRow *, ItemRow *> result = std::make_tuple(nullptr, nullptr);

		if (mRoot != nullptr)
		{
			entry *head = findMin(mRoot);
			entry *tail = reorder(mRoot);

			tail->mRow->mNext = nullptr;

			result = std::make_tuple(head->mRow, tail->mRow);
		}

		return result;
	}

	size_t size() const;
	//	bool isValid() const;

  private:
	struct entry
	{
		entry(ItemRow *r)
			: mRow(r)
			, mLeft(nullptr)
			, mRight(nullptr)
			, mRed(true)
		{
		}

		~entry()
		{
			delete mLeft;
			delete mRight;
		}

		ItemRow *mRow;
		entry *mLeft;
		entry *mRight;
		bool mRed;
	};

	entry *insert(entry *h, ItemRow *v);
	entry *erase(entry *h, ItemRow *k);

	//	void validate(entry* h, bool isParentRed, uint32_t blackDepth, uint32_t& minBlack, uint32_t& maxBlack) const;

	entry *rotateLeft(entry *h)
	{
		entry *x = h->mRight;
		h->mRight = x->mLeft;
		x->mLeft = h;
		x->mRed = h->mRed;
		h->mRed = true;
		return x;
	}

	entry *rotateRight(entry *h)
	{
		entry *x = h->mLeft;
		h->mLeft = x->mRight;
		x->mRight = h;
		x->mRed = h->mRed;
		h->mRed = true;
		return x;
	}

	void flipColour(entry *h)
	{
		h->mRed = not h->mRed;

		if (h->mLeft != nullptr)
			h->mLeft->mRed = not h->mLeft->mRed;

		if (h->mRight != nullptr)
			h->mRight->mRed = not h->mRight->mRed;
	}

	bool isRed(entry *h) const
	{
		return h != nullptr and h->mRed;
	}

	entry *moveRedLeft(entry *h)
	{
		flipColour(h);

		if (h->mRight != nullptr and isRed(h->mRight->mLeft))
		{
			h->mRight = rotateRight(h->mRight);
			h = rotateLeft(h);
			flipColour(h);
		}

		return h;
	}

	entry *moveRedRight(entry *h)
	{
		flipColour(h);

		if (h->mLeft != nullptr and isRed(h->mLeft->mLeft))
		{
			h = rotateRight(h);
			flipColour(h);
		}

		return h;
	}

	entry *fixUp(entry *h)
	{
		if (isRed(h->mRight))
			h = rotateLeft(h);

		if (isRed(h->mLeft) and isRed(h->mLeft->mLeft))
			h = rotateRight(h);

		if (isRed(h->mLeft) and isRed(h->mRight))
			flipColour(h);

		return h;
	}

	entry *findMin(entry *h)
	{
		while (h->mLeft != nullptr)
			h = h->mLeft;

		return h;
	}

	entry *eraseMin(entry *h)
	{
		if (h->mLeft == nullptr)
		{
			delete h;
			h = nullptr;
		}
		else
		{
			if (not isRed(h->mLeft) and not isRed(h->mLeft->mLeft))
				h = moveRedLeft(h);

			h->mLeft = eraseMin(h->mLeft);

			h = fixUp(h);
		}

		return h;
	}

	// Fix mNext fields for rows in order of this index
	entry *reorder(entry *e)
	{
		auto result = e;

		if (e->mLeft != nullptr)
		{
			auto l = reorder(e->mLeft);
			l->mRow->mNext = e->mRow;
		}

		if (e->mRight != nullptr)
		{
			auto mr = findMin(e->mRight);
			e->mRow->mNext = mr->mRow;

			result = reorder(e->mRight);
		}

		return result;
	}

	Category &mCat;
	RowComparator mComp;
	entry *mRoot;
};

CatIndex::CatIndex(Category *cat)
	: mCat(*cat)
	, mComp(cat)
	, mRoot(nullptr)
{
}

CatIndex::~CatIndex()
{
	delete mRoot;
}

ItemRow *CatIndex::find(ItemRow *k) const
{
	const entry *r = mRoot;
	while (r != nullptr)
	{
		int d = mComp(k, r->mRow);
		if (d < 0)
			r = r->mLeft;
		else if (d > 0)
			r = r->mRight;
		else
			break;
	}

	return r ? r->mRow : nullptr;
}

void CatIndex::insert(ItemRow *k)
{
	mRoot = insert(mRoot, k);
	mRoot->mRed = false;
}

CatIndex::entry *CatIndex::insert(entry *h, ItemRow *v)
{
	if (h == nullptr)
		return new entry(v);

	int d = mComp(v, h->mRow);
	if (d < 0)
		h->mLeft = insert(h->mLeft, v);
	else if (d > 0)
		h->mRight = insert(h->mRight, v);
	else
		throw std::runtime_error("Duplicate Key violation, cat: " + mCat.name() + " values: " + v->str());

	if (isRed(h->mRight) and not isRed(h->mLeft))
		h = rotateLeft(h);

	if (isRed(h->mLeft) and isRed(h->mLeft->mLeft))
		h = rotateRight(h);

	if (isRed(h->mLeft) and isRed(h->mRight))
		flipColour(h);

	return h;
}

void CatIndex::erase(ItemRow *k)
{
	mRoot = erase(mRoot, k);
	if (mRoot != nullptr)
		mRoot->mRed = false;
}

CatIndex::entry *CatIndex::erase(entry *h, ItemRow *k)
{
	if (mComp(k, h->mRow) < 0)
	{
		if (h->mLeft != nullptr)
		{
			if (not isRed(h->mLeft) and not isRed(h->mLeft->mLeft))
				h = moveRedLeft(h);

			h->mLeft = erase(h->mLeft, k);
		}
	}
	else
	{
		if (isRed(h->mLeft))
			h = rotateRight(h);

		if (mComp(k, h->mRow) == 0 and h->mRight == nullptr)
		{
			delete h;
			return nullptr;
		}

		if (h->mRight != nullptr)
		{
			if (not isRed(h->mRight) and not isRed(h->mRight->mLeft))
				h = moveRedRight(h);

			if (mComp(k, h->mRow) == 0)
			{
				h->mRow = findMin(h->mRight)->mRow;
				h->mRight = eraseMin(h->mRight);
			}
			else
				h->mRight = erase(h->mRight, k);
		}
	}

	return fixUp(h);
}

void CatIndex::reconstruct()
{
	delete mRoot;
	mRoot = nullptr;

	for (auto r : mCat)
		insert(r.mData);

	// maybe reconstruction can be done quicker by using the following commented code.
	// however, I've not had the time to think of a way to set the red/black flag correctly in that case.

	//	std::vector<ItemRow*> rows;
	//	transform(mCat.begin(), mCat.end(), backInserter(rows),
	//		[](Row r) -> ItemRow* { assert(r.mData); return r.mData; });
	//
	//	assert(std::find(rows.begin(), rows.end(), nullptr) == rows.end());
	//
	//	// don't use sort here, it will run out of the stack of something.
	//	// quicksort is notorious for using excessive recursion.
	//	// Besides, most of the time, the data is ordered already anyway.
	//
	//	stable_sort(rows.begin(), rows.end(), [this](ItemRow* a, ItemRow* b) -> bool { return this->mComp(a, b) < 0; });
	//
	//	for (size_t i = 0; i < rows.size() - 1; ++i)
	//		assert(mComp(rows[i], rows[i + 1]) < 0);
	//
	//	deque<entry*> e;
	//	transform(rows.begin(), rows.end(), back_inserter(e),
	//		[](ItemRow* r) -> entry* { return new entry(r); });
	//
	//	while (e.size() > 1)
	//	{
	//		deque<entry*> ne;
	//
	//		while (not e.empty())
	//		{
	//			entry* a = e.front();
	//			e.pop_front();
	//
	//			if (e.empty())
	//				ne.push_back(a);
	//			else
	//			{
	//				entry* b = e.front();
	//				b->mLeft = a;
	//
	//				assert(mComp(a->mRow, b->mRow) < 0);
	//
	//				e.pop_front();
	//
	//				if (not e.empty())
	//				{
	//					entry* c = e.front();
	//					e.pop_front();
	//
	//					assert(mComp(b->mRow, c->mRow) < 0);
	//
	//					b->mRight = c;
	//				}
	//
	//				ne.push_back(b);
	//
	//				if (not e.empty())
	//				{
	//					ne.push_back(e.front());
	//					e.pop_front();
	//				}
	//			}
	//		}
	//
	//		swap (e, ne);
	//	}
	//
	//	assert(e.size() == 1);
	//	mRoot = e.front();
}

size_t CatIndex::size() const
{
	std::stack<entry *> s;
	s.push(mRoot);

	size_t result = 0;

	while (not s.empty())
	{
		entry *e = s.top();
		s.pop();

		if (e == nullptr)
			continue;

		++result;

		s.push(e->mLeft);
		s.push(e->mRight);
	}

	return result;
}

// --------------------------------------------------------------------

RowSet::RowSet(Category &cat)
	: mCat(&cat)
{
}

RowSet::RowSet(Category &cat, Condition &&cond)
	: mCat(&cat)
{
	cond.prepare(cat);

	for (auto r : cat)
	{
		if (cond(cat, r))
			mItems.push_back(r.mData);
	}
}

RowSet::RowSet(const RowSet &rhs)
	: mCat(rhs.mCat)
	, mItems(rhs.mItems)
{
}

RowSet::RowSet(RowSet &&rhs)
	: mCat(rhs.mCat)
	, mItems(std::move(rhs.mItems))
{
}

RowSet::~RowSet()
{
}

RowSet &RowSet::operator=(const RowSet &rhs)
{
	if (this != &rhs)
	{
		mItems = rhs.mItems;
		mCat = rhs.mCat;
	}

	return *this;
}

RowSet &RowSet::operator=(RowSet &&rhs)
{
	if (this != &rhs)
	{
		std::swap(mItems, rhs.mItems);
		mCat = rhs.mCat;
	}

	return *this;
}

RowSet &RowSet::orderBy(std::initializer_list<std::string> items)
{
	RowComparator c(mCat, items.begin(), items.end());

	stable_sort(mItems.begin(), mItems.end(), c);

	return *this;
}

// --------------------------------------------------------------------

Category::Category(Datablock &db, const std::string_view name, const Validator *Validator)
	: mDb(db)
	, mName(name)
	, mValidator(Validator)
	, mHead(nullptr)
	, mTail(nullptr)
	, mIndex(nullptr)
{
	if (mName.empty())
		throw ValidationError("invalid empty name for Category");

	if (mValidator != nullptr)
	{
		mCatValidator = mValidator->getValidatorForCategory(mName);
		if (mCatValidator != nullptr)
		{
			// make sure all required columns are added

			for (auto &k : mCatValidator->mKeys)
				addColumn(k);

			for (auto &k : mCatValidator->mMandatoryFields)
				addColumn(k);

			mIndex = new CatIndex(this);
		}
	}
}

Category::~Category()
{
	delete mHead;
	delete mIndex;
}

void Category::setValidator(const Validator *v)
{
	mValidator = v;

	if (mIndex != nullptr)
	{
		delete mIndex;
		mIndex = nullptr;
	}

	if (mValidator != nullptr)
	{
		mCatValidator = mValidator->getValidatorForCategory(mName);
		if (mCatValidator != nullptr)
		{
			mIndex = new CatIndex(this);
			mIndex->reconstruct();
			//#if DEBUG
			//			assert(mIndex->size() == size());
			//			mIndex->validate();
			//#endif
		}
	}
	else
		mCatValidator = nullptr;

	updateLinks();
}

void Category::updateLinks()
{
	mChildLinks.clear();
	mParentLinks.clear();

	if (mValidator != nullptr)
	{
		for (auto link : mValidator->getLinksForParent(mName))
		{
			auto childCat = mDb.get(link->mChildCategory);
			if (childCat == nullptr)
				continue;
			mChildLinks.push_back({childCat, link});
		}

		for (auto link : mValidator->getLinksForChild(mName))
		{
			auto parentCat = mDb.get(link->mParentCategory);
			if (parentCat == nullptr)
				continue;
			mParentLinks.push_back({parentCat, link});
		}
	}
}

bool Category::hasColumn(std::string_view name) const
{
	return getColumnIndex(name) < mColumns.size();
}

size_t Category::getColumnIndex(std::string_view name) const
{
	size_t result;

	for (result = 0; result < mColumns.size(); ++result)
	{
		if (iequals(name, mColumns[result].mName))
			break;
	}

	if (VERBOSE > 0 and result == mColumns.size() and mCatValidator != nullptr) // validate the name, if it is known at all (since it was not found)
	{
		auto iv = mCatValidator->getValidatorForItem(name);
		if (iv == nullptr)
			std::cerr << "Invalid name used '" << name << "' is not a known column in " + mName << std::endl;
	}

	return result;
}

const std::string &Category::getColumnName(size_t columnIx) const
{
	return mColumns.at(columnIx).mName;
}

std::vector<std::string> Category::getColumnNames() const
{
	std::vector<std::string> result;
	for (auto &c : mColumns)
		result.push_back(c.mName);
	return result;
}

size_t Category::addColumn(std::string_view name)
{
	using namespace std::literals;

	size_t result = getColumnIndex(name);

	if (result == mColumns.size())
	{
		const ValidateItem *itemValidator = nullptr;

		if (mCatValidator != nullptr)
		{
			itemValidator = mCatValidator->getValidatorForItem(name);
			if (itemValidator == nullptr)
				mValidator->reportError("tag " + std::string(name) + " not allowed in Category " + mName, false);
		}

		mColumns.push_back(ItemColumn{std::string(name), itemValidator});
	}

	return result;
}

void Category::reorderByIndex()
{
	if (mIndex != nullptr)
		std::tie(mHead, mTail) = mIndex->reorder();
}

void Category::sort(std::function<int(const Row &, const Row &)> comparator)
{
	if (mHead == nullptr)
		return;

	std::vector<ItemRow *> rows;
	for (auto itemRow = mHead; itemRow != nullptr; itemRow = itemRow->mNext)
		rows.push_back(itemRow);

	std::stable_sort(rows.begin(), rows.end(),
		[&comparator](ItemRow *ia, ItemRow *ib)
		{
			Row ra(ia);
			Row rb(ib);
			return comparator(ra, rb) < 0;
		});

	mHead = rows.front();
	mTail = rows.back();

	auto r = mHead;
	for (size_t i = 1; i < rows.size(); ++i)
		r = r->mNext = rows[i];
	r->mNext = nullptr;

	assert(r == mTail);
	assert(size() == rows.size());
}

std::string Category::getUniqueID(std::function<std::string(int)> generator)
{
	using namespace cif::literals;

	std::string key = "id";
	if (mCatValidator != nullptr and mCatValidator->mKeys.size() == 1)
		key = mCatValidator->mKeys.front();

	// calling size() often is a waste of resources
	if (mLastUniqueNr == 0)
		mLastUniqueNr = size();

	for (;;)
	{
		std::string result = generator(static_cast<int>(mLastUniqueNr++));

		if (exists(Key(key) == result))
			continue;

		return result;
	}
}

size_t Category::size() const
{
	size_t result = 0;

	for (auto pi = mHead; pi != nullptr; pi = pi->mNext)
		++result;

	return result;
}

bool Category::empty() const
{
	return mHead == nullptr or mHead->mValues == nullptr;
}

void Category::drop(const std::string &field)
{
	using namespace std::placeholders;
	auto ci = find_if(mColumns.begin(), mColumns.end(),
		[field](ItemColumn &c) -> bool
		{ return iequals(c.mName, field); });

	if (ci != mColumns.end())
	{
		size_t columnIx = ci - mColumns.begin();

		for (auto pi = mHead; pi != nullptr; pi = pi->mNext)
			pi->drop(columnIx);

		mColumns.erase(ci);
	}
}

const Row Category::operator[](Condition &&cond) const
{
	return const_cast<Category*>(this)->operator[](std::forward<Condition>(cond));
}

Row Category::operator[](Condition &&cond)
{
	Row result;

	cond.prepare(*this);

	for (auto r : *this)
	{
		if (cond(*this, r))
		{
			result = r;
			break;
		}
	}

	return result;
}

bool Category::exists(Condition &&cond) const
{
	bool result = false;

	cond.prepare(*this);

	for (auto r : *this)
	{
		if (cond(*this, r))
		{
			result = true;
			break;
		}
	}

	return result;
}

RowSet Category::orderBy(std::initializer_list<std::string> items)
{
	RowSet result(*this);
	result.insert(result.begin(), begin(), end());

	return result.orderBy(items);
}

void Category::clear()
{
	delete mHead;
	mHead = mTail = nullptr;

	if (mIndex != nullptr)
	{
		delete mIndex;
		mIndex = new CatIndex(this);
	}
}

template <class Iter>
std::tuple<Row, bool> Category::emplace(Iter b, Iter e)
{
	// First, make sure all mandatory fields are supplied
	Row result;
	bool isNew = true;

	if (mCatValidator != nullptr and b != e)
	{
		for (auto &col : mColumns)
		{
			auto iv = mCatValidator->getValidatorForItem(col.mName);

			if (iv == nullptr)
				continue;

			bool seen = false;

			for (auto v = b; v != e; ++v)
			{
				if (iequals(v->name(), col.mName))
				{
					seen = true;
					break;
				}
			}

			if (not seen and iv->mMandatory)
				throw std::runtime_error("missing mandatory field " + col.mName + " for Category " + mName);
		}

		if (mIndex != nullptr)
		{
			std::unique_ptr<ItemRow> nr(new ItemRow{nullptr, this, nullptr});
			Row r(nr.get());
			auto keys = keyFields();

			for (auto v = b; v != e; ++v)
			{
				if (keys.count(v->name()))
					r.assign(v->name(), v->value(), true);
			}

			auto test = mIndex->find(nr.get());
			if (test != nullptr)
			{
				if (VERBOSE > 1)
					std::cerr << "Not inserting new record in " << mName << " (duplicate Key)" << std::endl;
				result = test;
				isNew = false;
			}
		}
	}

	if (isNew)
	{
		auto nr = new ItemRow{nullptr, this, nullptr};

		Row r(nr);

		for (auto v = b; v != e; ++v)
			r.assign(*v, true);

		// if (isOrphan(r))
		// 	throw std::runtime_error("Cannot insert row in category " + mName + " since it would be an orphan");

		if (mHead == nullptr)
		{
			assert(mTail == nullptr);
			mHead = mTail = nr;
		}
		else
		{
			assert(mTail != nullptr);
			assert(mHead != nullptr);
			mTail->mNext = nr;
			mTail = nr;
		}

		result = r;

		if (mIndex != nullptr)
			mIndex->insert(nr);
	}

	return {result, isNew};
}

std::tuple<Row, bool> Category::emplace(Row r)
{
	return emplace(r.begin(), r.end());
}

size_t Category::erase(Condition &&cond)
{
	size_t result = 0;

	cond.prepare(*this);

	auto ri = begin();
	while (ri != end())
	{
		if (cond(*this, *ri))
		{
			ri = erase(ri);
			++result;
		}
		else
			++ri;
	}

	return result;
}

size_t Category::erase(Condition &&cond, std::function<void(const Row &)> &&verbose)
{
	size_t result = 0;

	cond.prepare(*this);

	auto ri = begin();
	while (ri != end())
	{
		if (cond(*this, *ri))
		{
			verbose(*ri);
			ri = erase(ri);
			++result;
		}
		else
			++ri;
	}

	return result;
}

void Category::eraseOrphans(Condition &&cond)
{
	std::vector<ItemRow *> remove;

	cond.prepare(*this);

	for (auto r : *this)
	{
		if (cond(*this, r) and isOrphan(r))
		{
			if (VERBOSE > 1)
				std::cerr << "Removing orphaned record: " << std::endl
						  << r << std::endl
						  << std::endl;

			remove.push_back(r.mData);
		}
	}

	for (auto r : remove)
		erase(iterator(r));
}

void Category::erase(Row r)
{
	erase(iterator(r.mData));
}

auto Category::erase(iterator pos) -> iterator
{
	auto r = *pos;
	iterator result = ++pos;

	iset keys;
	if (mCatValidator)
		keys = iset(mCatValidator->mKeys.begin(), mCatValidator->mKeys.end());

	if (mHead == nullptr)
		throw std::runtime_error("erase");

	if (mIndex != nullptr)
		mIndex->erase(r.mData);

	if (r == mHead)
	{
		mHead = mHead->mNext;
		r.mData->mNext = nullptr;
	}
	else
	{
		for (auto pi = mHead; pi != nullptr; pi = pi->mNext)
		{
			if (pi->mNext == r.mData)
			{
				pi->mNext = r.mData->mNext;
				r.mData->mNext = nullptr;
				break;
			}
		}
	}

	// links are created based on the _pdbx_item_linked_group_list entries
	// in mmcif_pdbx_v50.dic dictionary.
	//
	// For each link group in _pdbx_item_linked_group_list
	// a std::set of keys from one category is mapped to another.
	// If all values in a child are the same as the specified parent ones
	// the child is removed as well, recursively of course.

	if (mValidator != nullptr)
	{
		for (auto &&[childCat, link] : mChildLinks)
		{
			Condition cond;

			for (size_t ix = 0; ix < link->mParentKeys.size(); ++ix)
			{
				const char *value = r[link->mParentKeys[ix]].c_str();

				cond = std::move(cond) && (Key(link->mChildKeys[ix]) == value);
			}

			childCat->eraseOrphans(std::move(cond));
		}
	}

	delete r.mData;

	// reset mTail, if needed
	if (r == mTail)
	{
		mTail = mHead;
		if (mTail != nullptr)
			while (mTail->mNext != nullptr)
				mTail = mTail->mNext;
	}

	return result;
}

Row Category::copyRow(const Row &row)
{
	// copy the values
	std::vector<Item> items;
	std::copy(row.begin(), row.end(), std::back_inserter(items));

	if (mCatValidator and mCatValidator->mKeys.size() == 1)
	{
		auto key = mCatValidator->mKeys.front();
		auto kv = mCatValidator->getValidatorForItem(key);

		for (auto &item : items)
		{
			if (item.name() != key)
				continue;

			if (kv->mType->mPrimitiveType == DDL_PrimitiveType::Numb)
				item.value(getUniqueID(""));
			else
				item.value(getUniqueID(mName + "_id_"));
			break;
		}
	}

	auto &&[result, inserted] = emplace(items.begin(), items.end());
	// assert(inserted);

	return result;
}

void Category::getTagOrder(std::vector<std::string> &tags) const
{
	for (auto &c : mColumns)
		tags.push_back("_" + mName + "." + c.mName);
}

Category::iterator Category::begin()
{
	return iterator(mHead);
}

Category::iterator Category::end()
{
	return iterator();
}

Category::const_iterator Category::cbegin() const
{
	return const_iterator(mHead);
}

Category::const_iterator Category::cend() const
{
	return const_iterator();
}

Category::const_iterator Category::begin() const
{
	return const_iterator(mHead);
}

Category::const_iterator Category::end() const
{
	return const_iterator();
}

bool Category::hasParent(Row r, const Category &parentCat, const ValidateLink &link) const
{
	assert(mValidator != nullptr);
	assert(mCatValidator != nullptr);

	bool result = true;

	Condition cond;
	for (size_t ix = 0; ix < link.mChildKeys.size(); ++ix)
	{
		auto &name = link.mChildKeys[ix];
		auto field = r[name];
		if (field.empty())
		{
			if (mCatValidator->mMandatoryFields.count(name) and field.is_null())
				cond = std::move(cond) and (Key(link.mParentKeys[ix]) == Empty());
		}
		else if (parentCat.mCatValidator->mMandatoryFields.count(link.mParentKeys[ix]))
		{
			const char *value = field.c_str();
			cond = std::move(cond) and (Key(link.mParentKeys[ix]) == value);
		}
		else
		{
			const char *value = field.c_str();
			cond = std::move(cond) and (Key(link.mParentKeys[ix]) == value or Key(link.mParentKeys[ix]) == Empty());
		}
	}

	if (result and not cond.empty())
	{
		result = parentCat.exists(std::move(cond));

		// if (VERBOSE > 3 or (result == false and VERBOSE > 2))
		// 	std::cerr << "result = " << std::boolalpha << result << " for: '" << cond << "' in parent category " << link.mParentCategory << " for child cat " << mName << std::endl;
	}
	// else if (VERBOSE > 3 and cond.empty())
	// 	std::cerr << "Condition is empty due to missing data in parent category " << link.mParentCategory << " for child cat " << mName << std::endl;

	return result;
}

bool Category::isOrphan(Row r)
{
	// be safe
	if (mCatValidator == nullptr)
		return false;

	bool isOrphan = true;
	for (auto &&[parentCat, link] : mParentLinks)
	{
		Condition cond;
		for (size_t ix = 0; ix < link->mChildKeys.size(); ++ix)
		{
			const char *value = r[link->mChildKeys[ix]].c_str();
			cond = std::move(cond) && (Key(link->mParentKeys[ix]) == value);
		}

		// if (VERBOSE > 2)
		// 	std::cerr << "Check condition '" << cond << "' in parent category " << link->mParentCategory << " for child cat " << mName << std::endl;

		if (parentCat->exists(std::move(cond)))
		{
			if (VERBOSE > 2)
				std::cerr << "Not removing because row has a parent in category " << link->mParentCategory << std::endl;

			isOrphan = false;
			break;
		}
	}

	return isOrphan;
}

bool Category::hasChildren(Row r) const
{
	assert(mValidator != nullptr);
	assert(mCatValidator != nullptr);

	bool result = false;

	for (auto &&[childCat, link] : mChildLinks)
	{
		Condition cond;

		for (size_t ix = 0; ix < link->mParentKeys.size(); ++ix)
		{
			const char *value = r[link->mParentKeys[ix]].c_str();

			cond = std::move(cond) && (Key(link->mChildKeys[ix]) == value);
		}

		result = not childCat->find(std::move(cond)).empty();

		if (result)
			break;
	}

	return result;
}

bool Category::hasParents(Row r) const
{
	assert(mValidator != nullptr);
	assert(mCatValidator != nullptr);

	bool result = false;

	for (auto &&[parentCat, link] : mParentLinks)
	{
		Condition cond;

		for (size_t ix = 0; ix < link->mChildKeys.size(); ++ix)
		{
			const char *value = r[link->mChildKeys[ix]].c_str();

			cond = std::move(cond) && (Key(link->mParentKeys[ix]) == value);
		}

		result = not parentCat->find(std::move(cond)).empty();

		if (result)
			break;
	}

	return result;
}

RowSet Category::getChildren(Row r, const char *childCat)
{
	return getChildren(r, mDb[childCat]);
}

RowSet Category::getChildren(Row r, Category &childCat)
{
	assert(mValidator != nullptr);
	assert(mCatValidator != nullptr);

	RowSet result(childCat);

	for (auto &link : mValidator->getLinksForParent(mName))
	{
		if (link->mChildCategory != childCat.mName)
			continue;

		Condition cond;

		for (size_t ix = 0; ix < link->mParentKeys.size(); ++ix)
		{
			const char *value = r[link->mParentKeys[ix]].c_str();

			cond = std::move(cond) && (Key(link->mChildKeys[ix]) == value);
		}

		auto children = childCat.find(std::move(cond));
		result.insert(result.end(), children.begin(), children.end());
	}

	// remove duplicates
	result.make_unique();

	return result;
}

RowSet Category::getParents(Row r, const char *parentCat)
{
	return getParents(r, mDb[parentCat]);
}

RowSet Category::getParents(Row r, Category &parentCat)
{
	assert(mValidator != nullptr);
	assert(mCatValidator != nullptr);

	RowSet result(parentCat);

	for (auto &link : mValidator->getLinksForChild(mName))
	{
		if (link->mParentCategory != parentCat.mName)
			continue;

		Condition cond;

		for (size_t ix = 0; ix < link->mChildKeys.size(); ++ix)
		{
			const char *value = r[link->mChildKeys[ix]].c_str();

			cond = std::move(cond) && (Key(link->mParentKeys[ix]) == value);
		}

		auto parents = parentCat.find(std::move(cond));
		result.insert(result.end(), parents.begin(), parents.end());
	}

	// remove duplicates
	result.make_unique();

	return result;
}

RowSet Category::getLinked(Row r, const char *cat)
{
	return getLinked(r, mDb[cat]);
}

RowSet Category::getLinked(Row r, Category &cat)
{
	RowSet result = getChildren(r, cat);
	if (result.empty())
		result = getParents(r, cat);
	return result;
}

bool Category::isValid()
{
	bool result = true;

	if (mValidator == nullptr)
		throw std::runtime_error("no Validator specified");

	if (empty())
	{
		if (VERBOSE > 2)
			std::cerr << "Skipping validation of empty Category " << mName << std::endl;
		return true;
	}

	if (mCatValidator == nullptr)
	{
		mValidator->reportError("undefined Category " + mName, false);
		return false;
	}

	auto mandatory = mCatValidator->mMandatoryFields;

	for (auto &col : mColumns)
	{
		auto iv = mCatValidator->getValidatorForItem(col.mName);
		if (iv == nullptr)
		{
			mValidator->reportError("Field " + col.mName + " is not valid in Category " + mName, false);
			result = false;
		}

		col.mValidator = iv;

		mandatory.erase(col.mName);
	}

	if (not mandatory.empty())
	{
		mValidator->reportError("In Category " + mName + " the following mandatory fields are missing: " + ba::join(mandatory, ", "), false);
		result = false;
	}

	//#if not defined(NDEBUG)
	//	// check index?
	//	if (mIndex)
	//	{
	//		mIndex->validate();
	//		for (auto r: *this)
	//		{
	//			if (mIndex->find(r.mData) != r.mData)
	//				mValidator->reportError("Key not found in index for Category " + mName);
	//		}
	//	}
	//#endif

	// validate all values
	mandatory = mCatValidator->mMandatoryFields;

	for (auto ri = mHead; ri != nullptr; ri = ri->mNext)
	{
		for (size_t cix = 0; cix < mColumns.size(); ++cix)
		{
			bool seen = false;
			auto iv = mColumns[cix].mValidator;

			if (iv == nullptr)
			{
				mValidator->reportError("invalid field " + mColumns[cix].mName + " for Category " + mName, false);
				result = false;
				continue;
			}

			for (auto vi = ri->mValues; vi != nullptr; vi = vi->mNext)
			{
				if (vi->mColumnIndex == cix)
				{
					seen = true;
					try
					{
						(*iv)(vi->mText);
					}
					catch (const std::exception &e)
					{
						mValidator->reportError("Error validating " + mColumns[cix].mName + ": " + e.what(), false);
						continue;
					}
				}
			}

			if (seen or ri != mHead)
				continue;

			if (iv != nullptr and iv->mMandatory)
			{
				mValidator->reportError("missing mandatory field " + mColumns[cix].mName + " for Category " + mName, false);
				result = false;
			}
		}
	}

	return result;
}

void Category::validateLinks() const
{
	for (auto &&[parentCat, link] : mParentLinks)
	{
		size_t missing = 0;
		for (auto r : *this)
			if (not hasParent(r, *parentCat, *link))
			{
				if (cif::VERBOSE > 1)
				{
					if (missing == 0)
					{
						std::cerr << "Links for " << link->mLinkGroupLabel << " are incomplete" << std::endl
								  << "  These are the items in " << mName << " that don't have matching parent items in " << parentCat->mName << std::endl
								  << std::endl;
					}

					for (auto k : link->mChildKeys)
						std::cerr << k << ": " << r[k].as<std::string>() << std::endl;
					
					std::cerr << std::endl;
				}

				++missing;

			}

		if (missing and VERBOSE == 1)
		{
			std::cerr << "Links for " << link->mLinkGroupLabel << " are incomplete" << std::endl
					  << "  There are " << missing << " items in " << mName << " that don't have matching parent items in " << parentCat->mName << std::endl;
		}
	}
}

const Validator &Category::getValidator() const
{
	if (mValidator == nullptr)
		throw std::runtime_error("no Validator defined yet");
	return *mValidator;
}

iset Category::fields() const
{
	if (mValidator == nullptr)
		throw std::runtime_error("No Validator specified");

	if (mCatValidator == nullptr)
		mValidator->reportError("undefined Category", true);

	iset result;
	for (auto &iv : mCatValidator->mItemValidators)
		result.insert(iv.mTag);
	return result;
}

iset Category::mandatoryFields() const
{
	if (mValidator == nullptr)
		throw std::runtime_error("No Validator specified");

	if (mCatValidator == nullptr)
		mValidator->reportError("undefined Category", true);

	return mCatValidator->mMandatoryFields;
}

iset Category::keyFields() const
{
	if (mValidator == nullptr)
		throw std::runtime_error("No Validator specified");

	if (mCatValidator == nullptr)
		mValidator->reportError("undefined Category", true);

	return iset{mCatValidator->mKeys.begin(), mCatValidator->mKeys.end()};
}

std::set<size_t> Category::keyFieldsByIndex() const
{
	if (mValidator == nullptr)
		throw std::runtime_error("No Validator specified");

	if (mCatValidator == nullptr)
		mValidator->reportError("undefined Category", true);

	std::set<size_t> result;
	for (auto &k : mCatValidator->mKeys)
		result.insert(getColumnIndex(k));

	return result;
}

bool operator==(const Category &a, const Category &b)
{
	using namespace std::placeholders;

	bool result = true;

	//	set<std::string> tagsA(a.fields()), tagsB(b.fields());
	//
	//	if (tagsA != tagsB)
	//		std::cout << "Unequal number of fields" << std::endl;

	auto &validator = a.getValidator();
	auto catValidator = validator.getValidatorForCategory(a.name());
	if (catValidator == nullptr)
		throw std::runtime_error("missing cat validator");

	typedef std::function<int(const char *, const char *)> compType;
	std::vector<std::tuple<std::string, compType>> tags;
	auto keys = catValidator->mKeys;
	std::vector<size_t> keyIx;

	for (auto &tag : a.fields())
	{
		auto iv = catValidator->getValidatorForItem(tag);
		if (iv == nullptr)
			throw std::runtime_error("missing item validator");
		auto tv = iv->mType;
		if (tv == nullptr)
			throw std::runtime_error("missing type validator");
		tags.push_back(std::make_tuple(tag, std::bind(&cif::ValidateType::compare, tv, std::placeholders::_1, std::placeholders::_2)));

		auto pred = [tag](const std::string &s) -> bool
		{ return cif::iequals(tag, s) == 0; };
		if (find_if(keys.begin(), keys.end(), pred) == keys.end())
			keyIx.push_back(tags.size() - 1);
	}

	// a.reorderByIndex();
	// b.reorderByIndex();

	auto rowEqual = [&](const cif::Row &ra, const cif::Row &rb)
	{
		int d = 0;

		for (auto kix : keyIx)
		{
			std::string tag;
			compType compare;

			std::tie(tag, compare) = tags[kix];

			d = compare(ra[tag].c_str(), rb[tag].c_str());

			if (d != 0)
			{
				if (cif::VERBOSE > 1)
					std::cerr << "Values in _" << a.name() << '.' << tag << " are not equal: '" << ra[tag].c_str() << "' != '" << rb[tag].c_str() << '\'' << std::endl;
				break;
			}
		}

		return d == 0;
	};

	auto ai = a.begin(), bi = b.begin();
	while (ai != a.end() or bi != b.end())
	{
		if (ai == a.end() or bi == b.end())
		{
			if (cif::VERBOSE > 1)
			{
				std::cerr << "Unequal number of rows in " << a.name() << std::endl;
				result = false;
				break;
			}
			else
				return false;
		}

		cif::Row ra = *ai, rb = *bi;

		if (not rowEqual(ra, rb))
		{
			if (cif::VERBOSE > 1)
				result = false;
			else
				return false;
		}

		std::vector<std::string> missingA, missingB, different;

		for (auto &tt : tags)
		{
			std::string tag;
			compType compare;

			std::tie(tag, compare) = tt;

			// make it an option to compare unapplicable to empty or something

			const char *ta = ra[tag].c_str();
			if (strcmp(ta, ".") == 0 or strcmp(ta, "?") == 0)
				ta = "";
			const char *tb = rb[tag].c_str();
			if (strcmp(tb, ".") == 0 or strcmp(tb, "?") == 0)
				tb = "";

			if (compare(ta, tb) != 0)
			{
				if (cif::VERBOSE > 1)
				{
					std::cerr << "Values in _" << a.name() << '.' << tag << " are not equal: '" << ta << "' != '" << tb << '\'' << std::endl;
					result = false;
				}
				else
					return false;
			}
		}

		++ai;
		++bi;
	}

	return result;
}

namespace detail
{

	size_t writeValue(std::ostream &os, std::string value, size_t offset, size_t width)
	{
		if (value.find('\n') != std::string::npos or width == 0 or value.length() > 132) // write as text field
		{
			ba::replace_all(value, "\n;", "\n\\;");

			if (offset > 0)
				os << std::endl;
			os << ';' << value;
			if (not ba::ends_with(value, "\n"))
				os << std::endl;
			os << ';' << std::endl;
			offset = 0;
		}
		else if (isUnquotedString(value.c_str()))
		{
			os << value;

			if (value.length() < width)
			{
				os << std::string(width - value.length(), ' ');
				offset += width;
			}
			else
			{
				os << ' ';
				offset += value.length() + 1;
			}
		}
		else
		{
			bool done = false;
			for (char q : {'\'', '"'})
			{
				auto p = value.find(q); // see if we can use the quote character
				while (p != std::string::npos and isNonBlank(value[p + 1]) and value[p + 1] != q)
					p = value.find(q, p + 1);

				if (p != std::string::npos)
					continue;

				os << q << value << q;

				if (value.length() + 2 < width)
				{
					os << std::string(width - value.length() - 2, ' ');
					offset += width;
				}
				else
				{
					os << ' ';
					offset += value.length() + 1;
				}

				done = true;
				break;
			}

			if (not done)
			{
				if (offset > 0)
					os << std::endl;
				os << ';' << value << std::endl
				   << ';' << std::endl;
				offset = 0;
			}
		}

		return offset;
	}

} // namespace detail

void Category::write(std::ostream &os, const std::vector<size_t> &order, bool includeEmptyColumns)
{
	if (empty())
		return;

	// If the first Row has a next, we need a loop_
	bool needLoop = (mHead->mNext != nullptr);

	if (needLoop)
	{
		os << "loop_" << std::endl;

		std::vector<size_t> columnWidths;

		for (auto cix : order)
		{
			auto &col = mColumns[cix];
			os << '_' << mName << '.' << col.mName << ' ' << std::endl;
			columnWidths.push_back(2);
		}

		for (auto Row = mHead; Row != nullptr; Row = Row->mNext)
		{
			for (auto v = Row->mValues; v != nullptr; v = v->mNext)
			{
				if (strchr(v->mText, '\n') == nullptr)
				{
					size_t l = strlen(v->mText);

					if (not isUnquotedString(v->mText))
						l += 2;

					if (l > 132)
						continue;

					if (columnWidths[v->mColumnIndex] < l + 1)
						columnWidths[v->mColumnIndex] = l + 1;
				}
			}
		}

		for (auto Row = mHead; Row != nullptr; Row = Row->mNext) // loop over rows
		{
			size_t offset = 0;

			for (size_t cix : order)
			{
				size_t w = columnWidths[cix];

				std::string s;
				for (auto iv = Row->mValues; iv != nullptr; iv = iv->mNext)
				{
					if (iv->mColumnIndex == cix)
					{
						s = iv->mText;
						break;
					}
				}

				if (s.empty())
					s = "?";

				size_t l = s.length();
				if (not isUnquotedString(s.c_str()))
					l += 2;
				if (l < w)
					l = w;

				if (offset + l > 132 and offset > 0)
				{
					os << std::endl;
					offset = 0;
				}

				offset = detail::writeValue(os, s, offset, w);

				if (offset > 132)
				{
					os << std::endl;
					offset = 0;
				}
			}

			if (offset > 0)
				os << std::endl;
		}
	}
	else
	{
		// first find the indent level
		size_t l = 0;

		for (auto &col : mColumns)
		{
			std::string tag = '_' + mName + '.' + col.mName;

			if (l < tag.length())
				l = tag.length();
		}

		l += 3;

		for (size_t cix : order)
		{
			auto &col = mColumns[cix];

			os << '_' << mName << '.' << col.mName << std::string(l - col.mName.length() - mName.length() - 2, ' ');

			std::string s;
			for (auto iv = mHead->mValues; iv != nullptr; iv = iv->mNext)
			{
				if (iv->mColumnIndex == cix)
				{
					s = iv->mText;
					break;
				}
			}

			if (s.empty())
				s = "?";

			size_t offset = l;
			if (s.length() + l >= kMaxLineLength)
			{
				os << std::endl;
				offset = 0;
			}

			if (detail::writeValue(os, s, offset, 1) != 0)
				os << std::endl;
		}
	}

	os << "# " << std::endl;
}

void Category::write(std::ostream &os)
{
	std::vector<size_t> order(mColumns.size());
	iota(order.begin(), order.end(), 0);
	write(os, order, false);
}

void Category::write(std::ostream &os, const std::vector<std::string> &columns)
{
	// make sure all columns are present
	for (auto &c : columns)
		addColumn(c);

	std::vector<size_t> order;
	order.reserve(mColumns.size());

	for (auto &c : columns)
		order.push_back(getColumnIndex(c));

	for (size_t i = 0; i < mColumns.size(); ++i)
	{
		if (std::find(order.begin(), order.end(), i) == order.end())
			order.push_back(i);
	}

	write(os, order, true);
}

// --------------------------------------------------------------------

void Category::update_value(RowSet &&rows, const std::string &tag, const std::string &value)
{
	if (rows.empty())
		return;

	auto colIx = getColumnIndex(tag);
	if (colIx >= mColumns.size())
		throw std::runtime_error("Invalid column " + value + " for " + mName);

	auto &col = mColumns[colIx];

	// check the value
	if (col.mValidator)
		(*col.mValidator)(value);

	// first some sanity checks, what was the old value and is it the same for all rows?
	std::string oldValue = rows.front()[tag].c_str();
	for (auto &row : rows)
	{
		if (oldValue != row[tag].c_str())
			throw std::runtime_error("Inconsistent old values in update_value");
	}

	if (oldValue == value) // no need to do anything
		return;

	// update rows, but do not cascade
	for (auto &row : rows)
		row.assign(colIx, value, true);

	// see if we need to update any child categories that depend on this value
	for (auto parent : rows)
	{
		for (auto &&[childCat, linked] : mChildLinks)
		{
			if (std::find(linked->mParentKeys.begin(), linked->mParentKeys.end(), tag) == linked->mParentKeys.end())
				continue;

			Condition cond;
			std::string childTag;

			for (size_t ix = 0; ix < linked->mParentKeys.size(); ++ix)
			{
				std::string pk = linked->mParentKeys[ix];
				std::string ck = linked->mChildKeys[ix];

				// TODO add code to *NOT* test mandatory fields for Empty

				if (pk == tag)
				{
					childTag = ck;
					cond = std::move(cond) && Key(ck) == oldValue;
				}
				else
					cond = std::move(cond) && Key(ck) == parent[pk].c_str();
			}

			auto children = RowSet{*childCat, std::move(cond)};
			if (children.empty())
				continue;

			// now be careful. If we search back from child to parent and still find a valid parent row
			// we cannot simply rename the child but will have to create a new child. Unless that new
			// child already exists of course.

			RowSet process(*childCat);

			for (auto child : children)
			{
				Condition cond_c;

				for (size_t ix = 0; ix < linked->mParentKeys.size(); ++ix)
				{
					std::string pk = linked->mParentKeys[ix];
					std::string ck = linked->mChildKeys[ix];

					// TODO add code to *NOT* test mandatory fields for Empty

					cond_c = std::move(cond_c) && Key(pk) == child[ck].c_str();
				}

				auto parents = find(std::move(cond_c));
				if (parents.empty())
				{
					process.push_back(child);
					continue;
				}

				// oops, we need to split this child, unless a row already exists for the new value
				Condition check;

				for (size_t ix = 0; ix < linked->mParentKeys.size(); ++ix)
				{
					std::string pk = linked->mParentKeys[ix];
					std::string ck = linked->mChildKeys[ix];

					// TODO add code to *NOT* test mandatory fields for Empty

					if (pk == tag)
						check = std::move(check) && Key(ck) == value;
					else
						check = std::move(check) && Key(ck) == parent[pk].c_str();
				}

				if (childCat->exists(std::move(check))) // phew..., narrow escape
					continue;

				// create the actual copy, if we can...
				if (childCat->mCatValidator != nullptr and childCat->mCatValidator->mKeys.size() == 1)
				{
					auto copy = childCat->copyRow(child);
					if (copy != child)
					{
						process.push_back(child);
						continue;
					}
				}

				// cannot update this...
				if (cif::VERBOSE > 0)
					std::cerr << "Cannot update child " << childCat->mName << "." << childTag << " with value " << value << std::endl;
			}

			// finally, update the children
			if (not process.empty())
				childCat->update_value(std::move(process), childTag, value);
		}
	}
}

// --------------------------------------------------------------------

Row::Row(const Row &rhs)
	: mData(rhs.mData)
	, mCascade(rhs.mCascade)
{
}

Row::Row(Row &&rhs)
	: mData(rhs.mData)
	, mCascade(rhs.mCascade)
{
	rhs.mData = nullptr;
}

Row::~Row()
{
}

void Row::next() const
{
	if (mData != nullptr)
		mData = mData->mNext;
}

Row &Row::operator=(Row &&rhs)
{
	mData = rhs.mData;
	rhs.mData = nullptr;
	mCascade = rhs.mCascade;
	return *this;
}

Row &Row::operator=(const Row &rhs)
{
	mData = rhs.mData;
	mCascade = rhs.mCascade;
	return *this;
}

void Row::assign(const std::vector<Item> &values)
{
	auto cat = mData->mCategory;

	std::map<std::string, std::tuple<size_t, std::string, std::string>> changed;

	for (auto &value : values)
	{
		auto columnIx = cat->addColumn(value.name());
		auto &col = cat->mColumns[columnIx];
		std::string tag = col.mValidator ? col.mValidator->mTag : std::to_string(columnIx);

		changed[tag] = std::make_tuple(columnIx, operator[](columnIx).c_str(), value.value());

		assign(columnIx, value.value(), true);
	}

	// see if we need to update any child categories that depend on these values
	// auto iv = col.mValidator;
	if (mCascade)
	{
		for (auto &&[childCat, linked] : cat->mChildLinks)
		{
			Condition cond;
			std::string childTag;

			std::vector<Item> newValues;

			for (size_t ix = 0; ix < linked->mParentKeys.size(); ++ix)
			{
				std::string pk = linked->mParentKeys[ix];
				std::string ck = linked->mChildKeys[ix];

				if (changed.count(pk) > 0)
				{
					childTag = ck;
					cond = std::move(cond) && (Key(ck) == std::get<1>(changed[pk]));
					newValues.emplace_back(ck, std::get<2>(changed[pk]));
				}
				else
				{
					const char *value = (*this)[pk].c_str();
					cond = std::move(cond) && (Key(ck) == value);
				}
			}

			auto rows = childCat->find(std::move(cond));
			for (auto &cr : rows)
				cr.assign(newValues);
		}
	}
}

void Row::assign(const Item &value, bool skipUpdateLinked)
{
	assign(value.name(), value.value(), skipUpdateLinked);
}

void Row::assign(std::string_view name, const std::string &value, bool skipUpdateLinked, bool validate)
{
	try
	{
		auto cat = mData->mCategory;
		assign(cat->addColumn(name), value, skipUpdateLinked, validate);
	}
	catch (const std::exception &ex)
	{
		if (cif::VERBOSE >= 0)
			std::cerr << "Could not assign value '" << value << "' to column _" << mData->mCategory->name() << '.' << name << std::endl;
		throw;
	}
}

void Row::assign(size_t column, const std::string &value, bool skipUpdateLinked, bool validate)
{
	if (mData == nullptr)
		throw std::logic_error("invalid Row, no data assigning value '" + value + "' to column with index " + std::to_string(column));

	auto cat = mData->mCategory;
	auto &col = cat->mColumns[column];

	const char *oldValue = nullptr;
	for (auto iv = mData->mValues; iv != nullptr; iv = iv->mNext)
	{
		assert(iv != iv->mNext and (iv->mNext == nullptr or iv != iv->mNext->mNext));

		if (iv->mColumnIndex == column)
		{
			oldValue = iv->mText;
			break;
		}
	}

	if (oldValue != nullptr and value == oldValue) // no need to update
		return;

	std::string oldStrValue = oldValue ? oldValue : "";

	// check the value
	if (col.mValidator and validate)
		(*col.mValidator)(value);

	// If the field is part of the Key for this Category, remove it from the index
	// before updating

	bool reinsert = false;

	if (not skipUpdateLinked and // an update of an Item's value
		cat->mIndex != nullptr and cat->keyFieldsByIndex().count(column))
	{
		reinsert = cat->mIndex->find(mData);
		if (reinsert)
			cat->mIndex->erase(mData);
	}

	// first remove old value with cix

	if (mData->mValues == nullptr)
		; // nothing to do
	else if (mData->mValues->mColumnIndex == column)
	{
		auto iv = mData->mValues;
		mData->mValues = iv->mNext;
		iv->mNext = nullptr;
		delete iv;
	}
	else
	{
		for (auto iv = mData->mValues; iv->mNext != nullptr; iv = iv->mNext)
		{
			if (iv->mNext->mColumnIndex == column)
			{
				auto nv = iv->mNext;
				iv->mNext = nv->mNext;
				nv->mNext = nullptr;
				delete nv;

				break;
			}
		}
	}

	if (not value.empty())
	{
		auto nv = new (value.length()) ItemValue(value.c_str(), column);

		if (mData->mValues == nullptr)
			mData->mValues = nv;
		else
		{
			auto iv = mData->mValues;
			while (iv->mNext != nullptr)
				iv = iv->mNext;
			iv->mNext = nv;
		}
	}

	if (reinsert)
		cat->mIndex->insert(mData);

	// see if we need to update any child categories that depend on this value
	auto iv = col.mValidator;
	if (not skipUpdateLinked and iv != nullptr and mCascade)
	{
		for (auto &&[childCat, linked] : cat->mChildLinks)
		{
			if (find(linked->mParentKeys.begin(), linked->mParentKeys.end(), iv->mTag) == linked->mParentKeys.end())
				continue;

			Condition cond;
			std::string childTag;

			for (size_t ix = 0; ix < linked->mParentKeys.size(); ++ix)
			{
				std::string pk = linked->mParentKeys[ix];
				std::string ck = linked->mChildKeys[ix];

				// TODO add code to *NOT* test mandatory fields for Empty

				if (pk == iv->mTag)
				{
					childTag = ck;
					cond = std::move(cond) && Key(ck) == oldStrValue;
				}
				else
				{
					const char *pk_value = (*this)[pk].c_str();
					if (*pk_value == 0)
						cond = std::move(cond) && Key(ck) == Empty();
					else
						cond = std::move(cond) && ((Key(ck) == pk_value) or Key(ck) == Empty());
				}
			}

			auto rows = childCat->find(std::move(cond));
			if (rows.empty())
				continue;

			// if (cif::VERBOSE > 2)
			// {
			// 	std::cerr << "Parent: " << linked->mParentCategory << " Child: " << linked->mChildCategory << std::endl
			// 			  << cond << std::endl;
			// }

			// Now, suppose there are already rows in child that conform to the new value,
			// we then skip this renam

			Condition cond_n;

			for (size_t ix = 0; ix < linked->mParentKeys.size(); ++ix)
			{
				std::string pk = linked->mParentKeys[ix];
				std::string ck = linked->mChildKeys[ix];

				// TODO add code to *NOT* test mandatory fields for Empty

				if (pk == iv->mTag)
					cond_n = std::move(cond_n) && Key(ck) == value;
				else
				{
					const char *pk_value = (*this)[pk].c_str();
					if (*pk_value == 0)
						cond_n = std::move(cond_n) && Key(ck) == Empty();
					else
						cond_n = std::move(cond_n) && ((Key(ck) == pk_value) or Key(ck) == Empty());
				}
			}

			auto rows_n = childCat->find(std::move(cond_n));
			if (not rows_n.empty())
			{
				if (cif::VERBOSE > 0)
					std::cerr << "Will not rename in child category since there are already rows that link to the parent" << std::endl;

				continue;
			}

			for (auto &cr : rows)
				cr.assign(childTag, value, false);
		}
	}
}

void Row::swap(size_t cix, ItemRow *a, ItemRow *b)
{
	if (a == nullptr or b == nullptr)
		throw std::logic_error("invalid Rows in swap");

	assert(a->mCategory == b->mCategory);
	if (a->mCategory != b->mCategory)
		throw std::logic_error("Categories not same in swap");

	auto cat = a->mCategory;

	// If the field is part of the Key for this Category, remove it from the index
	// before updating

	bool reinsert = false;

	if (cat->mIndex != nullptr and cat->keyFieldsByIndex().count(cix))
	{
		reinsert = true;
		cat->mIndex->erase(a);
		cat->mIndex->erase(b);
	}

	ItemValue *ap = nullptr; // parent of ai
	ItemValue *ai = nullptr;
	ItemValue *bp = nullptr; // parent of bi
	ItemValue *bi = nullptr;

	if (a->mValues == nullptr)
		;
	else if (a->mValues->mColumnIndex == cix)
		ai = a->mValues;
	else
	{
		ap = a->mValues;
		while (ap->mNext != nullptr)
		{
			if (ap->mNext->mColumnIndex == cix)
			{
				ai = ap->mNext;
				ap->mNext = ai->mNext;
				ai->mNext = nullptr;
				break;
			}
			ap = ap->mNext;
		}
	}

	if (b->mValues == nullptr)
		;
	else if (b->mValues->mColumnIndex == cix)
		bi = b->mValues;
	else
	{
		bp = b->mValues;
		while (bp->mNext != nullptr)
		{
			if (bp->mNext->mColumnIndex == cix)
			{
				bi = bp->mNext;
				bp->mNext = bi->mNext;
				bi->mNext = nullptr;
				break;
			}
			bp = bp->mNext;
		}
	}

	if (ai != nullptr)
	{
		if (bp == nullptr)
			b->mValues = ai;
		else
		{
			ai->mNext = bp->mNext;
			bp->mNext = ai;
		}
	}

	if (bi != nullptr)
	{
		if (ap == nullptr)
			a->mValues = bi;
		else
		{
			bi->mNext = ap->mNext;
			ap->mNext = bi;
		}
	}

	if (reinsert)
	{
		cat->mIndex->insert(a);
		cat->mIndex->insert(b);
	}

	if ((ai != nullptr or bi != nullptr))
	{
		auto parentColName = cat->getColumnName(cix);

		// see if we need to update any child categories that depend on these values
		auto parentCatValidator = cat->getCatValidator();

		for (auto &&[childCat, link] : cat->mChildLinks)
		{
			if (find(link->mParentKeys.begin(), link->mParentKeys.end(), parentColName) == link->mParentKeys.end())
				continue;

			auto childCatValidator = childCat->getCatValidator();
			if (childCatValidator == nullptr)
				continue;

			std::string linkChildColName;

			Condition cond[2];
			for (size_t ab = 0; ab < 2; ++ab)
			{
				auto i = ab == 0 ? ai : bi;
				auto r = ab == 0 ? a : b;

				for (size_t ix = 0; ix < link->mChildKeys.size(); ++ix)
				{
					assert(ix < link->mParentKeys.size());
					auto pcix = cat->getColumnIndex(link->mParentKeys[ix]);

					auto childColName = link->mChildKeys[ix];
					bool mandatory =
						find(childCatValidator->mMandatoryFields.begin(), childCatValidator->mMandatoryFields.end(), childColName) != childCatValidator->mMandatoryFields.end() or
						find(parentCatValidator->mMandatoryFields.begin(), parentCatValidator->mMandatoryFields.end(), link->mParentKeys[ix]) != parentCatValidator->mMandatoryFields.end();

					std::string childValue;

					if (pcix == cix)
					{
						linkChildColName = childColName;
						if (not(i == nullptr or strcmp(i->mText, ".") == 0 or strcmp(i->mText, "?") == 0))
							childValue = i->mText;
					}
					else
					{
						std::string ps = r->c_str(pcix);
						if (not(ps.empty() or ps == "." or ps == "?"))
							childValue = ps;
					}

					if (not childValue.empty())
					{
						if (mandatory or pcix == cix)
							cond[ab] = std::move(cond[ab]) and Key(childColName) == childValue;
						else
							cond[ab] = std::move(cond[ab]) and (Key(childColName) == childValue or Key(childColName) == Empty());
					}
					else
						cond[ab] = std::move(cond[ab]) and Key(childColName) == Empty();
				}
			}

			std::vector<conditional_iterator_proxy<Category>> rs;

			// first find the respective rows, then flip values, otherwise you won't find them anymore!
			for (size_t ab = 0; ab < 2; ++ab)
			{
				if (cond[ab].empty())
					continue;

				// if (VERBOSE > 1)
				// 	std::cerr << "Fixing link from " << cat->mName << " to " << childCat->mName << " with " << std::endl
				// 		 << cond[ab] << std::endl;

				rs.push_back(childCat->find(std::move(cond[ab])));
			}

			for (size_t ab = 0; ab < 2; ++ab)
			{
				auto i = ab == 0 ? bi : ai;

				for (auto r : rs[ab])
				{
					// now due to the way links are defined, we might have found a row
					// that contains an empty value for all child columns...
					// Now, that's not a real hit, is it?

					size_t n = 0;
					for (auto c : link->mChildKeys)
						if (r[c].empty())
							++n;

					if (n == link->mChildKeys.size())
					{
						if (VERBOSE > 1)
							std::cerr << "All empty columns, skipping" << std::endl;
					}
					else
					{
						if (VERBOSE > 0)
							std::cerr << "In " << childCat->mName << " changing " << linkChildColName << ": " << r[linkChildColName].as<std::string>() << " => " << (i ? i->mText : "") << std::endl;
						r[linkChildColName] = i ? i->mText : "";
					}
				}
			}
		}
	}
}

size_t Row::ColumnForItemTag(std::string_view itemTag) const
{
	size_t result = 0;
	if (mData != nullptr)
	{
		auto cat = mData->mCategory;
		result = cat->getColumnIndex(itemTag);
	}
	return result;
}

bool Row::empty() const
{
	return mData == nullptr or mData->mValues == nullptr;
}

auto Row::begin() const -> const_iterator
{
	return const_iterator(mData, mData ? mData->mValues : nullptr);
}

auto Row::end() const -> const_iterator
{
	return const_iterator(mData, nullptr);
}

uint32_t Row::lineNr() const
{
	return mData ? mData->mLineNr : 0;
}

void Row::lineNr(uint32_t l)
{
	if (mData)
		mData->mLineNr = l;
}

Row::const_iterator::const_iterator(ItemRow *data, ItemValue *ptr)
	: mData(data)
	, mPtr(ptr)
{
	if (mPtr != nullptr)
		fetch();
}

Row::const_iterator &Row::const_iterator::operator++()
{
	if (mPtr != nullptr)
		mPtr = mPtr->mNext;

	if (mPtr != nullptr)
		fetch();

	return *this;
}

void Row::const_iterator::fetch()
{
	mCurrent = Item(
		mData->mCategory->getColumnName(mPtr->mColumnIndex),
		mPtr->mText);
}

std::ostream &operator<<(std::ostream &os, const Row &row)
{
	auto category = row.mData->mCategory;
	std::string catName = category->name();
	for (auto item = row.mData->mValues; item != nullptr; item = item->mNext)
	{
		std::string tagName = category->getColumnName(item->mColumnIndex);
		os << '_' << catName << '.' << tagName << ' ' << item->mText << std::endl;
	}

	return os;
}

// --------------------------------------------------------------------

File::File()
	: mHead(nullptr)
	, mValidator(nullptr)
{
}

File::File(std::istream &is, bool validate)
	: File()
{
	load(is);
}

File::File(const std::filesystem::path &path, bool validate)
	: File()
{
	try
	{
		load(path);
	}
	catch (const std::exception &ex)
	{
		if (cif::VERBOSE >= 0)
			std::cerr << "Error while loading file " << path << std::endl;
		throw;
	}
}

File::File(const char *data, std::size_t length)
{
	load(data, length);
}

File::File(File &&rhs)
	: mHead(nullptr)
	, mValidator(nullptr)
{
	std::swap(mHead, rhs.mHead);
	std::swap(mValidator, rhs.mValidator);
}

File::~File()
{
	delete mHead;
}

void File::append(Datablock *e)
{
	e->setValidator(mValidator);

	if (mHead == nullptr)
		mHead = e;
	else
	{
		auto ie = mHead;
		for (;;)
		{
			if (iequals(ie->getName(), e->getName()))
				throw ValidationError("Datablock " + e->getName() + " already defined in File");

			if (ie->mNext == nullptr)
			{
				ie->mNext = e;
				break;
			}

			ie = ie->mNext;
		}
	}
}

void File::load(const std::filesystem::path &p)
{
	fs::path path(p);

	std::ifstream inFile(p, std::ios_base::in | std::ios_base::binary);
	if (not inFile.is_open())
		throw std::runtime_error("Could not open file: " + path.string());

	io::filtering_stream<io::input> in;
	std::string ext;

	if (path.extension() == ".gz")
	{
		in.push(io::gzip_decompressor());
		ext = path.stem().extension().string();
	}

	in.push(inFile);

	try
	{
		load(in);
	}
	catch (const std::exception &ex)
	{
		if (cif::VERBOSE >= 0)
			std::cerr << "Error loading file " << path << std::endl;
		throw;
	}
}

void File::save(const std::filesystem::path &p)
{
	fs::path path(p);

	std::ofstream outFile(p, std::ios_base::out | std::ios_base::binary);
	io::filtering_stream<io::output> out;

	if (path.extension() == ".gz")
	{
		out.push(io::gzip_compressor());
		path = path.stem();
	}

	out.push(outFile);
	save(out);
}

void File::load(std::istream &is)
{
	auto saved = mValidator;
	setValidator(nullptr);

	Parser p(is, *this);
	p.parseFile();

	if (saved != nullptr)
	{
		setValidator(saved);
		(void)isValid();
	}
}

void File::load(std::istream &is, const std::string &datablock)
{
	auto saved = mValidator;
	setValidator(nullptr);

	Parser p(is, *this);
	p.parseSingleDatablock(datablock);

	if (saved != nullptr)
	{
		setValidator(saved);
		(void)isValid();
	}
}

void File::load(const char *data, std::size_t length)
{
	bool gzipped = length > 2 and data[0] == static_cast<char>(0x1f) and data[1] == static_cast<char>(0x8b);

	struct membuf : public std::streambuf
	{
		membuf(char *data, size_t length) { this->setg(data, data, data + length); }
	} buffer(const_cast<char *>(data), length);

	std::istream is(&buffer);

	io::filtering_stream<io::input> in;
	if (gzipped)
		in.push(io::gzip_decompressor());
	in.push(is);

	load(is);
}

void File::save(std::ostream &os)
{
	Datablock *e = mHead;
	while (e != nullptr)
	{
		e->write(os);
		e = e->mNext;
	}
}

void File::write(std::ostream &os, const std::vector<std::string> &order)
{
	Datablock *e = mHead;
	while (e != nullptr)
	{
		e->write(os, order);
		e = e->mNext;
	}
}

Datablock *File::get(std::string_view name) const
{
	const Datablock *result = mHead;
	while (result != nullptr and not iequals(result->mName, name))
		result = result->mNext;
	return const_cast<Datablock *>(result);
}

Datablock &File::operator[](std::string_view name)
{
	using namespace std::literals;

	Datablock *result = mHead;
	while (result != nullptr and not iequals(result->mName, name))
		result = result->mNext;
	if (result == nullptr)
		throw std::runtime_error("Datablock " + std::string(name) + " does not exist");
	return *result;
}

Datablock &File::front()
{
	assert(mHead);
	return *mHead;
}

Datablock &File::back()
{
	assert(mHead);
	auto *block = mHead;
	while (block->mNext != nullptr)
		block = block->mNext;
	return *block;
}

bool File::isValid()
{
	if (mValidator == nullptr)
	{
		if (VERBOSE > 0)
			std::cerr << "No dictionary loaded explicitly, loading default" << std::endl;

		loadDictionary();
	}

	bool result = true;
	for (auto d = mHead; d != nullptr; d = d->mNext)
		result = d->isValid() and result;
	return result;
}

void File::validateLinks() const
{
	for (auto d = mHead; d != nullptr; d = d->mNext)
		d->validateLinks();
}

const Validator &File::getValidator() const
{
	if (mValidator == nullptr)
		throw std::runtime_error("no Validator defined yet");
	return *mValidator;
}

void File::loadDictionary()
{
	loadDictionary("mmcif_ddl");
}

void File::loadDictionary(const char *dict)
{
	setValidator(&ValidatorFactory::instance()[dict]);
}

void File::setValidator(const Validator *v)
{
	mValidator = v;

	for (auto d = mHead; d != nullptr; d = d->mNext)
		d->setValidator(mValidator);
}

void File::getTagOrder(std::vector<std::string> &tags) const
{
	for (auto d = mHead; d != nullptr; d = d->mNext)
		d->getTagOrder(tags);
}

auto File::iterator::operator++() -> iterator &
{
	mCurrent = mCurrent->mNext;
	return *this;
}

auto File::begin() const -> iterator
{
	return iterator(mHead);
}

auto File::end() const -> iterator
{
	return iterator(nullptr);
}

} // namespace cif
