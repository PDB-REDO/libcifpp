// cif parsing library

#include <cassert>

#include <stack>
#include <tuple>
#include <regex>
#include <set>
#include <unordered_map>
#include <numeric>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#if defined(USE_RSRC)
#include "cif++/mrsrc.h"
#endif

#include "cif++/Cif++.h"
#include "cif++/CifParser.h"
#include "cif++/CifValidator.h"
#include "cif++/CifUtils.h"

using namespace std;
namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;

extern int VERBOSE;

namespace cif
{

static const char* kEmptyResult = "";
	
// --------------------------------------------------------------------
// most internal data structures are stored as linked lists
// Item values are stored in a simple struct. They should be const anyway

struct ItemValue
{
	ItemValue*				mNext;
	uint32					mColumnIndex;
	char					mText[0];
	
	ItemValue(const char* v, uint32 columnIndex);
	~ItemValue();
	
	void* operator new(size_t size, size_t dataSize);
	void operator delete(void* p);
};

// --------------------------------------------------------------------

ItemValue::ItemValue(const char* value, uint32 columnIndex)
	: mNext(nullptr), mColumnIndex(columnIndex)
{
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

void* ItemValue::operator new(size_t size, size_t dataSize)
{
	return malloc(size + dataSize + 1);
}

void ItemValue::operator delete(void* p)
{
	free(p);
}

// --------------------------------------------------------------------

// itemColumn contains info about a column or field in a Category

struct ItemColumn
{
	string					mName;		// store lower-case, for optimization
	const ValidateItem*		mValidator;
};

// itemRow contains the actual values for a Row in a Category

struct ItemRow
{
	~ItemRow();

	void drop(uint32 columnIx);
	const char* c_str(uint32 columnIx) const;
	
	string str() const
	{
		stringstream s;

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
	
	ItemRow*				mNext;
	Category*				mCategory;
	ItemValue*				mValues;
	uint32					mLineNr = 0;
};

ostream& operator<<(ostream& os, const ItemRow& r)
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

void ItemRow::drop(uint32 columnIx)
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

const char* ItemRow::c_str(uint32 columnIx) const
{
	const char* result = kEmptyResult;
	
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

template<>
ItemReference& ItemReference::operator=(const string& value)
{
	Row(mRow).assign(mName, value, false);
	return *this;
}

const char* ItemReference::c_str() const
{
	const char* result = kEmptyResult;
	
	if (mRow != nullptr /* and mRow->mCategory != nullptr*/)
	{
//		assert(mRow->mCategory);
		
		auto cix = mRow->mCategory->getColumnIndex(mName);
		
		for (auto iv = mRow->mValues; iv != nullptr; iv = iv->mNext)
		{
			if (iv->mColumnIndex == cix)
			{
				if (iv->mText[0] != '.' or iv->mText[1] != 0)
					result = iv->mText;
	
				break;
			}
		}
	}
	
	return result;
}

const char* ItemReference::c_str(const char* defaultValue) const
{
	const char* result = defaultValue;
	
	if (mRow != nullptr and mRow->mCategory != nullptr)
	{
//		assert(mRow->mCategory);
		
		auto cix = mRow->mCategory->getColumnIndex(mName);
		
		if (cix < mRow->mCategory->mColumns.size())
		{
			auto iv = mRow->mCategory->mColumns[cix].mValidator;
			if (iv != nullptr and not iv->mDefault.empty())
				result = iv->mDefault.c_str();
			
			for (auto iv = mRow->mValues; iv != nullptr; iv = iv->mNext)
			{
				if (iv->mColumnIndex == cix)
				{
					// only really non-empty values
					if (iv->mText[0] != 0 and ((iv->mText[0] != '.' and iv->mText[0] != '?') or iv->mText[1] != 0))
						result = iv->mText;
	
					break;
				}
			}
		}
	}
	
	return result;
}

bool ItemReference::empty() const
{
	return c_str() == kEmptyResult;
}

void ItemReference::swap(ItemReference& b)
{
	Row::swap(mName, mRow, b.mRow);
}

}

// --------------------------------------------------------------------
// Datablock implementation

Datablock::Datablock(const string& name)
	: mName(name), mValidator(nullptr), mNext(nullptr)
{
}

Datablock::~Datablock()
{
	delete mNext;
}

string Datablock::firstItem(const string& tag) const
{
	string result;

	string catName, itemName;
	std::tie(catName, itemName) = splitTagName(tag);
	
	for (auto& cat: mCategories)
	{
		if (iequals(cat.name(), catName))
		{
			result = cat.getFirstItem(itemName.c_str()).as<string>();
			break;
		}
	}

	return result;
}

auto Datablock::emplace(const string& name) -> tuple<iterator,bool>
{
	bool isNew = false;
	iterator i = find_if(begin(), end(), [name](const Category& cat) -> bool
		{ return iequals(cat.name(), name); });
	
	if (i == end())
	{
		isNew = true;
		i = mCategories.emplace(end(), *this, name, mValidator);
	}
	
	return make_tuple(i, isNew);
}

Category& Datablock::operator[](const string& name)
{
	iterator i;
	std::tie(i, ignore) = emplace(name);
	return *i;
}

Category* Datablock::get(const string& name)
{
	auto i = find_if(begin(), end(), [name](const Category& cat) -> bool 
		{ return iequals(cat.name(), name); });
	
	return i == end() ? nullptr : &*i;
}

bool Datablock::isValid()
{
	if (mValidator == nullptr)
		throw runtime_error("Validator not specified");

	bool result = true;
	for (auto& cat: *this)
		result = cat.isValid() and result;
	return result;
}

void Datablock::setValidator(Validator* v)
{
	mValidator = v;

	for (auto& cat: *this)
		cat.setValidator(v);
}

void Datablock::getTagOrder(vector<string>& tags) const
{
	for (auto& cat: *this)
		cat.getTagOrder(tags);
}

void Datablock::write(ostream& os)
{
	os << "data_" << mName << endl
	   << "# " << endl;
	
	// mmcif support, sort of. First write the 'entry' Category
	// and if it exists, _AND_ we have a Validator, write out the
	// audit_conform record.

	for (auto& cat: mCategories)
	{
		if (cat.name() == "entry")
		{
			cat.write(os);
			
			if (mValidator != nullptr)
			{
				Category auditConform(*this, "audit_conform", nullptr);
				auditConform.emplace({
					{ "dict_name", mValidator->dictName() },
					{ "dict_version", mValidator->dictVersion() }
				});
				auditConform.write(os);
			}
			
			break;
		}
	}

	for (auto& cat: mCategories)
	{
		if (cat.name() != "entry" and cat.name() != "audit_conform")
			cat.write(os);
	}
}

void Datablock::write(ostream& os, const vector<string>& order)
{
	os << "data_" << mName << endl
	   << "# " << endl;
	
	vector<string> catOrder;
	for (auto& o: order)
	{
		string cat, Item;
		std::tie(cat, Item) = splitTagName(o);
		if (find_if(catOrder.rbegin(), catOrder.rend(), [cat](const string& s) -> bool { return iequals(cat, s); }) == catOrder.rend())
			catOrder.push_back(cat);
	}

	for (auto& c: catOrder)
	{
		auto cat = get(c);
		if (cat == nullptr)
			continue;
		
		vector<string> items;
		for (auto& o: order)
		{
			string catName, Item;
			std::tie(catName, Item) = splitTagName(o);
			
			if (catName == c)
				items.push_back(Item);
		}
		
		cat->write(os, items);
	}
	
	// for any Category we missed in the catOrder
	for (auto& cat: mCategories)
	{
		if (find_if(catOrder.begin(), catOrder.end(), [&](const string& s) -> bool { return iequals(cat.name(), s); }) != catOrder.end())
			continue;
		
		cat.write(os);
	}
	
	
//	// mmcif support, sort of. First write the 'entry' Category
//	// and if it exists, _AND_ we have a Validator, write out the
//	// auditConform record.
//
//	for (auto& cat: mCategories)
//	{
//		if (cat.name() == "entry")
//		{
//			cat.write(os);
//			
//			if (mValidator != nullptr)
//			{
//				Category auditConform(*this, "audit_conform", nullptr);
//				auditConform.emplace({
//					{ "dict_name", mValidator->dictName() },
//					{ "dict_version", mValidator->dictVersion() }
//				});
//				auditConform.write(os);
//			}
//			
//			break;
//		}
//	}
//
//	for (auto& cat: mCategories)
//	{
//		if (cat.name() != "entry" and cat.name() != "audit_conform")
//			cat.write(os);
//	}
}

// --------------------------------------------------------------------
//

namespace detail
{
	
void KeyIsEmptyConditionImpl::prepare(const Category& c)
{
	mItemIx = c.getColumnIndex(mItemTag);
}

void KeyMatchesConditionImpl::prepare(const Category& c)
{
	mItemIx = c.getColumnIndex(mItemTag);
}

}

// --------------------------------------------------------------------
//
//	class to compare two rows based on their keys.

class RowComparator
{
  public:

	RowComparator(Category* cat)
		: RowComparator(cat, cat->getCatValidator()->mKeys.begin(), cat->getCatValidator()->mKeys.end())
	{
	}
		
	template<typename KeyIter>
	RowComparator(Category* cat, KeyIter b, KeyIter e);
	
	int operator()(const ItemRow* a, const ItemRow* b) const;

	int operator()(const Row& a, const Row& b) const
	{
		return operator()(a.mData, b.mData);
	}
	
  private:
	typedef function<int(const char*,const char*)>	compareFunc;

	typedef tuple<size_t,compareFunc>	keyComp;

	vector<keyComp>	mComp;
};

template<typename KeyIter>
RowComparator::RowComparator(Category* cat, KeyIter b, KeyIter e)
{
	auto cv = cat->getCatValidator();
	
	for (auto ki = b; ki != e; ++ki)
	{
		string k = *ki;
		
		size_t ix = cat->getColumnIndex(k);

		auto iv = cv->getValidatorForItem(k);
		if (iv == nullptr)
			throw runtime_error("Incomplete dictionary, no Item Validator for Key " + k);
		
		auto tv = iv->mType;
		if (tv == nullptr)
			throw runtime_error("Incomplete dictionary, no type Validator for Item " + k);
		
		using namespace placeholders;
		
		mComp.emplace_back(ix, bind(&ValidateType::compare, tv, _1, _2));
	}
}

int RowComparator::operator()(const ItemRow* a, const ItemRow* b) const
{
	assert(a);
	assert(b);

	int d = 0;
	for (auto& c: mComp)
	{
		size_t k;
		compareFunc f;
		
		std::tie(k, f) = c;
		
		const char* ka = a->c_str(k);
		const char* kb = b->c_str(k);
		
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
	CatIndex(Category* cat);
	~CatIndex();
	
	ItemRow* find(ItemRow* k) const;

	void insert(ItemRow* r);
	void erase(ItemRow* r);

	// batch create
	void reconstruct();
	
	// reorder the ItemRow's and returns new head and tail
	tuple<ItemRow*,ItemRow*> reorder()
	{
		tuple<ItemRow*,ItemRow*> result = make_tuple(nullptr, nullptr);
		
		if (mRoot != nullptr)
		{
			entry* head = findMin(mRoot);
			entry* tail = reorder(mRoot);
			
			tail->mRow->mNext = nullptr;
			
			result = make_tuple(head->mRow, tail->mRow);
		}
			
		return result;
	}

	size_t size() const;
//	bool isValid() const;
	
  private:

	struct entry
	{
		entry(ItemRow* r)
			: mRow(r), mLeft(nullptr), mRight(nullptr), mRed(true) {}
		
		~entry()
		{
			delete mLeft;
			delete mRight;
		}

		ItemRow*		mRow;
		entry*			mLeft;
		entry*			mRight;
		bool			mRed;
	};

	entry* insert(entry* h, ItemRow* v);
	entry* erase(entry* h, ItemRow* k);

//	void validate(entry* h, bool isParentRed, uint32 blackDepth, uint32& minBlack, uint32& maxBlack) const;

	entry* rotateLeft(entry* h)
	{
		entry* x = h->mRight;
		h->mRight = x->mLeft;
		x->mLeft = h;
		x->mRed = h->mRed;
		h->mRed = true;
		return x;
	}
	
	entry* rotateRight(entry* h)
	{
		entry* x = h->mLeft;
		h->mLeft = x->mRight;
		x->mRight = h;
		x->mRed = h->mRed;
		h->mRed = true;
		return x;
	}
	
	void flipColour(entry* h)
	{
		h->mRed = not h->mRed;
		
		if (h->mLeft != nullptr)
			h->mLeft->mRed = not h->mLeft->mRed;
	
		if (h->mRight != nullptr)
			h->mRight->mRed = not h->mRight->mRed;
	}
	
	bool isRed(entry* h) const
	{
		return h != nullptr and h->mRed;
	}
	
	entry* moveRedLeft(entry* h)
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
	
	entry* moveRedRight(entry* h)
	{
		flipColour(h);
		
		if (h->mLeft != nullptr and isRed(h->mLeft->mLeft))
		{
			h = rotateRight(h);
			flipColour(h);
		}
		
		return h;
	}
	
	entry* fixUp(entry* h)
	{
		if (isRed(h->mRight))
			h = rotateLeft(h);
		
		if (isRed(h->mLeft) and isRed(h->mLeft->mLeft))
			h = rotateRight(h);
		
		if (isRed(h->mLeft) and isRed(h->mRight))
			flipColour(h);
		
		return h;
	}
	
	entry* findMin(entry* h)
	{
		while (h->mLeft != nullptr)
			h = h->mLeft;

		return h;
	}
	
	entry* eraseMin(entry* h)
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
	entry* reorder(entry* e)
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
	
	Category&			mCat;
	RowComparator		mComp;
	entry*				mRoot;
};

CatIndex::CatIndex(Category* cat)
	: mCat(*cat), mComp(cat), mRoot(nullptr)
{
}

CatIndex::~CatIndex()
{
	delete mRoot;
}

ItemRow* CatIndex::find(ItemRow* k) const
{
	const entry* r = mRoot;
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

void CatIndex::insert(ItemRow* k)
{
	mRoot = insert(mRoot, k);
	mRoot->mRed = false;
}

CatIndex::entry* CatIndex::insert(entry* h, ItemRow* v)
{
	if (h == nullptr)
		return new entry(v);
	
	int d = mComp(v, h->mRow);
	if (d < 0)		h->mLeft = insert(h->mLeft, v);
	else if (d > 0)	h->mRight = insert(h->mRight, v);
	else
		throw runtime_error("Duplicate Key violation, cat: " + mCat.name() + " values: " + v->str());

	if (isRed(h->mRight) and not isRed(h->mLeft))
		h = rotateLeft(h);

	if (isRed(h->mLeft) and isRed(h->mLeft->mLeft))
		h = rotateRight(h);
	
	if (isRed(h->mLeft) and isRed(h->mRight))	
		flipColour(h);
	
	return h;
}

void CatIndex::erase(ItemRow* k)
{
	mRoot = erase(mRoot, k);
	if (mRoot != nullptr)
		mRoot->mRed = false;
}

CatIndex::entry* CatIndex::erase(entry* h, ItemRow* k)
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
	
	for (auto r: mCat)
		insert(r.mData);

// maybe reconstruction can be done quicker by using the following commented code.
// however, I've not had the time to think of a way to set the red/black flag correctly in that case.
	
//	vector<ItemRow*> rows;
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
	stack<entry*> s;
	s.push(mRoot);
	
	size_t result = 0;
	
	while (not s.empty())
	{
		entry* e = s.top();
		s.pop();
		
		if (e == nullptr)
			continue;
		
		++result;
		
		s.push(e->mLeft);
		s.push(e->mRight);
	}
	
	return result;
}

//bool CatIndex::isValid() const
//{
//	bool result = true;
//	
//	if (mRoot != nullptr)
//	{
//		uint32 minBlack = numeric_limits<uint32>::max();
//		uint32 maxBlack = 0;
//		
//		assert(not mRoot->mRed);
//		
//		result = isValid(mRoot, false, 0, minBlack, maxBlack);
//		assert(minBlack == maxBlack);
//	}
//	
//	return result;
//}
//
//bool CatIndex::validate(entry* h, bool isParentRed, uint32 blackDepth, uint32& minBlack, uint32& maxBlack) const
//{
//	bool result = true;
//	
//	if (h->mRed)
//		assert(not isParentRed);
//	else
//		++blackDepth;
//	
//	if (isParentRed)
//		assert(not h->mRed);
//	
//	if (h->mLeft != nullptr and h->mRight != nullptr)
//	{
//		if (isRed(h->mLeft))
//			assert(not isRed(h->mRight));
//		if (isRed(h->mRight))
//			assert(not isRed(h->mLeft));
//	}
//	
//	if (h->mLeft != nullptr)
//	{
//		assert(mComp(h->mLeft->mRow, h->mRow) < 0);
//		validate(h->mLeft, h->mRed, blackDepth, minBlack, maxBlack);
//	}
//	else
//	{
//		if (minBlack > blackDepth)
//			minBlack = blackDepth;
//		if (maxBlack < blackDepth)
//			maxBlack = blackDepth;
//	}
//	
//	if (h->mRight != nullptr)
//	{
//		assert(mComp(h->mRight->mRow, h->mRow) > 0);
//		validate(h->mRight, h->mRight, blackDepth, minBlack, maxBlack);
//	}
//	else
//	{
//		if (minBlack > blackDepth)
//			minBlack = blackDepth;
//		if (maxBlack < blackDepth)
//			maxBlack = blackDepth;
//	}
//}

// --------------------------------------------------------------------

RowSet::RowSet(const RowSet& rhs)
	: base_type(rhs)
	, mCat(rhs.mCat)
{
}

RowSet::RowSet(RowSet&& rhs)
	: base_type(move(rhs))
	, mCat(rhs.mCat)
{
}

RowSet& RowSet::operator=(const RowSet& rhs)
{
	if (this != &rhs)
	{
		base_type::operator=(rhs);
		mCat = rhs.mCat;
	}
	
	return *this;
}

RowSet& RowSet::operator=(RowSet&& rhs)
{
	if (this != &rhs)
	{
		base_type::operator=(move(rhs));
		mCat = rhs.mCat;
	}
	
	return *this;
}

RowSet::RowSet(Category& cat)
	: mCat(&cat)
{
}

RowSet& RowSet::orderBy(initializer_list<string> items)
{
	RowComparator c(mCat, items.begin(), items.end());
	
	stable_sort(begin(), end(), c);
	
	return *this;
}

// --------------------------------------------------------------------

Category::Category(Datablock& db, const string& name, Validator* Validator)
	: mDb(db), mName(name), mValidator(Validator)
	, mHead(nullptr), mTail(nullptr), mIndex(nullptr)
{
	if (mName.empty())
		throw ValidationError("invalid empty name for Category");
	
	if (mValidator != nullptr)
	{
		mCatValidator = mValidator->getValidatorForCategory(mName);
		if (mCatValidator != nullptr)
		{
			// make sure all required columns are added
			
			for (auto& k: mCatValidator->mKeys)
				addColumn(k);

			for (auto& k: mCatValidator->mMandatoryFields)
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

void Category::setValidator(Validator* v)
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
}

size_t Category::getColumnIndex(const string& name) const
{
	size_t result;

	for (result = 0; result < mColumns.size(); ++result)
	{
		if (iequals(name, mColumns[result].mName))
			break;
	}
	
	return result;
}

const string& Category::getColumnName(size_t columnIx) const
{
	return mColumns.at(columnIx).mName;
}

vector<string> Category::getColumnNames() const
{
	vector<string> result;
	for (auto& c: mColumns)
		result.push_back(c.mName);
	return result;
}

size_t Category::addColumn(const string& name)
{
	size_t result = getColumnIndex(name);
	
	if (result == mColumns.size())
	{
		const ValidateItem* itemValidator = nullptr;
		
		if (mCatValidator != nullptr)
		{
			itemValidator = mCatValidator->getValidatorForItem(name);
			if (itemValidator == nullptr)
				mValidator->reportError("tag " + name + " not allowed in Category " + mName);
		}
		
		mColumns.push_back({name, itemValidator});
	}
	
	return result;
}

void Category::reorderByIndex()
{
	if (mIndex != nullptr)
		std::tie(mHead, mTail) = mIndex->reorder();
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

void Category::drop(const string& field)
{
	using namespace placeholders;
	auto ci = find_if(mColumns.begin(), mColumns.end(),
		[field](ItemColumn& c) -> bool { return iequals(c.mName, field); });

	if (ci != mColumns.end())
	{
		uint32 columnIx = ci - mColumns.begin();
		
		for (auto pi = mHead; pi != nullptr; pi = pi->mNext)
			pi->drop(columnIx);
		
		mColumns.erase(ci);
	}
}

Row Category::operator[](Condition&& cond)
{
	Row result;

	cond.prepare(*this);
	
	for (auto r: *this)
	{
		if (cond(*this, r))
		{
			result = r;
			break;
		}
	}

	return result;
}	

RowSet Category::find(Condition&& cond)
{
	RowSet result(*this);
	
	cond.prepare(*this);
	
	for (auto r: *this)
	{
		if (cond(*this, r))
			result.push_back(r);
	}
	return result;
}

bool Category::exists(Condition&& cond)
{
	bool result = false;

	cond.prepare(*this);
	
	for (auto r: *this)
	{
		if (cond(*this, r))
		{
			result = true;
			break;
		}
	}

	return result;
}

RowSet Category::orderBy(std::initializer_list<string> items)
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

template<class Iter>
tuple<Row,bool> Category::emplace(Iter b, Iter e)
{
	// First, make sure all mandatory fields are supplied
	tuple<Row,bool> result = make_tuple(Row(), true);

	if (mCatValidator != nullptr and b != e)
	{
		for (auto& col: mColumns)
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
				throw runtime_error("missing mandatory field " + col.mName + " for Category " + mName);
		}
		
		if (mIndex != nullptr)
		{
			unique_ptr<ItemRow> nr(new ItemRow{nullptr, this, nullptr});
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
					cerr << "Not inserting new record in " << mName << " (duplicate Key)" << endl;
				result = make_tuple(Row(test), false);
			}
		}
	}
	
	if (get<1>(result))
	{
		auto nr = new ItemRow{nullptr, this, nullptr};

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

		Row r(nr);

		for (auto v = b; v != e; ++v)
			r.assign(*v, true);
		
		get<0>(result) = r;

		if (mIndex != nullptr)
			mIndex->insert(nr);
	}
	
	return result;
}

tuple<Row,bool> Category::emplace(Row r)
{
	return emplace(r.begin(), r.end());
}

void Category::erase(Condition&& cond)
{
	RowSet remove(*this);
	
	cond.prepare(*this);

	for (auto r: *this)
	{
		if (cond(*this, r))
			remove.push_back(r);
	}

	for (auto r: remove)
		erase(r);
}

void Category::erase(iterator p)
{
	erase(*p);
}

void Category::erase(Row r)
{
	iset keys;
	if (mCatValidator)
		keys = iset(mCatValidator->mKeys.begin(), mCatValidator->mKeys.end());
	
	for (auto& col: mColumns)
	{
		auto iv = col.mValidator;
		if (iv == nullptr or iv->mChildren.empty())
			continue;
		
		if (not keys.count(col.mName))
			continue;
		
		const char* value = r[col.mName].c_str();
		
		for (auto child: iv->mChildren)
		{
			if (child->mCategory == nullptr)
				continue;
			
			auto childCat = mDb.get(child->mCategory->mName);
			if (childCat == nullptr)
				continue;
				
			auto rows = childCat->find(Key(child->mTag) == value);
			for (auto& cr: rows)
				childCat->erase(cr);
		}
	}

	if (mHead == nullptr)
		throw runtime_error("erase");

	if (mIndex != nullptr)
		mIndex->erase(r.mData);
	
	if (r == mHead)
	{
		mHead = mHead->mNext;
		r.mData->mNext = nullptr;
		delete r.mData;
	}
	else
	{
		for (auto pi = mHead; pi != nullptr; pi = pi->mNext)
		{
			if (pi->mNext == r.mData)
			{
				pi->mNext = r.mData->mNext;
				r.mData->mNext = nullptr;
				delete r.mData;
				break;
			}
		}
	}
}

void Category::getTagOrder(vector<string>& tags) const
{
	for (auto& c: mColumns)
		tags.push_back("_" + mName + "." + c.mName);
}

const detail::ItemReference Category::getFirstItem(const char* itemName) const
{
	return detail::ItemReference{itemName, mHead};
}

Category::iterator Category::begin()
{
	return iterator(mHead);
}

Category::iterator Category::end()
{
	return iterator(nullptr);
}

bool Category::isValid()
{
	bool result = true;
	
	if (mValidator == nullptr)
		throw runtime_error("no Validator specified");

	if (empty())
	{
		if (VERBOSE > 2)
			cerr << "Skipping validation of empty Category " << mName << endl;
		return true;
	}
	
	if (mCatValidator == nullptr)
	{
		mValidator->reportError("undefined Category " + mName);
		return false;
	}
	
	auto mandatory = mCatValidator->mMandatoryFields;

	for (auto& col: mColumns)
	{
		auto iv = mCatValidator->getValidatorForItem(col.mName);
		if (iv == nullptr)
		{
			mValidator->reportError("Field " + col.mName + " is not valid in Category " + mName);
			result = false;
		}
		
		col.mValidator = iv;
		
		mandatory.erase(col.mName);
	}
	
	if (not mandatory.empty())
	{
		mValidator->reportError("In Category " + mName + " the following mandatory fields are missing: " + ba::join(mandatory, ", "));
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
				mValidator->reportError("invalid field " + mColumns[cix].mName + " for Category " + mName);
				result = false;
				continue;
			}
			
			for (auto vi = ri->mValues; vi != nullptr; vi = vi->mNext)
			{
				if (vi->mColumnIndex == cix)
				{
					seen = true;
 					(*iv)(vi->mText);
				}
			}
			
			if (seen or ri != mHead)
				continue;
			
			if (iv != nullptr and iv->mMandatory)
			{
				mValidator->reportError("missing mandatory field " + mColumns[cix].mName + " for Category " + mName);
				result = false;
			}
		}
	}
	
	return result;
}

const Validator& Category::getValidator() const
{
	if (mValidator == nullptr)
		throw runtime_error("no Validator defined yet");
	return *mValidator;
}

iset Category::fields() const
{
	if (mValidator == nullptr)
		throw runtime_error("No Validator specified");
	
	if (mCatValidator == nullptr)
		mValidator->reportError("undefined Category");
	
	iset result;
	for (auto& iv: mCatValidator->mItemValidators)
		result.insert(iv.mTag);
	return result;
}

iset Category::mandatoryFields() const
{
	if (mValidator == nullptr)
		throw runtime_error("No Validator specified");
	
	if (mCatValidator == nullptr)
		mValidator->reportError("undefined Category");
	
	return mCatValidator->mMandatoryFields;
}

iset Category::keyFields() const
{
	if (mValidator == nullptr)
		throw runtime_error("No Validator specified");
	
	if (mCatValidator == nullptr)
		mValidator->reportError("undefined Category");
	
	return iset{ mCatValidator->mKeys.begin(), mCatValidator->mKeys.end() };
}

auto Category::iterator::operator++() -> iterator&
{
	mCurrent = Row(mCurrent.data()->mNext);
	return *this;
}

namespace detail
{

size_t writeValue(ostream& os, string value, size_t offset, size_t width)
{
	if (value.find('\n') != string::npos or width == 0 or value.length() >= 132)					// write as text field
	{
		ba::replace_all(value, "\n;", "\n\\;");

		if (offset > 0)
			os << endl;
		os << ';' << value;
		if (not ba::ends_with(value, "\n"))
			os << endl;
		os << ';' << endl;
		offset = 0;
	}
	else if (isUnquotedString(value.c_str()))
	{
		os << value;

		if (value.length() < width)
		{
			os << string(width - value.length(), ' ');
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
		for (char q: { '\'', '"'})
		{
			auto p = value.find(q);	// see if we can use the quote character
			while (p != string::npos and isNonBlank(value[p + 1]) and value[p + 1] != q)
				p = value.find(q, p + 1);
			
			if (p != string::npos)
				continue;
			
			os << q << value << q;

			if (value.length() + 2 < width)
			{
				os << string(width - value.length() - 2, ' ');
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
				os << endl;
			os << ';' << value << endl
			   << ';' << endl;
			offset = 0;
		}
	}

	return offset;
}
	
}

void Category::write(ostream& os, const vector<int>& order, bool includeEmptyColumns)
{
	if (empty())
		return;
	
	// If the first Row has a next, we need a loop_
	bool needLoop = (mHead->mNext != nullptr);
	
	if (needLoop)
	{
		os << "loop_" << endl;
		
		vector<size_t> columnWidths;
		
		for (auto cix: order)
		{
			auto& col = mColumns[cix];
			os << '_' << mName << '.' << col.mName << ' ' << endl;
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

					if (l >= 132)
						continue;

					if (columnWidths[v->mColumnIndex] < l + 1)
						columnWidths[v->mColumnIndex] = l + 1;
				}
			}
		}
		
		for (auto Row = mHead; Row != nullptr; Row = Row->mNext)	// loop over rows
		{
			size_t offset = 0;
		
			for (size_t cix: order)
			{
				size_t w = columnWidths[cix];
				
				string s;
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

				if (offset + l >= 132 and offset > 0)
				{
					os << endl;
					offset = 0;
				}
				
				offset = detail::writeValue(os, s, offset, w);
				
				if (offset >= 132)
				{
					os << endl;
					offset = 0;
				}
			}
			
			if (offset > 0)
				os << endl;
		}
	}
	else
	{
		// first find the indent level 
		size_t l = 0;
		
		for (auto& col: mColumns)
		{
			string tag = '_' + mName + '.' + col.mName;
			
			if (l < tag.length())
				l = tag.length();
		}
		
		l += 3;
		
		for (size_t cix: order)
		{
			auto& col = mColumns[cix];
			
			os << '_' << mName << '.' << col.mName << string(l - col.mName.length() - mName.length() - 2, ' ');
			
			string s;
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
				os << endl;
				offset = 0;
			}

			if (detail::writeValue(os, s, offset, 1) != 0)
				os << endl;
		}
	}

	os << "# " << endl;
}

void Category::write(ostream& os)
{
	vector<int> order(mColumns.size());
	iota(order.begin(), order.end(), 0);
	write(os, order, false);
}

void Category::write(ostream& os, const vector<string>& columns)
{
	// make sure all columns are present
	for (auto& c: columns)
		addColumn(c);
	
	vector<int> order;
	order.reserve(mColumns.size());

	for (auto& c: columns)
		order.push_back(getColumnIndex(c));

	for (size_t i = 0; i < mColumns.size(); ++i)
	{
		if (std::find(order.begin(), order.end(), i) == order.end())
			order.push_back(i);
	}

	write(os, order, true);
}

// --------------------------------------------------------------------

Row::Row(const Row& rhs)
	: mData(rhs.mData)
{
}

Row& Row::operator=(const Row& rhs)
{
	mData = rhs.mData;
	return *this;
}

void Row::assign(const string& name, const string& value, bool emplacing)
{
	if (mData == nullptr)
		throw logic_error("invalid Row, no data assigning value '" + value + "' to " + name);
	
	auto cat = mData->mCategory;
	auto cix = cat->addColumn(name);
	auto& col = cat->mColumns[cix];
	auto& db = cat->mDb;

	const char* oldValue = nullptr;
	for (auto iv = mData->mValues; iv != nullptr; iv = iv->mNext)
	{
		assert(iv != iv->mNext and (iv->mNext == nullptr or iv != iv->mNext->mNext));

		if (iv->mColumnIndex == cix)
		{
			oldValue = iv->mText;
			break;
		}
	}
	
	if (oldValue != nullptr and value == oldValue)	// no need to update
		return;

	string oldStrValue = oldValue ? oldValue : "";

	// check the value
	if (col.mValidator)
		(*col.mValidator)(value);

	// If the field is part of the Key for this Category, remove it from the index
	// before updating
	
	bool reinsert = false;
	
	if (not emplacing and	// an update of an Item's value
		cat->mIndex != nullptr and cat->keyFields().count(name))
	{
		reinsert = cat->mIndex->find(mData);
		if (reinsert)
			cat->mIndex->erase(mData);
	}

	// first remove old value with cix

	if (mData->mValues == nullptr)
		;	// nothing to do
	else if (mData->mValues->mColumnIndex == cix)
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
			if (iv->mNext->mColumnIndex == cix)
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
		auto nv = new(value.length()) ItemValue(value.c_str(), cix);
	
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
	if (not emplacing and iv != nullptr and not iv->mChildren.empty())
	{
		for (auto child: iv->mChildren)
		{
			if (child->mCategory == nullptr)
				continue;
			
			auto childCat = db.get(child->mCategory->mName);
			if (childCat == nullptr)
				continue;

#if DEBUG
cerr << "fixing linked item " << child->mCategory->mName << '.' << child->mTag << endl;
#endif
				
			auto rows = childCat->find(Key(child->mTag) == oldStrValue);
			for (auto& cr: rows)
				cr.assign(child->mTag, value, false);
		}
	}
}

void Row::swap(const string& name, ItemRow* a, ItemRow* b)
{
	if (a == nullptr or b == nullptr)
		throw logic_error("invalid Rows in swap");

	assert(a->mCategory == b->mCategory);
	if (a->mCategory != b->mCategory)
		throw logic_error("Categories not same in swap");
	
	auto cat = a->mCategory;
	auto cix = cat->addColumn(name);
	auto& col = cat->mColumns[cix];
	auto& db = cat->mDb;

	// If the field is part of the Key for this Category, remove it from the index
	// before updating
	
	bool reinsert = false;
	
	if (cat->mIndex != nullptr and cat->keyFields().count(name))
	{
		reinsert = true;
		cat->mIndex->erase(a);
		cat->mIndex->erase(b);
	}

	ItemValue* ap = nullptr;
	ItemValue* ai = nullptr;
	ItemValue* bp = nullptr;
	ItemValue* bi = nullptr;
	
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

	// see if we need to update any child categories that depend on these values
	auto iv = col.mValidator;
	if ((ai != nullptr or bi != nullptr) and
		iv != nullptr and not iv->mChildren.empty())
	{
		for (auto child: iv->mChildren)
		{
			if (child->mCategory == nullptr)
				continue;
			
			auto childCat = db.get(child->mCategory->mName);
			if (childCat == nullptr)
				continue;

#if DEBUG
cerr << "fixing linked item " << child->mCategory->mName << '.' << child->mTag << endl;
#endif
			if (ai != nullptr)
			{
				auto rows = childCat->find(Key(child->mTag) == string(ai->mText));
				for (auto& cr: rows)
					cr.assign(child->mTag, bi == nullptr ? "" : bi->mText, false);
			}

			if (bi != nullptr)
			{
				auto rows = childCat->find(Key(child->mTag) == string(bi->mText));
				for (auto& cr: rows)
					cr.assign(child->mTag, ai == nullptr ? "" : ai->mText, false);
			}
		}
	}
}


void Row::assign(const Item& value, bool emplacing)
{
	assign(value.name(), value.value(), emplacing);
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

uint32 Row::lineNr() const
{
	return mData ? mData->mLineNr : 0;
}

void Row::lineNr(uint32 l)
{
	if (mData)
		mData->mLineNr = l;
}

Row::const_iterator::const_iterator(ItemRow* data, ItemValue* ptr)
	: mData(data), mPtr(ptr)
{
	if (mPtr != nullptr)
		fetch();
}

Row::const_iterator& Row::const_iterator::operator++()
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

// --------------------------------------------------------------------

File::File()
	: mHead(nullptr)
	, mValidator(nullptr)
{
}

File::File(istream& is, bool validate)
	: File()
{
//	parser p(is, *this);
//	p.parseFile();
	load(is);
}

File::File(boost::filesystem::path p, bool validate)
	: File()
{
	load(p);
}

File::File(File&& rhs)
	: mHead(nullptr), mValidator(nullptr)
{
	std::swap(mHead, rhs.mHead);
	std::swap(mValidator, rhs.mValidator);
}

File::~File()
{
	delete mHead;
	delete mValidator;
}

void File::append(Datablock* e)
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

void File::load(fs::path p)
{
	fs::ifstream inFile(p, ios_base::in | ios_base::binary);
	if (not inFile.is_open())
		throw runtime_error("Could not open file: " + p.string());
	
	io::filtering_stream<io::input> in;
	string ext;
	
	if (p.extension() == ".bz2")
	{
		in.push(io::bzip2_decompressor());
		ext = p.stem().extension().string();
	}
	else if (p.extension() == ".gz")
	{
		in.push(io::gzip_decompressor());
		ext = p.stem().extension().string();
	}
	
	in.push(inFile);

	load(in);
}

void File::save(fs::path p)
{
	fs::ofstream outFile(p, ios_base::out | ios_base::binary);
	io::filtering_stream<io::output> out;
	
	if (p.extension() == ".gz")
	{
		out.push(io::gzip_compressor());
		p = p.stem();
	}
	else if (p.extension() == ".bz2")
	{
		out.push(io::bzip2_compressor());
		p = p.stem();
	}
	
	out.push(outFile);
	save(out);
}

void File::load(istream& is)
{
	Validator* saved = mValidator;
	setValidator(nullptr);

	Parser p(is, *this);
	p.parseFile();
	
	if (saved != nullptr)
	{
		setValidator(saved);
		(void)isValid();
	}
}

void File::save(ostream& os)
{
	Datablock* e = mHead;
	while (e != nullptr)
	{
		e->write(os);
		e = e->mNext;
	}
}

void File::write(ostream& os, const vector<string>& order)
{
	Datablock* e = mHead;
	while (e != nullptr)
	{
		e->write(os, order);
		e = e->mNext;
	}
}

Datablock* File::get(const string& name) const
{
	const Datablock* result = mHead;
	while (result != nullptr and not iequals(result->mName, name))
		result = result->mNext;
	return const_cast<Datablock*>(result);
}

Datablock& File::operator[](const string& name)
{
	Datablock* result = mHead;
	while (result != nullptr and not iequals(result->mName, name))
		result = result->mNext;
	if (result == nullptr)
		throw runtime_error("Datablock " + name + " does not exist");
	return *result;
}

bool File::isValid()
{
	if (mValidator == nullptr)
	{
		if (VERBOSE)
			cerr << "No dictionary loaded explicitly, loading default" << endl;
		
		loadDictionary();
	}

	bool result = true;
	for (auto d = mHead; d != nullptr; d = d->mNext)
		result = d->isValid() and result;
	return result;
}

const Validator& File::getValidator() const
{
	if (mValidator == nullptr)
		throw runtime_error("no Validator defined yet");
	return *mValidator;
}

void File::loadDictionary()
{
	loadDictionary("mmcif_ddl");
}

void File::loadDictionary(const char* dict)
{
	for (;;)
	{
		string name(dict);
		
		if (fs::exists(name))
		{
			fs::ifstream is(name);
			loadDictionary(is);
			break;
		}
		
		fs::path dictFile = string("dictionaries/") + dict + ".dic";
		if (fs::exists(dictFile))
		{
			fs::ifstream is(dictFile);
			loadDictionary(is);
			break;
		}
		
#if defined(USE_RSRC)
		mrsrc::rsrc dictData(dictFile.string());
	
		if (dictData)
		{
			struct membuf : public streambuf
			{
				membuf(char* dict, size_t length)
				{
					this->setg(dict, dict, dict + length);
				}
			} buffer(const_cast<char*>(dictData.data()), dictData.size());
			
			istream is(&buffer);
			
			loadDictionary(is);
			break;
		}
#endif
		
		throw runtime_error("Dictionary not found or defined (" + name + ")");
	}
}

void File::loadDictionary(istream& is)
{
	unique_ptr<Validator> v(new Validator());

	DictParser p(*v, is);
	p.loadDictionary();

	setValidator(v.release());
}

void File::setValidator(Validator* v)
{
	mValidator = v;

	for (auto d = mHead; d != nullptr; d = d->mNext)
		d->setValidator(mValidator);
}

void File::getTagOrder(vector<string>& tags) const
{
	for (auto d = mHead; d != nullptr; d = d->mNext)
		d->getTagOrder(tags);
}

auto File::iterator::operator++() -> iterator&
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

}
