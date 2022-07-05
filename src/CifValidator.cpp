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

#include <fstream>
#include <filesystem>

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "cif++/Cif++.hpp"
#include "cif++/CifParser.hpp"
#include "cif++/CifValidator.hpp"

namespace ba = boost::algorithm;
namespace fs = std::filesystem;
namespace io = boost::iostreams;

extern int VERBOSE;

namespace cif
{

ValidationError::ValidationError(const std::string &msg)
	: mMsg(msg)
{
}

ValidationError::ValidationError(const std::string &cat, const std::string &item, const std::string &msg)
	: mMsg("When validating _" + cat + '.' + item + ": " + msg)
{
}

// --------------------------------------------------------------------

DDL_PrimitiveType mapToPrimitiveType(std::string_view s)
{
	DDL_PrimitiveType result;
	if (iequals(s, "char"))
		result = DDL_PrimitiveType::Char;
	else if (iequals(s, "uchar"))
		result = DDL_PrimitiveType::UChar;
	else if (iequals(s, "numb"))
		result = DDL_PrimitiveType::Numb;
	else
		throw ValidationError("Not a known primitive type");
	return result;
}

// --------------------------------------------------------------------

int ValidateType::compare(const char *a, const char *b) const
{
	int result = 0;

	if (*a == 0)
		result = *b == 0 ? 0 : -1;
	else if (*b == 0)
		result = *a == 0 ? 0 : +1;
	else
	{
		try
		{
			switch (mPrimitiveType)
			{
				case DDL_PrimitiveType::Numb:
				{
					double da = strtod(a, nullptr);
					double db = strtod(b, nullptr);

					auto d = da - db;
					if (std::abs(d) > std::numeric_limits<double>::epsilon())
					{
						if (d > 0)
							result = 1;
						else if (d < 0)
							result = -1;
					}
					break;
				}

				case DDL_PrimitiveType::UChar:
				case DDL_PrimitiveType::Char:
				{
					// CIF is guaranteed to have ascii only, therefore this primitive code will do
					// also, we're collapsing spaces

					auto ai = a, bi = b;
					for (;;)
					{
						if (*ai == 0)
						{
							if (*bi != 0)
								result = -1;
							break;
						}
						else if (*bi == 0)
						{
							result = 1;
							break;
						}

						char ca = *ai;
						char cb = *bi;

						if (mPrimitiveType == DDL_PrimitiveType::UChar)
						{
							ca = tolower(ca);
							cb = tolower(cb);
						}

						result = ca - cb;

						if (result != 0)
							break;

						if (ca == ' ')
						{
							while (ai[1] == ' ')
								++ai;
							while (bi[1] == ' ')
								++bi;
						}

						++ai;
						++bi;
					}

					break;
				}
			}
		}
		catch (const std::invalid_argument &ex)
		{
			result = 1;
		}
	}

	return result;
}

// --------------------------------------------------------------------

//void ValidateItem::addLinked(ValidateItem* parent, const std::string& parentItem, const std::string& childItem)
//{
////	if (mParent != nullptr and VERBOSE)
////		cerr << "replacing parent in " << mCategory->mName << " from " << mParent->mCategory->mName << " to " << parent->mCategory->mName << endl;
////	mParent = parent;
//
//	if (mType == nullptr and parent != nullptr)
//		mType = parent->mType;
//
//	if (parent != nullptr)
//	{
//		mLinked.push_back({parent, parentItem, childItem});
//
//		parent->mChildren.insert(this);
////
////		if (mCategory->mKeys == std::vector<std::string>{mTag})
////			parent->mForeignKeys.insert(this);
//	}
//}

void ValidateItem::operator()(std::string value) const
{
	if (not value.empty() and value != "?" and value != ".")
	{
		if (mType != nullptr and not regex_match(value, mType->mRx))
			throw ValidationError(mCategory->mName, mTag, "Value '" + value + "' does not match type expression for type " + mType->mName);

		if (not mEnums.empty())
		{
			if (mEnums.count(value) == 0)
				throw ValidationError(mCategory->mName, mTag, "Value '" + value + "' is not in the list of allowed values");
		}
	}
}

// --------------------------------------------------------------------

void ValidateCategory::addItemValidator(ValidateItem &&v)
{
	if (v.mMandatory)
		mMandatoryFields.insert(v.mTag);

	v.mCategory = this;

	auto r = mItemValidators.insert(std::move(v));
	if (not r.second and VERBOSE >= 4)
		std::cout << "Could not add validator for item " << v.mTag << " to category " << mName << std::endl;
}

const ValidateItem *ValidateCategory::getValidatorForItem(std::string_view tag) const
{
	const ValidateItem *result = nullptr;
	auto i = mItemValidators.find(ValidateItem{std::string(tag)});
	if (i != mItemValidators.end())
		result = &*i;
	else if (VERBOSE > 4)
		std::cout << "No validator for tag " << tag << std::endl;
	return result;
}

// --------------------------------------------------------------------

Validator::Validator(std::string_view name, std::istream &is)
	: mName(name)
{
	DictParser p(*this, is);
	p.loadDictionary();
}

Validator::~Validator()
{
}

void Validator::addTypeValidator(ValidateType &&v)
{
	auto r = mTypeValidators.insert(std::move(v));
	if (not r.second and VERBOSE > 4)
		std::cout << "Could not add validator for type " << v.mName << std::endl;
}

const ValidateType *Validator::getValidatorForType(std::string_view typeCode) const
{
	const ValidateType *result = nullptr;

	auto i = mTypeValidators.find(ValidateType{std::string(typeCode), DDL_PrimitiveType::Char, boost::regex()});
	if (i != mTypeValidators.end())
		result = &*i;
	else if (VERBOSE > 4)
		std::cout << "No validator for type " << typeCode << std::endl;
	return result;
}

void Validator::addCategoryValidator(ValidateCategory &&v)
{
	auto r = mCategoryValidators.insert(std::move(v));
	if (not r.second and VERBOSE > 4)
		std::cout << "Could not add validator for category " << v.mName << std::endl;
}

const ValidateCategory *Validator::getValidatorForCategory(std::string_view category) const
{
	const ValidateCategory *result = nullptr;
	auto i = mCategoryValidators.find(ValidateCategory{std::string(category)});
	if (i != mCategoryValidators.end())
		result = &*i;
	else if (VERBOSE > 4)
		std::cout << "No validator for category " << category << std::endl;
	return result;
}

ValidateItem *Validator::getValidatorForItem(std::string_view tag) const
{
	ValidateItem *result = nullptr;

	std::string cat, item;
	std::tie(cat, item) = splitTagName(tag);

	auto *cv = getValidatorForCategory(cat);
	if (cv != nullptr)
		result = const_cast<ValidateItem *>(cv->getValidatorForItem(item));

	if (result == nullptr and VERBOSE > 4)
		std::cout << "No validator for item " << tag << std::endl;

	return result;
}

void Validator::addLinkValidator(ValidateLink &&v)
{
	assert(v.mParentKeys.size() == v.mChildKeys.size());
	if (v.mParentKeys.size() != v.mChildKeys.size())
		throw std::runtime_error("unequal number of keys for parent and child in link");

	auto pcv = getValidatorForCategory(v.mParentCategory);
	auto ccv = getValidatorForCategory(v.mChildCategory);

	if (pcv == nullptr)
		throw std::runtime_error("unknown parent category " + v.mParentCategory);

	if (ccv == nullptr)
		throw std::runtime_error("unknown child category " + v.mChildCategory);

	for (size_t i = 0; i < v.mParentKeys.size(); ++i)
	{
		auto piv = pcv->getValidatorForItem(v.mParentKeys[i]);

		if (piv == nullptr)
			throw std::runtime_error("unknown parent tag _" + v.mParentCategory + '.' + v.mParentKeys[i]);

		auto civ = ccv->getValidatorForItem(v.mChildKeys[i]);
		if (civ == nullptr)
			throw std::runtime_error("unknown child tag _" + v.mChildCategory + '.' + v.mChildKeys[i]);

		if (civ->mType == nullptr and piv->mType != nullptr)
			const_cast<ValidateItem *>(civ)->mType = piv->mType;
	}

	mLinkValidators.emplace_back(std::move(v));
}

std::vector<const ValidateLink *> Validator::getLinksForParent(std::string_view category) const
{
	std::vector<const ValidateLink *> result;

	for (auto &l : mLinkValidators)
	{
		if (l.mParentCategory == category)
			result.push_back(&l);
	}

	return result;
}

std::vector<const ValidateLink *> Validator::getLinksForChild(std::string_view category) const
{
	std::vector<const ValidateLink *> result;

	for (auto &l : mLinkValidators)
	{
		if (l.mChildCategory == category)
			result.push_back(&l);
	}

	return result;
}

void Validator::reportError(const std::string &msg, bool fatal) const
{
	if (mStrict or fatal)
		throw ValidationError(msg);
	else if (VERBOSE > 0)
		std::cerr << msg << std::endl;
}

// --------------------------------------------------------------------

ValidatorFactory ValidatorFactory::sInstance;

ValidatorFactory::ValidatorFactory()
{
}

const Validator &ValidatorFactory::operator[](std::string_view dictionary)
{
	std::lock_guard lock(mMutex);

	for (auto &validator : mValidators)
	{
		if (iequals(validator.mName, dictionary))
			return validator;
	}

	// not found, add it

	// too bad clang version 10 did not have a constructor for fs::path that accepts a std::string_view
	fs::path dict_name(dictionary.data(), dictionary.data() + dictionary.length());

	auto data = loadResource(dict_name);

	if (not data and dict_name.extension().string() != ".dic")
		data = loadResource(dict_name.parent_path() / (dict_name.filename().string() + ".dic"));

	if (data)
		mValidators.emplace_back(dictionary, *data);
	else
	{
		std::error_code ec;

		// might be a compressed dictionary on disk
		fs::path p = dict_name;
		if (p.extension() == ".dic")
			p = p.parent_path() / (p.filename().string() + ".gz");
		else
			p = p.parent_path() / (p.filename().string() + ".dic.gz");

#if defined(CACHE_DIR) and defined(DATA_DIR)
		if (not fs::exists(p, ec) or ec)
		{
			for (const char *dir : {CACHE_DIR, DATA_DIR})
			{
				auto p2 = fs::path(dir) / p;
				if (fs::exists(p2, ec) and not ec)
				{
					swap(p, p2);
					break;
				}
			}
		}
#endif

		if (fs::exists(p, ec) and not ec)
		{
			std::ifstream file(p, std::ios::binary);
			if (not file.is_open())
				throw std::runtime_error("Could not open dictionary (" + p.string() + ")");

			io::filtering_stream<io::input> in;
			in.push(io::gzip_decompressor());
			in.push(file);

			mValidators.emplace_back(dictionary, in);
		}
		else
			throw std::runtime_error("Dictionary not found or defined (" + dict_name.string() + ")");
	}

	assert(iequals(mValidators.back().mName, dictionary));

	return mValidators.back();
}

} // namespace cif
