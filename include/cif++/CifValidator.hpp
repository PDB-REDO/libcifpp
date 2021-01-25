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

#pragma once

#include "cif++/Cif++.hpp"

// duh.. https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86164
// #include <regex>
#include <boost/regex.hpp>

#include <set>

namespace cif
{
	
struct ValidateCategory;

// --------------------------------------------------------------------

class ValidationError : public std::exception
{
  public:
	ValidationError(const std::string& msg);
	ValidationError(const std::string& cat, const std::string& item,
		const std::string& msg);
	const char* what() const noexcept		{ return mMsg.c_str(); }
	std::string mMsg;
};

// --------------------------------------------------------------------

enum class DDL_PrimitiveType
{
	Char, UChar, Numb
};

DDL_PrimitiveType mapToPrimitiveType(const std::string& s);

struct ValidateType
{
	std::string				mName;
	DDL_PrimitiveType		mPrimitiveType;
	// std::regex				mRx;
	boost::regex			mRx;

	bool operator<(const ValidateType& rhs) const
	{
		return icompare(mName, rhs.mName) < 0;
	}

	// compare values based on type	
//	int compare(const std::string& a, const std::string& b) const
//	{
//		return compare(a.c_str(), b.c_str());
//	}
	
	int compare(const char* a, const char* b) const;
};

struct ValidateItem
{
	std::string				mTag;
	bool					mMandatory;
	const ValidateType*		mType;
	cif::iset				mEnums;
	std::string				mDefault;
	bool					mDefaultIsNull;
	ValidateCategory*		mCategory = nullptr;

	// ItemLinked is used for non-key links
	struct ItemLinked
	{
		ValidateItem*		mParent;
		std::string			mParentItem;
		std::string			mChildItem;
	};

	std::vector<ItemLinked>	mLinked;
	
	bool operator<(const ValidateItem& rhs) const
	{
		return icompare(mTag, rhs.mTag) < 0;
	}

	bool operator==(const ValidateItem& rhs) const
	{
		return iequals(mTag, rhs.mTag);
	}

	void operator()(std::string value) const;
};

struct ValidateCategory
{
	std::string					mName;
	std::vector<std::string>	mKeys;
	cif::iset					mGroups;
	cif::iset					mMandatoryFields;
	std::set<ValidateItem>		mItemValidators;

	bool operator<(const ValidateCategory& rhs) const
	{
		return icompare(mName, rhs.mName) < 0;
	}

	void addItemValidator(ValidateItem&& v);
	
	const ValidateItem* getValidatorForItem(std::string tag) const;
	
	const std::set<ValidateItem>& itemValidators() const
	{
		return mItemValidators;
	}
};

struct ValidateLink
{
	int							mLinkGroupID;
	std::string					mParentCategory;
	std::vector<std::string>	mParentKeys;
	std::string					mChildCategory;
	std::vector<std::string>	mChildKeys;
	std::string					mLinkGroupLabel;
};

// --------------------------------------------------------------------

class Validator
{
  public:
	friend class DictParser;

	Validator();
	~Validator();

	Validator(const Validator& rhs) = delete;
	Validator& operator=(const Validator& rhs) = delete;
	
	Validator(Validator&& rhs);
	Validator& operator=(Validator&& rhs);
	
	void addTypeValidator(ValidateType&& v);
	const ValidateType* getValidatorForType(std::string typeCode) const;

	void addCategoryValidator(ValidateCategory&& v);
	const ValidateCategory* getValidatorForCategory(std::string category) const;

	void addLinkValidator(ValidateLink&& v);
	std::vector<const ValidateLink*> getLinksForParent(const std::string& category) const;
	std::vector<const ValidateLink*> getLinksForChild(const std::string& category) const;

	void reportError(const std::string& msg, bool fatal);
	
	std::string dictName() const					{ return mName; }
	void dictName(const std::string& name)			{ mName = name; }

	std::string dictVersion() const				{ return mVersion; }
	void dictVersion(const std::string& version)	{ mVersion = version; }

  private:

	// name is fully qualified here:
	ValidateItem* getValidatorForItem(std::string name) const;

	std::string					mName;
	std::string					mVersion;
	bool						mStrict = false;
//	std::set<uint32_t>			mSubCategories;
	std::set<ValidateType>		mTypeValidators;
	std::set<ValidateCategory>	mCategoryValidators;
	std::vector<ValidateLink>	mLinkValidators;
};

}
