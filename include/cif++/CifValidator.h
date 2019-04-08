// cif parsing library

#include "cif++/Cif++.h"

#include <boost/filesystem/path.hpp>

//// the std regex of gcc is crashing....
//#include <boost/regex.hpp>

#include <regex>
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

enum DDL_PrimitiveType
{
	ptChar, ptUChar, ptNumb
};

DDL_PrimitiveType mapToPrimitiveType(const std::string& s);

struct ValidateType
{
	std::string				mName;
	DDL_PrimitiveType		mPrimitiveType;
	std::regex				mRx;

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
	std::string				mName;
	std::vector<string>		mKeys;
	cif::iset				mGroups;
	cif::iset				mMandatoryFields;
	std::set<ValidateItem>	mItemValidators;

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
	ValidateItem*			mParent;
	ValidateItem*			mChild;
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
	std::vector<ValidateLink> getLinksForParent(std::string category, std::string item) const;
	std::vector<ValidateLink> getLinksForChild(std::string category, std::string item) const;

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
//	std::set<uint32>			mSubCategories;
	std::set<ValidateType>		mTypeValidators;
	std::set<ValidateCategory>	mCategoryValidators;
	std::vector<ValidateLink>	mLinkValidators;
};

}
