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

#include "cif++/Config.hpp"

#include <string>

#include <regex>
#include <iostream>
#include <sstream>
#include <set>
#include <list>

#include "cif++/CifUtils.hpp"

/*
	Simple C++ interface to CIF files.
	
	Assumptions: a file contains one or more datablocks modelled by the class datablock.
	Each datablock contains categories. These map to the original tables used to fill
	the mmCIF file. Each Category can contain multiple Items, the columns in the table.
	
	Values are stored as character strings internally.
	
	Synopsis:
	
	// create a cif file
	
	cif::datablock e("1MVE");
	e.append(cif::Category{"_entry", { "id", "1MVE" } });
	
	cif::Category atomSite("atom_site");
	size_t nr{};
	for (auto& myAtom: atoms)
	{
		atomSite.push_back({
			{ "group_PDB", "ATOM" },
			{ "id", ++nr },
			{ "type_symbol", myAtom.type.str() },
			...
		});
	}
	
	e.append(move(atomSite));
	
	cif::File f;
	f.append(e);
	
	ofstream os("1mve.cif");
	f.write(os);

	// read
	f.read(ifstream{"1mve.cif"});
	
	auto& e = f.firstDatablock();
	
	cout << "ID of datablock: " << e.id() << endl;
	
	auto& atoms = e["atom_site"];
	for (auto& atom: atoms)
	{
		cout << atom["group_PDB"] << ", "
			 << atom["id"] << ", "
			 ...

		float x, y, z;
		cif::tie(x, y, z) = atom.get("Cartn_x", "Cartn_y", "Cartn_z");
		...
	}

	Another way of querying a Category is by using this construct:
	
	auto cat& = e["atom_site"];
	auto Rows = cat.find(Key("label_asym_id") == "A" and Key("label_seq_id") == 1);


*/

namespace cif
{

// flag for verbose output
extern int VERBOSE;

// mmCIF mapping
// A CIF data file in this case contains entries (data blocks) which can contain
// one or more Category objects. Each Category object contains arrays of Items.
// Better, you can consider the categories as tables containing columns which
// are the Items.

class File;
class Datablock;
class Category;
class Row;			// a flyweight class that references data in categories
class RowSet;
class Item;
class Validator;

struct ValidateItem;
struct ValidateCategory;
struct ValidateLink;

struct ItemColumn;
struct ItemRow;
struct ItemValue;

// --------------------------------------------------------------------
// class Item
//
//	This class is only transient, it is used to construct new Rows.
//	Access to already stored data is through an ItemReference object.

class Item
{
  public:
	Item() {}

	template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	Item(const std::string& name, const T& value)
		: mName(name), mValue(std::to_string(value)) {}

	Item(const std::string& name, const std::string& value)
		: mName(name), mValue(value) {}

	Item(const Item& rhs) : mName(rhs.mName), mValue(rhs.mValue) {}
	Item(Item&& rhs) : mName(std::move(rhs.mName)), mValue(std::move(rhs.mValue)) {}

	Item& operator=(const Item& rhs)
	{
		if (this != &rhs)
		{
			mName = rhs.mName;
			mValue = rhs.mValue;
		}
		
		return *this;
	}
	
	Item& operator=(Item&& rhs)
	{
		if (this != &rhs)
		{
			mName = std::move(rhs.mName);
			mValue = std::move(rhs.mValue);
		}
		
		return *this;
	}
	
	const std::string& name() const	{ return mName; }
	const std::string& value() const	{ return mValue; }

	void value(const std::string& v)	{ mValue = v; }
	
	// empty means either null or unknown
	bool empty() const;

	// is_null means the field contains '.'
	bool is_null() const;

	// is_unknown means the field contains '?'
	bool is_unknown() const;

	size_t length() const		{ return mValue.length(); }
	const char* c_str() const	{ return mValue.c_str(); }
	
  private:
	std::string	mName;
  	std::string	mValue;
};

// --------------------------------------------------------------------
// class datablock acts as an STL container for Category objects

class Datablock
{
  public:
	friend class File;
	
	using CategoryList = std::list<Category>;
	using iterator = CategoryList::iterator;
	using const_iterator = CategoryList::const_iterator;
	
	Datablock(const std::string& name);
	~Datablock();

	Datablock(const Datablock&) = delete;
	Datablock& operator=(const Datablock&) = delete;

	std::string getName() const							{ return mName; }
	void setName(const std::string& n)					{ mName = n; }
	
	std::string firstItem(const std::string& tag) const;

	iterator begin()		{ return mCategories.begin(); }
	iterator end()			{ return mCategories.end(); }

	const_iterator begin() const	{ return mCategories.begin(); }
	const_iterator end() const		{ return mCategories.end(); }

	Category& operator[](const std::string& name);

	std::tuple<iterator,bool> emplace(const std::string& name);
	
	bool isValid();
	void validateLinks() const;

	void setValidator(Validator* v);

	// this one only looks up a Category, returns nullptr if it does not exist
	Category* get(const std::string& name);

	void getTagOrder(std::vector<std::string>& tags) const;
	void write(std::ostream& os, const std::vector<std::string>& order);
	void write(std::ostream& os);
	
	// convenience function, add a line to the software category
	void add_software(const std::string& name, const std::string& classification,
		const std::string& versionNr, const std::string& versionDate);

  private:

	std::list<Category>	mCategories;
	std::string			mName;
	Validator*			mValidator;
	Datablock*			mNext;
};

// --------------------------------------------------------------------
// class Row acts as a container for Item objects, It has a more useful
// interface for accessing the contained columns. The get() method
// returns a RowResult object that can be used to access only a subset
// of column values by index or by name.

namespace detail
{
	// ItemReference is a helper class
	class ItemReference
	{
	  public:

		// conversion helper class
		template<typename T, typename = void>
		struct item_value_as;

		template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
		ItemReference& operator=(const T& value)
		{
			this->operator=(std::to_string(value));
			return *this;
		}

		ItemReference& operator=(const std::string& value);

		template<typename... Ts>
		void os(const Ts& ... v)
		{
			std::ostringstream ss;
			((ss << v), ...);
			this->operator=(ss.str());
		}

		void swap(ItemReference& b);
		
		template<typename T>
		T as() const
		{
			return item_value_as<T>::convert(*this);
		}

		template<typename T>
		int compare(const T& value, bool icase) const
		{
			return item_value_as<T>::compare(*this, value, icase);
		}

		// empty means either null or unknown
		bool empty() const;
		explicit operator bool() const		{ return not empty(); }

		// is_null means the field contains '.'
		bool is_null() const;

		// is_unknown means the field contains '?'
		bool is_unknown() const;

		const char* c_str() const;
		
		// the following returns the defaultValue from either the parameter
		// or, if specified, the value from _item_default.value in the dictionary
		const char* c_str(const char* defaultValue) const;
		
		bool operator!=(const std::string& s) const		{ return s != c_str(); }
		bool operator==(const std::string& s) const		{ return s == c_str(); }

	private:
		friend class ::cif::Row;

		ItemReference(const char* name, size_t column, Row& row)
			: mName(name), mColumn(column), mRow(row) {}

		ItemReference(const char* name, size_t column, const Row& row)
			: mName(name), mColumn(column), mRow(const_cast<Row&>(row)), mConst(true) {}

		const char*		mName;
		size_t			mColumn;
		Row&			mRow;
		bool			mConst = false;
	};

	template<typename T>
	struct ItemReference::item_value_as<T, std::enable_if_t<std::is_floating_point_v<T>>>
	{
		static T convert(const ItemReference& ref)
		{
			T result = {};
			if (not ref.empty())
				result = static_cast<T>(std::stod(ref.c_str()));
			return result;
		}

		static int compare(const ItemReference& ref, double value, bool icase)
		{
			int result = 0;
			try
			{
				double v = std::stod(ref.c_str());
				if (v < value)
					result = -1;
				else if (v > value)
					result = 1;
			}
			catch (...)
			{
				if (VERBOSE)
					std::cerr << "conversion error in compare for '" << ref.c_str() << '\'' << std::endl;
				result = 1;
			}
			
			return result;
		}
	};

	template<typename T>
	struct ItemReference::item_value_as<T, std::enable_if_t<std::is_integral_v<T> and std::is_unsigned_v<T>>>
	{
		static T convert(const ItemReference& ref)
		{
			T result = {};
			if (not ref.empty())
				result = static_cast<T>(std::stoul(ref.c_str()));
			return result;
		}

		static int compare(const ItemReference& ref, unsigned long value, bool icase)
		{
			int result = 0;
			try
			{
				auto v = std::stoul(ref.c_str());
				if (v < value)
					result = -1;
				else if (v > value)
					result = 1;
			}
			catch (...)
			{
				if (VERBOSE)
					std::cerr << "conversion error in compare for '" << ref.c_str() << '\'' << std::endl;
				result = 1;
			}
			
			return result;
		}
	};

	template<typename T>
	struct ItemReference::item_value_as<T, std::enable_if_t<std::is_integral_v<T> and std::is_signed_v<T>>>
	{
		static T convert(const ItemReference& ref)
		{
			T result = {};
			if (not ref.empty())
				result = static_cast<T>(std::stol(ref.c_str()));
			return result;
		}

		static int compare(const ItemReference& ref, long value, bool icase)
		{
			int result = 0;
			try
			{
				auto v = std::stol(ref.c_str());
				if (v < value)
					result = -1;
				else if (v > value)
					result = 1;
			}
			catch (...)
			{
				if (VERBOSE)
					std::cerr << "conversion error in compare for '" << ref.c_str() << '\'' << std::endl;
				result = 1;
			}
			
			return result;
		}
	};

	template<typename T>
	struct ItemReference::item_value_as<std::optional<T>>
	{
		static std::optional<T> convert(const ItemReference& ref)
		{
			std::optional<T> result;
			if (ref)
				result = ref.as<T>();
			return result;
		}

		static int compare(const ItemReference& ref, std::optional<T> value, bool icase)
		{
			if (ref.empty() and not value)
				return 0;
			
			if (ref.empty())
				return -1;
			else if (not value)
				return 1;
			else
				return ref.compare(*value, icase);
		}
	};

	template<size_t N>
	struct ItemReference::item_value_as<char[N]>
	{
		static int compare(const ItemReference& ref, const char (&value)[N], bool icase)
		{
			return icase ? cif::icompare(ref.c_str(), value) : std::strcmp(ref.c_str(), value);
		}
	};

	template<>
	struct ItemReference::item_value_as<const char*>
	{
		static const char* convert(const ItemReference& ref)
		{
			return ref.c_str();
		}

		static int compare(const ItemReference& ref, const char* value, bool icase)
		{
			return icase ? cif::icompare(ref.c_str(), value) : std::strcmp(ref.c_str(), value);
		}
	};

	template<>
	struct ItemReference::item_value_as<std::string>
	{
		static std::string convert(const ItemReference& ref)
		{
			return ref.c_str();
		}

		static int compare(const ItemReference& ref, const std::string& value, bool icase)
		{
			return icase ? cif::icompare(ref.c_str(), value) : std::strcmp(ref.c_str(), value.c_str());
		}
	};

	// some helper classes to help create tuple result types
	template<typename... C>
	struct getRowResult
	{
		static constexpr size_t N = sizeof...(C);
		
		getRowResult(const Row& r, std::array<size_t, N>&& columns)
			: mRow(r), mColumns(std::move(columns))
		{
		}
	
		const ItemReference operator[](size_t ix) const
		{
			return mRow[mColumns[ix]];
		}
		
		template<typename... Ts, std::enable_if_t<N == sizeof...(Ts), int> = 0>
		operator std::tuple<Ts...>() const
		{
			return get<Ts...>(std::index_sequence_for<Ts...>{});
		}
	
		template<typename... Ts, std::size_t... Is>
		std::tuple<Ts...> get(std::index_sequence<Is...>) const
		{
			return std::tuple<Ts...>{mRow[mColumns[Is]].template as<Ts>()...};
		}

		const Row& mRow;
		std::array<size_t, N> mColumns;
	};
	
	// we want to be able to tie some variables to a RowResult, for this we use tiewraps
	template<typename... Ts>
	struct tieWrap
	{
		tieWrap(Ts... args) : mVal(args...) {}

		template<typename RR>
		void operator=(const RR&& rr)
		{
			// getRowResult will do the conversion, but only if the types
			// are compatible. That means the number of parameters to the get()
			// of the row should be equal to the number of items in the tuple
			// you are trying to tie.

            using RType = std::tuple<typename std::remove_reference<Ts>::type...>;

			mVal = static_cast<RType>(rr);
		}

		std::tuple<Ts...>	mVal;
	};
}

template<typename... Ts>
auto tie(Ts&... v) -> detail::tieWrap<Ts&...>
{
	return detail::tieWrap<Ts&...>(std::forward<Ts&>(v)...);
}

class Row
{
  public:
	friend class Category;
	friend class CatIndex;
	friend class RowComparator;
	friend class detail::ItemReference;
	friend class RowSet;

	Row()
		: mData(nullptr) {}

	Row(ItemRow* data)
		: mData(data) {}

	Row(const ItemRow* data)
		: mData(const_cast<ItemRow*>(data)), mCascade(false) {}

	Row(const Row& rhs);
	Row& operator=(const Row& rhs);

	Row(Row&& rhs);
	Row& operator=(Row&& rhs);

	~Row();

	void setCascading(bool cascade)
	{
		mCascade = cascade;
	}

	void next();	///< make this row point to the next ItemRow

	struct const_iterator : public std::iterator<std::forward_iterator_tag, const Item>
	{
		typedef std::iterator<std::forward_iterator_tag, Item>	baseType;
		typedef typename baseType::pointer						pointer;
		typedef typename baseType::reference					reference;
		
		const_iterator(ItemRow* data, ItemValue* ptr);
		
		reference operator*()								{ return mCurrent; }
		pointer operator->()								{ return &mCurrent; }
		
		const_iterator& operator++();
		const_iterator operator++(int)						{ const_iterator result(*this); this->operator++(); return result; } 

		bool operator==(const const_iterator& rhs) const	{ return mPtr == rhs.mPtr; } 
		bool operator!=(const const_iterator& rhs) const	{ return mPtr != rhs.mPtr; } 
		
	  private:

		void fetch();

	  	ItemRow*	mData;
		ItemValue*	mPtr;
		Item		mCurrent;
	};
	
	// checks for an initialized Row:
	explicit operator bool() const							{ return mData != nullptr; }
	
	// for debugging
	uint32_t lineNr() const;
	void lineNr(uint32_t l);
	
	bool empty() const;
	const_iterator begin() const;
	const_iterator end() const;

	// TODO: implement real const version?
	
	friend class detail::ItemReference;

	const detail::ItemReference operator[](size_t column) const
	{
		return detail::ItemReference("<anonymous column>", column, *this);
	}

	const detail::ItemReference operator[](const char* itemTag) const
	{
		size_t column = ColumnForItemTag(itemTag);
		return detail::ItemReference(itemTag, column, *this);
	}

	detail::ItemReference operator[](const char* itemTag)
	{
		size_t column = ColumnForItemTag(itemTag);
		return detail::ItemReference(itemTag, column, *this);
	}

	const detail::ItemReference operator[](const std::string& itemTag) const
	{
		size_t column = ColumnForItemTag(itemTag.c_str());
		return detail::ItemReference(itemTag.c_str(), column, *this);
	}

	detail::ItemReference operator[](const std::string& itemTag)
	{
		size_t column = ColumnForItemTag(itemTag.c_str());
		return detail::ItemReference(itemTag.c_str(), column, *this);
	}

	template<typename... C>
	auto get(C... columns) const -> detail::getRowResult<C...>
	{
		std::array<size_t,sizeof...(C)> cix;
		
		auto c = cix.begin();
		for (auto col: { columns... })
			*c++ = ColumnForItemTag(col);
		
		return detail::getRowResult<C...>(*this, std::move(cix));
	}
	
	void assign(const std::vector<Item>& values);

	bool operator==(const Row& rhs) const
	{
		return mData == rhs.mData;
	}

	bool operator!=(const Row& rhs) const
	{
		return mData != rhs.mData;
	}

	ItemRow* data() const							{ return mData; }

	void swap(Row& rhs)
	{
		std::swap(mData, rhs.mData);
	}

	friend std::ostream& operator<<(std::ostream& os, const Row& row);

  private:

	void assign(const std::string& name, const std::string& value, bool updateLinked);
	void assign(size_t column, const std::string& value, bool updateLinked);
	void assign(const Item& i, bool updateLinked);
	
	static void swap(size_t column, ItemRow* a, ItemRow* b);

	size_t ColumnForItemTag(const char* itemTag) const;

	ItemRow*	mData;
	uint32_t	mLineNr = 0;
	bool		mCascade = true;
};

// --------------------------------------------------------------------
// some more templates to be able to do querying

namespace detail
{

struct ConditionImpl
{
	virtual ~ConditionImpl() {}
	
	virtual void prepare(const Category& c) {}
	virtual bool test(const Category& c, const Row& r) const = 0;
	virtual void str(std::ostream& os) const = 0;
};

struct AllConditionImpl : public ConditionImpl
{
	virtual bool test(const Category& c, const Row& r) const { return true; }
	virtual void str(std::ostream& os) const { os << "*"; }
};

struct orConditionImpl;
struct andConditionImpl;
struct notConditionImpl;

}

class Condition
{
  public:

	Condition() : mImpl(nullptr) {}
	Condition(detail::ConditionImpl* impl) : mImpl(impl) {}

	Condition(const Condition&) = delete;

	Condition(Condition&& rhs)
		: mImpl(nullptr)
	{
		std::swap(mImpl, rhs.mImpl);
	}

	Condition& operator=(const Condition&) = delete;

	Condition& operator=(Condition&& rhs)
	{
		std::swap(mImpl, rhs.mImpl);
		return *this;
	}

	~Condition()
	{
		delete mImpl;
		mImpl = nullptr;
	}
	
	void prepare(const Category& c)
	{
		if (mImpl)
			mImpl->prepare(c);
		mPrepared = true;
	}
	
	bool operator()(const Category& c, const Row& r) const
	{
		assert(mImpl);
		assert(mPrepared);
		return mImpl ? mImpl->test(c, r) : false;
	}
	
	bool empty() const		{ return mImpl == nullptr; }

	friend Condition operator||(Condition&& a, Condition&& b);
	friend Condition operator&&(Condition&& a, Condition&& b);

	friend struct detail::orConditionImpl;
	friend struct detail::andConditionImpl;
	friend struct detail::notConditionImpl;

	void swap(Condition& rhs)
	{
		std::swap(mImpl, rhs.mImpl);
		std::swap(mPrepared, rhs.mPrepared);
	}

	friend std::ostream& operator<<(std::ostream& os, const Condition& cond);

  private:
	detail::ConditionImpl*	mImpl;
	bool mPrepared = false;
};

inline std::ostream& operator<<(std::ostream& os, const Condition& cond)
{
	if (cond.mImpl)
		cond.mImpl->str(os);
	return os;
}

namespace detail
{

struct KeyIsEmptyConditionImpl : public ConditionImpl
{
	KeyIsEmptyConditionImpl(const std::string& ItemTag)
		: mItemTag(ItemTag) {}

	virtual void prepare(const Category& c);
	
	virtual bool test(const Category& c, const Row& r) const
	{
		return r[mItemIx].empty();
	}

	virtual void str(std::ostream& os) const
	{
		os << mItemTag << " IS NULL";
	}

	std::string mItemTag;
	size_t mItemIx;
};

struct KeyCompareConditionImpl : public ConditionImpl
{
	template<typename COMP>
	KeyCompareConditionImpl(const std::string& ItemTag, COMP&& comp, const std::string& s)
		: mItemTag(ItemTag), mComp(std::move(comp)), mStr(s) {}
	
	virtual void prepare(const Category& c);
	
	virtual bool test(const Category& c, const Row& r) const
	{
		return mComp(c, r, mCaseInsensitive);
	}
	
	virtual void str(std::ostream& os) const
	{
		os << mItemTag << (mCaseInsensitive ? "^ " : " ") << mStr;
	}

	std::string mItemTag;
	size_t mItemIx;
	bool mCaseInsensitive = false;
	std::function<bool(const Category&, const Row&, bool)> mComp;
	std::string mStr;
};

struct KeyMatchesConditionImpl : public ConditionImpl
{
	KeyMatchesConditionImpl(const std::string& ItemTag, const std::regex& rx)
		: mItemTag(ItemTag), mRx(rx) {}
	
	virtual void prepare(const Category& c);
	
	virtual bool test(const Category& c, const Row& r) const
	{
		return std::regex_match(r[mItemIx].as<std::string>(), mRx);
	}

	virtual void str(std::ostream& os) const
	{
		os << mItemTag << " =~ expression";
	}
	
	std::string mItemTag;
	size_t mItemIx;
	std::regex mRx;
};

template<typename T>
struct anyIsConditionImpl : public ConditionImpl
{
	typedef T valueType;
	
	anyIsConditionImpl(const valueType& value)
		: mValue(value) {}
	
	virtual bool test(const Category& c, const Row& r) const;
	virtual void str(std::ostream& os) const
	{
		os << "<any> == " << mValue;
	}

	valueType mValue;
};

struct anyMatchesConditionImpl : public ConditionImpl
{
	anyMatchesConditionImpl(const std::regex& rx)
		: mRx(rx) {}
	
	virtual bool test(const Category& c, const Row& r) const;
	virtual void str(std::ostream& os) const
	{
		os << "<any> =~ expression";
	}

	std::regex mRx;
};

struct allConditionImpl : public ConditionImpl
{
	allConditionImpl() {}
	
	virtual bool test(const Category& c, const Row& r) const
    {
        return true;
    }

	virtual void str(std::ostream& os) const
	{
		os << "*";
	}
};

struct andConditionImpl : public ConditionImpl
{
	andConditionImpl(Condition&& a, Condition&& b)
		: mA(nullptr), mB(nullptr)
	{
		std::swap(mA, a.mImpl);
		std::swap(mB, b.mImpl);
	}
	
	~andConditionImpl()
	{
		delete mA;
		delete mB;
	}
	
	virtual void prepare(const Category& c)
	{
		mA->prepare(c);
		mB->prepare(c);
	}
	
	virtual bool test(const Category& c, const Row& r) const
	{
		return mA->test(c, r) and mB->test(c, r);
	}

	virtual void str(std::ostream& os) const
	{
		os << '(';
		mA->str(os);
		os << ") AND (";
		mB->str(os);
		os << ')';
	}

	ConditionImpl* mA;
	ConditionImpl* mB;
};

struct orConditionImpl : public ConditionImpl
{
	orConditionImpl(Condition&& a, Condition&& b)
		: mA(nullptr), mB(nullptr)
	{
		std::swap(mA, a.mImpl);
		std::swap(mB, b.mImpl);
	}
	
	~orConditionImpl()
	{
		delete mA;
		delete mB;
	}
	
	virtual void prepare(const Category& c)
	{
		mA->prepare(c);
		mB->prepare(c);
	}
	
	virtual bool test(const Category& c, const Row& r) const
	{
		return mA->test(c, r) or mB->test(c, r);
	}

	virtual void str(std::ostream& os) const
	{
		os << '(';
		mA->str(os);
		os << ") OR (";
		mB->str(os);
		os << ')';
	}

	ConditionImpl* mA;
	ConditionImpl* mB;
};

struct notConditionImpl : public ConditionImpl
{
	notConditionImpl(Condition&& a)
		: mA(nullptr)
	{
		std::swap(mA, a.mImpl);
	}
	
	~notConditionImpl()
	{
		delete mA;
	}
	
	virtual void prepare(const Category& c)
	{
		mA->prepare(c);
	}
	
	virtual bool test(const Category& c, const Row& r) const
	{
		return not mA->test(c, r);
	}

	virtual void str(std::ostream& os) const
	{
		os << "NOT (";
		mA->str(os);
		os << ')';
	}

	ConditionImpl* mA;
};

}

inline Condition operator&&(Condition&& a, Condition&& b)
{
	if (a.mImpl and b.mImpl)
		return Condition(new detail::andConditionImpl(std::move(a), std::move(b)));
	if (a.mImpl)
		return Condition(std::move(a));
	return Condition(std::move(b));
}

inline Condition operator||(Condition&& a, Condition&& b)
{
	if (a.mImpl and b.mImpl)
		return Condition(new detail::orConditionImpl(std::move(a), std::move(b)));
	if (a.mImpl)
		return Condition(std::move(a));
	return Condition(std::move(b));
}

struct Empty {};

struct Key
{
	Key(const std::string& itemTag) : mItemTag(itemTag) {}
	Key(const char* itemTag) : mItemTag(itemTag) {}

	Key(const Key&) = delete;
	Key& operator=(const Key&) = delete;
	
	std::string mItemTag;
};

template<typename T>
Condition operator==(const Key& key, const T& v)
{
	std::ostringstream s;
	s << "== " << v;

	return Condition(new detail::KeyCompareConditionImpl(key.mItemTag, [tag = key.mItemTag, v](const Category& c, const Row& r, bool icase)
		{ return r[tag].template compare<T>(v, icase) == 0; }, s.str()));
}

// inline Condition operator==(const Key& key, const detail::ItemReference& v)
// {
// 	if (v.empty())
// 		return Condition(new detail::KeyIsEmptyConditionImpl(key.mItemTag));
// 	else
// 		return Condition(new detail::KeyCompareConditionImpl(key.mItemTag, [tag = key.mItemTag, v](const Category& c, const Row& r, bool icase)
// 			{ return r[tag].template compare<(v, icase) == 0; }));
// }

inline Condition operator==(const Key& key, const Empty&)
{
	return Condition(new detail::KeyIsEmptyConditionImpl(key.mItemTag));
}

template<typename T>
Condition operator!=(const Key& key, const T& v)
{
	return Condition(new detail::notConditionImpl(operator==(key, v)));
}

inline Condition operator!=(const Key& key, const char* v)
{
	std::string value(v ? v : "");
	return Condition(new detail::notConditionImpl(operator==(key, value)));
}

template<typename T>
Condition operator>(const Key& key, const T& v)
{
	std::ostringstream s;
	s << ">" << v;

	return Condition(new detail::KeyCompareConditionImpl(key.mItemTag, [tag = key.mItemTag, v](const Category& c, const Row& r, bool icase)
		{ return r[tag].template compare<T>(v, icase) > 0; }, s.str()));
}

template<typename T>
Condition operator>=(const Key& key, const T& v)
{
	std::ostringstream s;
	s << ">=" << v;

	return Condition(new detail::KeyCompareConditionImpl(key.mItemTag, [tag = key.mItemTag, v](const Category& c, const Row& r, bool icase)
		{ return r[tag].template compare<T>(v, icase) >= 0; }, s.str()));
}

template<typename T>
Condition operator<(const Key& key, const T& v)
{
	std::ostringstream s;
	s << "<" << v;

	return Condition(new detail::KeyCompareConditionImpl(key.mItemTag, [tag = key.mItemTag, v](const Category& c, const Row& r, bool icase)
		{ return r[tag].template compare<T>(v, icase) < 0; }, s.str()));
}

template<typename T>
Condition operator<=(const Key& key, const T& v)
{
	std::ostringstream s;
	s << "<=" << v;

	return Condition(new detail::KeyCompareConditionImpl(key.mItemTag, [tag = key.mItemTag, v](const Category& c, const Row& r, bool icase)
		{ return r[tag].template compare<T>(v, icase) <= 0; }, s.str()));
}

template<>
inline
Condition operator==(const Key& key, const std::regex& rx)
{
	return Condition(new detail::KeyMatchesConditionImpl(key.mItemTag, rx));
}

template<>
inline
Condition operator==(const Key& key, const Empty&)
{
	return Condition(new detail::KeyIsEmptyConditionImpl(key.mItemTag));
}

struct any
{
	template<typename T>
	Condition operator==(const T& v) const
	{
		return Condition(new detail::anyIsConditionImpl<T>(v));
	}
};

template<>
inline
Condition any::operator==(const std::regex& rx) const
{
	return Condition(new detail::anyMatchesConditionImpl(rx));
}

inline Condition All()
{
    return Condition(new detail::allConditionImpl());
}

inline Condition Not(Condition&& cond)
{
    return Condition(new detail::notConditionImpl(std::move(cond)));
}

// -----------------------------------------------------------------------
// iterators

template<typename RowType>
class iterator_impl
{
  public:
	using iterator_category = std::forward_iterator_tag;
	using value_type = RowType;
	using difference_type = std::ptrdiff_t;
	using pointer = RowType*;
	using reference = RowType&;

	friend class Category;

	iterator_impl(ItemRow* data) : mCurrent(data) {}
	iterator_impl(const iterator_impl& i) : mCurrent(i.mCurrent) {}

	template<typename IteratorType>
	iterator_impl(IteratorType& i)
		: mCurrent(*i) {}

	iterator_impl& operator=(const iterator_impl& i)
	{
		mCurrent = i.mCurrent;
		return *this;
	}

	virtual ~iterator_impl() = default;

	reference operator*()							{ return mCurrent; }
	pointer operator->()							{ return &mCurrent; }
	
	iterator_impl& operator++()
	{
		mCurrent.next();
		return *this;
	}

	iterator_impl operator++(int)
	{
		iterator_impl result(*this);
		this->operator++();
		return result;
	} 

	bool operator==(const iterator_impl& rhs) const	{ return mCurrent == rhs.mCurrent; } 
	bool operator!=(const iterator_impl& rhs) const	{ return mCurrent != rhs.mCurrent; } 

  private:
	Row mCurrent;
};

// --------------------------------------------------------------------
// Iterator proxy to iterate over a subset of rows selected by a Condition

template<typename RowType>
class conditional_iterator_proxy
{
  public:

	using base_iterator = iterator_impl<RowType>;

	class conditional_iterator_impl
	{
	  public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = RowType;
		using difference_type = std::ptrdiff_t;
		using pointer = RowType*;
		using reference = RowType&;

		conditional_iterator_impl(Category& cat, base_iterator pos, const Condition& cond);
		conditional_iterator_impl(const conditional_iterator_impl& i) = default;
		conditional_iterator_impl& operator=(const conditional_iterator_impl& i) = default;

		virtual ~conditional_iterator_impl() = default;

		reference operator*()							{ return *mBegin; }
		pointer operator->()							{ return &*mBegin; }
		
		conditional_iterator_impl& operator++()
		{
			while (mBegin != mEnd)
			{
				if (++mBegin == mEnd)
					break;
				
				if ((*mCondition)(*mCat, *mBegin))
					break;
			}

			return *this;
		}

		conditional_iterator_impl operator++(int)
		{
			conditional_iterator_impl result(*this);
			this->operator++();
			return result;
		} 

		bool operator==(const conditional_iterator_impl& rhs) const	{ return mBegin == rhs.mBegin; } 
		bool operator!=(const conditional_iterator_impl& rhs) const	{ return mBegin != rhs.mBegin; } 

	  private:
	  	Category* mCat;
		base_iterator mBegin, mEnd;
		const Condition* mCondition;
	};

	using iterator = conditional_iterator_impl;
	using reference = typename iterator::reference;

	conditional_iterator_proxy(Category& cat, base_iterator pos, Condition&& cond);

	conditional_iterator_proxy(conditional_iterator_proxy&& p);
	conditional_iterator_proxy& operator=(conditional_iterator_proxy&& p);

	conditional_iterator_proxy(const conditional_iterator_proxy&) = delete;
	conditional_iterator_proxy& operator=(const conditional_iterator_proxy&) = delete;

	iterator begin() const;
	iterator end() const;

	bool empty() const;
	size_t size() const				{ return std::distance(begin(), end()); }

	RowType front()					{ return *begin(); }

	Category& category() const		{ return *mCat;}

	void swap(conditional_iterator_proxy& rhs);

  private:
	Category* mCat;
	Condition mCondition;
	base_iterator mCBegin, mCEnd;
};

// --------------------------------------------------------------------
// class RowSet is used to return find results. Use it to re-order the results
// or to group them 

class RowSet
{
	typedef std::vector<Row>	base_type;

  public:

	using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

	RowSet(Category& cat);
	RowSet(Category& cat, Condition&& cond);
	RowSet(const RowSet& rhs);
	RowSet(RowSet&& rhs);
	virtual ~RowSet();

	RowSet& operator=(const RowSet& rhs);
	RowSet& operator=(RowSet&& rhs);

	RowSet& orderBy(const std::string& Item)
		{ return orderBy({ Item }); }
	
	RowSet& orderBy(std::initializer_list<std::string> Items);

	class iterator
	{
	  public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = Row;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type*;
		using reference = value_type&;

		iterator() {}

		iterator(const iterator& i)
			: mPos(i.mPos)
			, mEnd(i.mEnd)
			, mCurrentRow(i.mCurrentRow) {}

		iterator(const std::vector<ItemRow*>::iterator& p, const std::vector<ItemRow*>::iterator& e)
			: mPos(p)
			, mEnd(e)
			, mCurrentRow(p == e ? nullptr : *p) {}

		iterator& operator=(const iterator& i)
		{
			mPos = i.mPos;
			mEnd = i.mEnd;
			mCurrentRow = i.mCurrentRow;
			return *this;
		}

		reference operator*()		{ return mCurrentRow; }
		pointer operator->()		{ return &mCurrentRow; }

		iterator& operator++()
		{
			++mPos;
			if (mPos != mEnd)
				mCurrentRow = Row(*mPos);
			else
				mCurrentRow = Row();
			return *this;
		}

		iterator operator++(int)
		{
			iterator t(*this);
			operator++();
			return t;
		}

		iterator& operator+=(difference_type i)
		{
			while (i-- > 0) operator++();
			return *this;
		}

		iterator operator+(difference_type i) const
		{
			auto result = *this;
			result += i;
			return result;
		}

		friend iterator operator+(difference_type i, const iterator& iter)
		{
			auto result = iter;
			result += i;
			return result;
		}

		friend difference_type operator-(const iterator& a, const iterator& b)
		{
			return std::distance(a.mPos, b.mPos);
		}

		bool operator==(const iterator& i) const		{ return mPos == i.mPos; }
		bool operator!=(const iterator& i) const		{ return mPos != i.mPos; }

	  private:

		friend class RowSet;

		std::vector<ItemRow*>::iterator current() const	{ return mPos; }

		std::vector<ItemRow*>::iterator mPos, mEnd;
		Row mCurrentRow;
	};

	iterator begin()		{ return iterator(mItems.begin(), mItems.end()); }
	iterator end()			{ return iterator(mItems.end(), mItems.end()); }

	Row front() const		{ return Row(mItems.front()); }
	size_t size() const		{ return mItems.size(); }
	bool empty() const		{ return mItems.empty(); }

	template<typename InputIterator>
	iterator insert(iterator pos, InputIterator b, InputIterator e)
	{
		difference_type offset = pos - begin();
		for (auto i = b; i != e; ++i, ++offset)
			insert(begin() + offset, *i);
		return begin() + offset;
	}

	iterator insert(iterator pos, Row& row)
	{
		return insert(pos, row.mData);
	}

	iterator insert(iterator pos, ItemRow* item)
	{
		return iterator(mItems.insert(pos.current(), item), mItems.end());
	}

  private:
	Category*			mCat;
	std::vector<ItemRow*>	mItems;

	// Condition*	mCond;
};

// --------------------------------------------------------------------
// class Category acts as an STL container for Row objects 

class Category
{
  public:
	friend class Datablock;
	friend class Row;
	friend class detail::ItemReference;

	Category(Datablock& db, const std::string& name, Validator* Validator);
	Category(const Category&) = delete;
	Category& operator=(const Category&) = delete;
	~Category();

	const std::string name() const						{ return mName; }

	using iterator = iterator_impl<Row>;
	using const_iterator = iterator_impl<const Row>;
	
	iterator begin();
	iterator end();

	const_iterator cbegin();
	const_iterator cend();

	const_iterator begin() const;
	const_iterator end() const;

	bool empty() const;
	size_t size() const;
	
	void clear();
	
	Row front()										{ return Row(mHead); }
	Row back()										{ return Row(mTail); }
	
	Row operator[](Condition&& cond);

	conditional_iterator_proxy<Row> find(Condition&& cond)
	{
		return { *this, begin(), std::forward<Condition>(cond) };
	}

	conditional_iterator_proxy<Row> find(const_iterator pos, Condition&& cond)
	{
		return { *this, pos, std::forward<Condition>(cond) };
	}

	bool exists(Condition&& cond) const;
	
	RowSet orderBy(const std::string& Item)
		{ return orderBy({ Item }); }
	
	RowSet orderBy(std::initializer_list<std::string> Items);
	
	std::tuple<Row,bool> emplace(Item value)		{ return emplace({ value }); }

	std::tuple<Row,bool> emplace(std::initializer_list<Item> values)
		{ return emplace(values.begin(), values.end()); }

	std::tuple<Row,bool> emplace(Row r);
	
	template<class Iter>
	std::tuple<Row,bool> emplace(Iter b, Iter e);

	size_t erase(Condition&& cond);
	size_t erase(Condition&& cond, std::function<void(const Row&)>&& visit);

	void erase(Row r);
	iterator erase(iterator ri);

	// erase without cascade, should only be used when speed is needed

	size_t erase_nocascade(Condition&& cond)
	{
		return erase_nocascade(std::forward<Condition>(cond), [](auto r){});
	}

	size_t erase_nocascade(Condition&& cond, std::function<void(const Row&)>&& visit)
	{
		auto savedValidator = mValidator;
		mValidator = nullptr;
		auto result = erase(std::forward<Condition>(cond), std::forward<std::function<void(const Row&)>>(visit));
		mValidator = savedValidator;
		return result;
	}

	void eraseOrphans(Condition&& cond);

	/// an orphan is a row that is the child side of one or more
	/// links and for which there is no single parent left.
	bool isOrphan(Row r);
	bool hasParent(Row r, const Category& parentCat, const ValidateLink& link) const;

	bool hasChildren(Row r) const;

	bool isValid();
	void validateLinks() const;

	const Validator& getValidator() const;
	const ValidateCategory* getCatValidator() const		{ return mCatValidator; }

	Datablock& db()										{ return mDb; }
	
	void setValidator(Validator* v);

	iset fields() const;
	iset mandatoryFields() const;
	iset keyFields() const;
	
	std::set<size_t> keyFieldsByIndex() const;
	
	void drop(const std::string& field);

	void getTagOrder(std::vector<std::string>& tags) const;

	// return index for known column, or the next available column index
	size_t getColumnIndex(const std::string& name) const;
    bool hasColumn(const std::string& name) const;
	const std::string& getColumnName(size_t columnIndex) const;
	std::vector<std::string> getColumnNames() const;

	void reorderByIndex();
	void sort(std::function<int(const Row&, const Row&)> comparator);

  private:

	void write(std::ostream& os);
	void write(std::ostream& os, const std::vector<std::string>& order);
	void write(std::ostream& os, const std::vector<int>& order, bool includeEmptyColumns);

	size_t addColumn(const std::string& name);
	
	Datablock&			mDb;
	std::string			mName;
	Validator*			mValidator;
	const ValidateCategory*	mCatValidator = nullptr;
	std::vector<ItemColumn>	mColumns;
	ItemRow*			mHead;
	ItemRow*			mTail;
	class CatIndex*		mIndex;
};

// --------------------------------------------------------------------

class File
{
  public:
	friend class parser;
	friend class Validator;

	File();
	File(std::istream& is, bool validate = false);
	File(const std::string& path, bool validate = false);
	File(File&& rhs);
	File(const File& rhs) = delete;
	File& operator=(const File& rhs) = delete;
	
	~File();

	void load(const std::string& p);
	void save(const std::string& p);

	void load(std::istream& is);
	void save(std::ostream& os);

	void save(std::ostream& os, const std::vector<std::string>& order)	{ write(os, order); }
	void write(std::ostream& os, const std::vector<std::string>& order);

	void loadDictionary();						// load the default dictionary, that is mmcifDdl in this case
	void loadDictionary(const char* dict);		// load one of the compiled in dictionaries 
	void loadDictionary(std::istream& is);		// load dictionary from input stream

	bool isValid();
	void validateLinks() const;
	
	Datablock& firstDatablock()			{ return *mHead; }
	void append(Datablock* e);
	
	Datablock* get(const std::string& name) const;
	Datablock& operator[](const std::string& name);

	struct iterator : public std::iterator<std::forward_iterator_tag, Datablock>
	{
		typedef std::iterator<std::forward_iterator_tag, Datablock>	baseType;
		typedef typename baseType::pointer							pointer;
		typedef typename baseType::reference						reference;
		
		iterator(Datablock* db) : mCurrent(db) {}
		
		reference operator*()						{ return *mCurrent; }
		pointer operator->()						{ return mCurrent; }
		
		iterator& operator++();
		iterator operator++(int)					{ iterator result(*this); this->operator++(); return result; } 

		bool operator==(const iterator& rhs) const	{ return mCurrent == rhs.mCurrent; } 
		bool operator!=(const iterator& rhs) const	{ return not (mCurrent == rhs.mCurrent); } 
		
	  private:
		Datablock*		mCurrent;
	};
	
	iterator begin() const;
	iterator end() const;
	
	bool empty() const								{ return mHead == nullptr; }

	const Validator& getValidator() const;
	void getTagOrder(std::vector<std::string>& tags) const;
	
  private:

	void setValidator(Validator* v);

	Datablock*	mHead;
	Validator*	mValidator;
};

// --------------------------------------------------------------------
// some postponed inlines

namespace detail
{

template<typename T>
inline
bool anyIsConditionImpl<T>::test(const Category& c, const Row& r) const
{
	bool result = false;
	for (auto& f: c.fields())
	{
		try
		{
			if (r[f].as<valueType>() == mValue)
			{
				result = true;
				break;
			}
		}
		catch (...) {}
	}
	
	return result;
}

inline bool anyMatchesConditionImpl::test(const Category& c, const Row& r) const
{
	bool result = false;
	for (auto& f: c.fields())
	{
		try
		{
			if (std::regex_match(r[f].as<std::string>(), mRx))
			{
				result = true;
				break;
			}
		}
		catch (...) {}
	}
	
	return result;
}
	
}

// these should be here, as I learned today

inline void swap(cif::Row& a, cif::Row& b)
{
	a.swap(b);
}

inline void swap(cif::detail::ItemReference& a, cif::detail::ItemReference& b)
{
	a.swap(b);
}

// --------------------------------------------------------------------

template<typename RowType>
conditional_iterator_proxy<RowType>::conditional_iterator_impl::conditional_iterator_impl(Category& cat, base_iterator pos, const Condition& cond)
	: mCat(&cat), mBegin(pos), mEnd(cat.end()), mCondition(&cond)
{
	// skip until the first row matching cond

	while (mBegin != mEnd and not (*mCondition)(*mCat, *mBegin))
		++mBegin;
}

template<typename RowType>
conditional_iterator_proxy<RowType>::conditional_iterator_proxy(conditional_iterator_proxy&& p)
	: mCat(nullptr), mCBegin(p.mCBegin), mCEnd(p.mCEnd)
{
	std::swap(mCat, p.mCat);
	mCondition.swap(p.mCondition);
}

template<typename RowType>
conditional_iterator_proxy<RowType>::conditional_iterator_proxy(Category& cat, base_iterator pos, Condition&& cond)
	: mCat(&cat), mCondition(std::move(cond))
	, mCBegin(pos)
	, mCEnd(cat.end())
{
	mCondition.prepare(cat);

	while (mCBegin != mCEnd and not mCondition(*mCat, *mCBegin))
		++mCBegin;
}

template<typename RowType>
conditional_iterator_proxy<RowType>& conditional_iterator_proxy<RowType>::operator=(conditional_iterator_proxy&& p)
{
	swap(p);
	return *this;
}

template<typename RowType>
typename conditional_iterator_proxy<RowType>::iterator conditional_iterator_proxy<RowType>::begin() const
{
	return iterator(*mCat, mCBegin, mCondition);
}

template<typename RowType>
typename conditional_iterator_proxy<RowType>::iterator conditional_iterator_proxy<RowType>::end() const
{
	return iterator(*mCat, mCEnd, mCondition);
}

template<typename RowType>
bool conditional_iterator_proxy<RowType>::empty() const
{
	return mCBegin == mCEnd;
}

template<typename RowType>
void conditional_iterator_proxy<RowType>::swap(conditional_iterator_proxy& rhs)
{
	std::swap(mCat, rhs.mCat);
	mCondition.swap(rhs.mCondition);
	std::swap(mCBegin, rhs.mCBegin);
	std::swap(mCEnd, rhs.mCEnd);
}

}

