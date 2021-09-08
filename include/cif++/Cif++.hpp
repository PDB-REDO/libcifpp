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

#include <string>

#include <array>
#include <functional>
#include <iostream>
#include <list>
#include <optional>
#include <regex>
#include <set>
#include <sstream>

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
CIFPP_EXPORT extern int VERBOSE;

// mmCIF mapping
// A CIF data file in this case contains entries (data blocks) which can contain
// one or more Category objects. Each Category object contains arrays of Items.
// Better, you can consider the categories as tables containing columns which
// are the Items.

class File;
class Datablock;
class Category;
class Row; // a flyweight class that references data in categories
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

	template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	Item(const std::string &name, const T &value)
		: mName(name)
		, mValue(std::to_string(value))
	{
	}

	Item(const std::string &name, const std::string &value)
		: mName(name)
		, mValue(value)
	{
	}

	Item(const Item &rhs)
		: mName(rhs.mName)
		, mValue(rhs.mValue)
	{
	}
	Item(Item &&rhs) noexcept
		: mName(std::move(rhs.mName))
		, mValue(std::move(rhs.mValue))
	{
	}

	Item &operator=(const Item &rhs)
	{
		if (this != &rhs)
		{
			mName = rhs.mName;
			mValue = rhs.mValue;
		}

		return *this;
	}

	Item &operator=(Item &&rhs) noexcept
	{
		if (this != &rhs)
		{
			mName = std::move(rhs.mName);
			mValue = std::move(rhs.mValue);
		}

		return *this;
	}

	const std::string &name() const { return mName; }
	const std::string &value() const { return mValue; }

	void value(const std::string &v) { mValue = v; }

	// empty means either null or unknown
	bool empty() const;

	// is_null means the field contains '.'
	bool is_null() const;

	// is_unknown means the field contains '?'
	bool is_unknown() const;

	size_t length() const { return mValue.length(); }
	const char *c_str() const { return mValue.c_str(); }

  private:
	std::string mName;
	std::string mValue;
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

	Datablock(const std::string &name);
	~Datablock();

	Datablock(const Datablock &) = delete;
	Datablock &operator=(const Datablock &) = delete;

	std::string getName() const { return mName; }
	void setName(const std::string &n) { mName = n; }

	std::string firstItem(const std::string &tag) const;

	iterator begin() { return mCategories.begin(); }
	iterator end() { return mCategories.end(); }

	const_iterator begin() const { return mCategories.begin(); }
	const_iterator end() const { return mCategories.end(); }

	Category &operator[](const std::string &name);

	std::tuple<iterator, bool> emplace(const std::string &name);

	bool isValid();
	void validateLinks() const;

	void setValidator(Validator *v);

	// this one only looks up a Category, returns nullptr if it does not exist
	const Category *get(const std::string &name) const;
	Category *get(const std::string &name);

	void getTagOrder(std::vector<std::string> &tags) const;
	void write(std::ostream &os, const std::vector<std::string> &order);
	void write(std::ostream &os);

	// convenience function, add a line to the software category
	void add_software(const std::string &name, const std::string &classification,
		const std::string &versionNr, const std::string &versionDate);

	friend bool operator==(const Datablock &lhs, const Datablock &rhs);

	friend std::ostream& operator<<(std::ostream &os, const Datablock &data);

  private:
	std::list<Category> mCategories;
	std::string mName;
	Validator *mValidator;
	Datablock *mNext;
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
		template <typename T, typename = void>
		struct item_value_as;

		template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
		ItemReference &operator=(const T &value)
		{
			this->operator=(std::to_string(value));
			return *this;
		}

		template <typename T>
		ItemReference &operator=(const std::optional<T> &value)
		{
			if (value)
				this->operator=(*value);
			else
				this->operator=("?");
			return *this;
		}

		ItemReference &operator=(const std::string &value);

		template <typename... Ts>
		void os(const Ts &...v)
		{
			std::ostringstream ss;
			((ss << v), ...);
			this->operator=(ss.str());
		}

		void swap(ItemReference &b);

		template <typename T = std::string>
		auto as() const
		{
			using value_type = std::remove_cv_t<std::remove_reference_t<T>>;
			return item_value_as<value_type>::convert(*this);
		}

		template <typename T>
		int compare(const T &value, bool icase) const
		{
			return item_value_as<T>::compare(*this, value, icase);
		}

		// empty means either null or unknown
		bool empty() const;
		explicit operator bool() const { return not empty(); }

		// is_null means the field contains '.'
		bool is_null() const;

		// is_unknown means the field contains '?'
		bool is_unknown() const;

		const char *c_str() const;

		// the following returns the defaultValue from either the parameter
		// or, if specified, the value from _item_default.value in the dictionary
		const char *c_str(const char *defaultValue) const;

		bool operator!=(const std::string &s) const { return s != c_str(); }
		bool operator==(const std::string &s) const { return s == c_str(); }

	  private:
		friend class ::cif::Row;

		ItemReference(const char *name, size_t column, Row &row)
			: mName(name)
			, mColumn(column)
			, mRow(row)
		{
		}

		ItemReference(const char *name, size_t column, const Row &row)
			: mName(name)
			, mColumn(column)
			, mRow(const_cast<Row &>(row))
			, mConst(true)
		{
		}

		const char *mName;
		size_t mColumn;
		Row &mRow;
		bool mConst = false;
	};

	template <typename T>
	struct ItemReference::item_value_as<T, std::enable_if_t<std::is_floating_point_v<T>>>
	{
		using value_type = std::remove_reference_t<std::remove_cv_t<T>>;

		static value_type convert(const ItemReference &ref)
		{
			value_type result = {};
			if (not ref.empty())
				result = static_cast<value_type>(std::stod(ref.c_str()));
			return result;
		}

		static int compare(const ItemReference &ref, double value, bool icase)
		{
			int result = 0;

			const char *s = ref.c_str();

			if (s == nullptr or *s == 0)
				result = 1;
			else
			{
				try
				{
					auto v = std::strtod(s, nullptr);
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
			}

			return result;
		}
	};

	template <typename T>
	struct ItemReference::item_value_as<T, std::enable_if_t<std::is_integral_v<T> and std::is_unsigned_v<T> and not std::is_same_v<T, bool>>>
	{
		static T convert(const ItemReference &ref)
		{
			T result = {};
			if (not ref.empty())
				result = static_cast<T>(std::stoul(ref.c_str()));
			return result;
		}

		static int compare(const ItemReference &ref, unsigned long value, bool icase)
		{
			int result = 0;

			const char *s = ref.c_str();

			if (s == nullptr or *s == 0)
				result = 1;
			else
			{
				try
				{
					auto v = std::strtoul(s, nullptr, 10);
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
			}

			return result;
		}
	};

	template <typename T>
	struct ItemReference::item_value_as<T, std::enable_if_t<std::is_integral_v<T> and std::is_signed_v<T> and not std::is_same_v<T, bool>>>
	{
		static T convert(const ItemReference &ref)
		{
			T result = {};
			if (not ref.empty())
				result = static_cast<T>(std::stol(ref.c_str()));
			return result;
		}

		static int compare(const ItemReference &ref, long value, bool icase)
		{
			int result = 0;

			const char *s = ref.c_str();

			if (s == nullptr or *s == 0)
				result = 1;
			else
			{
				try
				{
					auto v = std::strtol(s, nullptr, 10);
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
			}

			return result;
		}
	};

	template <typename T>
	struct ItemReference::item_value_as<std::optional<T>>
	{
		static std::optional<T> convert(const ItemReference &ref)
		{
			std::optional<T> result;
			if (ref)
				result = ref.as<T>();
			return result;
		}

		static int compare(const ItemReference &ref, std::optional<T> value, bool icase)
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

	template <typename T>
	struct ItemReference::item_value_as<T, std::enable_if_t<std::is_same_v<T, bool>>>
	{
		static bool convert(const ItemReference &ref)
		{
			bool result = false;
			if (not ref.empty())
				result = iequals(ref.c_str(), "y");
			return result;
		}

		static int compare(const ItemReference &ref, bool value, bool icase)
		{
			bool rv = convert(ref);
			return value && rv ? 0
			                   : (rv < value ? -1 : 1);
		}
	};

	template <size_t N>
	struct ItemReference::item_value_as<char[N]>
	{
		static int compare(const ItemReference &ref, const char (&value)[N], bool icase)
		{
			return icase ? cif::icompare(ref.c_str(), value) : std::strcmp(ref.c_str(), value);
		}
	};

	template <>
	struct ItemReference::item_value_as<const char *>
	{
		static const char *convert(const ItemReference &ref)
		{
			return ref.c_str();
		}

		static int compare(const ItemReference &ref, const char *value, bool icase)
		{
			return icase ? cif::icompare(ref.c_str(), value) : std::strcmp(ref.c_str(), value);
		}
	};

	template <>
	struct ItemReference::item_value_as<std::string>
	{
		static std::string convert(const ItemReference &ref)
		{
			return ref.c_str();
		}

		static int compare(const ItemReference &ref, const std::string &value, bool icase)
		{
			return icase ? cif::icompare(ref.c_str(), value) : std::strcmp(ref.c_str(), value.c_str());
		}
	};

	// some helper classes to help create tuple result types
	template <typename... C>
	struct getRowResult
	{
		static constexpr size_t N = sizeof...(C);

		getRowResult(const Row &r, std::array<size_t, N> &&columns)
			: mRow(r)
			, mColumns(std::move(columns))
		{
		}

		const ItemReference operator[](size_t ix) const
		{
			return mRow[mColumns[ix]];
		}

		template <typename... Ts, std::enable_if_t<N == sizeof...(Ts), int> = 0>
		operator std::tuple<Ts...>() const
		{
			return get<Ts...>(std::index_sequence_for<Ts...>{});
		}

		template <typename... Ts, std::size_t... Is>
		std::tuple<Ts...> get(std::index_sequence<Is...>) const
		{
			return std::tuple<Ts...>{mRow[mColumns[Is]].template as<Ts>()...};
		}

		const Row &mRow;
		std::array<size_t, N> mColumns;
	};

	// we want to be able to tie some variables to a RowResult, for this we use tiewraps
	template <typename... Ts>
	struct tieWrap
	{
		tieWrap(Ts... args)
			: mVal(args...)
		{
		}

		template <typename RR>
		void operator=(const RR &&rr)
		{
			// getRowResult will do the conversion, but only if the types
			// are compatible. That means the number of parameters to the get()
			// of the row should be equal to the number of items in the tuple
			// you are trying to tie.

			using RType = std::tuple<typename std::remove_reference<Ts>::type...>;

			mVal = static_cast<RType>(rr);
		}

		std::tuple<Ts...> mVal;
	};
} // namespace detail

template <typename... Ts>
auto tie(Ts &...v)
{
	return detail::tieWrap<Ts &...>(std::forward<Ts &>(v)...);
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
		: mData(nullptr)
	{
	}

	Row(ItemRow *data)
		: mData(data)
	{
	}

	Row(const ItemRow *data)
		: mData(const_cast<ItemRow *>(data))
		, mCascade(false)
	{
	}

	Row(const Row &rhs);
	Row &operator=(const Row &rhs);

	Row(Row &&rhs);
	Row &operator=(Row &&rhs);

	~Row();

	void setCascading(bool cascade)
	{
		mCascade = cascade;
	}

	void next(); ///< make this row point to the next ItemRow

	struct const_iterator
	{
		using iterator_category = std::forward_iterator_tag;
		using value_type = const Item;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type *;
		using reference = value_type &;

		const_iterator(ItemRow *data, ItemValue *ptr);

		reference operator*() { return mCurrent; }
		pointer operator->() { return &mCurrent; }

		const_iterator &operator++();
		const_iterator operator++(int)
		{
			const_iterator result(*this);
			this->operator++();
			return result;
		}

		bool operator==(const const_iterator &rhs) const { return mPtr == rhs.mPtr; }
		bool operator!=(const const_iterator &rhs) const { return mPtr != rhs.mPtr; }

	  private:
		void fetch();

		ItemRow *mData;
		ItemValue *mPtr;
		Item mCurrent;
	};

	// checks for an initialized Row:
	explicit operator bool() const { return mData != nullptr; }

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

	const detail::ItemReference operator[](const char *itemTag) const
	{
		size_t column = ColumnForItemTag(itemTag);
		return detail::ItemReference(itemTag, column, *this);
	}

	detail::ItemReference operator[](const char *itemTag)
	{
		size_t column = ColumnForItemTag(itemTag);
		return detail::ItemReference(itemTag, column, *this);
	}

	const detail::ItemReference operator[](const std::string &itemTag) const
	{
		size_t column = ColumnForItemTag(itemTag.c_str());
		return detail::ItemReference(itemTag.c_str(), column, *this);
	}

	detail::ItemReference operator[](const std::string &itemTag)
	{
		size_t column = ColumnForItemTag(itemTag.c_str());
		return detail::ItemReference(itemTag.c_str(), column, *this);
	}

	template <typename... Ts, size_t N>
	std::tuple<Ts...> get(char const *const (&columns)[N]) const
	{
		static_assert(sizeof...(Ts) == N, "Number of columns should be equal to number of types to return");

		std::array<size_t, N> cix;
		for (size_t i = 0; i < N; ++i)
			cix[i] = ColumnForItemTag(columns[i]);
		return detail::getRowResult<Ts...>(*this, std::move(cix));
	}

	template <typename... C>
	auto get(C... columns) const
	{
		return detail::getRowResult<C...>(*this, {ColumnForItemTag(columns)...});
	}

	void assign(const std::vector<Item> &values);
	void assign(const std::string &name, const std::string &value, bool updateLinked);

	bool operator==(const Row &rhs) const
	{
		return mData == rhs.mData;
	}

	bool operator!=(const Row &rhs) const
	{
		return mData != rhs.mData;
	}

	ItemRow *data() const { return mData; }

	void swap(Row &rhs)
	{
		std::swap(mData, rhs.mData);
	}

	friend std::ostream &operator<<(std::ostream &os, const Row &row);

  private:
	void assign(size_t column, const std::string &value, bool updateLinked);
	void assign(const Item &i, bool updateLinked);

	static void swap(size_t column, ItemRow *a, ItemRow *b);

	size_t ColumnForItemTag(const char *itemTag) const;

	ItemRow *mData;
	uint32_t mLineNr = 0;
	bool mCascade = true;
};

// --------------------------------------------------------------------
// some more templates to be able to do querying

namespace detail
{

	struct ConditionImpl
	{
		virtual ~ConditionImpl() {}

		virtual void prepare(const Category &c) {}
		virtual bool test(const Category &c, const Row &r) const = 0;
		virtual void str(std::ostream &os) const = 0;
	};

	struct AllConditionImpl : public ConditionImpl
	{
		virtual bool test(const Category &c, const Row &r) const { return true; }
		virtual void str(std::ostream &os) const { os << "*"; }
	};

	struct OrConditionImpl;
	struct AndConditionImpl;
	struct NotConditionImpl;

} // namespace detail

class Condition
{
  public:
	Condition()
		: mImpl(nullptr)
	{
	}
	Condition(detail::ConditionImpl *impl)
		: mImpl(impl)
	{
	}

	Condition(const Condition &) = delete;

	Condition(Condition &&rhs) noexcept
		: mImpl(nullptr)
	{
		std::swap(mImpl, rhs.mImpl);
	}

	Condition &operator=(const Condition &) = delete;

	Condition &operator=(Condition &&rhs) noexcept
	{
		std::swap(mImpl, rhs.mImpl);
		return *this;
	}

	~Condition()
	{
		delete mImpl;
		mImpl = nullptr;
	}

	void prepare(const Category &c)
	{
		if (mImpl)
			mImpl->prepare(c);
		mPrepared = true;
	}

	bool operator()(const Category &c, const Row &r) const
	{
		assert(mImpl);
		assert(mPrepared);
		return mImpl ? mImpl->test(c, r) : false;
	}

	bool empty() const { return mImpl == nullptr; }

	friend Condition operator||(Condition &&a, Condition &&b);
	friend Condition operator&&(Condition &&a, Condition &&b);

	friend struct detail::OrConditionImpl;
	friend struct detail::AndConditionImpl;
	friend struct detail::NotConditionImpl;

	void swap(Condition &rhs)
	{
		std::swap(mImpl, rhs.mImpl);
		std::swap(mPrepared, rhs.mPrepared);
	}

	friend std::ostream &operator<<(std::ostream &os, const Condition &cond);

  private:
	detail::ConditionImpl *mImpl;
	bool mPrepared = false;
};

inline std::ostream &operator<<(std::ostream &os, const Condition &cond)
{
	if (cond.mImpl)
		cond.mImpl->str(os);
	return os;
}

namespace detail
{

	struct KeyIsEmptyConditionImpl : public ConditionImpl
	{
		KeyIsEmptyConditionImpl(const std::string &ItemTag)
			: mItemTag(ItemTag)
		{
		}

		virtual void prepare(const Category &c);

		virtual bool test(const Category &c, const Row &r) const
		{
			return r[mItemIx].empty();
		}

		virtual void str(std::ostream &os) const
		{
			os << mItemTag << " IS NULL";
		}

		std::string mItemTag;
		size_t mItemIx = 0;
	};

	struct KeyCompareConditionImpl : public ConditionImpl
	{
		template <typename COMP>
		KeyCompareConditionImpl(const std::string &ItemTag, COMP &&comp, const std::string &s)
			: mItemTag(ItemTag)
			, mComp(std::move(comp))
			, mStr(s)
		{
		}

		virtual void prepare(const Category &c);

		virtual bool test(const Category &c, const Row &r) const
		{
			return mComp(c, r, mCaseInsensitive);
		}

		virtual void str(std::ostream &os) const
		{
			os << mItemTag << (mCaseInsensitive ? "^ " : " ") << mStr;
		}

		std::string mItemTag;
		size_t mItemIx = 0;
		bool mCaseInsensitive = false;
		std::function<bool(const Category &, const Row &, bool)> mComp;
		std::string mStr;
	};

	struct KeyMatchesConditionImpl : public ConditionImpl
	{
		KeyMatchesConditionImpl(const std::string &ItemTag, const std::regex &rx)
			: mItemTag(ItemTag)
			, mItemIx(0)
			, mRx(rx)
		{
		}

		virtual void prepare(const Category &c);

		virtual bool test(const Category &c, const Row &r) const
		{
			return std::regex_match(r[mItemIx].as<std::string>(), mRx);
		}

		virtual void str(std::ostream &os) const
		{
			os << mItemTag << " =~ expression";
		}

		std::string mItemTag;
		size_t mItemIx;
		std::regex mRx;
	};

	template <typename T>
	struct AnyIsConditionImpl : public ConditionImpl
	{
		typedef T valueType;

		AnyIsConditionImpl(const valueType &value)
			: mValue(value)
		{
		}

		virtual bool test(const Category &c, const Row &r) const;
		virtual void str(std::ostream &os) const
		{
			os << "<any> == " << mValue;
		}

		valueType mValue;
	};

	struct AnyMatchesConditionImpl : public ConditionImpl
	{
		AnyMatchesConditionImpl(const std::regex &rx)
			: mRx(rx)
		{
		}

		virtual bool test(const Category &c, const Row &r) const;
		virtual void str(std::ostream &os) const
		{
			os << "<any> =~ expression";
		}

		std::regex mRx;
	};

	struct AndConditionImpl : public ConditionImpl
	{
		AndConditionImpl(Condition &&a, Condition &&b)
			: mA(nullptr)
			, mB(nullptr)
		{
			std::swap(mA, a.mImpl);
			std::swap(mB, b.mImpl);
		}

		~AndConditionImpl()
		{
			delete mA;
			delete mB;
		}

		virtual void prepare(const Category &c)
		{
			mA->prepare(c);
			mB->prepare(c);
		}

		virtual bool test(const Category &c, const Row &r) const
		{
			return mA->test(c, r) and mB->test(c, r);
		}

		virtual void str(std::ostream &os) const
		{
			os << '(';
			mA->str(os);
			os << ") AND (";
			mB->str(os);
			os << ')';
		}

		ConditionImpl *mA;
		ConditionImpl *mB;
	};

	struct OrConditionImpl : public ConditionImpl
	{
		OrConditionImpl(Condition &&a, Condition &&b)
			: mA(nullptr)
			, mB(nullptr)
		{
			std::swap(mA, a.mImpl);
			std::swap(mB, b.mImpl);
		}

		~OrConditionImpl()
		{
			delete mA;
			delete mB;
		}

		virtual void prepare(const Category &c)
		{
			mA->prepare(c);
			mB->prepare(c);
		}

		virtual bool test(const Category &c, const Row &r) const
		{
			return mA->test(c, r) or mB->test(c, r);
		}

		virtual void str(std::ostream &os) const
		{
			os << '(';
			mA->str(os);
			os << ") OR (";
			mB->str(os);
			os << ')';
		}

		ConditionImpl *mA;
		ConditionImpl *mB;
	};

	struct NotConditionImpl : public ConditionImpl
	{
		NotConditionImpl(Condition &&a)
			: mA(nullptr)
		{
			std::swap(mA, a.mImpl);
		}

		~NotConditionImpl()
		{
			delete mA;
		}

		virtual void prepare(const Category &c)
		{
			mA->prepare(c);
		}

		virtual bool test(const Category &c, const Row &r) const
		{
			return not mA->test(c, r);
		}

		virtual void str(std::ostream &os) const
		{
			os << "NOT (";
			mA->str(os);
			os << ')';
		}

		ConditionImpl *mA;
	};

} // namespace detail

inline Condition operator&&(Condition &&a, Condition &&b)
{
	if (a.mImpl and b.mImpl)
		return Condition(new detail::AndConditionImpl(std::move(a), std::move(b)));
	if (a.mImpl)
		return Condition(std::move(a));
	return Condition(std::move(b));
}

inline Condition operator||(Condition &&a, Condition &&b)
{
	if (a.mImpl and b.mImpl)
		return Condition(new detail::OrConditionImpl(std::move(a), std::move(b)));
	if (a.mImpl)
		return Condition(std::move(a));
	return Condition(std::move(b));
}

struct Empty
{
};

struct Key
{
	Key(const std::string &itemTag)
		: mItemTag(itemTag)
	{
	}
	Key(const char *itemTag)
		: mItemTag(itemTag)
	{
	}

	Key(const Key &) = delete;
	Key &operator=(const Key &) = delete;

	std::string mItemTag;
};

template <typename T>
Condition operator==(const Key &key, const T &v)
{
	std::ostringstream s;
	s << " == " << v;

	return Condition(new detail::KeyCompareConditionImpl(
		key.mItemTag, [tag = key.mItemTag, v](const Category &c, const Row &r, bool icase)
		{ return r[tag].template compare<T>(v, icase) == 0; },
		s.str()));
}

inline Condition operator==(const Key &key, const char *value)
{
	if (value != nullptr and *value != 0)
	{
		std::ostringstream s;
		s << " == " << value;

		return Condition(new detail::KeyCompareConditionImpl(
			key.mItemTag, [tag = key.mItemTag, value](const Category &c, const Row &r, bool icase)
			{ return r[tag].compare(value, icase) == 0; },
			s.str()));
	}
	else
		return Condition(new detail::KeyIsEmptyConditionImpl(key.mItemTag));
}

// inline Condition operator==(const Key& key, const detail::ItemReference& v)
// {
// 	if (v.empty())
// 		return Condition(new detail::KeyIsEmptyConditionImpl(key.mItemTag));
// 	else
// 		return Condition(new detail::KeyCompareConditionImpl(key.mItemTag, [tag = key.mItemTag, v](const Category& c, const Row& r, bool icase)
// 			{ return r[tag].template compare<(v, icase) == 0; }));
// }

inline Condition operator==(const Key &key, const Empty &)
{
	return Condition(new detail::KeyIsEmptyConditionImpl(key.mItemTag));
}

template <typename T>
Condition operator!=(const Key &key, const T &v)
{
	return Condition(new detail::NotConditionImpl(operator==(key, v)));
}

inline Condition operator!=(const Key &key, const char *v)
{
	std::string value(v ? v : "");
	return Condition(new detail::NotConditionImpl(operator==(key, value)));
}

template <typename T>
Condition operator>(const Key &key, const T &v)
{
	std::ostringstream s;
	s << " > " << v;

	return Condition(new detail::KeyCompareConditionImpl(
		key.mItemTag, [tag = key.mItemTag, v](const Category &c, const Row &r, bool icase)
		{ return r[tag].template compare<T>(v, icase) > 0; },
		s.str()));
}

template <typename T>
Condition operator>=(const Key &key, const T &v)
{
	std::ostringstream s;
	s << " >= " << v;

	return Condition(new detail::KeyCompareConditionImpl(
		key.mItemTag, [tag = key.mItemTag, v](const Category &c, const Row &r, bool icase)
		{ return r[tag].template compare<T>(v, icase) >= 0; },
		s.str()));
}

template <typename T>
Condition operator<(const Key &key, const T &v)
{
	std::ostringstream s;
	s << " < " << v;

	return Condition(new detail::KeyCompareConditionImpl(
		key.mItemTag, [tag = key.mItemTag, v](const Category &c, const Row &r, bool icase)
		{ return r[tag].template compare<T>(v, icase) < 0; },
		s.str()));
}

template <typename T>
Condition operator<=(const Key &key, const T &v)
{
	std::ostringstream s;
	s << " <= " << v;

	return Condition(new detail::KeyCompareConditionImpl(
		key.mItemTag, [tag = key.mItemTag, v](const Category &c, const Row &r, bool icase)
		{ return r[tag].template compare<T>(v, icase) <= 0; },
		s.str()));
}

template <>
inline Condition operator==(const Key &key, const std::regex &rx)
{
	return Condition(new detail::KeyMatchesConditionImpl(key.mItemTag, rx));
}

template <>
inline Condition operator==(const Key &key, const Empty &)
{
	return Condition(new detail::KeyIsEmptyConditionImpl(key.mItemTag));
}

struct any
{
	template <typename T>
	Condition operator==(const T &v) const
	{
		return Condition(new detail::AnyIsConditionImpl<T>(v));
	}
};

template <>
inline Condition any::operator==(const std::regex &rx) const
{
	return Condition(new detail::AnyMatchesConditionImpl(rx));
}

inline Condition All()
{
	return Condition(new detail::AllConditionImpl());
}

inline Condition Not(Condition &&cond)
{
	return Condition(new detail::NotConditionImpl(std::move(cond)));
}

namespace literals
{

	inline Key operator""_key(const char *text, size_t length)
	{
		return Key(std::string(text, length));
	}

	inline constexpr Empty Null = Empty();

} // namespace literals

// -----------------------------------------------------------------------
// iterators

template <typename RowType, typename... Ts>
class iterator_impl
{
  public:
	template <typename, typename...>
	friend class iterator_impl;

	static constexpr size_t N = sizeof...(Ts);

	using iterator_category = std::forward_iterator_tag;
	using value_type = std::conditional_t<N == 0, RowType, std::tuple<Ts...>>;
	using difference_type = std::ptrdiff_t;
	using pointer = value_type *;
	using reference = value_type &;

	friend class Category;

	// default constructor, equal to end()
	iterator_impl() {}

	iterator_impl(const iterator_impl &rhs)
		: mCurrent(rhs.mCurrent)
		, mValue(rhs.mValue)
		, mCix(rhs.mCix)
	{
	}

	iterator_impl(ItemRow *data)
		: mCurrent(data)
	{
		static_assert(N == 0, "Only valid if this is a row iterator, not a row<xxx> iterator");
	}

	iterator_impl(ItemRow *data, const std::array<size_t, N> &cix)
		: mCurrent(data)
		, mCix(cix)
	{
	}

	template <typename IRowType>
	iterator_impl(iterator_impl<IRowType, Ts...> &rhs)
		: mCurrent(rhs.mCurrent)
		, mCix(rhs.mCix)
	{
		if constexpr (N > 0)
			mValue = get(mCurrent, std::make_index_sequence<N>());
	}

	template <typename IRowType>
	iterator_impl(const iterator_impl<IRowType> &rhs, const std::array<size_t, N> &cix)
		: mCurrent(rhs.mCurrent)
		, mCix(cix)
	{
		if constexpr (N > 0)
			mValue = get(mCurrent, std::make_index_sequence<N>());
	}

	iterator_impl &operator=(const iterator_impl &i)
	{
		mCurrent = i.mCurrent;
		if constexpr (N != 0)
		{
			mCix = i.mCix;
			mValue = i.mValue;
		}
		return *this;
	}

	virtual ~iterator_impl() = default;

	reference operator*()
	{
		if constexpr (N == 0)
			return mCurrent;
		else
			return mValue;
	}

	pointer operator->()
	{
		if constexpr (N == 0)
			return &mCurrent;
		else
			return &mValue;
	}

	RowType row() const
	{
		return mCurrent;
	}

	iterator_impl &operator++()
	{
		mCurrent.next();

		if constexpr (N != 0)
			mValue = get(mCurrent, std::make_index_sequence<N>());

		return *this;
	}

	iterator_impl operator++(int)
	{
		iterator_impl result(*this);
		this->operator++();
		return result;
	}

	bool operator==(const iterator_impl &rhs) const { return mCurrent == rhs.mCurrent; }
	bool operator!=(const iterator_impl &rhs) const { return mCurrent != rhs.mCurrent; }

	template <typename IRowType, typename... ITs>
	bool operator==(const iterator_impl<IRowType, ITs...> &rhs) const
	{
		return mCurrent == rhs.mCurrent;
	}

	template <typename IRowType, typename... ITs>
	bool operator!=(const iterator_impl<IRowType, ITs...> &rhs) const
	{
		return mCurrent != rhs.mCurrent;
	}

  private:
	template <std::size_t... Is>
	std::tuple<Ts...> get(Row row, std::index_sequence<Is...>) const
	{
		if (row)
			return std::tuple<Ts...>{row[mCix[Is]].template as<Ts>()...};
		return {};
	}

	Row mCurrent;
	value_type mValue;
	std::array<size_t, N> mCix;
};

// --------------------------------------------------------------------
// iterator proxy

template <typename RowType, typename... Ts>
class iterator_proxy
{
  public:
	static constexpr const size_t N = sizeof...(Ts);

	using iterator = iterator_impl<RowType, Ts...>;
	using row_iterator = iterator_impl<RowType>;

	iterator_proxy(Category &cat, row_iterator pos, char const *const columns[N]);
	iterator_proxy(Category &cat, row_iterator pos, std::initializer_list<char const *> columns);

	iterator_proxy(iterator_proxy &&p);
	iterator_proxy &operator=(iterator_proxy &&p);

	iterator_proxy(const iterator_proxy &) = delete;
	iterator_proxy &operator=(const iterator_proxy &) = delete;

	iterator begin() const { return iterator(mCBegin, mCix); }
	iterator end() const { return iterator(mCEnd, mCix); }

	bool empty() const { return mCBegin == mCEnd; }

	explicit operator bool() const { return not empty(); }

	size_t size() const { return std::distance(begin(), end()); }

	RowType front() { return *begin(); }
	RowType back() { return *(std::prev(end())); }

	Category &category() const { return *mCat; }

	void swap(iterator_proxy &rhs)
	{
		std::swap(mCat, rhs.mCat);
		std::swap(mCBegin, rhs.mCBegin);
		std::swap(mCEnd, rhs.mCEnd);
		std::swap(mCix, rhs.mCix);
	}

  private:
	Category *mCat;
	row_iterator mCBegin, mCEnd;
	std::array<size_t, N> mCix;
};

// --------------------------------------------------------------------
// conditional iterator proxy

template <typename RowType, typename... Ts>
class conditional_iterator_proxy
{
  public:
	static constexpr const size_t N = sizeof...(Ts);

	using base_iterator = iterator_impl<RowType, Ts...>;
	using value_type = typename base_iterator::value_type;

	using row_type = std::remove_cv_t<RowType>;
	using row_iterator = iterator_impl<row_type>;

	class conditional_iterator_impl
	{
	  public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = conditional_iterator_proxy::value_type;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type *;
		using reference = value_type &;

		conditional_iterator_impl(Category &cat, row_iterator pos, const Condition &cond, const std::array<size_t, N> &cix);
		conditional_iterator_impl(const conditional_iterator_impl &i) = default;
		conditional_iterator_impl &operator=(const conditional_iterator_impl &i) = default;

		virtual ~conditional_iterator_impl() = default;

		reference operator*()
		{
			return *mBegin;
		}

		pointer operator->()
		{
			return &*mBegin;
		}

		conditional_iterator_impl &operator++()
		{
			while (mBegin != mEnd)
			{
				if (++mBegin == mEnd)
					break;

				if ((*mCondition)(*mCat, mBegin.row()))
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

		bool operator==(const conditional_iterator_impl &rhs) const { return mBegin == rhs.mBegin; }
		bool operator!=(const conditional_iterator_impl &rhs) const { return mBegin != rhs.mBegin; }

		template <typename IRowType, typename... ITs>
		bool operator==(const iterator_impl<IRowType, ITs...> &rhs) const { return mBegin == rhs; }

		template <typename IRowType, typename... ITs>
		bool operator!=(const iterator_impl<IRowType, ITs...> &rhs) const { return mBegin != rhs; }

	  private:
		Category *mCat;
		base_iterator mBegin, mEnd;
		const Condition *mCondition;
	};

	using iterator = conditional_iterator_impl;
	using reference = typename iterator::reference;

	conditional_iterator_proxy(Category &cat, row_iterator pos, Condition &&cond);

	template <std::size_t TN = N, std::enable_if_t<TN != 0, bool> = true>
	conditional_iterator_proxy(Category &cat, row_iterator pos, Condition &&cond, char const *const columns[N]);

	conditional_iterator_proxy(conditional_iterator_proxy &&p);
	conditional_iterator_proxy &operator=(conditional_iterator_proxy &&p);

	conditional_iterator_proxy(const conditional_iterator_proxy &) = delete;
	conditional_iterator_proxy &operator=(const conditional_iterator_proxy &) = delete;

	iterator begin() const;
	iterator end() const;

	bool empty() const;

	explicit operator bool() const { return not empty(); }

	size_t size() const { return std::distance(begin(), end()); }

	RowType front() { return *begin(); }

	Category &category() const { return *mCat; }

	void swap(conditional_iterator_proxy &rhs);

  private:
	Category *mCat;
	Condition mCondition;
	row_iterator mCBegin, mCEnd;
	std::array<size_t, N> mCix;
};

// --------------------------------------------------------------------
// class RowSet is used to return find results. Use it to re-order the results
// or to group them

class RowSet
{
	typedef std::vector<Row> base_type;

  public:
	using size_type = std::size_t;
	using difference_type = std::ptrdiff_t;

	RowSet(Category &cat);
	RowSet(Category &cat, Condition &&cond);
	RowSet(const RowSet &rhs);
	RowSet(RowSet &&rhs);
	virtual ~RowSet();

	RowSet &operator=(const RowSet &rhs);
	RowSet &operator=(RowSet &&rhs);

	RowSet &orderBy(const std::string &Item)
	{
		return orderBy({Item});
	}

	RowSet &orderBy(std::initializer_list<std::string> Items);

	class iterator
	{
	  public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = Row;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type *;
		using reference = value_type &;

		iterator() {}

		iterator(const iterator &i)
			: mPos(i.mPos)
			, mEnd(i.mEnd)
			, mCurrentRow(i.mCurrentRow)
		{
		}

		iterator(const std::vector<ItemRow *>::iterator &p, const std::vector<ItemRow *>::iterator &e)
			: mPos(p)
			, mEnd(e)
			, mCurrentRow(p == e ? nullptr : *p)
		{
		}

		iterator &operator=(const iterator &i)
		{
			mPos = i.mPos;
			mEnd = i.mEnd;
			mCurrentRow = i.mCurrentRow;
			return *this;
		}

		reference operator*() { return mCurrentRow; }
		pointer operator->() { return &mCurrentRow; }

		iterator &operator++()
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

		iterator &operator+=(difference_type i)
		{
			while (i-- > 0)
				operator++();
			return *this;
		}

		iterator operator+(difference_type i) const
		{
			auto result = *this;
			result += i;
			return result;
		}

		friend iterator operator+(difference_type i, const iterator &iter)
		{
			auto result = iter;
			result += i;
			return result;
		}

		friend difference_type operator-(const iterator &a, const iterator &b)
		{
			return std::distance(a.mPos, b.mPos);
		}

		bool operator==(const iterator &i) const { return mPos == i.mPos; }
		bool operator!=(const iterator &i) const { return mPos != i.mPos; }

	  private:
		friend class RowSet;

		std::vector<ItemRow *>::iterator current() const { return mPos; }

		std::vector<ItemRow *>::iterator mPos, mEnd;
		Row mCurrentRow;
	};

	iterator begin() { return iterator(mItems.begin(), mItems.end()); }
	iterator end() { return iterator(mItems.end(), mItems.end()); }

	Row front() const { return Row(mItems.front()); }
	size_t size() const { return mItems.size(); }
	bool empty() const { return mItems.empty(); }

	template <typename InputIterator>
	iterator insert(iterator pos, InputIterator b, InputIterator e)
	{
		difference_type offset = pos - begin();
		for (auto i = b; i != e; ++i, ++offset)
			insert(begin() + offset, *i);
		return begin() + offset;
	}

	iterator insert(iterator pos, Row &row)
	{
		return insert(pos, row.mData);
	}

	iterator insert(iterator pos, ItemRow *item)
	{
		auto p = mItems.insert(pos.current(), item);
		return iterator(p, mItems.end());
	}

	iterator push_back(ItemRow *item)
	{
		return insert(end(), item);
	}

	iterator push_back(Row &row)
	{
		return insert(end(), row.mData);
	}

	void make_unique()
	{
		std::sort(mItems.begin(), mItems.end());
		mItems.erase(std::unique(mItems.begin(), mItems.end()), mItems.end());
	}

  private:
	Category *mCat;
	std::vector<ItemRow *> mItems;

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

	Category(Datablock &db, const std::string &name, Validator *Validator);
	Category(const Category &) = delete;
	Category &operator=(const Category &) = delete;
	~Category();

	const std::string name() const { return mName; }

	using iterator = iterator_impl<Row>;
	using const_iterator = iterator_impl<const Row>;

	iterator begin();
	iterator end();

	const_iterator cbegin() const;
	const_iterator cend() const;

	const_iterator begin() const;
	const_iterator end() const;

	bool empty() const;
	size_t size() const;

	void clear();

	Row front() { return Row(mHead); }
	Row back() { return Row(mTail); }

	Row operator[](Condition &&cond);

	template <typename... Ts, size_t N>
	iterator_proxy<Row, Ts...> rows(char const *const (&columns)[N])
	{
		static_assert(sizeof...(Ts) == N, "The number of column titles should be equal to the number of types to return");
		return iterator_proxy<Row, Ts...>(*this, begin(), columns);
	}

	template <typename... Ts, typename... Ns>
	iterator_proxy<Row, Ts...> rows(Ns... names)
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of column titles should be equal to the number of types to return");
		return iterator_proxy<Row, Ts...>(*this, begin(), {names...});
	}

	conditional_iterator_proxy<Row> find(Condition &&cond)
	{
		return find(cbegin(), std::forward<Condition>(cond));
	}

	conditional_iterator_proxy<Row> find(const_iterator pos, Condition &&cond)
	{
		return {*this, pos, std::forward<Condition>(cond)};
	}

	template <typename... Ts, size_t N>
	conditional_iterator_proxy<Row, Ts...> find(Condition &&cond, char const *const (&columns)[N])
	{
		static_assert(sizeof...(Ts) == N, "The number of column titles should be equal to the number of types to return");
		return find<Ts...>(cbegin(), std::forward<Condition>(cond), std::forward<char const *const[N]>(columns));
	}

	template <typename... Ts, size_t N>
	conditional_iterator_proxy<Row, Ts...> find(const_iterator pos, Condition &&cond, char const *const (&columns)[N])
	{
		static_assert(sizeof...(Ts) == N, "The number of column titles should be equal to the number of types to return");
		return {*this, pos, std::forward<Condition>(cond), columns};
	}

	// --------------------------------------------------------------------
	// if you only expect a single row

	Row find1(Condition &&cond)
	{
		return find1(cbegin(), std::forward<Condition>(cond));
	}

	Row find1(const_iterator pos, Condition &&cond)
	{
		auto h = find(pos, std::forward<Condition>(cond));

		if (h.empty())
			throw std::runtime_error("No hits found");

		if (h.size() != 1)
			throw std::runtime_error("Hit not unique");

		return *h.begin();
	}

	template <typename... Ts, size_t N>
	std::tuple<Ts...> find1(Condition &&cond, char const *const (&columns)[N])
	{
		static_assert(sizeof...(Ts) == N, "The number of column titles should be equal to the number of types to return");
		return find1<Ts...>(cbegin(), std::forward<Condition>(cond), std::forward<char const *const[N]>(columns));
	}

	template <typename... Ts, size_t N>
	std::tuple<Ts...> find1(const_iterator pos, Condition &&cond, char const *const (&columns)[N])
	{
		static_assert(sizeof...(Ts) == N, "The number of column titles should be equal to the number of types to return");

		auto h = find<Ts...>(pos, std::forward<Condition>(cond), std::forward<char const *const[N]>(columns));

		if (h.empty())
			throw std::runtime_error("No hits found");

		if (h.size() != 1)
			throw std::runtime_error("Hit not unique");

		return *h.begin();
	}

	bool exists(Condition &&cond) const;

	RowSet orderBy(const std::string &Item)
	{
		return orderBy({Item});
	}

	RowSet orderBy(std::initializer_list<std::string> Items);

	std::tuple<Row, bool> emplace(Item value) { return emplace({value}); }

	std::tuple<Row, bool> emplace(std::initializer_list<Item> values)
	{
		return emplace(values.begin(), values.end());
	}

	std::tuple<Row, bool> emplace(Row r);

	template <class Iter>
	std::tuple<Row, bool> emplace(Iter b, Iter e);

	size_t erase(Condition &&cond);
	size_t erase(Condition &&cond, std::function<void(const Row &)> &&visit);

	void erase(Row r);
	iterator erase(iterator ri);

	/// \brief Create a copy of Row \a r and return the copy. If this row has
	/// a single key field, this will be updated with a new unique value.
	Row copyRow(const Row &r);

	// erase without cascade, should only be used when speed is needed

	size_t erase_nocascade(Condition &&cond)
	{
		return erase_nocascade(std::forward<Condition>(cond), [](auto r) {});
	}

	size_t erase_nocascade(Condition &&cond, std::function<void(const Row &)> &&visit)
	{
		auto savedValidator = mValidator;
		mValidator = nullptr;
		auto result = erase(std::forward<Condition>(cond), std::forward<std::function<void(const Row &)>>(visit));
		mValidator = savedValidator;
		return result;
	}

	void eraseOrphans(Condition &&cond);

	/// an orphan is a row that is the child side of one or more
	/// links and for which there is no single parent left.
	bool isOrphan(Row r);
	bool hasParent(Row r, const Category &parentCat, const ValidateLink &link) const;

	bool hasChildren(Row r) const;
	bool hasParents(Row r) const;

	RowSet getChildren(Row r, Category &childCat);
	RowSet getChildren(Row r, const char *childCat);

	RowSet getParents(Row r, Category &parentCat);
	RowSet getParents(Row r, const char *parentCat);

	RowSet getLinked(Row r, Category &cat);
	RowSet getLinked(Row r, const char *cat);

	bool isValid();
	void validateLinks() const;

	const Validator &getValidator() const;
	const ValidateCategory *getCatValidator() const { return mCatValidator; }

	Datablock &db() { return mDb; }

	void setValidator(Validator *v);

	iset fields() const;
	iset mandatoryFields() const;
	iset keyFields() const;

	std::set<size_t> keyFieldsByIndex() const;

	void drop(const std::string &field);

	void getTagOrder(std::vector<std::string> &tags) const;

	// return index for known column, or the next available column index
	size_t getColumnIndex(const std::string &name) const;
	bool hasColumn(const std::string &name) const;
	const std::string &getColumnName(size_t columnIndex) const;
	std::vector<std::string> getColumnNames() const;

	void reorderByIndex();
	void sort(std::function<int(const Row &, const Row &)> comparator);

	// --------------------------------------------------------------------
	/// Rename a single column in the rows that match \a cond to value \a value
	/// making sure the linked categories are updated according to the link.
	/// That means, child categories are updated if the links are absolute
	/// and unique. If they are not, the child category rows are split.

	void update_value(Condition &&cond, const std::string &tag, const std::string &value)
	{
		update_value(RowSet{*this, std::move(cond)}, tag, value);
	}

	void update_value(RowSet &&rows, const std::string &tag, const std::string &value);

	// --------------------------------------------------------------------
	// generate a new, unique ID. Pass it an ID generating function based on
	// a sequence number. This function will be called until the result is
	// unique in the context of this category
	std::string getUniqueID(std::function<std::string(int)> generator = cif::cifIdForNumber);
	std::string getUniqueID(const std::string &prefix)
	{
		return getUniqueID([prefix](int nr)
			{ return prefix + std::to_string(nr); });
	}

	// --------------------------------------------------------------------
	// for debugging

	friend bool operator==(const Category &lhs, const Category &rhs);

  private:
	void write(std::ostream &os);
	void write(std::ostream &os, const std::vector<std::string> &order);
	void write(std::ostream &os, const std::vector<size_t> &order, bool includeEmptyColumns);

	size_t addColumn(const std::string &name);

	Datablock &mDb;
	std::string mName;
	Validator *mValidator;
	const ValidateCategory *mCatValidator = nullptr;
	std::vector<ItemColumn> mColumns;
	ItemRow *mHead;
	ItemRow *mTail;
	class CatIndex *mIndex;
};

// --------------------------------------------------------------------

class File
{
  public:
	friend class parser;
	friend class Validator;

	File();
	File(std::istream &is, bool validate = false);
	File(const std::string &path, bool validate = false);
	File(File &&rhs);
	File(const File &rhs) = delete;
	File &operator=(const File &rhs) = delete;

	~File();

	void load(const std::string &p);
	void save(const std::string &p);

	void load(std::istream &is);
	void save(std::ostream &os);

	/// \brief Load only the data block \a datablock from the mmCIF file
	void load(std::istream &is, const std::string &datablock);

	void save(std::ostream &os, const std::vector<std::string> &order) { write(os, order); }
	void write(std::ostream &os, const std::vector<std::string> &order);

	void loadDictionary();                 // load the default dictionary, that is mmcifDdl in this case
	void loadDictionary(const char *dict); // load one of the compiled in dictionaries
	void loadDictionary(std::istream &is); // load dictionary from input stream

	bool isValid();
	void validateLinks() const;

	const Datablock &firstDatablock() const
	{
		if (mHead == nullptr)
			throw std::runtime_error("No datablocks in file");
		return *mHead;
	}

	Datablock &firstDatablock()
	{
		if (mHead == nullptr)
			throw std::runtime_error("No datablocks in file");
		return *mHead;
	}

	void append(Datablock *e);

	Datablock *get(const std::string &name) const;
	Datablock &operator[](const std::string &name);

	struct iterator
	{
		using iterator_category = std::forward_iterator_tag;
		using value_type = Datablock;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type *;
		using reference = value_type &;

		iterator(Datablock *db)
			: mCurrent(db)
		{
		}

		reference operator*() { return *mCurrent; }
		pointer operator->() { return mCurrent; }

		iterator &operator++();
		iterator operator++(int)
		{
			iterator result(*this);
			this->operator++();
			return result;
		}

		bool operator==(const iterator &rhs) const { return mCurrent == rhs.mCurrent; }
		bool operator!=(const iterator &rhs) const { return not(mCurrent == rhs.mCurrent); }

	  private:
		Datablock *mCurrent;
	};

	iterator begin() const;
	iterator end() const;

	bool empty() const { return mHead == nullptr; }

	const Validator &getValidator() const;
	void getTagOrder(std::vector<std::string> &tags) const;

  private:
	void setValidator(Validator *v);

	Datablock *mHead;
	Validator *mValidator;
};

// --------------------------------------------------------------------
// some postponed inlines

namespace detail
{

	template <typename T>
	inline bool AnyIsConditionImpl<T>::test(const Category &c, const Row &r) const
	{
		bool result = false;
		for (auto &f : c.fields())
		{
			try
			{
				if (r[f].as<valueType>() == mValue)
				{
					result = true;
					break;
				}
			}
			catch (...)
			{
			}
		}

		return result;
	}

	inline bool AnyMatchesConditionImpl::test(const Category &c, const Row &r) const
	{
		bool result = false;
		for (auto &f : c.fields())
		{
			try
			{
				if (std::regex_match(r[f].as<std::string>(), mRx))
				{
					result = true;
					break;
				}
			}
			catch (...)
			{
			}
		}

		return result;
	}

} // namespace detail

// these should be here, as I learned today

inline void swap(cif::Row &a, cif::Row &b)
{
	a.swap(b);
}

inline void swap(cif::detail::ItemReference &a, cif::detail::ItemReference &b)
{
	a.swap(b);
}

// --------------------------------------------------------------------

template <typename RowType, typename... Ts>
iterator_proxy<RowType, Ts...>::iterator_proxy(Category &cat, row_iterator pos, char const *const columns[N])
	: mCat(&cat)
	, mCBegin(pos)
	, mCEnd(cat.end())
{
	for (size_t i = 0; i < N; ++i)
		mCix[i] = mCat->getColumnIndex(columns[i]);
}

template <typename RowType, typename... Ts>
iterator_proxy<RowType, Ts...>::iterator_proxy(Category &cat, row_iterator pos, std::initializer_list<char const *> columns)
	: mCat(&cat)
	, mCBegin(pos)
	, mCEnd(cat.end())
{
	// static_assert(columns.size() == N, "The list of column names should be exactly the same as the list of requested columns");

	std::size_t i = 0;
	for (auto column : columns)
		mCix[i++] = mCat->getColumnIndex(column);
}

// --------------------------------------------------------------------

template <typename RowType, typename... Ts>
conditional_iterator_proxy<RowType, Ts...>::conditional_iterator_impl::conditional_iterator_impl(
	Category &cat, row_iterator pos, const Condition &cond, const std::array<size_t, N> &cix)
	: mCat(&cat)
	, mBegin(pos, cix)
	, mEnd(cat.end(), cix)
	, mCondition(&cond)
{
}

template <typename RowType, typename... Ts>
conditional_iterator_proxy<RowType, Ts...>::conditional_iterator_proxy(conditional_iterator_proxy &&p)
	: mCat(nullptr)
	, mCBegin(p.mCBegin)
	, mCEnd(p.mCEnd)
	, mCix(p.mCix)
{
	std::swap(mCat, p.mCat);
	std::swap(mCix, p.mCix);
	mCondition.swap(p.mCondition);
}

template <typename RowType, typename... Ts>
conditional_iterator_proxy<RowType, Ts...>::conditional_iterator_proxy(Category &cat, row_iterator pos, Condition &&cond)
	: mCat(&cat)
	, mCondition(std::move(cond))
	, mCBegin(pos)
	, mCEnd(cat.end())
{
	mCondition.prepare(cat);

	while (mCBegin != mCEnd and not mCondition(*mCat, mCBegin.row()))
		++mCBegin;
}

template <typename RowType, typename... Ts>
template <std::size_t TN, std::enable_if_t<TN != 0, bool>>
conditional_iterator_proxy<RowType, Ts...>::conditional_iterator_proxy(Category &cat, row_iterator pos, Condition &&cond, const char *const columns[N])
	: mCat(&cat)
	, mCondition(std::move(cond))
	, mCBegin(pos)
	, mCEnd(cat.end())
{
	mCondition.prepare(cat);

	while (mCBegin != mCEnd and not mCondition(*mCat, mCBegin.row()))
		++mCBegin;

	for (size_t i = 0; i < N; ++i)
		mCix[i] = mCat->getColumnIndex(columns[i]);
}

template <typename RowType, typename... Ts>
conditional_iterator_proxy<RowType, Ts...> &conditional_iterator_proxy<RowType, Ts...>::operator=(conditional_iterator_proxy &&p)
{
	swap(p);
	return *this;
}

template <typename RowType, typename... Ts>
typename conditional_iterator_proxy<RowType, Ts...>::iterator conditional_iterator_proxy<RowType, Ts...>::begin() const
{
	return iterator(*mCat, mCBegin, mCondition, mCix);
}

template <typename RowType, typename... Ts>
typename conditional_iterator_proxy<RowType, Ts...>::iterator conditional_iterator_proxy<RowType, Ts...>::end() const
{
	return iterator(*mCat, mCEnd, mCondition, mCix);
}

template <typename RowType, typename... Ts>
bool conditional_iterator_proxy<RowType, Ts...>::empty() const
{
	return mCBegin == mCEnd;
}

template <typename RowType, typename... Ts>
void conditional_iterator_proxy<RowType, Ts...>::swap(conditional_iterator_proxy &rhs)
{
	std::swap(mCat, rhs.mCat);
	mCondition.swap(rhs.mCondition);
	std::swap(mCBegin, rhs.mCBegin);
	std::swap(mCEnd, rhs.mCEnd);
	std::swap(mCix, rhs.mCix);
}

} // namespace cif
