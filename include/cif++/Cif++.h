// cif parsing library

#pragma once

#include "cif++/Config.h"

#include <regex>
#include <iostream>
#include <set>
#include <list>

#include <boost/lexical_cast.hpp>
#include <boost/any.hpp>
#include <boost/filesystem/path.hpp>

#include "cif++/CifUtils.h"

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

using std::string;
using std::vector;

// mmCIF mapping
// A CIF data file in this case contains entries (data blocks) which can contain
// one or more Category objects. Each Category object contains arrays of Items.
// Better, you can consider the categories as tables containing columns which
// are the Items.

class File;
class Datablock;
class Category;
class Row;			// a flyweight class that references data in categories
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
	template<typename T>
	Item(const string& name, const T& value);
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
	
	const string& name() const	{ return mName; }
	const string& value() const	{ return mValue; }

	void value(const string& v)	{ mValue = v; }
	
	// empty means either null or unknown
	bool empty() const;

	// is_null means the field contains '.'
	bool is_null() const;

	// is_unknown means the field contains '?'
	bool is_unknown() const;

	size_t length() const		{ return mValue.length(); }
	const char* c_str() const	{ return mValue.c_str(); }
	
  private:
	string	mName;
  	string	mValue;
};

template<typename T>
inline
Item::Item(const string& name, const T& value)
	: mName(name), mValue(boost::lexical_cast<string>(value))
{	
}

template<>
inline
Item::Item(const string& name, const string& value)
	: mName(name), mValue(value)
{
}

// --------------------------------------------------------------------
// class datablock acts as an STL container for Category objects

class Datablock
{
  public:
	friend class File;
	
	using CategoryList = std::list<Category>;
	using iterator = CategoryList::iterator;
	using const_iterator = CategoryList::const_iterator;
	
	Datablock(const string& name);
	~Datablock();

	Datablock(const Datablock&) = delete;
	Datablock& operator=(const Datablock&) = delete;

	string getName() const							{ return mName; }
	void setName(const string& n)					{ mName = n; }
	
	string firstItem(const string& tag) const;

	iterator begin()		{ return mCategories.begin(); }
	iterator end()			{ return mCategories.end(); }

	const_iterator begin() const	{ return mCategories.begin(); }
	const_iterator end() const		{ return mCategories.end(); }

	Category& operator[](const string& name);

	std::tuple<iterator,bool> emplace(const std::string& name);
	
	bool isValid();
	void validateLinks() const;

	void setValidator(Validator* v);

	// this one only looks up a Category, returns nullptr if it does not exist
	Category* get(const string& name);

	void getTagOrder(vector<string>& tags) const;
	void write(std::ostream& os, const vector<string>& order);
	void write(std::ostream& os);
	
	// convenience function, add a line to the software category
	void add_software(const std::string& name, const std::string& classification,
		const std::string& versionNr, const std::string& versionDate);

  private:

	std::list<Category>	mCategories;
	string				mName;
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

		template<typename T>
		ItemReference& operator=(const T& value)
		{
			this->operator=(boost::lexical_cast<string>(value));
			return *this;
		}
		
		void swap(ItemReference& b);
		
//		operator string() const	{ return c_str(); }
		
		template<typename T>
		T as() const
		{
			return boost::lexical_cast<T>(c_str("0"));
		}
		
		template<typename T>
		int compare(const T& value) const
		{
			int result = 0;
			try
			{
				double v = boost::lexical_cast<T>(c_str());
				if (v < value)
					result = -1;
				else if (v > value)
					result = 1;
			}
			catch (...)
			{
				if (VERBOSE)
					std::cerr << "conversion error in compare for '" << c_str() << '\'' << std::endl;
				result = 1;
			}
			
			return result;
		}
		
		// empty means either null or unknown
		bool empty() const;

		// is_null means the field contains '.'
		bool is_null() const;

		// is_unknown means the field contains '?'
		bool is_unknown() const;

		const char* c_str() const;
		
		// the following returns the defaultValue from either the parameter
		// or, if specified, the value from _item_default.value in the dictionary
		const char* c_str(const char* defaultValue) const;
		
		bool operator!=(const string& s) const		{ return s != c_str(); }
		bool operator==(const string& s) const		{ return s == c_str(); }

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

	template<>
	inline
	string ItemReference::as<string>() const
	{
		return string(c_str(""));
	}
	
	template<>
	inline
	const char* ItemReference::as<const char*>() const
	{
		return c_str("");
	}
	
	template<>
	inline
	int ItemReference::compare<string>(const string& value) const
	{
		return icompare(c_str(), value.c_str());
	}

	template<>
	inline
	int ItemReference::compare(const char* const& value) const
	{
		return cif::icompare(c_str(), value);
	}
	
	inline std::ostream& operator<<(std::ostream& os, const ItemReference& rhs)
	{
		os << rhs.c_str();
		return os;
	}

	template<>
	ItemReference& ItemReference::operator=(const string& value);

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

	Row(ItemRow* data = nullptr, bool cascadeUpdate = true)
		: mData(data), mCascadeUpdate(cascadeUpdate) {}

	Row(const ItemRow* data)
		: Row(const_cast<ItemRow*>(data), false)
	{}

	Row(const Row& rhs);
	Row& operator=(const Row& rhs);

	void next();	///< make this row point to the next ItemRow

	/// When updating a value, you might want to change linked records as well
	/// But not always.
	void setCascadeUpdate(bool cascadeUpdate)
	{
		mCascadeUpdate = cascadeUpdate;
	}
	
	void setCascadeDelete(bool cascadeDelete)
	{
		mCascadeDelete = cascadeDelete;
	}

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
	operator bool() const									{ return mData != nullptr; }
	
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

	const detail::ItemReference operator[](const string& itemTag) const
	{
		size_t column = ColumnForItemTag(itemTag.c_str());
		return detail::ItemReference(itemTag.c_str(), column, *this);
	}

	detail::ItemReference operator[](const string& itemTag)
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

	ItemRow* data() const							{ return mData; }

	void swap(Row& rhs)
	{
		std::swap(mData, rhs.mData);
	}

	friend std::ostream& operator<<(std::ostream& os, const Row& row);

  private:

	void assign(const string& name, const string& value, bool updateLinked);
	void assign(size_t column, const string& value, bool updateLinked);
	void assign(const Item& i, bool updateLinked);
	
	static void swap(size_t column, ItemRow* a, ItemRow* b);

	size_t ColumnForItemTag(const char* itemTag) const;

	ItemRow*	mData;
	uint32_t	mLineNr = 0;
	bool		mCascadeUpdate = true;
	bool		mCascadeDelete = true;
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
	virtual std::string str() const = 0;
};

struct AllConditionImpl : public ConditionImpl
{
	virtual bool test(const Category& c, const Row& r) const { return true; }
	virtual std::string str() const { return "ALL"; }
};

}

struct Condition
{
	Condition() : mImpl(nullptr) {}
	Condition(detail::ConditionImpl* impl) : mImpl(impl) {}

	Condition(Condition&& rhs)
		: mImpl(nullptr)
	{
		std::swap(mImpl, rhs.mImpl);
	}
	
	Condition& operator=(Condition&& rhs)
	{
		std::swap(mImpl, rhs.mImpl);
		return *this;
	}

	~Condition()
	{
		delete mImpl;
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
	
	std::string str() const
	{
		return mImpl ? mImpl->str() : "";
	}

	bool empty() const		{ return mImpl == nullptr; }

	detail::ConditionImpl*	mImpl;
	bool mPrepared = false;
};

namespace detail
{

struct KeyIsEmptyConditionImpl : public ConditionImpl
{
	KeyIsEmptyConditionImpl(const string& ItemTag)
		: mItemTag(ItemTag) {}

	virtual void prepare(const Category& c);
	
	virtual bool test(const Category& c, const Row& r) const
	{
		return r[mItemIx].empty();
	}
	
	virtual std::string str() const
	{
		return mItemTag + " == <empty>";
	}
	
	string mItemTag;
	size_t mItemIx;
};

template<typename T>
struct KeyIsConditionImpl : public ConditionImpl
{
	typedef T valueType;
	
	KeyIsConditionImpl(const string& ItemTag, const valueType& value)
		: mItemTag(ItemTag), mValue(value) {}
	
	virtual void prepare(const Category& c);
	
	virtual bool test(const Category& c, const Row& r) const
	{
		return r[mItemIx].template compare<valueType>(mValue) == 0;
	}
	
	virtual std::string str() const
	{
		return mItemTag + " == " + boost::lexical_cast<std::string>(mValue);
	}
	
	string mItemTag;
	size_t mItemIx;
	valueType mValue;
};

template<typename T>
struct KeyIsNotConditionImpl : public ConditionImpl
{
	typedef T valueType;
	
	KeyIsNotConditionImpl(const string& ItemTag, const valueType& value)
		: mItemTag(ItemTag), mValue(value) {}
	
	virtual void prepare(const Category& c);
	
	virtual bool test(const Category& c, const Row& r) const
	{
		return r[mItemIx].template compare<valueType>(mValue) != 0;
	}
	
	virtual std::string str() const
	{
		return mItemTag + " != " + boost::lexical_cast<std::string>(mValue);
	}
	
	string mItemTag;
	size_t mItemIx;
	valueType mValue;
};

template<typename COMP>
struct KeyCompareConditionImpl : public ConditionImpl
{
	KeyCompareConditionImpl(const string& ItemTag, COMP&& comp)
		: mItemTag(ItemTag), mComp(std::move(comp)) {}
	
	virtual void prepare(const Category& c);
	
	virtual bool test(const Category& c, const Row& r) const
	{
		return mComp(c, r);
	}
	
	virtual std::string str() const
	{
		return mItemTag + " compare " /*+ boost::lexical_cast<std::string>(mValue)*/;
	}
	
	string mItemTag;
	size_t mItemIx;
	COMP mComp;
};

struct KeyMatchesConditionImpl : public ConditionImpl
{
	KeyMatchesConditionImpl(const string& ItemTag, const std::regex& rx)
		: mItemTag(ItemTag), mRx(rx) {}
	
	virtual void prepare(const Category& c);
	
	virtual bool test(const Category& c, const Row& r) const
	{
		return std::regex_match(r[mItemIx].as<string>(), mRx);
	}
	
	virtual std::string str() const
	{
		return mItemTag + " ~= " + "<rx>";
	}
	
	string mItemTag;
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

	virtual std::string str() const
	{
		return "any == " + boost::lexical_cast<std::string>(mValue);
	}
	
	valueType mValue;
};

struct anyMatchesConditionImpl : public ConditionImpl
{
	anyMatchesConditionImpl(const std::regex& rx)
		: mRx(rx) {}
	
	virtual bool test(const Category& c, const Row& r) const;

	virtual std::string str() const
	{
		return "any ~= <rx>";
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

	virtual std::string str() const
	{
		return "ALL";
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

	virtual std::string str() const
	{
		return "(" + mA->str() + ") and (" + mB->str() + ")";
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
		
	virtual std::string str() const
	{
		return "(" + mA->str() + ") or (" + mB->str() + ")";
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
		
	virtual std::string str() const
	{
		return "NOT (" + mA->str() + ")";
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

inline
std::ostream& operator<<(std::ostream& os, const Condition& rhs)
{
	os << rhs.str();
	return os;
}

struct Empty {};

struct Key
{
	Key(const string& ItemTag) : mItemTag(ItemTag) {}
	Key(const char* ItemTag) : mItemTag(ItemTag) {}
	
	template<typename T>
	Condition operator==(const T& v) const
	{
		return Condition(new detail::KeyIsConditionImpl<T>(mItemTag, v));
	}

	Condition operator==(const char* v) const
	{
		string value(v ? v : "");
		return Condition(new detail::KeyIsConditionImpl<std::string>(mItemTag, value));
	}

	Condition operator==(const detail::ItemReference& v) const
	{
		if (v.empty())
			return Condition(new detail::KeyIsEmptyConditionImpl(mItemTag));
		else
			return Condition(new detail::KeyIsConditionImpl<std::string>(mItemTag, v.c_str()));
	}
	
	Condition operator==(const Empty&) const
	{
		return Condition(new detail::KeyIsEmptyConditionImpl(mItemTag));
	}
	
	template<typename T>
	Condition operator!=(const T& v) const
	{
		return Condition(new detail::KeyIsNotConditionImpl<T>(mItemTag, v));
	}

	Condition operator!=(const char* v) const
	{
		string value(v ? v : "");
		return Condition(new detail::KeyIsNotConditionImpl<std::string>(mItemTag, value));
	}

	template<typename T>
	Condition operator>(const T& v) const
	{
		auto comp = [this, v](const Category& c, const Row& r) -> bool { return r[this->mItemTag].as<T>() > v; };
		return Condition(new detail::KeyCompareConditionImpl<decltype(comp)>(mItemTag, std::move(comp)));
	}

	template<typename T>
	Condition operator>=(const T& v) const
	{
		auto comp = [this, v](const Category& c, const Row& r) -> bool { return r[this->mItemTag].as<T>() >= v; };
		return Condition(new detail::KeyCompareConditionImpl<decltype(comp)>(mItemTag, std::move(comp)));
	}

	template<typename T>
	Condition operator<(const T& v) const
	{
		auto comp = [this, v](const Category& c, const Row& r) -> bool { return r[this->mItemTag].as<T>() < v; };
		return Condition(new detail::KeyCompareConditionImpl<decltype(comp)>(mItemTag, std::move(comp)));
	}

	template<typename T>
	Condition operator<=(const T& v) const
	{
		auto comp = [this, v](const Category& c, const Row& r) -> bool { return r[this->mItemTag].as<T>() <= v; };
		return Condition(new detail::KeyCompareConditionImpl<decltype(comp)>(mItemTag, std::move(comp)));
	}
	
	string mItemTag;
};

template<>
inline
Condition Key::operator==(const std::regex& rx) const
{
	return Condition(new detail::KeyMatchesConditionImpl(mItemTag, rx));
}

template<>
inline
Condition Key::operator==(const Empty&) const
{
	return Condition(new detail::KeyIsEmptyConditionImpl(mItemTag));
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

template<typename RowType, typename ConditionType = void>
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

	template<typename T, std::enable_if_t<std::is_same_v<ConditionType,T>, int> = 0>
	iterator_impl(ItemRow* data, Category& cat, const T& cond)
		: mCurrent(data), mCat(&cat), mCondition(&cond)
	{
		skip();
	}

	virtual ~iterator_impl() = default;

	reference operator*()							{ return mCurrent; }
	pointer operator->()							{ return &mCurrent; }
	
	iterator_impl& operator++()
	{
		mCurrent.next();
		skip();

		return *this;
	}

	iterator_impl operator++(int)					{ iterator_impl result(*this); this->operator++(); return result; } 

	bool operator==(const iterator_impl& rhs) const	{ return mCurrent == rhs.mCurrent; } 
	bool operator!=(const iterator_impl& rhs) const	{ return not (mCurrent == rhs.mCurrent); } 

  private:

	void skip()
	{
		if constexpr (not std::is_void_v<ConditionType>)
		{
			while (mCurrent and not (*mCondition)(*mCat, mCurrent))
				mCurrent.next();
		}
	}

	Row mCurrent;
	Category* mCat = nullptr;
	const Condition* mCondition = nullptr;
};

template<typename RowType>
class conditional_iterator_proxy
{
  public:
	using iterator = iterator_impl<RowType,Condition>;
	using reference = typename iterator::reference;

	conditional_iterator_proxy(ItemRow* head, Category& cat, Condition&& cond)
		: mHead(head), mCat(cat), mCondition(std::forward<Condition>(cond))
	{
		mCondition.prepare(cat);
	}

	iterator begin() const			{ return iterator(mHead, mCat, mCondition); }
	iterator end() const			{ return iterator(nullptr, mCat, mCondition); }

	bool empty() { return begin() == end(); }
	size_t size() const { return std::distance(begin(), end()); }

	reference front() { return *begin(); }

	Category& category() const		{ return mCat;}

  private:
	ItemRow* mHead;
	Category& mCat;
	Condition mCondition;
};

// --------------------------------------------------------------------
// class RowSet is used to return find results. Use it to re-order the results
// or to group them 

class RowSet : public vector<Row>
{
	typedef vector<Row>	base_type;

  public:
	RowSet(const RowSet& rhs);
	RowSet(RowSet&& rhs);
	RowSet(conditional_iterator_proxy<Row>&& results);

	RowSet& operator=(const RowSet& rhs);
	RowSet& operator=(RowSet&& rhs);

	RowSet& operator=(conditional_iterator_proxy<Row>& results);
	RowSet& operator=(conditional_iterator_proxy<const Row>& results);

	RowSet(Category& cat);
	
	RowSet& orderBy(const string& Item)
		{ return orderBy({ Item }); }
	
	RowSet& orderBy(std::initializer_list<string> Items);

  private:
	Category*	mCat;
};

// --------------------------------------------------------------------
// class Category acts as an STL container for Row objects 

class Category
{
  public:
	friend class Datablock;
	friend class Row;
	friend class detail::ItemReference;

	Category(Datablock& db, const string& name, Validator* Validator);
	Category(const Category&) = delete;
	Category& operator=(const Category&) = delete;
	~Category();

	const string name() const						{ return mName; }

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
		return conditional_iterator_proxy<Row>(mHead, *this, std::forward<Condition>(cond));
	}

	bool exists(Condition&& cond) const;
	
	RowSet orderBy(const string& Item)
		{ return orderBy({ Item }); }
	
	RowSet orderBy(std::initializer_list<string> Items);
	
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
	
	void drop(const string& field);

	void getTagOrder(vector<string>& tags) const;

	// return index for known column, or the next available column index
	size_t getColumnIndex(const string& name) const;
    bool hasColumn(const string& name) const;
	const string& getColumnName(size_t columnIndex) const;
	vector<string> getColumnNames() const;

	void reorderByIndex();
	void sort(std::function<int(const Row&, const Row&)> comparator);

  private:

	void write(std::ostream& os);
	void write(std::ostream& os, const vector<string>& order);
	void write(std::ostream& os, const vector<int>& order, bool includeEmptyColumns);

	size_t addColumn(const string& name);
	
	Datablock&			mDb;
	string				mName;
	Validator*			mValidator;
	const ValidateCategory*	mCatValidator = nullptr;
	vector<ItemColumn>	mColumns;
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
	File(boost::filesystem::path p, bool validate = false);
	File(File&& rhs);
	File(const File& rhs) = delete;
	File& operator=(const File& rhs) = delete;
	
	~File();

	void load(boost::filesystem::path p);
	void save(boost::filesystem::path p);

	void load(std::istream& is);
	void save(std::ostream& os);

	void save(std::ostream& os, const vector<string>& order)	{ write(os, order); }
	void write(std::ostream& os, const vector<string>& order);

	void loadDictionary();						// load the default dictionary, that is mmcifDdl in this case
	void loadDictionary(const char* dict);		// load one of the compiled in dictionaries 
	void loadDictionary(std::istream& is);		// load dictionary from input stream

	bool isValid();
	void validateLinks() const;
	
	Datablock& firstDatablock()			{ return *mHead; }
	void append(Datablock* e);
	
	Datablock* get(const string& name) const;
	Datablock& operator[](const string& name);

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
	
	const Validator& getValidator() const;
	void getTagOrder(vector<string>& tags) const;
	
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
			if (std::regex_match(r[f].as<string>(), mRx))
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

namespace detail
{

template<typename T>
void KeyIsConditionImpl<T>::prepare(const Category& c)
{
	mItemIx = c.getColumnIndex(mItemTag);
}

template<typename T>
void KeyIsNotConditionImpl<T>::prepare(const Category& c)
{
	mItemIx = c.getColumnIndex(mItemTag);
}

template<typename T>
void KeyCompareConditionImpl<T>::prepare(const Category& c)
{
	mItemIx = c.getColumnIndex(mItemTag);
}

}

}

