// cif parsing library

#pragma once

#include "libcif/config.h"

#include <regex>
#include <iostream>
#include <set>

#include <boost/lexical_cast.hpp>
#include <boost/any.hpp>

#include "cif-utils.h"

extern int VERBOSE;

/*
	Simple C++ interface to CIF files.
	
	Assumptions: a file contains one or more datablocks modelled by the class datablock.
	Each datablock contains categories. These map to the original tables used to fill
	the mmCIF file. Each category can contain multiple items, the columns in the table.
	
	Values are stored as character strings internally.
	
	Synopsis:
	
	// create a cif file
	
	cif::datablock e("1MVE");
	e.append(cif::category{"_entry", { "id", "1MVE" } });
	
	cif::category atom_site("atom_site");
	size_t nr{};
	for (my_atom: atoms)
	{
		atom_site.push_back({
			{ "group_PDB", "ATOM" },
			{ "id", ++nr },
			{ "type_symbol", my_atom.type.str() },
			...
		});
	}
	
	e.append(move(atom_site));
	
	cif::file f;
	f.append(e);
	
	ofstream os("1mve.cif");
	f.write(os);

	// read
	f.read(ifstream{"1mve.cif"});
	
	auto& e = f.first_datablock();
	
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

	Another way of querying a category is by using this construct:
	
	auto cat& = e["atom_site"];
	auto rows = cat.find(key("label_asym_id") == "A" and key("label_seq_id") == 1);


*/

namespace cif
{

using std::string;
using std::vector;

// mmCIF mapping
// A CIF data file in this case contains entries (data blocks) which can contain
// one or more category objects. Each category object contains arrays of items.
// Better, you can consider the categories as tables containing columns which
// are the items.

class file;
class datablock;
class category;
class row;			// a flyweight class that references data in categories
class item;
class validator;

struct validate_item;
struct validate_category;

struct item_column;
struct item_row;
struct item_value;

// --------------------------------------------------------------------
// class item
//
//	This class is only transient, it is used to construct new rows.
//	Access to already stored data is through an item_reference object.

class item
{
  public:
	typedef enum { not_applicable, not_defined, text, number } item_content_type;
	
	item() {}
	template<typename T>
	item(const string& name, const T& value);
	item(const item& rhs) : m_name(rhs.m_name), m_value(rhs.m_value) {}
	item(item&& rhs) : m_name(std::move(rhs.m_name)), m_value(std::move(rhs.m_value)) {}

	item& operator=(const item& rhs)
	{
		if (this != &rhs)
		{
			m_name = rhs.m_name;
			m_value = rhs.m_value;
		}
		
		return *this;
	}
	
	item& operator=(item&& rhs)
	{
		if (this != &rhs)
		{
			m_name = std::move(rhs.m_name);
			m_value = std::move(rhs.m_value);
		}
		
		return *this;
	}
	
	const string& name() const	{ return m_name; }
	const string& value() const	{ return m_value; }

	void value(const string& v)	{ m_value = v; }
	
	bool empty() const			{ return m_value.empty(); }
	size_t length() const		{ return m_value.length(); }
	const char* c_str() const	{ return m_value.c_str(); }
	
  private:
	string	m_name;
  	string	m_value;
};

template<typename T>
inline
item::item(const string& name, const T& value)
	: m_name(name), m_value(boost::lexical_cast<string>(value))
{	
}

template<>
inline
item::item(const string& name, const string& value)
	: m_name(name), m_value(value)
{
}

// --------------------------------------------------------------------
// class datablock acts as an STL container for category objects

class datablock
{
  public:
	friend class file;
	
	typedef std::list<category> category_list;
	typedef category_list::iterator iterator;
	typedef category_list::const_iterator const_iterator;
	
	datablock(const string& name);
	~datablock();

	datablock(const datablock&) = delete;
	datablock& operator=(const datablock&) = delete;

	string name() const								{ return m_name; }
	void set_name(const string& n)					{ m_name = n; }
	
	string first_item(const string& tag) const;

	iterator begin()		{ return m_categories.begin(); }
	iterator end()			{ return m_categories.end(); }

	const_iterator begin() const	{ return m_categories.begin(); }
	const_iterator end() const		{ return m_categories.end(); }

	category& operator[](const string& name);

	std::tuple<iterator,bool> emplace(const std::string& name);
	
	void validate();
	void set_validator(validator* v);

	// this one only looks up a category, returns nullptr if it does not exist
	category* get(const string& name);

	void get_tag_order(vector<string>& tags) const;

  private:

	void write(std::ostream& os);
	void write(std::ostream& os, const vector<string>& order);

	std::list<category>	m_categories;
	string				m_name;
	validator*			m_validator;
	datablock*			m_next;
};

// --------------------------------------------------------------------
// class row acts as a container for item objects, It has a more useful
// interface for accessing the contained columns. The get() method
// returns a row_result object that can be used to access only a subset
// of column values by index or by name.

namespace detail
{
	// item_reference is a helper class
	struct item_reference
	{
		const char*		m_name;
		item_row*		m_row;

		template<typename T>
		item_reference& operator=(const T& value)
		{
			this->operator=(boost::lexical_cast<string>(value));
			return *this;
		}
		
//		operator string() const	{ return c_str(); }
		
		template<typename T>
		T as() const
		{
			T result = 0;
			if (not empty())
				result = boost::lexical_cast<T>(c_str());
			return result;
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
		
		bool empty() const;
//		bool unapplicable() const;
		
		const char* c_str() const;
		
		bool operator!=(const string& s) const		{ return s != c_str(); }
		bool operator==(const string& s) const		{ return s == c_str(); }
	};

	template<>
	inline
	string item_reference::as<string>() const
	{
		return string(c_str());
	}
	
	template<>
	inline
	const char* item_reference::as<const char*>() const
	{
		return c_str();
	}
	
	template<>
	inline
	int item_reference::compare<string>(const string& value) const
	{
		return icompare(c_str(), value.c_str());
	}

	template<>
	inline
	int item_reference::compare(const char* const& value) const
	{
		return cif::icompare(c_str(), value);
	}
	
	inline std::ostream& operator<<(std::ostream& os, const item_reference& rhs)
	{
		os << rhs.c_str();
		return os;
	}

	template<>
	item_reference& item_reference::operator=(const string& value);

	// some helper classes to help create tuple result types
	
	template<typename...> struct tuple_catter;
	
	template<typename... Ts>
	struct tuple_catter<std::tuple<Ts...>>
	{
		typedef std::tuple<Ts...> type;
	};
	
	template<typename... T1s, typename... T2s, typename... Rem>
	struct tuple_catter<std::tuple<T1s...>, std::tuple<T2s...>, Rem...>
	{
		typedef typename tuple_catter<std::tuple<T1s..., T2s...>, Rem...>::type type;
	};
	
	template<typename...> struct col_getter;
	
	template<typename T>
	struct col_getter<T>
	{
		typedef std::tuple<const item_reference>	type;
		
		template<typename Res>
		static type get(Res& rs)
		{
			size_t index = Res::N - 1;
			return std::tuple<const item_reference>{ rs[index] };
		}
	};
	
	template<typename T, typename... Ts>
	struct col_getter<T, Ts...>
	{
		typedef col_getter<Ts...> next;
		typedef typename tuple_catter<std::tuple<const item_reference>, typename next::type>::type type;
		
		template<typename Res>
		static type get(Res& rs)
		{
			typedef col_getter<Ts...> next;
			size_t index = Res::N - 1 - sizeof...(Ts);
			return std::tuple_cat(std::tuple<const item_reference>{ rs[index]}, next::get(rs));
		}
	};

	template<typename... C>
	struct get_row_result
	{
		enum { N = sizeof...(C) };
		typedef typename col_getter<C...>::type tuple_type;
	
//		const item_reference operator[](const string& col) const
//		{
//			return m_row[col];
//		}
		
		const item_reference operator[](size_t ix) const
		{
			return m_row[m_columns[ix]];
		}
		
		get_row_result(row& r, C... columns)
			: m_row(r), m_columns({{columns...}}) {}
	
		row& m_row;
		std::array<const char*, N> m_columns;
	};
	
	// we want to be able to tie some variables to a row_result, for this we use tiewraps

	template<int IX, typename... Ts>
	struct tie_wrap;
	
	template<int IX, typename T>
	struct tie_wrap<IX,T>
	{
		tie_wrap(T& t)
			: m_val(t) {}
	
		template<typename Res>
		void operator=(const Res& rr)
		{
			typedef typename std::remove_reference<T>::type basic_type;

			const item_reference v = rr[IX];
			basic_type tv = v.as<basic_type>();
			m_val = tv;
		}
		
		T& 		m_val;
	};
	
	template<int IX, typename T, typename... Ts>
	struct tie_wrap<IX, T, Ts...>
	{
		typedef tie_wrap<IX + 1, Ts...> next;
	
		tie_wrap(T& t, Ts&... ts)
			: m_val(t), m_next(ts...) {}
	
		template<typename Res>
		void operator=(const Res& rr)
		{
			typedef typename std::remove_reference<T>::type basic_type;
			
			const item_reference v = rr[IX];
			basic_type tv = v.as<basic_type>();
			m_val = tv;

			m_next.operator=(rr);
		}
		
		T& 		m_val;
		next	m_next;
	};
}

template<typename... Ts>
auto tie(Ts&... v) -> detail::tie_wrap<0, Ts...>
{
	return detail::tie_wrap<0, Ts...>(v...);
}

class row
{
  public:
	friend class category;
	friend class cat_index;
	friend class row_comparator;
	friend struct detail::item_reference;

	row(item_row* data = nullptr) : m_data(data) {}
	row(const row& rhs);
	row& operator=(const row& rhs);
	
	struct const_iterator : public std::iterator<std::forward_iterator_tag, const item>
	{
		typedef std::iterator<std::forward_iterator_tag, item>	base_type;
		typedef typename base_type::pointer						pointer;
		typedef typename base_type::reference					reference;
		
		const_iterator(item_row* data, item_value* ptr);
		
		reference operator*()								{ return m_current; }
		pointer operator->()								{ return &m_current; }
		
		const_iterator& operator++();
		const_iterator operator++(int)						{ const_iterator result(*this); this->operator++(); return result; } 

		bool operator==(const const_iterator& rhs) const	{ return m_ptr == rhs.m_ptr; } 
		bool operator!=(const const_iterator& rhs) const	{ return m_ptr != rhs.m_ptr; } 
		
	  private:

		void fetch();

	  	item_row*	m_data;
		item_value*	m_ptr;
		item		m_current;
	};
	
	// checks for an initialized row:
	operator bool() const									{ return m_data != nullptr; }
	
	bool empty() const;
	const_iterator begin() const;
	const_iterator end() const;

// TODO: implement real const version?
	
	const detail::item_reference operator[](const char* item_tag) const
	{
		return detail::item_reference{item_tag, m_data};
	}

	detail::item_reference operator[](const char* item_tag)
	{
		return detail::item_reference{item_tag, m_data};
	}

	const detail::item_reference operator[](const string& item_tag) const
	{
		return detail::item_reference{item_tag.c_str(), m_data};
	}

	detail::item_reference operator[](const string& item_tag)
	{
		return detail::item_reference{item_tag.c_str(), m_data};
	}

	template<typename... C>
	auto get(C... columns) -> detail::get_row_result<C...>
	{
		return detail::get_row_result<C...>(*this, columns...);
	}
	
	bool operator==(const row& rhs) const
	{
		return m_data == rhs.m_data;
	}

	item_row* data() const							{ return m_data; }

	void swap(row& rhs)
	{
		std::swap(m_data, rhs.m_data);
	}
	
  private:

	void assign(const string& name, const string& value, bool emplacing);
	void assign(const item& i, bool emplacing);

	item_row*	m_data;
};

// swap for rows is defined below

// --------------------------------------------------------------------
// some more templates to be able to do querying

namespace detail
{

struct condition_impl
{
	virtual ~condition_impl() {}
	
	virtual bool test(const category& c, const row& r) const = 0;
	virtual std::string str() const = 0;
};

}

struct condition
{
	condition(detail::condition_impl* impl) : m_impl(impl) {}

	condition(condition&& rhs)
		: m_impl(nullptr)
	{
		std::swap(m_impl, rhs.m_impl);
	}
	
	condition& operator=(condition&& rhs)
	{
		std::swap(m_impl, rhs.m_impl);
		return *this;
	}

	~condition()
	{
		delete m_impl;
	}
	
	bool operator()(const category& c, const row& r) const
	{
		assert(m_impl);
		return m_impl->test(c, r);
	}
	
	std::string str() const
	{
		return m_impl->str();
	}

	detail::condition_impl*	m_impl;
};

namespace detail
{

template<typename T>
struct key_is_condition_impl : public condition_impl
{
	typedef T value_type;
	
	key_is_condition_impl(const string& item_tag, const value_type& value)
		: m_item_tag(item_tag), m_value(value) {}
	
	virtual bool test(const category& c, const row& r) const
	{
		return r[m_item_tag].template compare<value_type>(m_value) == 0;
	}
	
	virtual std::string str() const
	{
		return m_item_tag + " == " + boost::lexical_cast<std::string>(m_value);
	}
	
	string m_item_tag;
	value_type m_value;
};

template<typename T>
struct key_is_not_condition_impl : public condition_impl
{
	typedef T value_type;
	
	key_is_not_condition_impl(const string& item_tag, const value_type& value)
		: m_item_tag(item_tag), m_value(value) {}
	
	virtual bool test(const category& c, const row& r) const
	{
		return r[m_item_tag].template compare<value_type>(m_value) != 0;
	}
	
	virtual std::string str() const
	{
		return m_item_tag + " != " + boost::lexical_cast<std::string>(m_value);
	}
	
	string m_item_tag;
	value_type m_value;
};

template<typename COMP>
struct key_compare_condition_impl : public condition_impl
{
	key_compare_condition_impl(const string& item_tag, COMP&& comp)
		: m_item_tag(item_tag), m_comp(std::move(comp)) {}
	
	virtual bool test(const category& c, const row& r) const
	{
		return m_comp(c, r);
	}
	
	virtual std::string str() const
	{
		return m_item_tag + " compare " /*+ boost::lexical_cast<std::string>(m_value)*/;
	}
	
	string m_item_tag;
	COMP m_comp;
};

struct key_matches_condition_impl : public condition_impl
{
	key_matches_condition_impl(const string& item_tag, const std::regex& rx)
		: m_item_tag(item_tag), m_rx(rx) {}
	
	virtual bool test(const category& c, const row& r) const
	{
		return std::regex_match(r[m_item_tag].as<string>(), m_rx);
	}
	
	virtual std::string str() const
	{
		return m_item_tag + " ~= " + "<rx>";
	}
	
	string m_item_tag;
	std::regex m_rx;
};

template<typename T>
struct any_is_condition_impl : public condition_impl
{
	typedef T value_type;
	
	any_is_condition_impl(const value_type& value)
		: m_value(value) {}
	
	virtual bool test(const category& c, const row& r) const;

	virtual std::string str() const
	{
		return "any == " + boost::lexical_cast<std::string>(m_value);
	}
	
	value_type m_value;
};

struct any_matches_condition_impl : public condition_impl
{
	any_matches_condition_impl(const std::regex& rx)
		: m_rx(rx) {}
	
	virtual bool test(const category& c, const row& r) const;

	virtual std::string str() const
	{
		return "any ~= <rx>";
	}
	
	std::regex m_rx;
};

struct and_condition_impl : public condition_impl
{
	and_condition_impl(condition&& a, condition&& b)
		: m_a(nullptr), m_b(nullptr)
	{
		std::swap(m_a, a.m_impl);
		std::swap(m_b, b.m_impl);
	}
	
	~and_condition_impl()
	{
		delete m_a;
		delete m_b;
	}
	
	virtual bool test(const category& c, const row& r) const
	{
		return m_a->test(c, r) and m_b->test(c, r);
	}

	virtual std::string str() const
	{
		return "(" + m_a->str() + ") and (" + m_b->str() + ")";
	}
		
	condition_impl* m_a;
	condition_impl* m_b;
};

struct or_condition_impl : public condition_impl
{
	or_condition_impl(condition&& a, condition&& b)
		: m_a(nullptr), m_b(nullptr)
	{
		std::swap(m_a, a.m_impl);
		std::swap(m_b, b.m_impl);
	}
	
	~or_condition_impl()
	{
		delete m_a;
		delete m_b;
	}
	
	virtual bool test(const category& c, const row& r) const
	{
		return m_a->test(c, r) or m_b->test(c, r);
	}
		
	virtual std::string str() const
	{
		return "(" + m_a->str() + ") or (" + m_b->str() + ")";
	}
		
	condition_impl* m_a;
	condition_impl* m_b;
};

}

inline condition operator&&(condition&& a, condition&& b)
{
	return condition(new detail::and_condition_impl(std::move(a), std::move(b)));
}

inline condition operator||(condition&& a, condition&& b)
{
	return condition(new detail::or_condition_impl(std::move(a), std::move(b)));
}
	
struct key
{
	key(const string& item_tag) : m_item_tag(item_tag) {}
	key(const char* item_tag) : m_item_tag(item_tag) {}
	
	template<typename T>
	condition operator==(const T& v) const
	{
		return condition(new detail::key_is_condition_impl<T>(m_item_tag, v));
	}

	condition operator==(const char* v) const
	{
		string value(v ? v : "");
		return condition(new detail::key_is_condition_impl<std::string>(m_item_tag, value));
	}
	
	template<typename T>
	condition operator!=(const T& v) const
	{
		return condition(new detail::key_is_not_condition_impl<T>(m_item_tag, v));
	}

	condition operator!=(const char* v) const
	{
		string value(v ? v : "");
		return condition(new detail::key_is_not_condition_impl<std::string>(m_item_tag, value));
	}

	template<typename T>
	condition operator>(const T& v) const
	{
		auto comp = [this, v](const category& c, const row& r) -> bool { return r[this->m_item_tag].as<T>() > v; };
		return condition(new detail::key_compare_condition_impl<decltype(comp)>(m_item_tag, std::move(comp)));
	}

	template<typename T>
	condition operator>=(const T& v) const
	{
		auto comp = [this, v](const category& c, const row& r) -> bool { return r[this->m_item_tag].as<T>() >= v; };
		return condition(new detail::key_compare_condition_impl<decltype(comp)>(m_item_tag, std::move(comp)));
	}

	template<typename T>
	condition operator<(const T& v) const
	{
		auto comp = [this, v](const category& c, const row& r) -> bool { return r[this->m_item_tag].as<T>() < v; };
		return condition(new detail::key_compare_condition_impl<decltype(comp)>(m_item_tag, std::move(comp)));
	}

	template<typename T>
	condition operator<=(const T& v) const
	{
		auto comp = [this, v](const category& c, const row& r) -> bool { return r[this->m_item_tag].as<T>() <= v; };
		return condition(new detail::key_compare_condition_impl<decltype(comp)>(m_item_tag, std::move(comp)));
	}
	
	string m_item_tag;
};

template<>
inline
condition key::operator==(const std::regex& rx) const
{
	return condition(new detail::key_matches_condition_impl(m_item_tag, rx));
}

struct any
{
	template<typename T>
	condition operator==(const T& v) const
	{
		return condition(new detail::any_is_condition_impl<T>(v));
	}
};

template<>
inline
condition any::operator==(const std::regex& rx) const
{
	return condition(new detail::any_matches_condition_impl(rx));
}

// --------------------------------------------------------------------
// class rowset is used to return find results. Use it to re-order the results
// or to group them 

class rowset : public vector<row>
{
  public:
	rowset(category& cat);
	
	rowset& orderBy(const string& item)
		{ return orderBy({ item }); }
	
	rowset& orderBy(std::initializer_list<string> items);

  private:
	category&	m_cat;
};

// --------------------------------------------------------------------
// class category acts as an STL container for row objects 

class category
{
  public:
	friend class datablock;
	friend class row;
	friend struct detail::item_reference;

	category(datablock& db, const string& name, validator* validator);
	category(const category&) = delete;
	category& operator=(const category&) = delete;
	~category();

	const string name() const						{ return m_name; }
	
	const detail::item_reference get_first_item(const char* item_name) const;

	struct iterator : public std::iterator<std::forward_iterator_tag, row>
	{
		friend class category;
		
		typedef std::iterator<std::forward_iterator_tag, row>	base_type;
		typedef typename base_type::pointer						pointer;
		typedef typename base_type::reference					reference;
		
		iterator(item_row* data) : m_current(data) {}
		
		reference operator*()						{ return m_current; }
		pointer operator->()						{ return &m_current; }
		
		iterator& operator++();
		iterator operator++(int)					{ iterator result(*this); this->operator++(); return result; } 

		bool operator==(const iterator& rhs) const	{ return m_current == rhs.m_current; } 
		bool operator!=(const iterator& rhs) const	{ return not (m_current == rhs.m_current); } 
		
	  private:
		row		m_current;
	};
	
	iterator begin();
	iterator end();

	bool empty() const;
	size_t size() const;
	
	void clear();
	
	row front()										{ return row(m_head); }
	row back()										{ return row(m_tail); }
	
	row operator[](condition&& cond);
	rowset find(condition&& cond);
	bool exists(condition&& cond);
	
	rowset orderBy(const string& item)
		{ return orderBy({ item }); }
	
	rowset orderBy(std::initializer_list<string> items);
	
	std::tuple<row,bool> emplace(item value)		{ return emplace({ value }); }

	std::tuple<row,bool> emplace(std::initializer_list<item> values)
		{ return emplace(values.begin(), values.end()); }

	std::tuple<row,bool> emplace(row r);
	
	template<class Iter>
	std::tuple<row,bool> emplace(Iter b, Iter e);

	void erase(condition&& cond);
	void erase(row r);
	void erase(iterator ri);

	void validate();

	const validator& get_validator() const;
	const validate_category* get_cat_validator() const		{ return m_cat_validator; }
	
	void set_validator(validator* v);

	iset fields() const;
	iset mandatory_fields() const;
	iset key_fields() const;
	
	void drop(const string& field);

	void get_tag_order(vector<string>& tags) const;

	// return index for known column, or the next available column index
	size_t get_column_index(const string& name) const;
	const string& get_column_name(size_t column_index) const;

	void reorderByIndex();

  private:

	void write(std::ostream& os);
	void write(std::ostream& os, const vector<string>& order);
	void write(std::ostream& os, const vector<int>& order, bool includeEmptyColumns);

	size_t add_column(const string& name);
	
	datablock&			m_db;
	string				m_name;
	validator*			m_validator;
	const validate_category*	m_cat_validator = nullptr;
	vector<item_column>	m_columns;
	item_row*			m_head;
	item_row*			m_tail;
	class cat_index*	m_index;
};

// --------------------------------------------------------------------

class file
{
  public:
	friend class parser;
	friend class validator;

	file();
	file(std::istream& is, bool validate = false);
	file(file&& rhs);
	file(const file& rhs) = delete;
	file& operator=(const file& rhs) = delete;
	
	~file();

	void load(std::istream& is);
	void save(std::ostream& os);

	void save(std::ostream& os, const vector<string>& order)	{ write(os, order); }
	void write(std::ostream& os, const vector<string>& order);

	void load_dictionary();						// load the default dictionary, that is mmcif_ddl in this case
	void load_dictionary(const char* dict);		// load one of the compiled in dictionaries 
	void load_dictionary(std::istream& is);		// load dictionary from input stream

	void validate();
	
	datablock& first_datablock()			{ return *m_head; }
	void append(datablock* e);
	
	datablock& operator[](const string& name);

	struct iterator : public std::iterator<std::forward_iterator_tag, datablock>
	{
		typedef std::iterator<std::forward_iterator_tag, datablock>	base_type;
		typedef typename base_type::pointer							pointer;
		typedef typename base_type::reference						reference;
		
		iterator(datablock* db) : m_current(db) {}
		
		reference operator*()						{ return *m_current; }
		pointer operator->()						{ return m_current; }
		
		iterator& operator++();
		iterator operator++(int)					{ iterator result(*this); this->operator++(); return result; } 

		bool operator==(const iterator& rhs) const	{ return m_current == rhs.m_current; } 
		bool operator!=(const iterator& rhs) const	{ return not (m_current == rhs.m_current); } 
		
	  private:
		datablock*		m_current;
	};
	
	iterator begin() const;
	iterator end() const;
	
	const validator& get_validator() const;
	void get_tag_order(vector<string>& tags) const;
	
  private:

	void set_validator(validator* v);

	datablock*	m_head;
	validator*	m_validator;
};

// --------------------------------------------------------------------
// some postponed inlines

namespace detail
{

template<typename T>
inline
bool any_is_condition_impl<T>::test(const category& c, const row& r) const
{
	bool result = false;
	for (auto& f: c.fields())
	{
		try
		{
			if (r[f].as<value_type>() == m_value)
			{
				result = true;
				break;
			}
		}
		catch (...) {}
	}
	
	return result;
}

inline bool any_matches_condition_impl::test(const category& c, const row& r) const
{
	bool result = false;
	for (auto& f: c.fields())
	{
		try
		{
			if (std::regex_match(r[f].as<string>(), m_rx))
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

}

namespace std
{

template<>
inline void swap(cif::row& a, cif::row& b)
{
	a.swap(b);
}

}

