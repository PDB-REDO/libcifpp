// cif parsing library

#include "libcif/cif++.h"

#include <boost/filesystem/path.hpp>

// the std regex of gcc is crashing....
#include <boost/regex.hpp>
#include <set>

namespace cif
{
	
struct validate_category;

// --------------------------------------------------------------------

class validation_error : public std::exception
{
  public:
	validation_error(const std::string& msg) : m_msg(msg) {}
	const char* what() const noexcept		{ return m_msg.c_str(); }
	std::string m_msg;
};

// --------------------------------------------------------------------

enum DDL_PrimitiveType
{
	ptChar, ptUChar, ptNumb
};

DDL_PrimitiveType map_to_primitive_type(const std::string& s);

struct validate_type
{
	std::string				m_name;
	DDL_PrimitiveType		m_primitive_type;
	boost::regex			m_rx;

	bool operator<(const validate_type& rhs) const
	{
		return icompare(m_name, rhs.m_name) < 0;
	}

	// compare values based on type	
//	int compare(const std::string& a, const std::string& b) const
//	{
//		return compare(a.c_str(), b.c_str());
//	}
	
	int compare(const char* a, const char* b) const;
};

struct validate_item
{
	std::string				m_tag;
	bool					m_mandatory;
	const validate_type*	m_type;
	cif::iset				m_enums;
	validate_item*			m_parent = nullptr;
	std::set<validate_item*>
							m_children;
	validate_category*		m_category = nullptr;
	std::set<validate_item*>
							m_foreign_keys;
	
	void set_parent(validate_item* parent);

	bool operator<(const validate_item& rhs) const
	{
		return icompare(m_tag, rhs.m_tag) < 0;
	}

	bool operator==(const validate_item& rhs) const
	{
		return iequals(m_tag, rhs.m_tag);
	}

	void operator()(std::string value) const;
};

struct validate_category
{
	std::string				m_name;
	std::vector<string>		m_keys;
	cif::iset				m_groups;
	cif::iset				m_mandatory_fields;
	std::set<validate_item>	m_item_validators;

	bool operator<(const validate_category& rhs) const
	{
		return icompare(m_name, rhs.m_name) < 0;
	}

	void add_item_validator(validate_item&& v);
	
	const validate_item* get_validator_for_item(std::string tag) const;
	
	const std::set<validate_item>& item_validators() const
	{
		return m_item_validators;
	}
};

// --------------------------------------------------------------------

class validator
{
  public:
	friend class dict_parser;

	validator();
	~validator();

	validator(const validator& rhs) = delete;
	validator& operator=(const validator& rhs) = delete;
	
	validator(validator&& rhs);
	validator& operator=(validator&& rhs);
	
	void add_type_validator(validate_type&& v);
	const validate_type* get_validator_for_type(std::string type_code) const;

	void add_category_validator(validate_category&& v);
	const validate_category* get_validator_for_category(std::string category) const;

	void report_error(const std::string& msg);
	
	std::string dict_name() const					{ return m_name; }
	void dict_name(const std::string& name)			{ m_name = name; }

	std::string dict_version() const				{ return m_version; }
	void dict_version(const std::string& version)	{ m_version = version; }

  private:

	// name is fully qualified here:
	validate_item* get_validator_for_item(std::string name) const;

	std::string					m_name;
	std::string					m_version;
	bool						m_strict = false;
//	std::set<uint32>			m_sub_categories;
	std::set<validate_type>		m_type_validators;
	std::set<validate_category>	m_category_validators;
};

}
