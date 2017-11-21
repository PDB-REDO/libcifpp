// cif parsing library

#include <boost/algorithm/string.hpp>

// since gcc's regex is crashing....
#include <boost/regex.hpp>

#include "libcif/cif++.h"
#include "libcif/cif-parser.h"
#include "libcif/cif-validator.h"

using namespace std;
namespace ba = boost::algorithm;

extern int VERBOSE;

namespace cif
{

DDL_PrimitiveType map_to_primitive_type(const string& s)
{
	DDL_PrimitiveType result;
	if (iequals(s, "char"))
		result = ptChar;
	else if (iequals(s, "uchar"))
		result = ptUChar;
	else if (iequals(s, "numb"))
		result = ptNumb;
	else
		throw validation_error("Not a known primitive type");
	return result;
}

// --------------------------------------------------------------------

int validate_type::compare(const char* a, const char* b) const
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
			switch (m_primitive_type)
			{
				case ptNumb:
				{
					double da = strtod(a, nullptr);
					double db = strtod(b, nullptr);
					
					auto d = da - db;
					if (abs(d) > numeric_limits<double>::epsilon())
					{
						if (d > 0)
							result = 1;
						else if (d < 0)
							result = -1;
					}
					break;
				}
				
				case ptUChar:
				case ptChar:
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
						
						char ca = toupper(*ai);
						char cb = toupper(*bi);
						
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
		catch (const std::invalid_argument& ex)
		{
			result = 1;
		}
	}
	
	return result;
}

// --------------------------------------------------------------------

void validate_item::set_parent(validate_item* parent)
{
	m_parent = parent;

	if (m_type == nullptr and m_parent != nullptr)
		m_type = m_parent->m_type;
		
	if (m_parent != nullptr)
	{
		m_parent->m_children.insert(this);
	
		if (m_category->m_keys == vector<string>{m_tag})
			m_parent->m_foreign_keys.insert(this);
	}
}

void validate_item::operator()(string value) const
{
	if (VERBOSE >= 4)
		cout << "validating '" << value << "' for '" << m_tag << "'" << endl;

	if (not value.empty() and value != "?" and value != ".")
	{
		if (m_type != nullptr and not boost::regex_match(value, m_type->m_rx))
			throw validation_error("Value '" + value + "' does not match type expression for type " + m_type->m_name + " in item " + m_tag);

		if (not m_enums.empty())
		{
			if (m_enums.count(value) == 0)
				throw validation_error("Value '" + value + "' is not in the list of allowed values for item " + m_tag);
		}
	}
}

// --------------------------------------------------------------------

void validate_category::add_item_validator(validate_item&& v)
{
	if (v.m_mandatory)
		m_mandatory_fields.insert(v.m_tag);

	v.m_category = this;

	auto r = m_item_validators.insert(move(v));
	if (not r.second and VERBOSE >= 4)
		cout << "Could not add validator for item " << v.m_tag << " to category " << m_name << endl;
}

const validate_item* validate_category::get_validator_for_item(string tag) const
{
	const validate_item* result = nullptr;
	auto i = m_item_validators.find(validate_item{tag});
	if (i != m_item_validators.end())
		result = &*i;
	else if (VERBOSE > 4)
		cout << "No validator for tag " << tag << endl;
	return result;
}

// --------------------------------------------------------------------

validator::validator()
{
}

validator::~validator()
{
}

void validator::add_type_validator(validate_type&& v)
{
	auto r = m_type_validators.insert(move(v));
	if (not r.second and VERBOSE > 4)
		cout << "Could not add validator for type " << v.m_name << endl;
}

const validate_type* validator::get_validator_for_type(string type_code) const
{
	const validate_type* result = nullptr;
	
	auto i = m_type_validators.find(validate_type{ type_code, ptChar, boost::regex() });
	if (i != m_type_validators.end())
		result = &*i;
	else if (VERBOSE > 4)
		cout << "No validator for type " << type_code << endl;
	return result;
}

void validator::add_category_validator(validate_category&& v)
{
	auto r = m_category_validators.insert(move(v));
	if (not r.second and VERBOSE > 4)
		cout << "Could not add validator for category " << v.m_name << endl;
}

const validate_category* validator::get_validator_for_category(string category) const
{
	const validate_category* result = nullptr;
	auto i = m_category_validators.find(validate_category{category});
	if (i != m_category_validators.end())
		result = &*i;
	else if (VERBOSE > 4)
		cout << "No validator for category " << category << endl;
	return result;
}

validate_item* validator::get_validator_for_item(string tag) const
{
	validate_item* result = nullptr;
	
	string cat, item;
	std::tie(cat, item) = split_tag_name(tag);

	auto* cv = get_validator_for_category(cat);
	if (cv != nullptr)
		result = const_cast<validate_item*>(cv->get_validator_for_item(item));

	if (result == nullptr and VERBOSE > 4)
		cout << "No validator for item " << tag << endl;

	return result;
}

void validator::report_error(const string& msg)
{
	if (m_strict)
		throw validation_error(msg);
	else if (VERBOSE)
		cerr << msg << endl;
}


}
