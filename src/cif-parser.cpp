// cif parsing library

#include <set>

#include <boost/algorithm/string.hpp>

#include "libcif/cif++.h"
#include "libcif/cif-parser.h"
#include "libcif/cif-validator.h"

using namespace std;
namespace ba = boost::algorithm;

extern int VERBOSE;

namespace cif
{

const uint32 kMaxLineLength = 132;

const uint8 kCharTraitsTable[128] = {
	//	0	1	2	3	4	5	6	7	8	9	a	b	c	d	e	f
		14,	15,	14,	14,	14,	15,	15,	14,	15,	15,	15,	15,	15,	15,	15,	15,	//	2
		15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	10,	15,	15,	15,	15,	//	3
		15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	//	4
		15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	14,	15,	14,	15,	14,	//	5
		15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	//	6
		15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	0,	//	7
};

// --------------------------------------------------------------------

cif_parser_error::cif_parser_error(uint32 line_nr, const string& message)
	: runtime_error("parse error at line " + to_string(line_nr) + ": " + message)
{
}

// --------------------------------------------------------------------

const char* sac_parser::kTokenName[] = {
	"unknown",
	"EOF",
	"DATA",
	"LOOP",
	"GLOBAL",
	"SAVE",
	"STOP",
	"Tag",
	"Value"
};

const char* sac_parser::kValueName[] = {
	"Int",
	"Float",
	"Numeric",
	"String",
	"TextField",
	"Inapplicable",
	"Unknown"
};

// --------------------------------------------------------------------

sac_parser::sac_parser(std::istream& is)
	: m_data(is)
{
	m_validate = true;
	m_line_nr = 1;
	m_bol = true;
	m_lookahead = get_next_token();
}

void sac_parser::error(const string& msg)
{
	throw cif_parser_error(m_line_nr, msg);
}

// get_next_char takes a char from the buffer, or if it is empty
// from the istream. This function also does carriage/linefeed
// translation.
int sac_parser::get_next_char()
{
	int result;

	if (m_buffer.empty())
		result = m_data.get();
	else
	{
		result = m_buffer.top();
		m_buffer.pop();
	}
	
	// very simple CR/LF translation into LF
	if (result == '\r')
	{
		int lookahead = m_data.get();
		if (lookahead != '\n')
			m_buffer.push(lookahead);
		result = '\n';
	}
	
	m_token_value += static_cast<char>(result);
	
	if (result == '\n')
		++m_line_nr;
	
	if (VERBOSE >= 6)
	{
		cerr << "get_next_char => ";
		if (iscntrl(result) or not isprint(result))
			cerr << int(result) << endl;
		else
			cerr << char(result) << endl;
	}
	
	return result;
}

void sac_parser::retract()
{
	assert(not m_token_value.empty());

	char ch = m_token_value.back();
	if (ch == '\n')
		--m_line_nr;
	
	m_buffer.push(ch);
	m_token_value.pop_back();
}

void sac_parser::restart()
{
	while (not m_token_value.empty())
		retract();
	
	switch (m_start)
	{
		case eStateStart:
			m_state = m_start = eStateFloat;
			break;
		
		case eStateFloat:
			m_state = m_start = eStateInt;
			break;
		
		case eStateInt:
			m_state = m_start = eStateValue;
			break;
		
		default:
			error("Invalid state in sac_parser");
	}
	
	m_bol = false;
}

void sac_parser::match(sac_parser::CIFToken t)
{
	if (m_lookahead != t)
		error(string("Unexpected token, expected ") + kTokenName[t] + " but found " + kTokenName[m_lookahead]);
	
	m_lookahead = get_next_token();
}

sac_parser::CIFToken sac_parser::get_next_token()
{
	const auto kEOF = char_traits<char>::eof();
	
	CIFToken result = eCIFTokenUnknown;
	int quoteChar = 0;
	m_state = m_start = eStateStart;
	m_bol = false;
	
	m_token_value.clear();
	m_token_type = eCIFValueUnknown;
	
	while (result == eCIFTokenUnknown)
	{
		auto ch = get_next_char();
		
		switch (m_state)
		{
			case eStateStart:
				if (ch == kEOF)
					result = eCIFTokenEOF;
				else if (ch == '\n')
				{
					m_bol = true;
					m_state = eStateWhite;
				}
				else if (ch == ' ' or ch == '\t')
					m_state = eStateWhite;
				else if (ch == '#')
					m_state = eStateComment;
				else if (ch == '.')
					m_state = eStateDot;
				else if (ch == '_')
					m_state = eStateTag;
				else if (ch == ';' and m_bol)
					m_state = eStateTextField;
				else if (ch == '\'' or ch == '"')
				{
					quoteChar = ch;
					m_state = eStateQuotedString;
				}
				else if (ch == '?')
					m_state = eStateQuestionMark;
				else
					restart();
				break;
			
			case eStateWhite:
				if (ch == kEOF)
					result = eCIFTokenEOF;
				else if (not isspace(ch))
				{
					m_state = eStateStart;
					retract();
					m_token_value.clear();
				}
				else
					m_bol = (ch == '\n');
				break;
			
			case eStateComment:
				if (ch == '\n')
				{
					m_state = eStateStart;
					m_bol = true;
					m_token_value.clear();
				}
				else if (ch == kEOF)
					result = eCIFTokenEOF;
				else if (not is_any_print(ch))
					error("invalid character in comment");
				break;
			
			case eStateQuestionMark:
				if (is_non_blank(ch))
					m_state = eStateValue;
				else
				{
					retract();
					result = eCIFTokenValue;
					m_token_value.clear();
					m_token_type = eCIFValueUnknown;
				}
				break;

			case eStateDot:
				if (isdigit(ch))
					m_state = eStateFloat + 2;
				else if (isspace(ch))
				{
					retract();
					result = eCIFTokenValue;
					m_token_type = eCIFValueInapplicable;
				}
				else
					m_state = eStateValue;
				break;

			case eStateTextField:
				if (ch == '\n')
					m_state = eStateTextField + 1;
				else if (ch == kEOF)
					error("unterminated textfield");
				else if (not is_any_print(ch))
//					error("invalid character in text field '" + string({ static_cast<char>(ch) }) + "' (" + to_string((int)ch) + ")");
					cerr << "invalid character in text field '" << string({ static_cast<char>(ch) }) << "' (" << ch << ") line: " << m_line_nr << endl;
				break;
			
			case eStateTextField + 1:
				if (is_text_lead(ch) or ch == ' ' or ch == '\t')
					m_state = eStateTextField;
				else if (ch == ';')
				{
					assert(m_token_value.length() >= 2);
					m_token_value = m_token_value.substr(1, m_token_value.length() - 3);
					m_token_type = eCIFValueTextField;
					result = eCIFTokenValue;
				}
				else if (ch == kEOF)
					error("unterminated textfield");
				else if (ch != '\n')
					error("invalid character in text field");
				break;
			
			case eStateQuotedString:
				if (ch == kEOF)
					error("unterminated quoted string");
				else if (ch == quoteChar)
					m_state = eStateQuotedStringQuote;
				else if (not is_any_print(ch))
					error("invalid character in quoted string");
				break;
			
			case eStateQuotedStringQuote:
				if (is_white(ch))
				{
					retract();
					result = eCIFTokenValue;
					m_token_type = eCIFValueString;
					
					assert(m_token_value.length() >= 3);
					m_token_value = m_token_value.substr(1, m_token_value.length() - 2);
				}
				else if (ch == quoteChar)
					;
				else if (is_any_print(ch))
					m_state = eStateQuotedString;
				else if (ch == kEOF)
					error("unterminated quoted string");
				else
					error("invalid character in quoted string");
				break;
			
			case eStateTag:
				if (not is_non_blank(ch))
				{
					retract();
					result = eCIFTokenTag;
				}
				break;
			
			case eStateFloat:
				if (ch == '+' or ch == '-')
				{
					m_state = eStateFloat + 1;
				}
				else if (isdigit(ch))
					m_state = eStateFloat + 1;
				else
					restart();
				break;
			
			case eStateFloat + 1:
//				if (ch == '(')	// numeric???
//					m_state = eStateNumericSuffix;
//				else
				if (ch == '.')
					m_state = eStateFloat + 2;
				else if (tolower(ch) == 'e')
					m_state = eStateFloat + 3;
				else if (is_white(ch) or ch == kEOF)
				{
					retract();
					result = eCIFTokenValue;
					m_token_type = eCIFValueInt;
				}
				else
					restart();
				break;
			
			// parsed '.'
			case eStateFloat + 2:
//				if (ch == '(')	// numeric???
//					m_state = eStateNumericSuffix;
//				else
				if (tolower(ch) == 'e')
					m_state = eStateFloat + 3;
				else if (is_white(ch) or ch == kEOF)
				{
					retract();
					result = eCIFTokenValue;
					m_token_type = eCIFValueFloat;
				}
				else
					restart();
				break;
			
			// parsed 'e'
			case eStateFloat + 3:
				if (ch == '-' or ch == '+')
					m_state = eStateFloat + 4;
				else if (isdigit(ch))
					m_state = eStateFloat + 5;
				else
					restart();
				break;

			case eStateFloat + 4:
				if (isdigit(ch))
					m_state = eStateFloat + 5;
				else
					restart();
				break;
			
			case eStateFloat + 5:
//				if (ch == '(')
//					m_state = eStateNumericSuffix;
//				else
				if (is_white(ch) or ch == kEOF)
				{
					retract();
					result = eCIFTokenValue;
					m_token_type = eCIFValueFloat;
				}
				else
					restart();
				break;
			
			case eStateInt:
				if (isdigit(ch) or ch == '+' or ch == '-')
					m_state = eStateInt + 1;
				else
					restart();
				break;
			
			case eStateInt + 1:
				if (is_white(ch) or ch == kEOF)
				{
					retract();
					result = eCIFTokenValue;
					m_token_type = eCIFValueInt;
				}
				else
					restart();
				break;
			
//			case eStateNumericSuffix:
//				if (isdigit(ch))
//					m_state = eStateNumericSuffix + 1;
//				else
//					restart();
//				break;
//			
//			case eStateNumericSuffix + 1:
//				if (ch == ')')
//				{
//					result = eCIFTokenValue;
//					m_token_type = eCIFValueNumeric;
//				}
//				else if (not isdigit(ch))
//					restart();
//				break;
			
			case eStateValue:
				if (is_non_blank(ch))
					m_state = eStateValue + 1;
				else
					error("invalid character at this position");
				break;
			
			case eStateValue + 1:
				if (ch == '_')		// first _, check for keywords
				{
					string s = to_lower_copy(m_token_value);
					
					if (s == "global_")
						result = eCIFTokenGLOBAL;
					else if (s == "stop_")
						result = eCIFTokenSTOP;
					else if (s == "loop_")
						result = eCIFTokenLOOP;
					else if (s == "data_" or s == "save_")
						m_state = eStateValue + 2;
				}
				else if (not is_non_blank(ch))
				{
					retract();
					result = eCIFTokenValue;
					m_token_type = eCIFValueString;
				}
				break;

			case eStateValue + 2:
				if (not is_non_blank(ch))
				{
					retract();
					
					if (tolower(m_token_value[0]) == 'd')
						result = eCIFTokenDATA;
					else
						result = eCIFTokenSAVE;
					
					m_token_value.erase(m_token_value.begin(), m_token_value.begin() + 5); 
				}
				break;
			
			default:
				assert(false);
				error("Invalid state in get_next_token");
				break;
		}
	}

	if (VERBOSE >= 5)
	{
		cerr << kTokenName[result];
		if (m_token_type != eCIFValueUnknown)
			cerr << ' ' << kValueName[m_token_type];
		if (result != eCIFTokenEOF)
			cerr << " '" << m_token_value << '\'';
		cerr << endl;
	}
	
	return result;
}

void sac_parser::parse_file()
{
	try
	{
		while (m_lookahead != eCIFTokenEOF)
		{
			switch (m_lookahead)
			{
				case eCIFTokenGLOBAL:
					parse_global();
					break;
				
				case eCIFTokenDATA:
					produce_datablock(m_token_value);

					match(eCIFTokenDATA);
					parse_data_block();
					break;
				
				default:
					error("This file does not seem to be an mmCIF file");
					break;
			}
		}
	}
	catch (const exception& ex)
	{
		error(string("Error parsing file: '") + ex.what() + "'");
	}
}

void sac_parser::parse_global()
{
	match(eCIFTokenGLOBAL);
	while (m_lookahead == eCIFTokenTag)
	{
		match(eCIFTokenTag);
		match(eCIFTokenValue);
	}
}

void sac_parser::parse_data_block()
{
	string cat;
	
	while (m_lookahead == eCIFTokenLOOP or m_lookahead == eCIFTokenTag or m_lookahead == eCIFTokenSAVE)
	{
		switch (m_lookahead)
		{
			case eCIFTokenLOOP:
			{
				cat.clear();	// should start a new category
				
				match(eCIFTokenLOOP);
				
				vector<string> tags;
				
				while (m_lookahead == eCIFTokenTag)
				{
					string cat_name, item_name;
					std::tie(cat_name, item_name) = split_tag_name(m_token_value);
					
					if (cat.empty())
					{
						produce_category(cat_name);
						cat = cat_name;
					}
					else if (not iequals(cat, cat_name))
						error("inconsistent categories in loop_");
					
					tags.push_back(item_name);

					match(eCIFTokenTag);
				}
				
				while (m_lookahead == eCIFTokenValue)
				{
					produce_row();
					
					for (auto tag: tags)
					{
						produce_item(cat, tag, m_token_value);
						match(eCIFTokenValue);
					}
				}
				
				cat.clear();
				break;
			}
		
			case eCIFTokenTag:
			{
				string cat_name, item_name;
				std::tie(cat_name, item_name) = split_tag_name(m_token_value);

				if (not iequals(cat, cat_name))
				{
					produce_category(cat_name);
					cat = cat_name;
					produce_row();
				}

				match(eCIFTokenTag);
				
				produce_item(cat, item_name, m_token_value);

				match(eCIFTokenValue);
				break;
			}
			
			case eCIFTokenSAVE:
				parse_save_frame();
				break;
			
			default:
				assert(false);
				break;
		}
	}
}

void sac_parser::parse_save_frame()
{
	error("A regular CIF file should not contain a save frame");
}

// --------------------------------------------------------------------

parser::parser(std::istream& is, file& f)
	: sac_parser(is), m_file(f), m_data_block(nullptr)
{
}

void parser::produce_datablock(const string& name)
{
	m_data_block = new datablock(name);
	m_file.append(m_data_block);
}

void parser::produce_category(const string& name)
{
	if (VERBOSE >= 4)
		cerr << "producing category " << name << endl;

	std::tie(m_cat, ignore) = m_data_block->emplace(name);
}

void parser::produce_row()
{
	if (VERBOSE >= 4)
		cerr << "producing row for category " << m_cat->name() << endl;

	m_cat->emplace({});
	m_row = m_cat->back();
}

void parser::produce_item(const string& category, const string& item, const string& value)
{
	if (VERBOSE >= 4)
		cerr << "producing _" << category << '.' << item << " -> " << value << endl;

	if (not iequals(category, m_cat->name()))
		error("inconsistent categories in loop_");

	m_row[item] = m_token_value;
}

// --------------------------------------------------------------------

struct dict_parser_data_impl
{
	// temporary values for constructing dictionaries
	vector<validate_category>			m_category_validators;
	map<string,vector<validate_item>>	m_item_validators;
};

dict_parser::dict_parser(validator& validator, std::istream& is)
	: parser(is, m_file), m_validator(validator), m_impl(new dict_parser_data_impl)
{
}

dict_parser::~dict_parser()
{
	delete m_impl;
}

void dict_parser::parse_save_frame()
{
	if (not m_collected_item_types)
		m_collected_item_types = collect_item_types();

	string saveFrameName = m_token_value;

	if (saveFrameName.empty())
		error("Invalid save frame, should contain more than just 'save_' here");
	
	bool isCategorySaveFrame = m_token_value[0] != '_';
	
	datablock dict(m_token_value);
	datablock::iterator cat = dict.end();

	match(eCIFTokenSAVE);
	while (m_lookahead == eCIFTokenLOOP or m_lookahead == eCIFTokenTag)
	{
		if (m_lookahead == eCIFTokenLOOP)
		{
			cat = dict.end();	// should start a new category
				
			match(eCIFTokenLOOP);
			
			vector<string> tags;
			while (m_lookahead == eCIFTokenTag)
			{
				string cat_name, item_name;
				std::tie(cat_name, item_name) = split_tag_name(m_token_value);
					
				if (cat == dict.end())
					std::tie(cat, ignore) = dict.emplace(cat_name);
				else if (not iequals(cat->name(), cat_name))
					error("inconsistent categories in loop_");
				
				tags.push_back(item_name);
				match(eCIFTokenTag);
			}
			
			while (m_lookahead == eCIFTokenValue)
			{
				cat->emplace({});
				auto row = cat->back();
				
				for (auto tag: tags)
				{
					row[tag] = m_token_value;
					match(eCIFTokenValue);
				}
			}
			
			cat = dict.end();
		}
		else
		{
			string cat_name, item_name;
			std::tie(cat_name, item_name) = split_tag_name(m_token_value);

			if (cat == dict.end() or not iequals(cat->name(), cat_name))
				std::tie(cat, ignore) = dict.emplace(cat_name);

			match(eCIFTokenTag);
			
			if (cat->empty())
				cat->emplace({});
			cat->back()[item_name] = m_token_value;
			
			match(eCIFTokenValue);
		}
	}

	match(eCIFTokenSAVE);
	
	if (isCategorySaveFrame)
	{
		string category = dict.first_item("_category.id");

		vector<string> keys;
		for (auto k: dict["category_key"])
			keys.push_back(get<1>(split_tag_name(k["name"].as<string>())));
		
		iset groups;
		for (auto g: dict["category_group"])
			groups.insert(g["id"].as<string>());
			
		m_impl->m_category_validators.push_back(validate_category{category, keys, groups});
	}
	else
	{
		// if the type code is missing, this must be a pointer, just skip it
		string type_code = dict.first_item("_item_type.code");

		const validate_type* tv = nullptr;
		if (not (type_code.empty() or type_code == "?"))
			tv = m_validator.get_validator_for_type(type_code);

		iset ess;
		for (auto e: dict["item_enumeration"])
			ess.insert(e["value"].as<string>());
		
		// collect the dict from our data_block and construct validators
		for (auto i: dict["item"])
		{
			string tag_name, category, mandatory;
			
			cif::tie(tag_name, category, mandatory) = i.get("name", "category_id", "mandatory_code");
			
			string cat_name, item_name;
			std::tie(cat_name, item_name) = split_tag_name(tag_name);
			
			if (cat_name.empty() or item_name.empty())
				error("Invalid tag name in _item.name " + tag_name);

			if (not iequals(category, cat_name) and not (category.empty() or category == "?"))
				error("specified category id does match the implicit category name for tag '" + tag_name + '\'');
			else
				category = cat_name;
			
			auto& ivs = m_impl->m_item_validators[category];
			
			auto vi = find(ivs.begin(), ivs.end(), validate_item{item_name});
			if (vi == ivs.end())
				ivs.push_back(validate_item{item_name, iequals(mandatory, "yes"), tv, ess});
			else
			{
				// need to update the item_validator?
				if (vi->m_mandatory != (iequals(mandatory, "yes")))
				{
					if (VERBOSE > 2)
					{
						cerr << "inconsistent mandatory value for " << tag_name << " in dictionary" << endl;
						
						if (iequals(tag_name, saveFrameName))
							cerr << "choosing " << mandatory << endl;
						else
							cerr << "choosing " << (vi->m_mandatory ? "Y" : "N") << endl;
					}

					if (iequals(tag_name, saveFrameName))
						vi->m_mandatory = (iequals(mandatory, "yes"));
				}

				if (vi->m_type != nullptr and tv != nullptr and vi->m_type != tv)
				{
					if (VERBOSE > 1)
						cerr << "inconsistent type for " << tag_name << " in dictionary" << endl;
				}

//				vi->m_mandatory = (iequals(mandatory, "yes"));
				if (vi->m_type == nullptr)
					vi->m_type = tv;

				vi->m_enums.insert(ess.begin(), ess.end());

				// anything else yet?
				// ...
			}
		}
	}
}

void dict_parser::link_items()
{
	if (not m_data_block)
		error("no datablock");
	
	auto& dict = *m_data_block;
	
	for (auto gl: dict["pdbx_item_linked_group_list"])
	{
		string child, parent;
		cif::tie(child, parent) = gl.get("child_name", "parent_name");
		
		auto civ = m_validator.get_validator_for_item(child);
		if (civ == nullptr)
			error("in pdbx_item_linked_group_list, item '" + child + "' is not specified");
		
		auto piv = m_validator.get_validator_for_item(parent);
		if (piv == nullptr)
			error("in pdbx_item_linked_group_list, item '" + parent + "' is not specified");
		
		civ->set_parent(piv);
	}
	
	// now make sure the item_type is specified for all item_validators
	
	for (auto& cv: m_validator.m_category_validators)
	{
		for (auto& iv: cv.m_item_validators)
		{
			if (iv.m_type == nullptr)
				cerr << "Missing item_type for " << iv.m_tag << endl;
		}
	}	
}

void dict_parser::load_dictionary()
{
	unique_ptr<datablock> dict;
	datablock* saved_datablock = m_data_block;
	
	try
	{
		while (m_lookahead != eCIFTokenEOF)
		{
			switch (m_lookahead)
			{
				case eCIFTokenGLOBAL:
					parse_global();
					break;
				
				default:
				{
					dict.reset(new datablock(m_token_value));	// dummy datablock, for constructing the validator only
					m_data_block = dict.get();
					
					match(eCIFTokenDATA);
					parse_data_block();
					break;
				}
			}
		}
	}
	catch (const exception& ex)
	{
		if (VERBOSE)
			cerr << "Error parsing dictionary: '" << ex.what() << "'" << endl;
	}

	// store all validators
	for (auto& ic: m_impl->m_category_validators)
		m_validator.add_category_validator(move(ic));
	m_impl->m_category_validators.clear();
	
	for (auto& iv: m_impl->m_item_validators)
	{
		auto cv = m_validator.get_validator_for_category(iv.first);
		if (cv == nullptr)
			error("Undefined category '" + iv.first);

		for (auto& v: iv.second)
			const_cast<validate_category*>(cv)->add_item_validator(move(v));
	}
		
	// check all item validators for having a type_validator
	
	if (dict)
		link_items();

	// store meta information
	datablock::iterator info;
	bool n;
	std::tie(info, n) = m_data_block->emplace("dictionary");
	if (n)
	{
		auto r = info->front();
		m_validator.dict_name(r["title"].as<string>());
		m_validator.dict_version(r["version"].as<string>());
	}

	m_data_block = saved_datablock;

	m_impl->m_item_validators.clear();
}

bool dict_parser::collect_item_types()
{
	bool result = false;
	
	if (not m_data_block)
		error("no datablock");
	
	auto& dict = *m_data_block;
	
	for (auto& t: dict["item_type_list"])
	{
		auto ts = t.get("code", "primitive_code", "construct");

		string code, primitive_code, construct;
		cif::tie(code, primitive_code, construct) = ts;
		
		ba::replace_all(construct, "\\n", "\n");
		ba::replace_all(construct, "\\t", "\t");
		ba::replace_all(construct, "\\\n", "");
		
		validate_type v = {
			code, map_to_primitive_type(primitive_code), boost::regex(construct, boost::regex::egrep)
		};

// Do not replace an already defined type validator, this won't work with pdbx_v40
// as it has a name that is too strict for its own names :-)
//		if (m_file_impl.m_type_validators.count(v))
//			m_file_impl.m_type_validators.erase(v);
		
		m_validator.add_type_validator(move(v));

		if (VERBOSE >= 5)
			cerr << "Added type " << code << " (" << primitive_code << ") => " << construct << endl;
		
		result = true;
	}
	
	return result;
}


}
