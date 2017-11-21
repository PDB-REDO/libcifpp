// CIF parser

#include "libcif/cif++.h"

#include <stack>

namespace cif
{

// --------------------------------------------------------------------

class cif_parser_error : public std::runtime_error
{
  public:
	cif_parser_error(uint32 line_nr, const std::string& message);
};

// --------------------------------------------------------------------

extern const uint32 kMaxLineLength;

extern const uint8 kCharTraitsTable[128];

enum CharTraitsMask: uint8 {
	kOrdinaryMask = 1 << 0,
	kNonBlankMask = 1 << 1,
	kTextLeadMask = 1 << 2,
	kAnyPrintMask = 1 << 3
};

inline bool is_white(int ch)
{
	return std::isspace(ch) or ch == '#';
}

inline bool is_ordinary(int ch)
{
	return ch >= 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kOrdinaryMask) != 0;
}

inline bool is_non_blank(int ch)
{
	return ch > 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kNonBlankMask) != 0;
}

inline bool is_text_lead(int ch)
{
	return ch >= 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kTextLeadMask) != 0;
}

inline bool is_any_print(int ch)	
{
	return ch == '\t' or 
		(ch >= 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kAnyPrintMask) != 0);
}

inline bool is_unquoted_string(const char* s)
{
	bool result = is_ordinary(*s++);
	while (result and *s != 0)
	{
		result = is_non_blank(*s);
		++s;
	}
	return result;
}

// --------------------------------------------------------------------

std::tuple<std::string,std::string> split_tag_name(const std::string& tag);

// --------------------------------------------------------------------
// sac parser, analogous to SAX parser (simple api for xml)

class sac_parser
{
  public:
	sac_parser(std::istream& is);
	virtual ~sac_parser() {}

	enum CIFToken
	{
		eCIFTokenUnknown,
		
		eCIFTokenEOF,
	
		eCIFTokenDATA,
		eCIFTokenLOOP,
		eCIFTokenGLOBAL,
		eCIFTokenSAVE,
		eCIFTokenSTOP,
		eCIFTokenTag,
		eCIFTokenValue,
	};

	static const char* kTokenName[];

	enum CIFValueType
	{
		eCIFValueInt,
		eCIFValueFloat,
		eCIFValueNumeric,
		eCIFValueString,
		eCIFValueTextField,
		eCIFValueInapplicable,
		eCIFValueUnknown
	};

	static const char* kValueName[];
	
	int get_next_char();

	void retract();
	void restart();
	
	CIFToken get_next_token();
	void match(CIFToken token);

	void parse_file();
	void parse_global();
	void parse_data_block();

	virtual void parse_save_frame();
	
	void parse_dictionary();
	
	void error(const std::string& msg);
	
	// production methods, these are pure virtual here
	
	virtual void produce_datablock(const std::string& name) = 0;
	virtual void produce_category(const std::string& name) = 0;
	virtual void produce_row() = 0;
	virtual void produce_item(const std::string& category, const std::string& item, const string& value) = 0;

  protected:

	enum State
	{
		eStateStart,
		eStateWhite,
		eStateComment,
		eStateQuestionMark,
		eStateDot,
		eStateQuotedString,
		eStateQuotedStringQuote,
		eStateUnquotedString,
		eStateTag,
		eStateTextField,
		eStateFloat = 100,
		eStateInt = 110,
//		eStateNumericSuffix = 200,
		eStateValue = 300
	};

	std::istream&			m_data;

	// parser state
	bool					m_validate;
	uint32					m_line_nr;
	bool					m_bol;
	int						m_state, m_start;
	CIFToken				m_lookahead;
	std::string				m_token_value;
	CIFValueType			m_token_type;
	std::stack<int>			m_buffer;
};

// --------------------------------------------------------------------

class parser : public sac_parser
{
  public:
	parser(std::istream& is, file& f);

	virtual void produce_datablock(const std::string& name);
	virtual void produce_category(const std::string& name);
	virtual void produce_row();
	virtual void produce_item(const std::string& category, const std::string& item, const std::string& value);

  protected:
	file&					m_file;
	datablock*				m_data_block;
	datablock::iterator		m_cat;
	row						m_row;
};

// --------------------------------------------------------------------

class dict_parser : public parser
{
  public:

	dict_parser(validator& validator, std::istream& is);
	~dict_parser();
	
	void load_dictionary();
	
  private:

	virtual void parse_save_frame();
	
	bool collect_item_types();
	void link_items();

	validator&						m_validator;
	file							m_file;
	struct dict_parser_data_impl*	m_impl;
	bool							m_collected_item_types = false;
};

}
