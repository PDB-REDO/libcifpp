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

#include "cif++/row.hpp"

#include <map>

namespace cif
{

// --------------------------------------------------------------------

class parse_error : public std::runtime_error
{
  public:
	parse_error(uint32_t line_nr, const std::string &message)
		: std::runtime_error("parse error at line " + std::to_string(line_nr) + ": " + message)
	{
	}
};

// --------------------------------------------------------------------

// TODO: Need to implement support for transformed long lines

class sac_parser
{
  public:
	using datablock_index = std::map<std::string, std::size_t>;

	virtual ~sac_parser() = default;

	enum CharTraitsMask : uint8_t
	{
		kOrdinaryMask = 1 << 0,
		kNonBlankMask = 1 << 1,
		kTextLeadMask = 1 << 2,
		kAnyPrintMask = 1 << 3
	};

	static constexpr bool is_space(int ch)
	{
		return ch == ' ' or ch == '\t' or ch == '\r' or ch == '\n';
	}

	static constexpr bool is_white(int ch)
	{
		return is_space(ch) or ch == '#';
	}

	static constexpr bool is_ordinary(int ch)
	{
		return ch >= 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kOrdinaryMask) != 0;
	}

	static constexpr bool is_non_blank(int ch)
	{
		return ch > 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kNonBlankMask) != 0;
	}

	static constexpr bool is_text_lead(int ch)
	{
		return ch >= 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kTextLeadMask) != 0;
	}

	static constexpr bool is_any_print(int ch)
	{
		return ch == '\t' or
		       (ch >= 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kAnyPrintMask) != 0);
	}

	static bool is_unquoted_string(std::string_view text);

  protected:
	static constexpr uint8_t kCharTraitsTable[128] = {
		//	0	1	2	3	4	5	6	7	8	9	a	b	c	d	e	f
		14, 15, 14, 14, 14, 15, 15, 14, 15, 15, 15, 15, 15, 15, 15, 15, //	2
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 10, 15, 15, 15, 15, //	3
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, //	4
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 14, 15, 14, 15, 14, //	5
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, //	6
		15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 0,  //	7
	};

	enum class CIFToken
	{
		Unknown,

		Eof,

		DATA,
		LOOP,
		GLOBAL,
		SAVE_,
		SAVE_NAME,
		STOP,
		Tag,
		Value
	};

	static constexpr const char *get_token_name(CIFToken token)
	{
		switch (token)
		{
			case CIFToken::Unknown: return "Unknown";
			case CIFToken::Eof: return "Eof";
			case CIFToken::DATA: return "DATA";
			case CIFToken::LOOP: return "LOOP";
			case CIFToken::GLOBAL: return "GLOBAL";
			case CIFToken::SAVE_: return "SAVE";
			case CIFToken::SAVE_NAME: return "SAVE+name";
			case CIFToken::STOP: return "STOP";
			case CIFToken::Tag: return "Tag";
			case CIFToken::Value: return "Value";
			default: return "Invalid token parameter";
		}
	}

	// get_next_char takes the next character from the istream.
	// This function also does carriage/linefeed translation.
	int get_next_char();

	// Put the last read character back into the istream
	void retract();

	CIFToken get_next_token();

	void match(CIFToken token);

  public:
	bool parse_single_datablock(const std::string &datablock);

	datablock_index index_datablocks();

	bool parse_single_datablock(const std::string &datablock, const datablock_index &index);

	void parse_file();

  protected:

	sac_parser(std::istream &is, bool init = true);

	void parse_global();

	void parse_datablock();

	virtual void parse_save_frame();

	void error(const std::string &msg)
	{
		if (cif::VERBOSE > 0)
			std::cerr << "Error parsing mmCIF: " << msg << std::endl;

		throw parse_error(m_line_nr, msg);
	}

	void warning(const std::string &msg)
	{
		if (cif::VERBOSE > 0)
			std::cerr << "parser warning at line " << m_line_nr << ": " << msg << std::endl;
	}

	// production methods, these are pure virtual here

	virtual void produce_datablock(std::string_view name) = 0;
	virtual void produce_category(std::string_view name) = 0;
	virtual void produce_row() = 0;
	virtual void produce_item(std::string_view category, std::string_view item, std::string_view value) = 0;

  protected:

	enum class State
	{
		Start,
		White,
		Esc,
		Comment,
		QuestionMark,
		Dot,
		QuotedString,
		QuotedStringQuote,
		UnquotedString,
		Tag,
		TextField,
		TextFieldNL,
		Reserved,
		Value
	};

	std::streambuf &m_source;

	// Parser state
	uint32_t m_line_nr;
	bool m_bol;
	CIFToken m_lookahead;

	// token buffer
	std::vector<char> m_token_buffer;
	std::string_view m_token_value;
};

// --------------------------------------------------------------------

class parser : public sac_parser
{
  public:
	parser(std::istream &is, file &file)
		: sac_parser(is)
		, m_file(file)
	{
	}

	void produce_datablock(std::string_view name) override;

	void produce_category(std::string_view name) override;

	void produce_row() override;

	void produce_item(std::string_view category, std::string_view item, std::string_view value) override;

  protected:
	file &m_file;
	datablock *m_datablock = nullptr;
	category *m_category = nullptr;
	row_handle m_row;
};

} // namespace cif
