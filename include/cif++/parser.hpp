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

/**
 * @file parser.hpp
 * 
 * This file contains the declaration of an mmCIF parser
 */

namespace cif
{

// --------------------------------------------------------------------

/** Exception that is thrown when the mmCIF file contains a parsing error */
class parse_error : public std::runtime_error
{
  public:
	/// \brief constructor
	parse_error(uint32_t line_nr, const std::string &message)
		: std::runtime_error("parse error at line " + std::to_string(line_nr) + ": " + message)
	{
	}
};

// --------------------------------------------------------------------

/**
 * @brief The sac_parser is a similar to SAX parsers (Simple API for XML, 
 * in our case it is Simple API for CIF)
 * 
 * This is a hand crafted, optimised parser for reading cif files,
 * both cif 1.0 and cif 1.1 is supported. But version 2.0 is not.
 * That means that the content of files strictly contains only
 * ASCII characters. Anything else will generate an error.
 * 
 * This class is an abstract base class. Derived classes should
 * implement the produce_ methods.
 */

// TODO: Need to implement support for transformed long lines

class sac_parser
{
  public:
	/** @cond */
	using datablock_index = std::map<std::string, std::size_t>;

	virtual ~sac_parser() = default;
	/** @endcond */

	/// \brief The parser only supports ASCII so we can
	/// create a table with character properties.
	enum CharTraitsMask : uint8_t
	{
		kOrdinaryMask = 1 << 0,	///< The character is in the Ordinary class
		kNonBlankMask = 1 << 1,	///< The character is in the NonBlank class
		kTextLeadMask = 1 << 2,	///< The character is in the TextLead class
		kAnyPrintMask = 1 << 3	///< The character is in the AnyPrint class
	};

	/// \brief Return true if the character @a ch is a *space* character
	static constexpr bool is_space(int ch)
	{
		return ch == ' ' or ch == '\t' or ch == '\r' or ch == '\n';
	}

	/// \brief Return true if the character @a ch is a *white* character
	static constexpr bool is_white(int ch)
	{
		return is_space(ch) or ch == '#';
	}

	/// \brief Return true if the character @a ch is a *ordinary* character
	static constexpr bool is_ordinary(int ch)
	{
		return ch >= 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kOrdinaryMask) != 0;
	}

	/// \brief Return true if the character @a ch is a *non_blank* character
	static constexpr bool is_non_blank(int ch)
	{
		return ch > 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kNonBlankMask) != 0;
	}

	/// \brief Return true if the character @a ch is a *text_lead* character
	static constexpr bool is_text_lead(int ch)
	{
		return ch >= 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kTextLeadMask) != 0;
	}

	/// \brief Return true if the character @a ch is a *any_print* character
	static constexpr bool is_any_print(int ch)
	{
		return ch == '\t' or
		       (ch >= 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kAnyPrintMask) != 0);
	}

	/// \brief Return true if the string in @a text can safely be written without quotation
	static bool is_unquoted_string(std::string_view text);

  protected:
	/** @cond */

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
		UNKNOWN,

		END_OF_FILE,

		DATA,
		LOOP,
		GLOBAL,
		SAVE_,
		SAVE_NAME,
		STOP,
		ITEM_NAME,
		VALUE
	};

	static constexpr const char *get_token_name(CIFToken token)
	{
		switch (token)
		{
			case CIFToken::UNKNOWN: return "Unknown";
			case CIFToken::END_OF_FILE: return "Eof";
			case CIFToken::DATA: return "DATA";
			case CIFToken::LOOP: return "LOOP";
			case CIFToken::GLOBAL: return "GLOBAL";
			case CIFToken::SAVE_: return "SAVE";
			case CIFToken::SAVE_NAME: return "SAVE+name";
			case CIFToken::STOP: return "STOP";
			case CIFToken::ITEM_NAME: return "Tag";
			case CIFToken::VALUE: return "Value";
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

	/** @endcond */

  public:

	/** \brief Parse only a single datablock in the string @a datablock
	 * The start of the datablock is first located and then data
	 * is parsed up until the next start of a datablock or the end of
	 * the data.
	 * */
	bool parse_single_datablock(const std::string &datablock);

	/** \brief Return an index for all the datablocks found, that is
	 * the index will contain the names and offsets for each.
	 */
	datablock_index index_datablocks();

	/**
	 * @brief Parse the datablock named @a datablock
	 * 
	 * This will first lookup the datablock's offset in the index @a index
	 * and then start parsing from that location until the next datablock.
	 * 
	 * @param datablock Name of the datablock to parse
	 * @param index The index created using index_datablocks
	 * @return true If the datablock was found
	 * @return false If the datablock was not found
	 */
	bool parse_single_datablock(const std::string &datablock, const datablock_index &index);

	/**
	 * @brief Parse the file
	 * 
	 */
	void parse_file();

  protected:

	/** @cond */

	sac_parser(std::istream &is, bool init = true);

	void parse_global();

	void parse_datablock();

	virtual void parse_save_frame();

	void error(const std::string &msg)
	{
		if (cif::VERBOSE > 0)
			std::cerr << "Error parsing mmCIF: " << msg << '\n';

		throw parse_error(m_line_nr, msg);
	}

	void warning(const std::string &msg)
	{
		if (cif::VERBOSE > 0)
			std::cerr << "parser warning at line " << m_line_nr << ": " << msg << '\n';
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
		ItemName,
		TextItem,
		TextItemNL,
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

	/** @endcond */
};

// --------------------------------------------------------------------

/**
 * @brief An actual implementation of a sac_parser generating data in a file
 * 
 * This parser will create the cif::file, cif::datablock and cif::category
 * objects required to contain all data
 */
class parser : public sac_parser
{
  public:
	/// \brief constructor, generates data into @a file from @a is
	parser(std::istream &is, file &file)
		: sac_parser(is)
		, m_file(file)
	{
	}

	/** @cond */
	void produce_datablock(std::string_view name) override;

	void produce_category(std::string_view name) override;

	void produce_row() override;

	void produce_item(std::string_view category, std::string_view item, std::string_view value) override;

  protected:
	file &m_file;
	datablock *m_datablock = nullptr;
	category *m_category = nullptr;
	row_handle m_row;

	/** @endcond */
};

} // namespace cif
