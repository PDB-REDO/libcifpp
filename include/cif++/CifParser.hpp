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

#include "cif++/Cif++.hpp"

#include <stack>
#include <map>

namespace cif
{

// --------------------------------------------------------------------

class CifParserError : public std::runtime_error
{
  public:
	CifParserError(uint32_t lineNr, const std::string& message);
};

// --------------------------------------------------------------------

extern const uint32_t kMaxLineLength;

extern const uint8_t kCharTraitsTable[128];

enum CharTraitsMask: uint8_t {
	kOrdinaryMask = 1 << 0,
	kNonBlankMask = 1 << 1,
	kTextLeadMask = 1 << 2,
	kAnyPrintMask = 1 << 3
};

inline bool isWhite(int ch)
{
	return std::isspace(ch) or ch == '#';
}

inline bool isOrdinary(int ch)
{
	return ch >= 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kOrdinaryMask) != 0;
}

inline bool isNonBlank(int ch)
{
	return ch > 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kNonBlankMask) != 0;
}

inline bool isTextLead(int ch)
{
	return ch >= 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kTextLeadMask) != 0;
}

inline bool isAnyPrint(int ch)	
{
	return ch == '\t' or 
		(ch >= 0x20 and ch <= 0x7f and (kCharTraitsTable[ch - 0x20] & kAnyPrintMask) != 0);
}

inline bool isUnquotedString(const char* s)
{
	bool result = isOrdinary(*s++);
	while (result and *s != 0)
	{
		result = isNonBlank(*s);
		++s;
	}
	return result;
}

// --------------------------------------------------------------------

using DatablockIndex = std::map<std::string,std::size_t>;

// --------------------------------------------------------------------
// sac Parser, analogous to SAX Parser (simple api for xml)

class SacParser
{
  public:
	SacParser(std::istream& is, bool init = true);
	virtual ~SacParser() {}

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
	
	int getNextChar();

	void retract();
	void restart();
	
	CIFToken getNextToken();
	void match(CIFToken token);

	bool parseSingleDatablock(const std::string& datablock);

	DatablockIndex indexDatablocks();
	bool parseSingleDatablock(const std::string& datablock, const DatablockIndex &index);

	void parseFile();
	void parseGlobal();
	void parseDataBlock();

	virtual void parseSaveFrame();
	
	void parseDictionary();
	
	void error(const std::string& msg);
	
	// production methods, these are pure virtual here
	
	virtual void produceDatablock(const std::string& name) = 0;
	virtual void produceCategory(const std::string& name) = 0;
	virtual void produceRow() = 0;
	virtual void produceItem(const std::string& category, const std::string& item, const std::string& value) = 0;

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

	std::istream&			mData;

	// Parser state
	bool					mValidate;
	uint32_t					mLineNr;
	bool					mBol;
	int						mState, mStart;
	CIFToken				mLookahead;
	std::string				mTokenValue;
	CIFValueType			mTokenType;
	std::stack<int>			mBuffer;
};

// --------------------------------------------------------------------

class Parser : public SacParser
{
  public:
	Parser(std::istream& is, File& f, bool init = true);

	virtual void produceDatablock(const std::string& name);
	virtual void produceCategory(const std::string& name);
	virtual void produceRow();
	virtual void produceItem(const std::string& category, const std::string& item, const std::string& value);

  protected:
	File&					mFile;
	Datablock*				mDataBlock;
	Datablock::iterator		mCat;
	Row						mRow;
};

// --------------------------------------------------------------------

class DictParser : public Parser
{
  public:

	DictParser(Validator& validator, std::istream& is);
	~DictParser();
	
	void loadDictionary();
	
  private:

	virtual void parseSaveFrame();
	
	bool collectItemTypes();
	void linkItems();

	Validator&						mValidator;
	File							mFile;
	struct DictParserDataImpl*		mImpl;
	bool							mCollectedItemTypes = false;
};

}
