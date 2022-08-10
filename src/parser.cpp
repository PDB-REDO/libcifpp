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

#include <set>

#include <cif++/v2/parser.hpp>

// extern int VERBOSE;

namespace cif::v2
{

const uint32_t kMaxLineLength = 132;

const uint8_t kCharTraitsTable[128] = {
	//	0	1	2	3	4	5	6	7	8	9	a	b	c	d	e	f
	14, 15, 14, 14, 14, 15, 15, 14, 15, 15, 15, 15, 15, 15, 15, 15, //	2
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 10, 15, 15, 15, 15, //	3
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, //	4
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 14, 15, 14, 15, 14, //	5
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, //	6
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 0,  //	7
};

// --------------------------------------------------------------------

parse_error::parse_error(uint32_t lineNr, const std::string &message)
	: std::runtime_error("parse error at line " + std::to_string(lineNr) + ": " + message)
{
}

// --------------------------------------------------------------------

const char *SacParser::kTokenName[] = {
	"unknown",
	"EOF",
	"DATA",
	"LOOP",
	"GLOBAL",
	"SAVE",
	"STOP",
	"Tag",
	"Value"};

const char *SacParser::kValueName[] = {
	"Int",
	"Float",
	"Numeric",
	"String",
	"TextField",
	"Inapplicable",
	"Unknown"};

// --------------------------------------------------------------------

bool isUnquotedString(const char *s)
{
	auto ss = s;

	bool result = isOrdinary(*s++);
	while (result and *s != 0)
	{
		result = isNonBlank(*s);
		++s;
	}

	// but be careful it does not contain e.g. stop_
	if (result)
	{
		static const std::regex reservedRx(R"((^(?:data|save)|.*(?:loop|stop|global))_.+)", std::regex_constants::icase);
		result = not std::regex_match(ss, reservedRx);
	}

	return result;
}

// --------------------------------------------------------------------

SacParser::SacParser(std::istream &is, bool init)
	: mData(is)
{
	mValidate = true;
	mLineNr = 1;
	mBol = true;

	if (init)
		mLookahead = getNextToken();
}

void SacParser::error(const std::string &msg)
{
	throw parse_error(mLineNr, msg);
}

// getNextChar takes a char from the buffer, or if it is empty
// from the istream. This function also does carriage/linefeed
// translation.
int SacParser::getNextChar()
{
	int result;

	if (mBuffer.empty())
		result = mData.get();
	else
	{
		result = mBuffer.top();
		mBuffer.pop();
	}

	// very simple CR/LF translation into LF
	if (result == '\r')
	{
		int lookahead = mData.get();
		if (lookahead != '\n')
			mBuffer.push(lookahead);
		result = '\n';
	}

	mTokenValue += static_cast<char>(result);

	if (result == '\n')
		++mLineNr;

	if (VERBOSE >= 6)
	{
		std::cerr << "getNextChar => ";
		if (iscntrl(result) or not isprint(result))
			std::cerr << int(result) << std::endl;
		else
			std::cerr << char(result) << std::endl;
	}

	return result;
}

void SacParser::retract()
{
	assert(not mTokenValue.empty());

	char ch = mTokenValue.back();
	if (ch == '\n')
		--mLineNr;

	mBuffer.push(ch);
	mTokenValue.pop_back();
}


int SacParser::restart(int start)
{
	int result = 0;

	while (not mTokenValue.empty())
		retract();

	switch (start)
	{
		case eStateStart:
			result = eStateFloat;
			break;

		case eStateFloat:
			result = eStateInt;
			break;

		case eStateInt:
			result = eStateValue;
			break;

		default:
			error("Invalid state in SacParser");
	}

	mBol = false;

	return result;
}

void SacParser::match(SacParser::CIFToken t)
{
	if (mLookahead != t)
		error(std::string("Unexpected token, expected ") + kTokenName[t] + " but found " + kTokenName[mLookahead]);

	mLookahead = getNextToken();
}

SacParser::CIFToken SacParser::getNextToken()
{
	const auto kEOF = std::char_traits<char>::eof();

	CIFToken result = eCIFTokenUnknown;
	int quoteChar = 0;
	int state = eStateStart, start = eStateStart;
	mBol = false;

	mTokenValue.clear();
	mTokenType = eCIFValueUnknown;

	while (result == eCIFTokenUnknown)
	{
		auto ch = getNextChar();

		switch (state)
		{
			case eStateStart:
				if (ch == kEOF)
					result = eCIFTokenEOF;
				else if (ch == '\n')
				{
					mBol = true;
					state = eStateWhite;
				}
				else if (ch == ' ' or ch == '\t')
					state = eStateWhite;
				else if (ch == '#')
					state = eStateComment;
				else if (ch == '_')
					state = eStateTag;
				else if (ch == ';' and mBol)
					state = eStateTextField;
				else if (ch == '\'' or ch == '"')
				{
					quoteChar = ch;
					state = eStateQuotedString;
				}
				else
					state = start = restart(start);
				break;

			case eStateWhite:
				if (ch == kEOF)
					result = eCIFTokenEOF;
				else if (not isspace(ch))
				{
					state = eStateStart;
					retract();
					mTokenValue.clear();
				}
				else
					mBol = (ch == '\n');
				break;

			case eStateComment:
				if (ch == '\n')
				{
					state = eStateStart;
					mBol = true;
					mTokenValue.clear();
				}
				else if (ch == kEOF)
					result = eCIFTokenEOF;
				else if (not isAnyPrint(ch))
					error("invalid character in comment");
				break;

			case eStateTextField:
				if (ch == '\n')
					state = eStateTextField + 1;
				else if (ch == kEOF)
					error("unterminated textfield");
				else if (not isAnyPrint(ch))
					//					error("invalid character in text field '" + string({ static_cast<char>(ch) }) + "' (" + to_string((int)ch) + ")");
					std::cerr << "invalid character in text field '" << std::string({static_cast<char>(ch)}) << "' (" << ch << ") line: " << mLineNr << std::endl;
				break;

			case eStateTextField + 1:
				if (isTextLead(ch) or ch == ' ' or ch == '\t')
					state = eStateTextField;
				else if (ch == ';')
				{
					assert(mTokenValue.length() >= 2);
					mTokenValue = mTokenValue.substr(1, mTokenValue.length() - 3);
					mTokenType = eCIFValueTextField;
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
					state = eStateQuotedStringQuote;
				else if (not isAnyPrint(ch))
					std::cerr << "invalid character in quoted string '" << std::string({static_cast<char>(ch)}) << "' (" << ch << ") line: " << mLineNr << std::endl;
					// error("invalid character in quoted string");
				break;

			case eStateQuotedStringQuote:
				if (isWhite(ch))
				{
					retract();
					result = eCIFTokenValue;
					mTokenType = eCIFValueString;

					if (mTokenValue.length() < 2)
						error("Invalid quoted string token");
						
					mTokenValue = mTokenValue.substr(1, mTokenValue.length() - 2);
				}
				else if (ch == quoteChar)
					;
				else if (isAnyPrint(ch))
					state = eStateQuotedString;
				else if (ch == kEOF)
					error("unterminated quoted string");
				else
					error("invalid character in quoted string");
				break;

			case eStateTag:
				if (not isNonBlank(ch))
				{
					retract();
					result = eCIFTokenTag;
				}
				break;

			case eStateFloat:
				if (ch == '+' or ch == '-')
				{
					state = eStateFloat + 1;
				}
				else if (isdigit(ch))
					state = eStateFloat + 1;
				else
					state = start = restart(start);
				break;

			case eStateFloat + 1:
				//				if (ch == '(')	// numeric???
				//					mState = eStateNumericSuffix;
				//				else
				if (ch == '.')
					state = eStateFloat + 2;
				else if (tolower(ch) == 'e')
					state = eStateFloat + 3;
				else if (isWhite(ch) or ch == kEOF)
				{
					retract();
					result = eCIFTokenValue;
					mTokenType = eCIFValueInt;
				}
				else
					state = start = restart(start);
				break;

			// parsed '.'
			case eStateFloat + 2:
				if (tolower(ch) == 'e')
					state = eStateFloat + 3;
				else if (isWhite(ch) or ch == kEOF)
				{
					retract();
					result = eCIFTokenValue;
					mTokenType = eCIFValueFloat;
				}
				else
					state = start = restart(start);
				break;

			// parsed 'e'
			case eStateFloat + 3:
				if (ch == '-' or ch == '+')
					state = eStateFloat + 4;
				else if (isdigit(ch))
					state = eStateFloat + 5;
				else
					state = start = restart(start);
				break;

			case eStateFloat + 4:
				if (isdigit(ch))
					state = eStateFloat + 5;
				else
					state = start = restart(start);
				break;

			case eStateFloat + 5:
				if (isWhite(ch) or ch == kEOF)
				{
					retract();
					result = eCIFTokenValue;
					mTokenType = eCIFValueFloat;
				}
				else
					state = start = restart(start);
				break;

			case eStateInt:
				if (isdigit(ch) or ch == '+' or ch == '-')
					state = eStateInt + 1;
				else
					state = start = restart(start);
				break;

			case eStateInt + 1:
				if (isWhite(ch) or ch == kEOF)
				{
					retract();
					result = eCIFTokenValue;
					mTokenType = eCIFValueInt;
				}
				else
					state = start = restart(start);
				break;

			case eStateValue:
				if (ch == '_')
				{
					std::string s = toLowerCopy(mTokenValue);

					if (s == "global_")
						result = eCIFTokenGLOBAL;
					else if (s == "stop_")
						result = eCIFTokenSTOP;
					else if (s == "loop_")
						result = eCIFTokenLOOP;
					else if (s == "data_")
					{
						state = eStateDATA;
						continue;
					}
					else if (s == "save_")
					{
						state = eStateSAVE;
						continue;
					}
				}

				if (result == eCIFTokenUnknown and not isNonBlank(ch))
				{
					retract();
					result = eCIFTokenValue;

					if (mTokenValue == ".")
						mTokenType = eCIFValueInapplicable;
					else if (mTokenValue == "?")
					{
						mTokenType = eCIFValueUnknown;
						mTokenValue.clear();
					}
				}
				break;

			case eStateDATA:
			case eStateSAVE:
				if (not isNonBlank(ch))
				{
					retract();

					if (state == eStateDATA)
						result = eCIFTokenDATA;
					else
						result = eCIFTokenSAVE;

					mTokenValue.erase(mTokenValue.begin(), mTokenValue.begin() + 5);
				}
				break;

			default:
				assert(false);
				error("Invalid state in getNextToken");
				break;
		}
	}

	if (VERBOSE >= 5)
	{
		std::cerr << kTokenName[result];
		if (mTokenType != eCIFValueUnknown)
			std::cerr << ' ' << kValueName[mTokenType];
		if (result != eCIFTokenEOF)
			std::cerr << " '" << mTokenValue << '\'';
		std::cerr << std::endl;
	}

	return result;
}


DatablockIndex SacParser::indexDatablocks()
{
	DatablockIndex index;

	// first locate the start, as fast as we can
	auto &sb = *mData.rdbuf();

	enum
	{
		start,
		comment,
		string,
		string_quote,
		qstring,
		data,
		data_name
	} state = start;

	int quote = 0;
	bool bol = true;
	const char dblk[] = "data_";
	std::string::size_type si = 0;
	std::string datablock;

	for (auto ch = sb.sbumpc(); ch != std::streambuf::traits_type::eof(); ch = sb.sbumpc())
	{
		switch (state)
		{
			case start:
				switch (ch)
				{
					case '#': state = comment; break;
					case 'd':
					case 'D':
						state = data;
						si = 1;
						break;
					case '\'':
					case '"':
						state = string;
						quote = ch;
						break;
					case ';':
						if (bol)
							state = qstring;
						break;
				}
				break;

			case comment:
				if (ch == '\n')
					state = start;
				break;

			case string:
				if (ch == quote)
					state = string_quote;
				break;

			case string_quote:
				if (std::isspace(ch))
					state = start;
				else
					state = string;
				break;

			case qstring:
				if (ch == ';' and bol)
					state = start;
				break;

			case data:
				if (dblk[si] == 0 and isNonBlank(ch))
				{
					datablock = {static_cast<char>(ch)};
					state = data_name;
				}
				else if (dblk[si++] != ch)
					state = start;
				break;

			case data_name:
				if (isNonBlank(ch))
					datablock.insert(datablock.end(), char(ch));
				else if (isspace(ch))
				{
					if (not datablock.empty())
						index[datablock] = mData.tellg();

					state = start;
				}
				else
					state = start;
				break;
		}

		bol = (ch == '\n');
	}

	return index;
}

bool SacParser::parseSingleDatablock(const std::string &datablock)
{
	// first locate the start, as fast as we can
	auto &sb = *mData.rdbuf();

	enum
	{
		start,
		comment,
		string,
		string_quote,
		qstring,
		data
	} state = start;

	int quote = 0;
	bool bol = true;
	std::string dblk = "data_" + datablock;
	std::string::size_type si = 0;
	bool found = false;

	for (auto ch = sb.sbumpc(); not found and ch != std::streambuf::traits_type::eof(); ch = sb.sbumpc())
	{
		switch (state)
		{
			case start:
				switch (ch)
				{
					case '#': state = comment; break;
					case 'd':
					case 'D':
						state = data;
						si = 1;
						break;
					case '\'':
					case '"':
						state = string;
						quote = ch;
						break;
					case ';':
						if (bol)
							state = qstring;
						break;
				}
				break;

			case comment:
				if (ch == '\n')
					state = start;
				break;

			case string:
				if (ch == quote)
					state = string_quote;
				break;

			case string_quote:
				if (std::isspace(ch))
					state = start;
				else
					state = string;
				break;

			case qstring:
				if (ch == ';' and bol)
					state = start;
				break;

			case data:
				if (isspace(ch) and dblk[si] == 0)
					found = true;
				else if (dblk[si++] != ch)
					state = start;
				break;
		}

		bol = (ch == '\n');
	}

	if (found)
	{
		produceDatablock(datablock);
		mLookahead = getNextToken();
		parseDataBlock();
	}

	return found;
}

bool SacParser::parseSingleDatablock(const std::string &datablock, const DatablockIndex &index)
{
	bool result = false;

	auto i = index.find(datablock);
	if (i != index.end())
	{
		mData.seekg(i->second);

		produceDatablock(datablock);
		mLookahead = getNextToken();
		parseDataBlock();

		result = true;
	}

	return result;
}

void SacParser::parseFile()
{
	while (mLookahead != eCIFTokenEOF)
	{
		switch (mLookahead)
		{
			case eCIFTokenGLOBAL:
				parseGlobal();
				break;

			case eCIFTokenDATA:
				produceDatablock(mTokenValue);

				match(eCIFTokenDATA);
				parseDataBlock();
				break;

			default:
				error("This file does not seem to be an mmCIF file");
				break;
		}
	}
}

void SacParser::parseGlobal()
{
	match(eCIFTokenGLOBAL);
	while (mLookahead == eCIFTokenTag)
	{
		match(eCIFTokenTag);
		match(eCIFTokenValue);
	}
}

void SacParser::parseDataBlock()
{
	std::string cat;

	while (mLookahead == eCIFTokenLOOP or mLookahead == eCIFTokenTag or mLookahead == eCIFTokenSAVE)
	{
		switch (mLookahead)
		{
			case eCIFTokenLOOP:
			{
				cat.clear(); // should start a new category

				match(eCIFTokenLOOP);

				std::vector<std::string> tags;

				while (mLookahead == eCIFTokenTag)
				{
					std::string catName, itemName;
					std::tie(catName, itemName) = splitTagName(mTokenValue);

					if (cat.empty())
					{
						produceCategory(catName);
						cat = catName;
					}
					else if (not iequals(cat, catName))
						error("inconsistent categories in loop_");

					tags.push_back(itemName);

					match(eCIFTokenTag);
				}

				while (mLookahead == eCIFTokenValue)
				{
					produceRow();

					for (auto tag : tags)
					{
						produceItem(cat, tag, mTokenValue);
						match(eCIFTokenValue);
					}
				}

				cat.clear();
				break;
			}

			case eCIFTokenTag:
			{
				std::string catName, itemName;
				std::tie(catName, itemName) = splitTagName(mTokenValue);

				if (not iequals(cat, catName))
				{
					produceCategory(catName);
					cat = catName;
					produceRow();
				}

				match(eCIFTokenTag);

				produceItem(cat, itemName, mTokenValue);

				match(eCIFTokenValue);
				break;
			}

			case eCIFTokenSAVE:
				parseSaveFrame();
				break;

			default:
				assert(false);
				break;
		}
	}
}

void SacParser::parseSaveFrame()
{
	error("A regular CIF file should not contain a save frame");
}

// --------------------------------------------------------------------

Parser::Parser(std::istream &is, File &f, bool init)
	: SacParser(is, init)
	, mFile(f)
	, mDataBlock(nullptr)
{
}

void Parser::produceDatablock(const std::string &name)
{
	mDataBlock = new Datablock(name);
	mFile.append(mDataBlock);
}

void Parser::produceCategory(const std::string &name)
{
	if (VERBOSE >= 4)
		std::cerr << "producing category " << name << std::endl;

	std::tie(mCat, std::ignore) = mDataBlock->emplace(name);
}

void Parser::produceRow()
{
	if (VERBOSE >= 4)
		std::cerr << "producing row for category " << mCat->name() << std::endl;

	mCat->emplace({});
	mRow = mCat->back();
	mRow.lineNr(mLineNr);
}

void Parser::produceItem(const std::string &category, const std::string &item, const std::string &value)
{
	if (VERBOSE >= 4)
		std::cerr << "producing _" << category << '.' << item << " -> " << value << std::endl;

	if (not iequals(category, mCat->name()))
		error("inconsistent categories in loop_");

	mRow[item] = mTokenValue;
}

// --------------------------------------------------------------------

struct DictParserDataImpl
{
	// temporary values for constructing dictionaries
	std::vector<ValidateCategory> mCategoryValidators;
	std::map<std::string, std::vector<ValidateItem>> mItemValidators;
	std::set<std::tuple<std::string, std::string>> mLinkedItems;
};

DictParser::DictParser(Validator &validator, std::istream &is)
	: Parser(is, mFile)
	, mValidator(validator)
	, mImpl(new DictParserDataImpl)
{
}

DictParser::~DictParser()
{
	delete mImpl;
}

void DictParser::parseSaveFrame()
{
	if (not mCollectedItemTypes)
		mCollectedItemTypes = collectItemTypes();

	std::string saveFrameName = mTokenValue;

	if (saveFrameName.empty())
		error("Invalid save frame, should contain more than just 'save_' here");

	bool isCategorySaveFrame = mTokenValue[0] != '_';

	Datablock dict(mTokenValue);
	Datablock::iterator cat = dict.end();

	match(eCIFTokenSAVE);
	while (mLookahead == eCIFTokenLOOP or mLookahead == eCIFTokenTag)
	{
		if (mLookahead == eCIFTokenLOOP)
		{
			cat = dict.end(); // should start a new category

			match(eCIFTokenLOOP);

			std::vector<std::string> tags;
			while (mLookahead == eCIFTokenTag)
			{
				std::string catName, itemName;
				std::tie(catName, itemName) = splitTagName(mTokenValue);

				if (cat == dict.end())
					std::tie(cat, std::ignore) = dict.emplace(catName);
				else if (not iequals(cat->name(), catName))
					error("inconsistent categories in loop_");

				tags.push_back(itemName);
				match(eCIFTokenTag);
			}

			while (mLookahead == eCIFTokenValue)
			{
				cat->emplace({});
				auto row = cat->back();

				for (auto tag : tags)
				{
					row[tag] = mTokenValue;
					match(eCIFTokenValue);
				}
			}

			cat = dict.end();
		}
		else
		{
			std::string catName, itemName;
			std::tie(catName, itemName) = splitTagName(mTokenValue);

			if (cat == dict.end() or not iequals(cat->name(), catName))
				std::tie(cat, std::ignore) = dict.emplace(catName);

			match(eCIFTokenTag);

			if (cat->empty())
				cat->emplace({});
			cat->back()[itemName] = mTokenValue;

			match(eCIFTokenValue);
		}
	}

	match(eCIFTokenSAVE);

	if (isCategorySaveFrame)
	{
		std::string category;
		cif::tie(category) = dict["category"].front().get("id");

		std::vector<std::string> keys;
		for (auto k : dict["category_key"])
			keys.push_back(std::get<1>(splitTagName(k["name"].as<std::string>())));

		iset groups;
		for (auto g : dict["category_group"])
			groups.insert(g["id"].as<std::string>());

		mImpl->mCategoryValidators.push_back(ValidateCategory{category, keys, groups});
	}
	else
	{
		// if the type code is missing, this must be a pointer, just skip it
		std::string typeCode;
		cif::tie(typeCode) = dict["item_type"].front().get("code");

		const ValidateType *tv = nullptr;
		if (not(typeCode.empty() or typeCode == "?"))
			tv = mValidator.getValidatorForType(typeCode);

		iset ess;
		for (auto e : dict["item_enumeration"])
			ess.insert(e["value"].as<std::string>());

		std::string defaultValue;
		cif::tie(defaultValue) = dict["item_default"].front().get("value");
		bool defaultIsNull = false;
		if (defaultValue.empty())
		{
			for (auto &r : dict["_item_default"])
			{
				defaultIsNull = r["value"].is_null();
				break;
			}
		}

		// collect the dict from our dataBlock and construct validators
		for (auto i : dict["item"])
		{
			std::string tagName, category, mandatory;

			cif::tie(tagName, category, mandatory) = i.get("name", "category_id", "mandatory_code");

			std::string catName, itemName;
			std::tie(catName, itemName) = splitTagName(tagName);

			if (catName.empty() or itemName.empty())
				error("Invalid tag name in _item.name " + tagName);

			if (not iequals(category, catName) and not(category.empty() or category == "?"))
				error("specified category id does match the implicit category name for tag '" + tagName + '\'');
			else
				category = catName;

			auto &ivs = mImpl->mItemValidators[category];

			auto vi = find(ivs.begin(), ivs.end(), ValidateItem{itemName});
			if (vi == ivs.end())
				ivs.push_back(ValidateItem{itemName, iequals(mandatory, "yes"), tv, ess, defaultValue, defaultIsNull});
			else
			{
				// need to update the itemValidator?
				if (vi->mMandatory != (iequals(mandatory, "yes")))
				{
					if (VERBOSE > 2)
					{
						std::cerr << "inconsistent mandatory value for " << tagName << " in dictionary" << std::endl;

						if (iequals(tagName, saveFrameName))
							std::cerr << "choosing " << mandatory << std::endl;
						else
							std::cerr << "choosing " << (vi->mMandatory ? "Y" : "N") << std::endl;
					}

					if (iequals(tagName, saveFrameName))
						vi->mMandatory = (iequals(mandatory, "yes"));
				}

				if (vi->mType != nullptr and tv != nullptr and vi->mType != tv)
				{
					if (VERBOSE > 1)
						std::cerr << "inconsistent type for " << tagName << " in dictionary" << std::endl;
				}

				//				vi->mMandatory = (iequals(mandatory, "yes"));
				if (vi->mType == nullptr)
					vi->mType = tv;

				vi->mEnums.insert(ess.begin(), ess.end());

				// anything else yet?
				// ...
			}
		}

		// collect the dict from our dataBlock and construct validators
		for (auto i : dict["item_linked"])
		{
			std::string childTagName, parentTagName;

			cif::tie(childTagName, parentTagName) = i.get("child_name", "parent_name");

			mImpl->mLinkedItems.emplace(childTagName, parentTagName);
		}
	}
}

void DictParser::linkItems()
{
	if (not mDataBlock)
		error("no datablock");

	auto &dict = *mDataBlock;

	// links are identified by a parent category, a child category and a group ID

	using key_type = std::tuple<std::string, std::string, int>;

	std::map<key_type, size_t> linkIndex;

	// Each link group consists of a set of keys
	std::vector<std::tuple<std::vector<std::string>, std::vector<std::string>>> linkKeys;

	auto addLink = [&](size_t ix, const std::string &pk, const std::string &ck)
	{
		auto &&[pkeys, ckeys] = linkKeys.at(ix);

		bool found = false;
		for (size_t i = 0; i < pkeys.size(); ++i)
		{
			if (pkeys[i] == pk and ckeys[i] == ck)
			{
				found = true;
				break;
			}
		}

		if (not found)
		{
			pkeys.push_back(pk);
			ckeys.push_back(ck);
		}
	};

	auto &linkedGroupList = dict["pdbx_item_linked_group_list"];

	for (auto gl : linkedGroupList)
	{
		std::string child, parent;
		int link_group_id;
		cif::tie(child, parent, link_group_id) = gl.get("child_name", "parent_name", "link_group_id");

		auto civ = mValidator.getValidatorForItem(child);
		if (civ == nullptr)
			error("in pdbx_item_linked_group_list, item '" + child + "' is not specified");

		auto piv = mValidator.getValidatorForItem(parent);
		if (piv == nullptr)
			error("in pdbx_item_linked_group_list, item '" + parent + "' is not specified");

		key_type key{piv->mCategory->mName, civ->mCategory->mName, link_group_id};
		if (not linkIndex.count(key))
		{
			linkIndex[key] = linkKeys.size();
			linkKeys.push_back({});
		}

		size_t ix = linkIndex.at(key);
		addLink(ix, piv->mTag, civ->mTag);
	}

	// Only process inline linked items if the linked group list is absent
	if (linkedGroupList.empty())
	{
		// for links recorded in categories but not in pdbx_item_linked_group_list
		for (auto li : mImpl->mLinkedItems)
		{
			std::string child, parent;
			std::tie(child, parent) = li;

			auto civ = mValidator.getValidatorForItem(child);
			if (civ == nullptr)
				error("in pdbx_item_linked_group_list, item '" + child + "' is not specified");

			auto piv = mValidator.getValidatorForItem(parent);
			if (piv == nullptr)
				error("in pdbx_item_linked_group_list, item '" + parent + "' is not specified");

			key_type key{piv->mCategory->mName, civ->mCategory->mName, 0};
			if (not linkIndex.count(key))
			{
				linkIndex[key] = linkKeys.size();
				linkKeys.push_back({});
			}

			size_t ix = linkIndex.at(key);
			addLink(ix, piv->mTag, civ->mTag);
		}
	}

	auto &linkedGroup = dict["pdbx_item_linked_group"];

	// now store the links in the validator
	for (auto &kv : linkIndex)
	{
		ValidateLink link = {};
		std::tie(link.mParentCategory, link.mChildCategory, link.mLinkGroupID) = kv.first;

		std::tie(link.mParentKeys, link.mChildKeys) = linkKeys[kv.second];

		// look up the label
		for (auto r : linkedGroup.find(cif::Key("category_id") == link.mChildCategory and cif::Key("link_group_id") == link.mLinkGroupID))
		{
			link.mLinkGroupLabel = r["label"].as<std::string>();
			break;
		}

		mValidator.addLinkValidator(std::move(link));
	}

	// now make sure the itemType is specified for all itemValidators

	for (auto &cv : mValidator.mCategoryValidators)
	{
		for (auto &iv : cv.mItemValidators)
		{
			if (iv.mType == nullptr and cif::VERBOSE >= 0)
				std::cerr << "Missing item_type for " << iv.mTag << std::endl;
		}
	}
}

void DictParser::loadDictionary()
{
	std::unique_ptr<Datablock> dict;
	Datablock *savedDatablock = mDataBlock;

	try
	{
		while (mLookahead != eCIFTokenEOF)
		{
			switch (mLookahead)
			{
				case eCIFTokenGLOBAL:
					parseGlobal();
					break;

				default:
				{
					dict.reset(new Datablock(mTokenValue)); // dummy datablock, for constructing the validator only
					mDataBlock = dict.get();

					match(eCIFTokenDATA);
					parseDataBlock();
					break;
				}
			}
		}
	}
	catch (const std::exception &)
	{
		if (cif::VERBOSE >= 0)
			std::cerr << "Error parsing dictionary" << std::endl;
		throw;
	}

	// store all validators
	for (auto &ic : mImpl->mCategoryValidators)
		mValidator.addCategoryValidator(std::move(ic));
	mImpl->mCategoryValidators.clear();

	for (auto &iv : mImpl->mItemValidators)
	{
		auto cv = mValidator.getValidatorForCategory(iv.first);
		if (cv == nullptr)
			error("Undefined category '" + iv.first);

		for (auto &v : iv.second)
			const_cast<ValidateCategory *>(cv)->addItemValidator(std::move(v));
	}

	// check all item validators for having a typeValidator

	if (dict)
		linkItems();

	// store meta information
	Datablock::iterator info;
	bool n;
	std::tie(info, n) = mDataBlock->emplace("dictionary");
	if (n)
	{
		auto r = info->front();
		mValidator.dictName(r["title"].as<std::string>());
		mValidator.dictVersion(r["version"].as<std::string>());
	}

	mDataBlock = savedDatablock;

	mImpl->mItemValidators.clear();
}

bool DictParser::collectItemTypes()
{
	bool result = false;

	if (not mDataBlock)
		error("no datablock");

	auto &dict = *mDataBlock;

	for (auto &t : dict["item_type_list"])
	{
		std::string code, primitiveCode, construct;
		cif::tie(code, primitiveCode, construct) = t.get("code", "primitive_code", "construct");

		cif::replace_all(construct, "\\n", "\n");
		cif::replace_all(construct, "\\t", "\t");
		cif::replace_all(construct, "\\\n", "");

		try
		{
			ValidateType v = {
				code, mapToPrimitiveType(primitiveCode), boost::regex(construct, boost::regex::extended | boost::regex::optimize)};

			mValidator.addTypeValidator(std::move(v));
		}
		catch (const std::exception &)
		{
			throw_with_nested(parse_error(t.lineNr(), "error in regular expression"));
		}

		// Do not replace an already defined type validator, this won't work with pdbx_v40
		// as it has a name that is too strict for its own names :-)
		//		if (mFileImpl.mTypeValidators.count(v))
		//			mFileImpl.mTypeValidators.erase(v);

		if (VERBOSE >= 5)
			std::cerr << "Added type " << code << " (" << primitiveCode << ") => " << construct << std::endl;

		result = true;
	}

	return result;
}

} // namespace cif
