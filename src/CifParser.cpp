// cif parsing library

#include <set>

#include <boost/algorithm/string.hpp>

#include "cif++/Cif++.h"
#include "cif++/CifParser.h"
#include "cif++/CifValidator.h"

using namespace std;
namespace ba = boost::algorithm;

extern int VERBOSE;

namespace cif
{

const uint32_t kMaxLineLength = 132;

const uint8_t kCharTraitsTable[128] = {
	//	0	1	2	3	4	5	6	7	8	9	a	b	c	d	e	f
		14,	15,	14,	14,	14,	15,	15,	14,	15,	15,	15,	15,	15,	15,	15,	15,	//	2
		15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	10,	15,	15,	15,	15,	//	3
		15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	//	4
		15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	14,	15,	14,	15,	14,	//	5
		15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	//	6
		15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	0,	//	7
};

// --------------------------------------------------------------------

CifParserError::CifParserError(uint32_t lineNr, const string& message)
	: runtime_error("parse error at line " + to_string(lineNr) + ": " + message)
{
}

// --------------------------------------------------------------------

const char* SacParser::kTokenName[] = {
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

const char* SacParser::kValueName[] = {
	"Int",
	"Float",
	"Numeric",
	"String",
	"TextField",
	"Inapplicable",
	"Unknown"
};

// --------------------------------------------------------------------

SacParser::SacParser(std::istream& is)
	: mData(is)
{
	mValidate = true;
	mLineNr = 1;
	mBol = true;
	mLookahead = getNextToken();
}

void SacParser::error(const string& msg)
{
	throw CifParserError(mLineNr, msg);
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
		cerr << "getNextChar => ";
		if (iscntrl(result) or not isprint(result))
			cerr << int(result) << endl;
		else
			cerr << char(result) << endl;
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

void SacParser::restart()
{
	while (not mTokenValue.empty())
		retract();
	
	switch (mStart)
	{
		case eStateStart:
			mState = mStart = eStateFloat;
			break;
		
		case eStateFloat:
			mState = mStart = eStateInt;
			break;
		
		case eStateInt:
			mState = mStart = eStateValue;
			break;
		
		default:
			error("Invalid state in SacParser");
	}
	
	mBol = false;
}

void SacParser::match(SacParser::CIFToken t)
{
	if (mLookahead != t)
		error(string("Unexpected token, expected ") + kTokenName[t] + " but found " + kTokenName[mLookahead]);
	
	mLookahead = getNextToken();
}

SacParser::CIFToken SacParser::getNextToken()
{
	const auto kEOF = char_traits<char>::eof();
	
	CIFToken result = eCIFTokenUnknown;
	int quoteChar = 0;
	mState = mStart = eStateStart;
	mBol = false;
	
	mTokenValue.clear();
	mTokenType = eCIFValueUnknown;
	
	while (result == eCIFTokenUnknown)
	{
		auto ch = getNextChar();
		
		switch (mState)
		{
			case eStateStart:
				if (ch == kEOF)
					result = eCIFTokenEOF;
				else if (ch == '\n')
				{
					mBol = true;
					mState = eStateWhite;
				}
				else if (ch == ' ' or ch == '\t')
					mState = eStateWhite;
				else if (ch == '#')
					mState = eStateComment;
				else if (ch == '.')
					mState = eStateDot;
				else if (ch == '_')
					mState = eStateTag;
				else if (ch == ';' and mBol)
					mState = eStateTextField;
				else if (ch == '\'' or ch == '"')
				{
					quoteChar = ch;
					mState = eStateQuotedString;
				}
				else if (ch == '?')
					mState = eStateQuestionMark;
				else
					restart();
				break;
			
			case eStateWhite:
				if (ch == kEOF)
					result = eCIFTokenEOF;
				else if (not isspace(ch))
				{
					mState = eStateStart;
					retract();
					mTokenValue.clear();
				}
				else
					mBol = (ch == '\n');
				break;
			
			case eStateComment:
				if (ch == '\n')
				{
					mState = eStateStart;
					mBol = true;
					mTokenValue.clear();
				}
				else if (ch == kEOF)
					result = eCIFTokenEOF;
				else if (not isAnyPrint(ch))
					error("invalid character in comment");
				break;
			
			case eStateQuestionMark:
				if (isNonBlank(ch))
					mState = eStateValue;
				else
				{
					retract();
					result = eCIFTokenValue;
					mTokenValue.clear();
					mTokenType = eCIFValueUnknown;
				}
				break;

			case eStateDot:
				if (isdigit(ch))
					mState = eStateFloat + 2;
				else if (isspace(ch))
				{
					retract();
					result = eCIFTokenValue;
					mTokenType = eCIFValueInapplicable;
				}
				else
					mState = eStateValue;
				break;

			case eStateTextField:
				if (ch == '\n')
					mState = eStateTextField + 1;
				else if (ch == kEOF)
					error("unterminated textfield");
				else if (not isAnyPrint(ch))
//					error("invalid character in text field '" + string({ static_cast<char>(ch) }) + "' (" + to_string((int)ch) + ")");
					cerr << "invalid character in text field '" << string({ static_cast<char>(ch) }) << "' (" << ch << ") line: " << mLineNr << endl;
				break;
			
			case eStateTextField + 1:
				if (isTextLead(ch) or ch == ' ' or ch == '\t')
					mState = eStateTextField;
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
					mState = eStateQuotedStringQuote;
				else if (not isAnyPrint(ch))
					error("invalid character in quoted string");
				break;
			
			case eStateQuotedStringQuote:
				if (isWhite(ch))
				{
					retract();
					result = eCIFTokenValue;
					mTokenType = eCIFValueString;
					
					assert(mTokenValue.length() >= 3);
					mTokenValue = mTokenValue.substr(1, mTokenValue.length() - 2);
				}
				else if (ch == quoteChar)
					;
				else if (isAnyPrint(ch))
					mState = eStateQuotedString;
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
					mState = eStateFloat + 1;
				}
				else if (isdigit(ch))
					mState = eStateFloat + 1;
				else
					restart();
				break;
			
			case eStateFloat + 1:
//				if (ch == '(')	// numeric???
//					mState = eStateNumericSuffix;
//				else
				if (ch == '.')
					mState = eStateFloat + 2;
				else if (tolower(ch) == 'e')
					mState = eStateFloat + 3;
				else if (isWhite(ch) or ch == kEOF)
				{
					retract();
					result = eCIFTokenValue;
					mTokenType = eCIFValueInt;
				}
				else
					restart();
				break;
			
			// parsed '.'
			case eStateFloat + 2:
//				if (ch == '(')	// numeric???
//					mState = eStateNumericSuffix;
//				else
				if (tolower(ch) == 'e')
					mState = eStateFloat + 3;
				else if (isWhite(ch) or ch == kEOF)
				{
					retract();
					result = eCIFTokenValue;
					mTokenType = eCIFValueFloat;
				}
				else
					restart();
				break;
			
			// parsed 'e'
			case eStateFloat + 3:
				if (ch == '-' or ch == '+')
					mState = eStateFloat + 4;
				else if (isdigit(ch))
					mState = eStateFloat + 5;
				else
					restart();
				break;

			case eStateFloat + 4:
				if (isdigit(ch))
					mState = eStateFloat + 5;
				else
					restart();
				break;
			
			case eStateFloat + 5:
//				if (ch == '(')
//					mState = eStateNumericSuffix;
//				else
				if (isWhite(ch) or ch == kEOF)
				{
					retract();
					result = eCIFTokenValue;
					mTokenType = eCIFValueFloat;
				}
				else
					restart();
				break;
			
			case eStateInt:
				if (isdigit(ch) or ch == '+' or ch == '-')
					mState = eStateInt + 1;
				else
					restart();
				break;
			
			case eStateInt + 1:
				if (isWhite(ch) or ch == kEOF)
				{
					retract();
					result = eCIFTokenValue;
					mTokenType = eCIFValueInt;
				}
				else
					restart();
				break;
			
//			case eStateNumericSuffix:
//				if (isdigit(ch))
//					mState = eStateNumericSuffix + 1;
//				else
//					restart();
//				break;
//			
//			case eStateNumericSuffix + 1:
//				if (ch == ')')
//				{
//					result = eCIFTokenValue;
//					mTokenType = eCIFValueNumeric;
//				}
//				else if (not isdigit(ch))
//					restart();
//				break;
			
			case eStateValue:
				if (isNonBlank(ch))
					mState = eStateValue + 1;
				else
					error("invalid character at this position");
				break;
			
			case eStateValue + 1:
				if (ch == '_')		// first _, check for keywords
				{
					string s = toLowerCopy(mTokenValue);
					
					if (s == "global_")
						result = eCIFTokenGLOBAL;
					else if (s == "stop_")
						result = eCIFTokenSTOP;
					else if (s == "loop_")
						result = eCIFTokenLOOP;
					else if (s == "data_" or s == "save_")
						mState = eStateValue + 2;
				}
				else if (not isNonBlank(ch))
				{
					retract();
					result = eCIFTokenValue;
					mTokenType = eCIFValueString;
				}
				break;

			case eStateValue + 2:
				if (not isNonBlank(ch))
				{
					retract();
					
					if (tolower(mTokenValue[0]) == 'd')
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
		cerr << kTokenName[result];
		if (mTokenType != eCIFValueUnknown)
			cerr << ' ' << kValueName[mTokenType];
		if (result != eCIFTokenEOF)
			cerr << " '" << mTokenValue << '\'';
		cerr << endl;
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
	string cat;
	
	while (mLookahead == eCIFTokenLOOP or mLookahead == eCIFTokenTag or mLookahead == eCIFTokenSAVE)
	{
		switch (mLookahead)
		{
			case eCIFTokenLOOP:
			{
				cat.clear();	// should start a new category
				
				match(eCIFTokenLOOP);
				
				vector<string> tags;
				
				while (mLookahead == eCIFTokenTag)
				{
					string catName, itemName;
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
					
					for (auto tag: tags)
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
				string catName, itemName;
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

Parser::Parser(std::istream& is, File& f)
	: SacParser(is), mFile(f), mDataBlock(nullptr)
{
}

void Parser::produceDatablock(const string& name)
{
	mDataBlock = new Datablock(name);
	mFile.append(mDataBlock);
}

void Parser::produceCategory(const string& name)
{
	if (VERBOSE >= 4)
		cerr << "producing category " << name << endl;

	std::tie(mCat, ignore) = mDataBlock->emplace(name);
}

void Parser::produceRow()
{
	if (VERBOSE >= 4)
		cerr << "producing row for category " << mCat->name() << endl;

	mCat->emplace({});
	mRow = mCat->back();
	mRow.lineNr(mLineNr);
}

void Parser::produceItem(const string& category, const string& item, const string& value)
{
	if (VERBOSE >= 4)
		cerr << "producing _" << category << '.' << item << " -> " << value << endl;

	if (not iequals(category, mCat->name()))
		error("inconsistent categories in loop_");

	mRow[item] = mTokenValue;
}

// --------------------------------------------------------------------

struct DictParserDataImpl
{
	// temporary values for constructing dictionaries
	vector<ValidateCategory>			mCategoryValidators;
	map<string,vector<ValidateItem>>	mItemValidators;
	set<tuple<string,string>>			mLinkedItems;
};

DictParser::DictParser(Validator& validator, std::istream& is)
	: Parser(is, mFile), mValidator(validator), mImpl(new DictParserDataImpl)
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

	string saveFrameName = mTokenValue;

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
			cat = dict.end();	// should start a new category
				
			match(eCIFTokenLOOP);
			
			vector<string> tags;
			while (mLookahead == eCIFTokenTag)
			{
				string catName, itemName;
				std::tie(catName, itemName) = splitTagName(mTokenValue);
					
				if (cat == dict.end())
					std::tie(cat, ignore) = dict.emplace(catName);
				else if (not iequals(cat->name(), catName))
					error("inconsistent categories in loop_");
				
				tags.push_back(itemName);
				match(eCIFTokenTag);
			}
			
			while (mLookahead == eCIFTokenValue)
			{
				cat->emplace({});
				auto row = cat->back();
				
				for (auto tag: tags)
				{
					row[tag] = mTokenValue;
					match(eCIFTokenValue);
				}
			}
			
			cat = dict.end();
		}
		else
		{
			string catName, itemName;
			std::tie(catName, itemName) = splitTagName(mTokenValue);

			if (cat == dict.end() or not iequals(cat->name(), catName))
				std::tie(cat, ignore) = dict.emplace(catName);

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
		string category = dict.firstItem("_category.id");

		vector<string> keys;
		for (auto k: dict["category_key"])
			keys.push_back(get<1>(splitTagName(k["name"].as<string>())));
		
		iset groups;
		for (auto g: dict["category_group"])
			groups.insert(g["id"].as<string>());
			
		mImpl->mCategoryValidators.push_back(ValidateCategory{category, keys, groups});
	}
	else
	{
		// if the type code is missing, this must be a pointer, just skip it
		string typeCode = dict.firstItem("_item_type.code");

		const ValidateType* tv = nullptr;
		if (not (typeCode.empty() or typeCode == "?"))
			tv = mValidator.getValidatorForType(typeCode);

		iset ess;
		for (auto e: dict["item_enumeration"])
			ess.insert(e["value"].as<string>());
		
		string defaultValue = dict.firstItem("_item_default.value");
		
		// collect the dict from our dataBlock and construct validators
		for (auto i: dict["item"])
		{
			string tagName, category, mandatory;
			
			cif::tie(tagName, category, mandatory) = i.get("name", "category_id", "mandatory_code");
			
			string catName, itemName;
			std::tie(catName, itemName) = splitTagName(tagName);
			
			if (catName.empty() or itemName.empty())
				error("Invalid tag name in _item.name " + tagName);

			if (not iequals(category, catName) and not (category.empty() or category == "?"))
				error("specified category id does match the implicit category name for tag '" + tagName + '\'');
			else
				category = catName;
			
			auto& ivs = mImpl->mItemValidators[category];
			
			auto vi = find(ivs.begin(), ivs.end(), ValidateItem{itemName});
			if (vi == ivs.end())
				ivs.push_back(ValidateItem{itemName, iequals(mandatory, "yes"), tv, ess, defaultValue});
			else
			{
				// need to update the itemValidator?
				if (vi->mMandatory != (iequals(mandatory, "yes")))
				{
					if (VERBOSE > 2)
					{
						cerr << "inconsistent mandatory value for " << tagName << " in dictionary" << endl;
						
						if (iequals(tagName, saveFrameName))
							cerr << "choosing " << mandatory << endl;
						else
							cerr << "choosing " << (vi->mMandatory ? "Y" : "N") << endl;
					}

					if (iequals(tagName, saveFrameName))
						vi->mMandatory = (iequals(mandatory, "yes"));
				}

				if (vi->mType != nullptr and tv != nullptr and vi->mType != tv)
				{
					if (VERBOSE > 1)
						cerr << "inconsistent type for " << tagName << " in dictionary" << endl;
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
		for (auto i: dict["item_linked"])
		{
			string childTagName, parentTagName;
			
			cif::tie(childTagName, parentTagName) = i.get("child_name", "parent_name");
			
			mImpl->mLinkedItems.emplace(childTagName, parentTagName);
		}
	}
}

void DictParser::linkItems()
{
	if (not mDataBlock)
		error("no datablock");
	
	auto& dict = *mDataBlock;

	map<tuple<string,string,int>,size_t> linkIndex;
	vector<tuple<vector<string>,vector<string>>> linkKeys;
	
	for (auto gl: dict["pdbx_item_linked_group_list"])
	{
		string child, parent;
		int link_group_id;
		cif::tie(child, parent, link_group_id) = gl.get("child_name", "parent_name", "link_group_id");
		
		auto civ = mValidator.getValidatorForItem(child);
		if (civ == nullptr)
			error("in pdbx_item_linked_group_list, item '" + child + "' is not specified");
		
		auto piv = mValidator.getValidatorForItem(parent);
		if (piv == nullptr)
			error("in pdbx_item_linked_group_list, item '" + parent + "' is not specified");
		
		auto key = make_tuple(piv->mCategory->mName, civ->mCategory->mName, link_group_id);
		if (not linkIndex.count(key))
		{
			linkIndex[key] = linkKeys.size();
			linkKeys.push_back({});
		}
		
		size_t ix = linkIndex.at(key);
		
		get<0>(linkKeys.at(ix)).push_back(piv->mTag);
		get<1>(linkKeys.at(ix)).push_back(civ->mTag);
	}

//	for (auto li: mImpl->mLinkedItems)
//	{
//		string child, parent;
//		std::tie(child, parent) = li;
//		
//		auto civ = mValidator.getValidatorForItem(child);
//		if (civ == nullptr)
//			error("in pdbx_item_linked_group_list, item '" + child + "' is not specified");
//		
//		auto piv = mValidator.getValidatorForItem(parent);
//		if (piv == nullptr)
//			error("in pdbx_item_linked_group_list, item '" + parent + "' is not specified");
//		
//		auto key = make_tuple(piv->mCategory->mName, civ->mCategory->mName, piv->mTag);
//		if (not linkIndex.count(key))
//		{
//			linkIndex[key] = linkKeys.size();
//			linkKeys.push_back({});
//		}
//		
//		size_t ix = linkIndex.at(key);
//		auto& keys = linkKeys.at(ix);
//		
//		keys.insert(civ->mTag);
//	}
	
	auto& linkedGroup = dict["pdbx_item_linked_group"];

	// now store the links in the validator
	for (auto& kv: linkIndex)
	{
		ValidateLink link;
		std::tie(link.mParentCategory, link.mChildCategory, link.mLinkGroupId) = kv.first;
		
		std::tie(link.mParentKeys, link.mChildKeys) = linkKeys[kv.second];

		// look up the label
		for (auto r:  linkedGroup.find(cif::Key("category_id") == link.mChildCategory and cif::Key("link_group_id") == link.mLinkGroupId))
		{
			link.mLinkGroupLabel = r["label"].as<string>();
			break;
		}

		mValidator.addLinkValidator(move(link));
	}
	
	// now make sure the itemType is specified for all itemValidators
	
	for (auto& cv: mValidator.mCategoryValidators)
	{
		for (auto& iv: cv.mItemValidators)
		{
			if (iv.mType == nullptr)
				cerr << "Missing item_type for " << iv.mTag << endl;
		}
	}	
}

void DictParser::loadDictionary()
{
	unique_ptr<Datablock> dict;
	Datablock* savedDatablock = mDataBlock;
	
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
					dict.reset(new Datablock(mTokenValue));	// dummy datablock, for constructing the validator only
					mDataBlock = dict.get();
					
					match(eCIFTokenDATA);
					parseDataBlock();
					break;
				}
			}
		}
	}
	catch (const exception& ex)
	{
		cerr << "Error parsing dictionary" << endl;
		throw;
	}

	// store all validators
	for (auto& ic: mImpl->mCategoryValidators)
		mValidator.addCategoryValidator(move(ic));
	mImpl->mCategoryValidators.clear();
	
	for (auto& iv: mImpl->mItemValidators)
	{
		auto cv = mValidator.getValidatorForCategory(iv.first);
		if (cv == nullptr)
			error("Undefined category '" + iv.first);

		for (auto& v: iv.second)
			const_cast<ValidateCategory*>(cv)->addItemValidator(move(v));
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
		mValidator.dictName(r["title"].as<string>());
		mValidator.dictVersion(r["version"].as<string>());
	}

	mDataBlock = savedDatablock;

	mImpl->mItemValidators.clear();
}

bool DictParser::collectItemTypes()
{
	bool result = false;
	
	if (not mDataBlock)
		error("no datablock");
	
	auto& dict = *mDataBlock;
	
	for (auto& t: dict["item_type_list"])
	{
		string code, primitiveCode, construct;
		cif::tie(code, primitiveCode, construct) = t.get("code", "primitive_code", "construct");
		
		ba::replace_all(construct, "\\n", "\n");
		ba::replace_all(construct, "\\t", "\t");
		ba::replace_all(construct, "\\\n", "");
		
		try
		{
			ValidateType v = {
				code, mapToPrimitiveType(primitiveCode), boost::regex(construct, boost::regex::extended | boost::regex::optimize)
			};
			
			mValidator.addTypeValidator(move(v));
		}
		catch (const exception& ex)
		{
			throw_with_nested(CifParserError(t.lineNr(), "error in regular expression"));
		}

// Do not replace an already defined type validator, this won't work with pdbx_v40
// as it has a name that is too strict for its own names :-)
//		if (mFileImpl.mTypeValidators.count(v))
//			mFileImpl.mTypeValidators.erase(v);

		if (VERBOSE >= 5)
			cerr << "Added type " << code << " (" << primitiveCode << ") => " << construct << endl;
		
		result = true;
	}
	
	return result;
}


}
