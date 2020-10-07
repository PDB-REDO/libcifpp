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

#include "cif++/Config.hpp"

#include <tuple>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <atomic>
#include <mutex>
#include <thread>
#include <filesystem>
#include <fstream>
#include <map>
#include <chrono>
#include <regex>
#include <iomanip>

#if defined(_MSC_VER)
#define TERM_WIDTH 80
#else
#include <termios.h>
#include <sys/ioctl.h>
#endif

#include <boost/algorithm/string.hpp>

#include "cif++/CifUtils.hpp"

namespace ba = boost::algorithm;
namespace fs = std::filesystem;

namespace cif
{

extern int VERBOSE;

// --------------------------------------------------------------------

std::string get_version_nr()
{
	const std::regex
		rxVersionNr(R"(build-(\d+)-g[0-9a-f]{7}(-dirty)?)");

#include "revision.hpp"

	struct membuf : public std::streambuf
	{
		membuf(char* data, size_t length)       { this->setg(data, data, data + length); }
	} buffer(const_cast<char*>(kRevision), sizeof(kRevision));

	std::istream is(&buffer);

	std::string line, result;

	while (getline(is, line))
	{
		std::smatch m;

		if (std::regex_match(line, m, rxVersionNr))
		{
			result = m[1];
			if (m[2].matched)
				result += '*';
			break;
		}
	}

	return result;
}

// --------------------------------------------------------------------
// This really makes a difference, having our own tolower routines

const uint8_t kCharToLowerMap[256] =
{
	0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 
	0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 
	0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27, 0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f, 
	0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f, 
	0x40, 0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69, 0x6a, 0x6b, 0x6c, 0x6d, 0x6e, 0x6f, 
	0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79, 0x7a, 0x5b, 0x5c, 0x5d, 0x5e, 0x5f, 
	0x60, 0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69, 0x6a, 0x6b, 0x6c, 0x6d, 0x6e, 0x6f, 
	0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79, 0x7a, 0x7b, 0x7c, 0x7d, 0x7e, 0x7f, 
	0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f, 
	0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f, 
	0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf, 
	0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf, 
	0xc0, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5, 0xc6, 0xc7, 0xc8, 0xc9, 0xca, 0xcb, 0xcc, 0xcd, 0xce, 0xcf, 
	0xd0, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7, 0xd8, 0xd9, 0xda, 0xdb, 0xdc, 0xdd, 0xde, 0xdf, 
	0xe0, 0xe1, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7, 0xe8, 0xe9, 0xea, 0xeb, 0xec, 0xed, 0xee, 0xef, 
	0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff
};

// --------------------------------------------------------------------

bool iequals(const std::string& a, const std::string& b)
{
	bool result = a.length() == b.length();
	for (auto ai = a.begin(), bi = b.begin(); result and ai != a.end() and bi != b.end(); ++ai, ++bi)
		result = tolower(*ai) == tolower(*bi);
	return result;
}

bool iequals(const char* a, const char* b)
{
	bool result = true;
	for (; result and *a and *b; ++a, ++b)
		result = tolower(*a) == tolower(*b);

	return result and *a == *b;
}

int icompare(const std::string& a, const std::string& b)
{
	int d = 0;
	auto ai = a.begin(), bi = b.begin();
	
	for (; d == 0 and ai != a.end() and bi != b.end(); ++ai, ++bi)
		d = tolower(*ai) - tolower(*bi);

	if (d == 0)
	{
		if (ai != a.end())
			d = 1;
		else if (bi != b.end())
			d = -1;
	}
	
	return d;
}

int icompare(const char* a, const char* b)
{
	int d = 0;
	
	for (; d == 0 and *a != 0 and *b != 0; ++a, ++b)
		d = tolower(*a) - tolower(*b);

	if (d == 0)
	{
		if (*a != 0)
			d = 1;
		else if (*b != 0)
			d = -1;
	}
	
	return d;
}

void toLower(std::string& s)
{
	for (auto& c: s)
		c = tolower(c);
}

std::string toLowerCopy(const std::string& s)
{
	std::string result(s);
	for (auto& c: result)
		c = tolower(c);
	return result;
}

// --------------------------------------------------------------------

std::tuple<std::string,std::string> splitTagName(const std::string& tag)
{
	if (tag.empty())
		throw std::runtime_error("empty tag");
	if (tag[0] != '_')
		throw std::runtime_error("tag does not start with underscore");

	auto s = tag.find('.');
	if (s == std::string::npos)
		throw std::runtime_error("tag does not contain dot");
	return std::tuple<std::string,std::string>{
		tag.substr(1, s - 1), tag.substr(s + 1)
	};
}	

// --------------------------------------------------------------------
// Simplified line breaking code taken from a decent text editor.
// In this case, simplified means it only supports ASCII.

enum LineBreakClass
{
	kLBC_OpenPunctuation,
	kLBC_ClosePunctuation,
	kLBC_CloseParenthesis,
	kLBC_Quotation,
	kLBC_NonBreaking,
	kLBC_Nonstarter,
	kLBC_Exlamation,
	kLBC_SymbolAllowingBreakAfter,
	kLBC_InfixNumericSeparator,
	kLBC_PrefixNumeric,
	kLBC_PostfixNumeric,
	kLBC_Numeric,
	kLBC_Alphabetic,
	kLBC_Ideographic,
	kLBC_Inseperable,
	kLBC_Hyphen,
	kLBC_BreakAfter,
	kLBC_BreakBefor,
	kLBC_BreakOpportunityBeforeAndAfter,
	kLBC_ZeroWidthSpace,
	kLBC_CombiningMark,
	kLBC_WordJoiner,
	kLBC_HangulLVSyllable,
	kLBC_HangulLVTSyllable,
	kLBC_HangulLJamo,
	kLBC_HangulVJamo,
	kLBC_HangulTJamo,

	kLBC_MandatoryBreak,
	kLBC_CarriageReturn,
	kLBC_LineFeed,
	kLBC_NextLine,
	kLBC_Surrogate,
	kLBC_Space,
	kLBC_ContigentBreakOpportunity,
	kLBC_Ambiguous,
	kLBC_ComplexContext,
	kLBC_Unknown
};

const LineBreakClass kASCII_LBTable[128] =
{
	kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark,
	kLBC_CombiningMark, kLBC_BreakAfter, kLBC_LineFeed, kLBC_MandatoryBreak, kLBC_MandatoryBreak, kLBC_CarriageReturn, kLBC_CombiningMark, kLBC_CombiningMark,
	kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark,
	kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark, kLBC_CombiningMark,
	kLBC_Space, kLBC_Exlamation, kLBC_Quotation, kLBC_Alphabetic, kLBC_PrefixNumeric, kLBC_PostfixNumeric, kLBC_Alphabetic, kLBC_Quotation,
	kLBC_OpenPunctuation, kLBC_CloseParenthesis, kLBC_Alphabetic, kLBC_PrefixNumeric, 
	
	// comma treated differently here, it is not a numeric separator in PDB
	kLBC_SymbolAllowingBreakAfter/*	kLBC_InfixNumericSeparator */,
	
	kLBC_Hyphen, kLBC_InfixNumericSeparator, kLBC_SymbolAllowingBreakAfter,
	kLBC_Numeric, kLBC_Numeric, kLBC_Numeric, kLBC_Numeric, kLBC_Numeric, kLBC_Numeric, kLBC_Numeric, kLBC_Numeric,
	kLBC_Numeric, kLBC_Numeric, kLBC_InfixNumericSeparator, kLBC_InfixNumericSeparator, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Exlamation,
	kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic,
	kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic,
	kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic,
	kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_OpenPunctuation, kLBC_PrefixNumeric, kLBC_CloseParenthesis, kLBC_Alphabetic, kLBC_Alphabetic,
	kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic,
	kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic,
	kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic,
	kLBC_Alphabetic, kLBC_Alphabetic, kLBC_Alphabetic, kLBC_OpenPunctuation, kLBC_BreakAfter, kLBC_ClosePunctuation, kLBC_Alphabetic, kLBC_CombiningMark
};

std::string::const_iterator nextLineBreak(std::string::const_iterator text, std::string::const_iterator end)
{
	if (text == end)
		return text;
	
	enum breakAction
	{ 
		DBK = 0, // direct break 	(blank in table)
		IBK, 	// indirect break	(% in table)
		PBK,	// prohibited break (^ in table)
		CIB,	// combining indirect break
		CPB		// combining prohibited break
	};

	const breakAction brkTable[27][27] = {
	//   	OP  	CL  	CP  	QU  	GL  	NS  	EX  	SY  	IS  	PR  	PO  	NU  	AL  	ID  	IN  	HY  	BA  	BB  	B2  	ZW  	CM  	WJ  	H2  	H3  	JL  	JV  	JT
/* OP */ { 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	CPB, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK, 	PBK },
/* CL */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	PBK, 	IBK, 	IBK, 	DBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* CP */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	IBK, 	DBK, 	DBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* QU */ { 	PBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	PBK, 	CIB, 	PBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK },
/* GL */ { 	IBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	PBK, 	CIB, 	PBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK },
/* NS */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* EX */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* SY */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	DBK, 	IBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* IS */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	DBK, 	IBK, 	IBK, 	DBK, 	DBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* PR */ { 	IBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	DBK, 	IBK, 	IBK, 	IBK, 	DBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK },
/* PO */ { 	IBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	DBK, 	IBK, 	IBK, 	DBK, 	DBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* NU */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	IBK, 	DBK, 	IBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* AL */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	DBK, 	IBK, 	IBK, 	DBK, 	IBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* ID */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	IBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* IN */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* HY */ { 	DBK, 	PBK, 	PBK, 	IBK, 	DBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	DBK, 	IBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* BA */ { 	DBK, 	PBK, 	PBK, 	IBK, 	DBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* BB */ { 	IBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	PBK, 	CIB, 	PBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK },
/* B2 */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK, 	DBK, 	PBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* ZW */ { 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* CM */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	DBK, 	IBK, 	IBK, 	DBK, 	IBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	DBK },
/* WJ */ { 	IBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK, 	PBK, 	CIB, 	PBK, 	IBK, 	IBK, 	IBK, 	IBK, 	IBK },
/* H2 */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	IBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK },
/* H3 */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	IBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	IBK },
/* JL */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	IBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	IBK, 	IBK, 	IBK, 	IBK, 	DBK },
/* JV */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	IBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK },
/* JT */ { 	DBK, 	PBK, 	PBK, 	IBK, 	IBK, 	IBK, 	PBK, 	PBK, 	PBK, 	DBK, 	IBK, 	DBK, 	DBK, 	DBK, 	IBK, 	IBK, 	IBK, 	DBK, 	DBK, 	PBK, 	CIB, 	PBK, 	DBK, 	DBK, 	DBK, 	DBK, 	IBK },
		};

	uint8_t ch = static_cast<uint8_t>(*text);

	LineBreakClass cls;
	
	if (ch == '\n')
		cls = kLBC_MandatoryBreak;
	else if (ch < 128)
	{
		cls = kASCII_LBTable[ch];
		if (cls > kLBC_MandatoryBreak and cls != kLBC_Space)	// duh...
			cls = kLBC_Alphabetic;
	}
	else
		cls = kLBC_Unknown;

	if (cls == kLBC_Space)
		cls = kLBC_WordJoiner;

	LineBreakClass ncls = cls;

	while (++text != end and cls != kLBC_MandatoryBreak)
	{
		ch = *text;
		
		LineBreakClass lcls = ncls;
		
		if (ch == '\n')
		{
			++text;
			break;
		}

		ncls = kASCII_LBTable[ch];
	
		if (ncls == kLBC_Space)
			continue;
		
		breakAction brk = brkTable[cls][ncls];
		
		if (brk == DBK or (brk == IBK and lcls == kLBC_Space))
			break;
		
		cls = ncls;
	}

	return text;
}

std::vector<std::string> wrapLine(const std::string& text, unsigned int width)
{
	std::vector<std::string> result;
	std::vector<size_t> offsets = { 0 };

	auto b = text.begin();
	while (b != text.end())
	{
		auto e = nextLineBreak(b, text.end());
		
		offsets.push_back(e - text.begin());
		
		b = e;
	}
	
	size_t count = offsets.size() - 1;
	
	std::vector<size_t> minima(count + 1, 1000000);
	minima[0] = 0;
	std::vector<size_t> breaks(count + 1, 0);
	
	for (size_t i = 0; i < count; ++i)
	{
		size_t j = i + 1;
		while (j <= count)
		{
			size_t w = offsets[j] - offsets[i];

			if (w > width)
				break;

			while (w > 0 and isspace(text[offsets[i] + w - 1]))
				--w;

			size_t cost = minima[i];
			if (j < count)	// last line may be shorter
				cost += (width - w) * (width - w);

			if (cost < minima[j])
			{
				minima[j] = cost;
				breaks[j] = i;
			}

			++j;
		}
	}
	
	size_t j = count;
	while (j > 0)
	{
		size_t i = breaks[j];
		result.push_back(text.substr(offsets[i], offsets[j] - offsets[i])); 
		j = i;
	}
	
	reverse(result.begin(), result.end());

	return result;
}

std::vector<std::string> wordWrap(const std::string& text, unsigned int width)
{
	std::vector<std::string> paragraphs;
	ba::split(paragraphs, text, ba::is_any_of("\n"));
	
	std::vector<std::string> result;
	for (auto& p: paragraphs)
	{
		if (p.empty())
		{
			result.push_back("");
			continue;
		}
		
		auto lines = wrapLine(p, width);
		result.insert(result.end(), lines.begin(), lines.end());
	}

	return result;
}

// --------------------------------------------------------------------

#ifdef _MSC_VER
uint32_t get_terminal_width()
{
	return TERM_WIDTH;
}

// I don't have a windows machine to test the following code, please accept my apologies in case it fails...
std::string GetExecutablePath()
{
	WCHAR buffer[4096];

	DWORD n = ::GetModuleFileNameW(nullptr, buffer, sizeof(buffer) / sizeof(WCHAR));
	if (n == 0)
		throw runtime_error("could not get exe path");

	wstring ws(buffer);

	return std::string(ws.begin(), ws.end());
}

#else
uint32_t get_terminal_width()
{
    uint32_t result = 80;

    if (isatty(STDOUT_FILENO))
    {
        struct winsize w;
        ioctl(0, TIOCGWINSZ, &w);
        result = w.ws_col;
    }
    return result;
}

std::string get_executable_path()
{
	using namespace std::literals;

	char path[PATH_MAX] = "";
	if (readlink("/proc/self/exe", path, sizeof(path)) == -1)
		throw std::runtime_error("could not get exe path "s + strerror(errno));
	return { path };
}

#endif

// --------------------------------------------------------------------

struct ProgressImpl
{
					ProgressImpl(int64_t inMax, const std::string& inAction)
						: mMax(inMax), mConsumed(0), mAction(inAction), mMessage(inAction)
						, mThread(std::bind(&ProgressImpl::Run, this)) {}

	void			Run();
	void			Stop()
	{
		mStop = true;
		if (mThread.joinable())
			mThread.join();
	}
	
	void			PrintProgress();
	void			PrintDone();

	int64_t				mMax;
	std::atomic<int64_t>mConsumed;
	int64_t				mLastConsumed = 0;
	int					mSpinnerIndex = 0;
	std::string			mAction, mMessage;
	std::mutex			mMutex;
	std::thread			mThread;
	std::chrono::time_point<std::chrono::system_clock>
						mStart = std::chrono::system_clock::now();
	bool				mStop = false;
};

void ProgressImpl::Run()
{
	bool printedAny = false;
	
	try
	{
		for (;;)
		{
			std::this_thread::sleep_for(std::chrono::milliseconds(100));

			std::unique_lock lock(mMutex);
			
			if (mStop or mConsumed == mMax)
				break;

			auto elapsed = std::chrono::system_clock::now() - mStart;

			if (elapsed < std::chrono::seconds(5))
				continue;

			PrintProgress();
			printedAny = true;

		}
	}
	catch (...) {}
	
	if (printedAny)
		PrintDone();
}

void ProgressImpl::PrintProgress()
{
//	const char* kBlocks[] = {
//		" ",				// 0
//		u8"\u258F",			// 1
//		u8"\u258E",			// 2
//		u8"\u258D",			// 3
//		u8"\u258C",			// 4
//		u8"\u258B",			// 5
//		u8"\u258A",			// 6
//		u8"\u2589",			// 7
//		u8"\u2588",			// 8
//	};

	const char* kBlocks[] = {
		" ",			// 0
		" ",			// 1
		" ",			// 2
		"-",			// 3
		"-",			// 4
		"-",			// 5
		"=",			// 6
		"=",			// 7
		"=",			// 8
	};

	uint32_t width = get_terminal_width();
	
	std::string msg;
	msg.reserve(width + 1);
	if (mMessage.length() <= 20)
	{
		msg = mMessage;
		if (msg.length() < 20)
			msg.append(20 - msg.length(), ' ');
	}
	else
		msg = mMessage.substr(0, 17) + "...";
	
	msg += " |";
	
	int64_t consumed = mConsumed;
	float progress = static_cast<float>(consumed) / mMax;
	int pi = static_cast<int>(std::ceil(progress * 33 * 8));
//	int tw = width - 28;
//	int twd = static_cast<int>(tw * progress + 0.5f);
//	msg.append(twd, '=');
//	msg.append(tw - twd, ' ');

	for (int i = 0; i < 33; ++i)
	{
		if (pi <= 0)
			msg += kBlocks[0];
		else if (pi >= 8)
			msg += kBlocks[8];
		else
			msg += kBlocks[pi];
		pi -= 8;
	}

	msg.append("| ");

//	const char		kSpinner[] = { '|', '/', '-', '\\' };
	const char		kSpinner[] = { ' ', '.', 'o', 'O', '0', 'O', 'o', '.' };
	const size_t	kSpinnerCount = sizeof(kSpinner);

	if (mLastConsumed < consumed)
	{
		mLastConsumed = consumed;
		mSpinnerIndex = (mSpinnerIndex + 1) % kSpinnerCount;
	}
	
	const char spinner[2] = { kSpinner[mSpinnerIndex], 0 };
	msg.append(spinner);

//	int perc = static_cast<int>(100 * progress);
//	if (perc < 100)
//		msg += ' ';
//	if (perc < 10)
//		msg += ' ';
//	msg += to_string(perc);
//	msg += '%';
	
	std::cout << '\r' << msg;
	std::cout.flush();
}

namespace
{

std::ostream& operator<<(std::ostream& os, const std::chrono::duration<double>& t)
{
	uint64_t s = static_cast<uint64_t>(std::trunc(t.count()));
	if (s > 24 * 60 * 60)
	{
		uint32_t days = s / (24 * 60 * 60);
		os << days << "d ";
		s %= 24 * 60 * 60;
	}
	
	if (s > 60 * 60)
	{
		uint32_t hours = s / (60 * 60);
		os << hours << "h ";
		s %= 60 * 60;
	}
	
	if (s > 60)
	{
		uint32_t minutes = s / 60;
		os << minutes << "m ";
		s %= 60;
	}
	
	double ss = s + 1e-6 * (t.count() - s);
	
	os << std::fixed << std::setprecision(1) << ss << 's';

	return os;
}

}

void ProgressImpl::PrintDone()
{
	std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - mStart;

	std::ostringstream msgstr;
	msgstr << mAction << " done in " << elapsed << " cpu / %ws wall";
	auto msg = msgstr.str();

	uint32_t width = get_terminal_width();

	if (msg.length() < width)
		msg += std::string(width - msg.length(), ' ');
	
	std::cout << '\r' << msg << std::endl;
}

Progress::Progress(int64_t inMax, const std::string& inAction)
	: mImpl(nullptr)
{
	if (isatty(STDOUT_FILENO))
		mImpl = new ProgressImpl(inMax, inAction);
}

Progress::~Progress()
{
	if (mImpl != nullptr)
		mImpl->Stop();

	delete mImpl;
}
	
void Progress::consumed(int64_t inConsumed)
{
	if (mImpl != nullptr and 
		(mImpl->mConsumed += inConsumed) >= mImpl->mMax)
	{
		mImpl->Stop();
	}
}

void Progress::progress(int64_t inProgress)
{
	if (mImpl != nullptr and 
		(mImpl->mConsumed = inProgress) >= mImpl->mMax)
	{
		mImpl->Stop();
	}
}

void Progress::message(const std::string& inMessage)
{
	if (mImpl != nullptr)
	{
		std::unique_lock lock(mImpl->mMutex);
		mImpl->mMessage = inMessage;
	}
}

// --------------------------------------------------------------------

std::unique_ptr<rsrc_loader> rsrc_loader::s_instance;

struct rsrc_loader_impl
{
	virtual ~rsrc_loader_impl() {}
	virtual rsrc load(const std::string& path) const = 0;
};

class rsrc_not_found_exception : public std::exception
{
  public:
    virtual const char* what() const throw()		{ return "resource not found"; }
};

struct mrsrc_rsrc_loader_impl : public rsrc_loader_impl
{
	mrsrc_rsrc_loader_impl(const rsrc_loader_info& info);

	virtual rsrc load(const std::string& path) const;

	const rsrc_imp*	m_index;
	const char*		m_data;
	const char*		m_name;
};

mrsrc_rsrc_loader_impl::mrsrc_rsrc_loader_impl(const rsrc_loader_info& info)
{
	m_index = static_cast<const rsrc_imp*>(info.ptrs[0]);
	m_data = reinterpret_cast<const char*>(info.ptrs[1]);
	m_name = reinterpret_cast<const char*>(info.ptrs[2]);
}

rsrc mrsrc_rsrc_loader_impl::load(const std::string& path) const
{
	auto node = m_index;
	
	std::string p(path);
	
	while (not p.empty())
	{
		if (node->m_child == 0)	// no children, this is an error
			return {};

		node = m_index + node->m_child;
		
		std::string::size_type s = p.find('/');
		std::string name;
		
		if (s != std::string::npos)
		{
			name = p.substr(0, s);
			p.erase(0, s + 1);
		}
		else
			std::swap(name, p);
		
		while (name != m_name + node->m_name)
		{
			if (node->m_next == 0)
				return {};
			node = m_index + node->m_next;
		}
	}

	return { m_data + node->m_data, node->m_size };
}

struct mrsrc_file_loader_impl : public rsrc_loader_impl
{
	mrsrc_file_loader_impl(const rsrc_loader_info& info);

	virtual rsrc load(const std::string& path) const;

	fs::path m_rsrc_dir;
	mutable std::mutex m_mutex;
	mutable std::map<fs::path,std::string> m_cache;
};

mrsrc_file_loader_impl::mrsrc_file_loader_impl(const rsrc_loader_info& info)
	: m_rsrc_dir(info.info)
{
	if (not fs::is_directory(m_rsrc_dir))
	{
		std::cerr << "The directory " << info.info << " does not exist" << std::endl
				  << "Loading dictionaries will most likely fail" << std::endl;
	}
}

rsrc mrsrc_file_loader_impl::load(const std::string& path) const
{
	std::unique_lock lock(m_mutex);
	
	fs::path p = m_rsrc_dir / path;

	if (fs::exists(p))
	{
		auto i = m_cache.find(p);

		if (i == m_cache.end())
		{
			std::ostringstream s;
			std::ifstream f(p);
			if (not f.is_open())
				throw std::runtime_error("Could not open " + p.string());

			std::string line;
			while (std::getline(f, line))
				s << line << std::endl;
			
			std::string data = s.str();
			std::tie(i, std::ignore) = m_cache.emplace(p, data);
		}

		return { i->second.c_str(), i->second.length() };
	}

	return {};
}

rsrc_loader::rsrc_loader(const std::initializer_list<rsrc_loader_info>& loaders)
{
	for (auto& l: loaders)
	{
		try
		{
			switch (l.type)
			{
				case rsrc_loader_type::file:
					m_loaders.push_back(new mrsrc_file_loader_impl(l));
					break;

				case rsrc_loader_type::mrsrc:
					m_loaders.push_back(new mrsrc_rsrc_loader_impl(l));
					break;
				
				default:
					break;
			}
		}
		catch (const std::exception& e)
		{
			if (VERBOSE)
				std::cerr << e.what() << std::endl;
		}
	}
}

rsrc_loader::~rsrc_loader()
{
	for (auto l: m_loaders)
		delete l;
}

rsrc rsrc_loader::do_load(const std::string& name)
{
	rsrc result;

	for (auto loader: m_loaders)
	{
		result = loader->load(name);
		if (result)
			break;
	}

	return result;
}

}
