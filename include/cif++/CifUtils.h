// cif parsing library

#pragma once

#include <vector>
#include <set>

#include "cif++/Config.h"

namespace cif
{

// some basic utilities: Since we're using ASCII input only, we define for optimisation
// our own case conversion routines.

bool iequals(const std::string& a, const std::string& b);
int icompare(const std::string& a, const std::string& b);

bool iequals(const char* a, const char* b);
int icompare(const char* a, const char* b);

void toLower(std::string& s);
std::string toLowerCopy(const std::string& s);

// To make life easier, we also define iless and iset using iequals

struct iless
{
	bool operator()(const std::string& a, const std::string& b) const
	{
		return icompare(a, b) < 0;
	}
};

typedef std::set<std::string, iless>	iset;

// --------------------------------------------------------------------
// This really makes a difference, having our own tolower routines

extern const uint8 kCharToLowerMap[256];

inline char tolower(char ch)
{
	return static_cast<char>(kCharToLowerMap[static_cast<uint8>(ch)]);
}

// --------------------------------------------------------------------

std::tuple<std::string,std::string> splitTagName(const std::string& tag);

// --------------------------------------------------------------------
//	custom wordwrapping routine

std::vector<std::string> wordWrap(const std::string& text, unsigned int width);

}
