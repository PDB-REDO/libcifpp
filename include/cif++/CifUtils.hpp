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

#include <charconv>
#include <cmath>
#include <filesystem>
#include <set>
#include <vector>

#ifndef STDOUT_FILENO
#define STDOUT_FILENO 1
#endif

#if _MSC_VER
#include <io.h>
#define isatty _isatty
#else
#include <unistd.h>
#endif

#include <cif++/Cif++Export.hpp>

#if _MSC_VER
#pragma warning(disable : 4996) // unsafe function or variable	(strcpy e.g.)
#pragma warning(disable : 4068) // unknown pragma
#pragma warning(disable : 4100) // unreferenced formal parameter
#pragma warning(disable : 4101) // unreferenced local variable
#define _SILENCE_CXX17_CODECVT_HEADER_DEPRECATION_WARNING 1
#endif

namespace cif
{

// the git 'build' number
std::string get_version_nr();
// std::string get_version_date();

// --------------------------------------------------------------------

// some basic utilities: Since we're using ASCII input only, we define for optimisation
// our own case conversion routines.

bool iequals(std::string_view a, std::string_view b);
int icompare(std::string_view a, std::string_view b);

bool iequals(const char *a, const char *b);
int icompare(const char *a, const char *b);

void toLower(std::string &s);
std::string toLowerCopy(std::string_view s);

void toUpper(std::string &s);
// std::string toUpperCopy(const std::string &s);

template <typename IterType>
std::string join(IterType b, IterType e, std::string_view sep)
{
	std::ostringstream s;

	if (b != e)
	{
		auto ai = b;
		auto ni = std::next(ai);

		for (;;)
		{
			s << *ai;

			if (ni == e)
				break;

			ai = ni;
			ni = std::next(ai);

			s << sep;
		}
	}

	return s.str();
}

template <typename V>
std::string join(const V &arr, std::string_view sep)
{
	return join(arr.begin(), arr.end(), sep);
}

template<typename StringType = std::string_view>
std::vector<StringType> split(std::string_view s, std::string_view separators, bool suppress_empty = false)
{
	std::vector<StringType> result;

	auto b = s.begin();
	auto e = b;

	while (e != s.end())
	{
		if (separators.find(*e) != std::string_view::npos)
		{
			if (e > b or not suppress_empty)
				result.emplace_back(b, e - b);
			b = e = e + 1;
			continue;
		}

		++e;
	}

	if (e > b or not suppress_empty)
		result.emplace_back(b, e - b);

	return result;
}

void replace_all(std::string &s, std::string_view what, std::string_view with = {});

#if defined(__cpp_lib_starts_ends_with)

inline bool starts_with(std::string s, std::string_view with)
{
	return s.starts_with(with);
}

inline bool ends_with(std::string_view s, std::string_view with)
{
	return s.ends_with(with);
}

#else

inline bool starts_with(std::string s, std::string_view with)
{
	return s.compare(0, with.length(), with) == 0;
}

inline bool ends_with(std::string_view s, std::string_view with)
{
	return s.length() >= with.length() and s.compare(s.length() - with.length(), with.length(), with) == 0;
}

#endif

#if defined(__cpp_lib_string_contains)

inline bool contains(std::string_view s, std::string_view q)
{
	return s.contains(q);
}

#else

inline bool contains(std::string_view s, std::string_view q)
{
	return s.find(q) != std::string_view::npos;
}

#endif

bool icontains(std::string_view s, std::string_view q);

void trim_left(std::string &s);
void trim_right(std::string &s);
void trim(std::string &s);

std::string trim_left_copy(std::string_view s);
std::string trim_right_copy(std::string_view s);
std::string trim_copy(std::string_view s);


// To make life easier, we also define iless and iset using iequals

struct iless
{
	bool operator()(const std::string &a, const std::string &b) const
	{
		return icompare(a, b) < 0;
	}
};

typedef std::set<std::string, iless> iset;

// --------------------------------------------------------------------
// This really makes a difference, having our own tolower routines

extern const uint8_t kCharToLowerMap[256];

inline char tolower(int ch)
{
	return static_cast<char>(kCharToLowerMap[static_cast<uint8_t>(ch)]);
}

// --------------------------------------------------------------------

std::tuple<std::string, std::string> splitTagName(std::string_view tag);

// --------------------------------------------------------------------
// generate a cif name, mainly used to generate asym_id's

std::string cifIdForNumber(int number);

// --------------------------------------------------------------------
//	custom wordwrapping routine

std::vector<std::string> wordWrap(const std::string &text, size_t width);

// --------------------------------------------------------------------
//	Code helping with terminal i/o

uint32_t get_terminal_width();

// --------------------------------------------------------------------
//	Path of the current executable

std::string get_executable_path();

// --------------------------------------------------------------------
//	some manipulators to write coloured text to terminals

enum StringColour
{
	scBLACK = 0,
	scRED,
	scGREEN,
	scYELLOW,
	scBLUE,
	scMAGENTA,
	scCYAN,
	scWHITE,
	scNONE = 9
};

template <typename String, typename CharT>
struct ColouredString
{
	static_assert(std::is_reference<String>::value or std::is_pointer<String>::value, "String type must be pointer or reference");

	ColouredString(String s, StringColour fore, StringColour back, bool bold = true)
		: m_s(s)
		, m_fore(fore)
		, m_back(back)
		, m_bold(bold)
	{
	}

	ColouredString &operator=(const ColouredString &) = delete;

	String m_s;
	StringColour m_fore, m_back;
	bool m_bold;
};

template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &os, const ColouredString<const CharT *, CharT> &s)
{
	if (isatty(STDOUT_FILENO))
	{
		std::basic_ostringstream<CharT, Traits> ostr;
		ostr << "\033[" << (30 + s.m_fore) << ';' << (s.m_bold ? "1" : "22") << ';' << (40 + s.m_back) << 'm'
			 << s.m_s
			 << "\033[0m";

		return os << ostr.str();
	}
	else
		return os << s.m_s;
}

template <typename CharT, typename Traits, typename String>
std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &os, const ColouredString<String, CharT> &s)
{
	if (isatty(STDOUT_FILENO))
	{
		std::basic_ostringstream<CharT, Traits> ostr;
		ostr << "\033[" << (30 + s.m_fore) << ';' << (s.m_bold ? "1" : "22") << ';' << (40 + s.m_back) << 'm'
			 << s.m_s
			 << "\033[0m";

		return os << ostr.str();
	}
	else
		return os << s.m_s;
}

template <typename CharT>
inline auto coloured(const CharT *s, StringColour fore = scWHITE, StringColour back = scRED, bool bold = true)
{
	return ColouredString<const CharT *, CharT>(s, fore, back, bold);
}

template <typename CharT, typename Traits, typename Alloc>
inline auto coloured(const std::basic_string<CharT, Traits, Alloc> &s, StringColour fore = scWHITE, StringColour back = scRED, bool bold = true)
{
	return ColouredString<const std::basic_string<CharT, Traits, Alloc>, CharT>(s, fore, back, bold);
}

template <typename CharT, typename Traits, typename Alloc>
inline auto coloured(std::basic_string<CharT, Traits, Alloc> &s, StringColour fore = scWHITE, StringColour back = scRED, bool bold = true)
{
	return ColouredString<std::basic_string<CharT, Traits, Alloc>, CharT>(s, fore, back, bold);
}

// --------------------------------------------------------------------

/// std::from_chars for floating point types.
template <typename FloatType, std::enable_if_t<std::is_floating_point_v<FloatType>, int> = 0>
std::from_chars_result from_chars(const char *first, const char *last, FloatType &value)
{
	std::from_chars_result result{first, {}};

	enum State {
		IntegerSign,
		Integer,
		Fraction,
		ExponentSign,
		Exponent
	} state = IntegerSign;
	int sign = 1;
	unsigned long long vi = 0;
	long double f = 1;
	int exponent_sign = 1;
	int exponent = 0;
	bool done = false;

	while (not done and result.ec == std::errc())
	{
		char ch = result.ptr != last ? *result.ptr : 0;
		++result.ptr;

		switch (state)
		{
			case IntegerSign:
				if (ch == '-')
				{
					sign = -1;
					state = Integer;
				}
				else if (ch == '+')
					state = Integer;
				else if (ch >= '0' and ch <= '9')
				{
					vi = ch - '0';
					state = Integer;
				}
				else if (ch == '.')
					state = Fraction;
				else
					result.ec = std::errc::invalid_argument;
				break;
			
			case Integer:
				if (ch >= '0' and ch <= '9')
					vi = 10 * vi + (ch - '0');
				else if (ch == 'e' or ch == 'E')
					state = ExponentSign;
				else if (ch == '.')
					state = Fraction;
				else
				{
					done = true;
					--result.ptr;
				}
				break;
			
			case Fraction:
				if (ch >= '0' and ch <= '9')
				{
					vi = 10 * vi + (ch - '0');
					f /= 10;
				}
				else if (ch == 'e' or ch == 'E')
					state = ExponentSign;
				else
				{
					done = true;
					--result.ptr;
				}
				break;

			case ExponentSign:
				if (ch == '-')
				{
					exponent_sign = -1;
					state = Exponent;
				}
				else if (ch == '+')
					state = Exponent;
				else if (ch >= '0' and ch <= '9')
				{
					exponent = ch - '0';
					state = Exponent;
				}
				else
					result.ec = std::errc::invalid_argument;
				break;
			
			case Exponent:
				if (ch >= '0' and ch <= '9')
					exponent = 10 * exponent + (ch - '0');
				else
				{
					done = true;
					--result.ptr;
				}
				break;
		}
	}

	if (result.ec == std::errc())
	{
		long double v = f * vi * sign;
		if (exponent != 0)
			v *= std::pow(10, exponent * exponent_sign);

		if (std::isnan(v))
			result.ec = std::errc::invalid_argument;
		else if (std::abs(v) > std::numeric_limits<FloatType>::max())
			result.ec = std::errc::result_out_of_range;

		value = static_cast<FloatType>(v);
	}

	return result;
}

enum class chars_format
{
    scientific		= 1,
    fixed			= 2,
    // hex,
    general = fixed | scientific
};

template <typename FloatType, std::enable_if_t<std::is_floating_point_v<FloatType>, int> = 0>
std::to_chars_result to_chars(char *first, char *last, FloatType &value, chars_format fmt)
{
	int size = last - first;
	int r;

	switch (fmt)
	{
		case chars_format::scientific:
			if constexpr (std::is_same_v<FloatType, long double>)
				r = snprintf(first, last - first, "%le", value);
			else
				r = snprintf(first, last - first, "%e", value);
			break;

		case chars_format::fixed:
			if constexpr (std::is_same_v<FloatType, long double>)
				r = snprintf(first, last - first, "%lf", value);
			else
				r = snprintf(first, last - first, "%f", value);
			break;

		case chars_format::general:
			if constexpr (std::is_same_v<FloatType, long double>)
				r = snprintf(first, last - first, "%lg", value);
			else
				r = snprintf(first, last - first, "%g", value);
			break;
	}

	std::to_chars_result result;
	if (r < 0 or r >= size)
		result = { first, std::errc::value_too_large };
	else
		result = { first + r, std::errc() };

	return result;
}

template <typename FloatType, std::enable_if_t<std::is_floating_point_v<FloatType>, int> = 0>
std::to_chars_result to_chars(char *first, char *last, FloatType &value, chars_format fmt, int precision)
{
	int size = last - first;
	int r;

	switch (fmt)
	{
		case chars_format::scientific:
			if constexpr (std::is_same_v<FloatType, long double>)
				r = snprintf(first, last - first, "%.*le", precision, value);
			else
				r = snprintf(first, last - first, "%.*e", precision, value);
			break;

		case chars_format::fixed:
			if constexpr (std::is_same_v<FloatType, long double>)
				r = snprintf(first, last - first, "%.*lf", precision, value);
			else
				r = snprintf(first, last - first, "%.*f", precision, value);
			break;

		case chars_format::general:
			if constexpr (std::is_same_v<FloatType, long double>)
				r = snprintf(first, last - first, "%.*lg", precision, value);
			else
				r = snprintf(first, last - first, "%.*g", precision, value);
			break;
	}

	std::to_chars_result result;
	if (r < 0 or r >= size)
		result = { first, std::errc::value_too_large };
	else
		result = { first + r, std::errc() };

	return result;
}

// --------------------------------------------------------------------
//	A progress bar

class Progress
{
  public:
	Progress(int64_t inMax, const std::string &inAction);
	virtual ~Progress();

	void consumed(int64_t inConsumed); // consumed is relative
	void progress(int64_t inProgress); // progress is absolute

	void message(const std::string &inMessage);

  private:
	Progress(const Progress &) = delete;
	Progress &operator=(const Progress &) = delete;

	struct ProgressImpl *mImpl;
};

// --------------------------------------------------------------------
// Resources

std::unique_ptr<std::istream> loadResource(std::filesystem::path name);
void addFileResource(const std::string &name, std::filesystem::path dataFile);
void addDataDirectory(std::filesystem::path dataDir);

} // namespace cif
