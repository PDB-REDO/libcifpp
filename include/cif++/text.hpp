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

#include "cif++/exports.hpp"

#include <charconv>
#include <cmath>
#include <cstdint>
#include <limits>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>

#if __has_include(<experimental/type_traits>)
#include <experimental/type_traits>
#else
// sub optimal, but replicating the same code is worse
#include <zeep/type-traits.hpp>
#endif

namespace cif
{

// --------------------------------------------------------------------

// some basic utilities: Since we're using ASCII input only, we define for optimisation
// our own case conversion routines.

bool iequals(std::string_view a, std::string_view b);
int icompare(std::string_view a, std::string_view b);

bool iequals(const char *a, const char *b);
int icompare(const char *a, const char *b);

void to_lower(std::string &s);
std::string to_lower_copy(std::string_view s);

void to_upper(std::string &s);
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

template <typename StringType = std::string_view>
std::vector<StringType> split(std::string_view s, std::string_view separators, bool suppress_empty = false)
{
	std::vector<StringType> result;

	auto b = s.data();
	auto e = b;

	while (e != s.data() + s.length())
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

extern CIFPP_EXPORT const uint8_t kCharToLowerMap[256];

inline char tolower(int ch)
{
	return static_cast<char>(kCharToLowerMap[static_cast<uint8_t>(ch)]);
}

// --------------------------------------------------------------------

std::tuple<std::string, std::string> split_tag_name(std::string_view tag);

// --------------------------------------------------------------------
// generate a cif name, mainly used to generate asym_id's

std::string cif_id_for_number(int number);

// --------------------------------------------------------------------
//	custom wordwrapping routine

std::vector<std::string> word_wrap(const std::string &text, size_t width);

// --------------------------------------------------------------------
/// std::from_chars for floating point types.
/// These are optional, there's a selected_charconv class below that selects
/// the best option to used based on support by the stl library
/// I.e. that in case of GNU < 12 (or something) the cif implementation will
/// be used, all other cases will use the stl version.

template <typename FloatType, std::enable_if_t<std::is_floating_point_v<FloatType>, int> = 0>
std::from_chars_result from_chars(const char *first, const char *last, FloatType &value)
{
	std::from_chars_result result{ first, {} };

	enum State
	{
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
	scientific = 1,
	fixed = 2,
	// hex,
	general = fixed | scientific
};

template <typename FloatType, std::enable_if_t<std::is_floating_point_v<FloatType>, int> = 0>
std::to_chars_result to_chars(char *first, char *last, FloatType &value, chars_format fmt)
{
	int size = static_cast<int>(last - first);
	int r = 0;

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
	int size = static_cast<int>(last - first);
	int r = 0;

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

template <typename T>
struct my_charconv
{
	static std::from_chars_result from_chars(const char *a, const char *b, T &d)
	{
		return cif::from_chars(a, b, d);
	}

	static std::to_chars_result to_chars(char *first, char *last, T &value, chars_format fmt)
	{
		return cif::to_chars(first, last, value, fmt);
	}
};

template <typename T>
struct std_charconv
{
	static std::from_chars_result from_chars(const char *a, const char *b, T &d)
	{
		return std::from_chars(a, b, d);
	}

	static std::to_chars_result to_chars(char *first, char *last, T &value, chars_format fmt)
	{
		return std::to_chars(first, last, value, fmt);
	}
};

template <typename T>
using from_chars_function = decltype(std::from_chars(std::declval<const char *>(), std::declval<const char *>(), std::declval<T &>()));

template <typename T>
using selected_charconv = typename std::conditional_t<std::experimental::is_detected_v<from_chars_function, T>, std_charconv<T>, my_charconv<T>>;

} // namespace cif