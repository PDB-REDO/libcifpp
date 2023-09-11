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

#include <filesystem>
#include <iostream>

#ifndef STDOUT_FILENO
#define STDOUT_FILENO 1
#endif

#ifndef STDERR_FILENO
#define STDERR_FILENO 2
#endif

#if _WIN32
#include <io.h>
#define isatty _isatty
#else
#include <unistd.h>
#endif

#if _MSC_VER
#pragma warning(disable : 4996) // unsafe function or variable	(strcpy e.g.)
#pragma warning(disable : 4068) // unknown pragma
#pragma warning(disable : 4100) // unreferenced formal parameter
#pragma warning(disable : 4101) // unreferenced local variable
#define _SILENCE_CXX17_CODECVT_HEADER_DEPRECATION_WARNING 1
#endif

/** \file utilities.hpp
 *
 * This file contains code that is very generic in nature like a progress_bar
 * and classes you can use to colourise output text.
 */

namespace cif
{

/**
 * @brief The global variable VERBOSE contains the level of verbosity
 * requested. A value of 0 is normal, with some output on error conditions.
 * A value > 0 will result in more output, the higher the value, the more
 * output. A value < 0 will make the library silent, even in error
 * conditions.
 */
extern CIFPP_EXPORT int VERBOSE;

/// return the git 'build' number
std::string get_version_nr();

// --------------------------------------------------------------------

/**
 * When writing out text to the terminal it is often useful to have
 * some of the text colourised. But only if the output is really a
 * terminal since colouring text is done using escape sequences
 * an if output is redirected to a file, these escape sequences end up
 * in the file making the real text less easy to read.
 *
 * The code presented here is rather basic. It mimics the std::quoted
 * manipulator in that it will colour a string with optionally
 * requested colours and text style.
 *
 * Example:
 *
 * @code {.cpp}
 * using namespace cif::colour;
 * std::cout << cif::coloured("Hello, world!", white, red, bold) << '\n';
 * @endcode
 *
 */

namespace colour
{
	/// @brief The defined colours
	enum colour_type
	{
		black = 0,
		red,
		green,
		yellow,
		blue,
		magenta,
		cyan,
		white,
		none = 9
	};

	enum style_type
	{
		bold = 1,
		underlined = 4,
		blink = 5,
		inverse = 7,
		regular = 22,
	};

	namespace detail
	{
		/**
		 * @brief Struct for delimited strings.
		 */
		template <typename StringType>
		struct coloured_string_t
		{
			static_assert(std::is_reference_v<StringType> or std::is_pointer_v<StringType>,
				"String type must be pointer or reference");

			coloured_string_t(StringType s, colour_type fc, colour_type bc, style_type st)
				: m_str(s)
				, m_fore_colour(static_cast<int>(fc) + 30)
				, m_back_colour(static_cast<int>(bc) + 40)
				, m_style(static_cast<int>(st))
			{
			}

			coloured_string_t &operator=(coloured_string_t &) = delete;

			template <typename char_type, typename traits_type>
			friend std::basic_ostream<char_type, traits_type> &operator<<(
				std::basic_ostream<char_type, traits_type> &os, const coloured_string_t &cs)
			{
				bool use_colour = false;

				if (os.rdbuf() == std::cout.rdbuf() and isatty(STDOUT_FILENO))
					use_colour = true;
				else if (os.rdbuf() == std::cerr.rdbuf() and isatty(STDERR_FILENO))
					use_colour = true;

				if (use_colour)
				{
					os << "\033[" << cs.m_fore_colour << ';' << cs.m_style << ';' << cs.m_back_colour << 'm'
					   << cs.m_str
					   << "\033[0m";
				}

				return os;
			}

			StringType m_str;
			int m_fore_colour, m_back_colour;
			int m_style;
		};

	} // namespace detail
} // namespace colour

/**
 * @brief Manipulator for coloured strings.
 * @param str String to quote.
 * @param fg Foreground (=text) colour to use
 * @param bg Background colour to use
 * @param st Text style to use
 */
template <typename char_type>
inline auto coloured(const char_type *str,
	colour::colour_type fg, colour::colour_type bg = colour::colour_type::none,
	colour::style_type st = colour::style_type::regular)
{
	return colour::detail::coloured_string_t<const char_type *>(str, fg, bg, st);
}

template <typename char_type, typename traits_type, typename allocator_type>
inline auto coloured(const std::basic_string<char_type, traits_type, allocator_type> &str,
	colour::colour_type fg, colour::colour_type bg = colour::colour_type::none,
	colour::style_type st = colour::style_type::regular)
{
	return colour::detail::coloured_string_t<const std::basic_string<char_type, traits_type, allocator_type> &>(str, fg, bg, st);
}

template <typename char_type, typename traits_type, typename allocator_type>
inline auto coloured(std::basic_string<char_type, traits_type, allocator_type> &str,
	colour::colour_type fg, colour::colour_type bg = colour::colour_type::none,
	colour::style_type st = colour::style_type::regular)
{
	return colour::detail::coloured_string_t<std::basic_string<char_type, traits_type, allocator_type> &>(str, fg, bg, st);
}

template <typename char_type, typename traits_type>
inline auto coloured(std::basic_string_view<char_type, traits_type> &str,
	colour::colour_type fg, colour::colour_type bg = colour::colour_type::none,
	colour::style_type st = colour::style_type::regular)
{
	return colour::detail::coloured_string_t<std::basic_string_view<char_type, traits_type> &>(str, fg, bg, st);
}

// --------------------------------------------------------------------
//	A progress bar

class progress_bar
{
  public:
	progress_bar(int64_t inMax, const std::string &inAction);
	~progress_bar();

	void consumed(int64_t inConsumed); // consumed is relative
	void progress(int64_t inProgress); // progress is absolute

	void message(const std::string &inMessage);

  private:
	progress_bar(const progress_bar &) = delete;
	progress_bar &operator=(const progress_bar &) = delete;

	struct progress_bar_impl *m_impl;
};

// --------------------------------------------------------------------
// Resources

std::unique_ptr<std::istream> load_resource(std::filesystem::path name);
void add_file_resource(const std::string &name, std::filesystem::path dataFile);
void add_data_directory(std::filesystem::path dataDir);

} // namespace cif
