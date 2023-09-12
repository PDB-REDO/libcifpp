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
/// @brief For systems that lack this value
#define STDOUT_FILENO 1
#endif

#ifndef STDERR_FILENO
/// @brief For systems that lack this value
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

	/// @brief The defined styles
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

			/**
			 * @brief Construct a new coloured string t object
			 */
			coloured_string_t(StringType s, colour_type fc, colour_type bc, style_type st)
				: m_str(s)
				, m_fore_colour(static_cast<int>(fc) + 30)
				, m_back_colour(static_cast<int>(bc) + 40)
				, m_style(static_cast<int>(st))
			{
			}

			coloured_string_t &operator=(coloured_string_t &) = delete;

			/**
			 * @brief Write out the string, either coloured or not
			 */
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

			/// @cond
			StringType m_str;
			int m_fore_colour, m_back_colour;
			int m_style;
			/// @endcond
		};

	} // namespace detail
} // namespace colour

/**
 * @brief Manipulator for coloured strings.
 * 
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

/// @brief Manipulator for coloured strings.
template <typename char_type, typename traits_type, typename allocator_type>
inline auto coloured(const std::basic_string<char_type, traits_type, allocator_type> &str,
	colour::colour_type fg, colour::colour_type bg = colour::colour_type::none,
	colour::style_type st = colour::style_type::regular)
{
	return colour::detail::coloured_string_t<const std::basic_string<char_type, traits_type, allocator_type> &>(str, fg, bg, st);
}

/// @brief Manipulator for coloured strings.
template <typename char_type, typename traits_type, typename allocator_type>
inline auto coloured(std::basic_string<char_type, traits_type, allocator_type> &str,
	colour::colour_type fg, colour::colour_type bg = colour::colour_type::none,
	colour::style_type st = colour::style_type::regular)
{
	return colour::detail::coloured_string_t<std::basic_string<char_type, traits_type, allocator_type> &>(str, fg, bg, st);
}

/// @brief Manipulator for coloured strings.
template <typename char_type, typename traits_type>
inline auto coloured(std::basic_string_view<char_type, traits_type> &str,
	colour::colour_type fg, colour::colour_type bg = colour::colour_type::none,
	colour::style_type st = colour::style_type::regular)
{
	return colour::detail::coloured_string_t<std::basic_string_view<char_type, traits_type> &>(str, fg, bg, st);
}

// --------------------------------------------------------------------
//	A progress bar

/**
 * @brief A simple progress bar class for terminal based output
 * 
 * Using a progress bar is very convenient for the end user when
 * you have long running code. It gives feed back on how fast an
 * operation is performed and may give an indication how long it
 * will take before it is finished.
 * 
 * Using this cif::progress_bar implementation is straightforward:
 * 
 * @code {.cpp}
 * using namespace std::chrono_literals;
 * 
 * cif::progress_bar pb(10, "counting to ten");
 * 
 * for (int i = 1; i <= 10; ++i)
 * {
 *   pb.consumed(1);
 *   std::this_thread::sleep_for(1s);
 * }
 * 
 * @endcode
 * 
 * When the progress_bar is created, it first checks
 * to see if stdout is to a real TTY and if the VERBOSE
 * flag is not less than zero (quiet mode). If this passes
 * a thread is started that waits for updates.
 * 
 * The first two seconds, nothing is written to the screen
 * so if the work is finished within those two seconds
 * the screen stays clean.
 * 
 * After this time, a progress bar is printed that may look
 * like this:
 * 
 * @code
 * step 3           ========================--------------------------------  40% ⢁
 * @endcode
 * 
 * The first characters contain the initial action name or
 * the message text if it was used afterwards.
 * 
 * The thermometer is made up with '=' and '-' characters.
 * 
 * A percentage is also shown and at the end there is a spinner
 * that gives feedback that the program is really still working.
 * 
 * The progress bar is removed if the max has been reached
 * or if the progress bar is destructed. If any output has
 * been generated, the initial action is printed out along
 * with the total time spent.
 */

class progress_bar
{
  public:
	/**
	 * @brief Construct a new progress bar object
	 * 
	 * Progress ranges from 0 (zero) to @a inMax
	 * 
	 * The action in @a inAction is used for display
	 * 
	 * @param inMax The maximum value
	 * @param inAction The description of what is
	 * going on
	 */

	progress_bar(int64_t inMax, const std::string &inAction);

	/**
	 * @brief Destroy the progress bar object
	 * 
	 */
	~progress_bar();

	/**
	 * @brief Notify the progress bar that @a inConsumed
	 * should be added to the internal progress counter
	 */
	void consumed(int64_t inConsumed); // consumed is relative

	/**
	 * @brief Notify the progress bar that the internal
	 * progress counter should be updated to @a inProgress
	 */
	void progress(int64_t inProgress); // progress is absolute

	/**
	 * @brief Replace the action string in the progress bar
	 * with @a inMessage
	 */
	void message(const std::string &inMessage);

  private:
	progress_bar(const progress_bar &) = delete;
	progress_bar &operator=(const progress_bar &) = delete;

	struct progress_bar_impl *m_impl;
};

// --------------------------------------------------------------------
// Resources

/**
 * @brief Resources are files required to perform some action, e.g.
 * dictionary files or the entire CCD file.
 * 
 * Resources can be compiled into the executable so that the resulting
 * application can be made portable to other machines. For this you
 * need to use https://github.com/mhekkel/mrc.git which only works
 * on Un*x like systems using the ELF executable format or on MS Windows
 * 
 * But resources may also be located as files on the filesytem at
 * specific locations. And you can specify your own location for
 * files (a directory) or even override named resources with your
 * own data.
 * 
 * The order in which resources are search for is:
 * 
 * * Use the resource that was defined by calling add_file_resource
 *   for this name.
 * 
 * * Search the paths specified by add_data_directory, last one
 *   added is searched first
 * 
 * * Search the so-called CACHE_DIR. This location is defined
 *   at compile time and based on the installation directory of
 *   libcifpp. Usually it is /var/cache/libcifpp.
 *   It is in this directory where the cron job for libcifpp will
 *   put the updated files weekly.
 * 
 * * If the CCP4 environment is available, the
 *   $ENV{CCP4}/share/libcifpp is searched.
 * 
 * * If the environment variable LIBCIFPP_DATA_DIR is set it
 *   is searched
 * 
 * * The DATA_DIR is searched, this is also a variable defined
 *   at compile time, also based on the installation directory
 *   of libcifpp. It usually is /usr/share/libcifpp
 * 
 * * As a last resort an attempt is made to load the data from
 *   resources compiled by mrc.
 * 
 * @param name 
 * @return std::unique_ptr<std::istream> 
 */

std::unique_ptr<std::istream> load_resource(std::filesystem::path name);
void add_file_resource(const std::string &name, std::filesystem::path dataFile);
void add_data_directory(std::filesystem::path dataDir);

} // namespace cif
