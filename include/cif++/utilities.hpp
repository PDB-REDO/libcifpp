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

/// return the width of the current output terminal, or 80 if it cannot be determined
uint32_t get_terminal_width();

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
		struct coloured_string_t
		{
			/**
			 * @brief Construct a new coloured string t object
			 */
			coloured_string_t(std::string_view s, colour_type fc, colour_type bc, style_type st)
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
				else
					os << cs.m_str;

				return os;
			}

			/// @cond
			std::string_view m_str;
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

template <typename T>
	requires std::is_assignable_v<std::string_view, T>
inline auto coloured(T str,
	colour::colour_type fg, colour::colour_type bg = colour::colour_type::none,
	colour::style_type st = colour::style_type::regular)
{
	return colour::detail::coloured_string_t(str, fg, bg, st);
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
 * @brief Load a resource from disk or the compiled in resources
 * 
 * @verbatim embed:rst
.. note::

   See the :doc:`documentation on resources </resources>` for more information.

   @endverbatim
 * 
 * @param name The named resource to load
 * @return std::unique_ptr<std::istream> A pointer to the std::istream or empty if not found
 */

std::unique_ptr<std::istream> load_resource(std::filesystem::path name);

/**
 * @brief Add a file specified by @a dataFile as the data for resource @a name
 * 
 * @verbatim embed:rst
.. note::

   See the :doc:`documentation on resources </resources>` for more information.

   @endverbatim
 * 
 * @param name The name of the resource to specify
 * @param dataFile Path to a file containing the data
 */

void add_file_resource(const std::string &name, std::filesystem::path dataFile);

/**
 * @brief List all the file resources added with cif::add_file_resource.
 * 
 * @param os The std::ostream to write the directories to
 */

void list_file_resources(std::ostream &os);

/**
 * @brief Add a directory to the list of search directories. This list is
 * searched in a last-in-first-out order.
 * 
 * @verbatim embed:rst
.. note::

   See the :doc:`documentation on resources </resources>` for more information.

   @endverbatim
 */

void add_data_directory(std::filesystem::path dataDir);

/**
 * @brief List all the data directories, for error reporting on missing resources.
 * 
 * @param os The std::ostream to write the directories to
 */

void list_data_directories(std::ostream &os);

} // namespace cif
