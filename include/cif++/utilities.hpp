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
 */

namespace cif
{

extern CIFPP_EXPORT int VERBOSE;

/// return the git 'build' number
std::string get_version_nr();

// --------------------------------------------------------------------

///	Return the width of the current output terminal, or the default 80 in case output is not to a terminal
uint32_t get_terminal_width();

// --------------------------------------------------------------------

///	Return the path of the current executable
std::string get_executable_path();

// --------------------------------------------------------------------
//	some manipulators to write coloured text to terminals

/// @brief The defined string colours 
enum class StringColour
{
	BLACK = 0,
	RED,
	GREEN,
	YELLOW,
	BLUE,
	MAGENTA,
	CYAN,
	WHITE,
	NONE = 9
};

template<StringColour C>
struct ColourDefinition
{
	static constexpr StringColour value = C;
};

enum class TextStyle
{
	REGULAR = 22,
	BOLD = 1
};

template<TextStyle S>
struct StyleDefinition
{
	static constexpr TextStyle value = S;
};

template<typename ForeColour, typename BackColour, typename Style>
struct ColourAndStyle
{
	static constexpr StringColour fore_colour = ForeColour::value;
	static constexpr StringColour back_colour = BackColour::value;
	static constexpr TextStyle text_style = Style::value;
	
	static constexpr int fore_colour_number = static_cast<int>(fore_colour) + 30;
	static constexpr int back_colour_number = static_cast<int>(back_colour) + 40;
	static constexpr int style_number = static_cast<int>(text_style);
	
	friend std::ostream &operator<<(std::ostream &os, ColourAndStyle clr)
	{
		bool use_colour = false;

		if (os.rdbuf() == std::cout.rdbuf() and isatty(STDOUT_FILENO))
			use_colour = true;
		else if (os.rdbuf() == std::cerr.rdbuf() and isatty(STDERR_FILENO))
			use_colour = true;

		if (use_colour)
		{
			if (fore_colour == StringColour::NONE and back_colour == StringColour::NONE)
				os << "\033[0m";
			else
				os << "\033[" << fore_colour_number << ';' << style_number << ';' << back_colour_number << 'm';
		}

		return os;
	}
};

template<typename ForeColour, typename BackColour>
constexpr auto coloured(const ForeColour fore, const BackColour back)
{
	return ColourAndStyle<ForeColour, BackColour, StyleDefinition<TextStyle::REGULAR>>{};
}

template<typename ForeColour, typename BackColour, typename Style>
constexpr auto coloured(const ForeColour fore, const BackColour back, Style style)
{
	return ColourAndStyle<ForeColour, BackColour, Style>{};
}

namespace colour
{
	constexpr ColourDefinition<StringColour::BLACK> black = ColourDefinition<StringColour::BLACK>();
	constexpr ColourDefinition<StringColour::RED> red = ColourDefinition<StringColour::RED>();
	constexpr ColourDefinition<StringColour::GREEN> green = ColourDefinition<StringColour::GREEN>();
	constexpr ColourDefinition<StringColour::YELLOW> yellow = ColourDefinition<StringColour::YELLOW>();
	constexpr ColourDefinition<StringColour::BLUE> blue = ColourDefinition<StringColour::BLUE>();
	constexpr ColourDefinition<StringColour::MAGENTA> magenta = ColourDefinition<StringColour::MAGENTA>();
	constexpr ColourDefinition<StringColour::CYAN> cyan = ColourDefinition<StringColour::CYAN>();
	constexpr ColourDefinition<StringColour::WHITE> white = ColourDefinition<StringColour::WHITE>();
	constexpr ColourDefinition<StringColour::NONE> none = ColourDefinition<StringColour::NONE>();

	constexpr StyleDefinition<TextStyle::REGULAR> regular = StyleDefinition<TextStyle::REGULAR>();
	constexpr StyleDefinition<TextStyle::BOLD> bold = StyleDefinition<TextStyle::BOLD>();

	constexpr auto reset = cif::coloured(none, none, regular);
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
