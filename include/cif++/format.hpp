/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2022 NKI/AVL, Netherlands Cancer Institute
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

#include <ostream>
#include <streambuf>


/**  \file format.hpp
 * 
 * Now using std::format instead of a home grown rip off
 */

namespace cif
{

// --------------------------------------------------------------------
/// A streambuf that fills out lines with spaces up until a specified width

class fill_out_streambuf : public std::streambuf
{
  public:

	/** @cond */

	using base_type = std::streambuf;
	using int_type = base_type::int_type;
	using char_type = base_type::char_type;
	using traits_type = base_type::traits_type;

	/** @endcond */

	/**
	 * @brief Construct a new fill out streambuf object based on ostream @a os and a
	 * width to fill out to of @a width
	 */
	fill_out_streambuf(std::ostream &os, int width = 80)
		: m_os(os)
		, m_upstream(os.rdbuf())
		, m_width(width)
	{
	}

	/** @cond */

	~fill_out_streambuf() override
	{
		m_os.rdbuf(m_upstream);
	}

	/** @endcond */

	/**
	 * @brief The magic happens here. Write out a couple of spaces when
	 * the last character to write is a newline to make the line as
	 * wide as the requested width.
	 */
	
	int_type overflow(int_type ic = traits_type::eof()) override
	{
		char ch = traits_type::to_char_type(ic);

		int_type result = ic;

		if (ch == '\n')
		{
			for (int i = m_column_count; result != traits_type::eof() and i < m_width; ++i)
				result = m_upstream->sputc(' ');
		}

		if (result != traits_type::eof())
			result = m_upstream->sputc(ch);

		if (result != traits_type::eof())
		{
			if (ch == '\n')
			{
				m_column_count = 0;
				++m_line_count;
			}
			else
				++m_column_count;
		}

		return result;
	}

	/** Return the upstream streambuf */
	[[nodiscard]] std::streambuf *get_upstream() const { return m_upstream; }

	/** Return how many lines have been written */
	[[nodiscard]] int get_line_count() const { return m_line_count; }

  private:
	std::ostream &m_os;
	std::streambuf *m_upstream;
	int m_width;
	int m_line_count = 0;
	int m_column_count = 0;
};

} // namespace pdbx
