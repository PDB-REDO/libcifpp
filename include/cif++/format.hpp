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

#include <string>

/**  \file format.hpp
 * 
 * File containing a basic reimplementation of boost::format
 * but then a bit more simplistic. Still this allowed me to move my code
 * from using boost::format to something without external dependency easily.
 */

namespace cif
{

namespace detail
{
	template <typename T>
	struct to_varg
	{
		using type = T;

		to_varg(const T &v)
			: m_value(v)
		{
		}

		type operator*() { return m_value; }

		T m_value;
	};

	template <>
	struct to_varg<const char *>
	{
		using type = const char *;

		to_varg(const char *v)
			: m_value(v)
		{
		}

		type operator*() { return m_value.c_str(); }

		std::string m_value;
	};

	template <>
	struct to_varg<std::string>
	{
		using type = const char *;

		to_varg(const std::string &v)
			: m_value(v)
		{
		}

		type operator*() { return m_value.c_str(); }

		std::string m_value;
	};

} // namespace

/** @cond */

template <typename... Args>
class format_plus_arg
{
  public:
	using args_vector_type = std::tuple<detail::to_varg<Args>...>;
	using vargs_vector_type = std::tuple<typename detail::to_varg<Args>::type...>;

	format_plus_arg(const format_plus_arg &) = delete;
	format_plus_arg &operator=(const format_plus_arg &) = delete;


	format_plus_arg(std::string_view fmt, Args... args)
		: m_fmt(fmt)
		, m_args(std::forward<Args>(args)...)
	{
		auto ix = std::make_index_sequence<sizeof...(Args)>();
		copy_vargs(ix);
	}

	std::string str()
	{
		char buffer[1024];
		std::string::size_type r = std::apply(snprintf, std::tuple_cat(std::make_tuple(buffer, sizeof(buffer), m_fmt.c_str()), m_vargs));
		return { buffer, r };
	}

	friend std::ostream &operator<<(std::ostream &os, const format_plus_arg &f)
	{
		char buffer[1024];
		std::string::size_type r = std::apply(snprintf, std::tuple_cat(std::make_tuple(buffer, sizeof(buffer), f.m_fmt.c_str()), f.m_vargs));
		os.write(buffer, r);
		return os;
	}

  private:

	template <std::size_t... I>
	void copy_vargs(std::index_sequence<I...>)
	{
		((std::get<I>(m_vargs) = *std::get<I>(m_args)), ...);
	}

	std::string m_fmt;
	args_vector_type m_args;
	vargs_vector_type m_vargs;
};

/** @endcond */

/**
 * @brief A simplistic reimplementation of boost::format, in fact it is
 * actually a way to call the C function snprintf to format the arguments
 * in @a args into the format string @a fmt
 * 
 * The string in @a fmt should thus be a C style format string.
 * 
 * TODO: Move to C++23 style of printing.
 * 
 * @tparam Args The types of the arguments
 * @param fmt The format string
 * @param args The arguments
 * @return An object that can be written out to a std::ostream using operator<<
 */

template <typename... Args>
constexpr auto format(std::string_view fmt, Args... args)
{
	return format_plus_arg(fmt, std::forward<Args>(args)...);
}

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

	~fill_out_streambuf()
	{
		m_os.rdbuf(m_upstream);
	}

	/** @endcond */

	/**
	 * @brief The magic happens here. Write out a couple of spaces when
	 * the last character to write is a newline to make the line as
	 * wide as the requested width.
	 */
	
	virtual int_type
	overflow(int_type ic = traits_type::eof())
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
	std::streambuf *get_upstream() const { return m_upstream; }

	/** Return how many lines have been written */
	int get_line_count() const { return m_line_count; }

  private:
	std::ostream &m_os;
	std::streambuf *m_upstream;
	int m_width;
	int m_line_count = 0;
	int m_column_count = 0;
};

} // namespace pdbx
