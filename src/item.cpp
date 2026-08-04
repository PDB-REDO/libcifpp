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

#include "cif++/item.hpp"

#include "cif++/row.hpp"
#include "cif++/text.hpp"

#include <algorithm>
#include <cassert>
#include <charconv>
#include <cmath>
#include <compare>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace cif
{

bool item_handle::empty() const
{
	return m_item_ix >= m_row.size() or m_row[m_item_ix].empty();
}

item_value &item_handle::value()
{
#ifndef NDEBUG
	assert(m_item_ix < m_row.size());
#else
	if (m_item_ix >= m_row.size())
		throw std::runtime_error("Invalid item handle");
#endif
	return m_row.operator[](m_item_ix);
}

const item_value &item_handle::value() const
{
#ifndef NDEBUG
	assert(m_item_ix < m_row.size());
#else
	if (m_item_ix >= m_row.size())
		throw std::runtime_error("Invalid item handle");
#endif
	return m_row.operator[](m_item_ix);
}

void swap(item_handle a, item_handle b) noexcept
{
	item_value v(std::move(a.value()));
	a.value() = std::move(b.value());
	b.value() = std::move(v);
}

void item_handle::set(item_value value, bool updateLinked)
{
	row_handle rh{ m_category, m_row };
	rh.assign(m_item_ix, std::move(value), updateLinked);
}

int item_value::compare(const item_value &b, bool ignore_case) const noexcept
{
	int d = static_cast<int>(m_data.m_type) - static_cast<int>(b.m_data.m_type);

	if (d == 0)
	{
		switch (m_data.m_type)
		{
			using enum item_value_type;

			case INT:
				d = m_data.m_value.m_integer - b.m_data.m_value.m_integer;
				break;
			case FLOAT:
				// stupid comparison based on chopped textual representation
				if (m_data.m_len > 0 or b.m_data.m_len > 0)
				{
					double fa = m_data.m_value.m_float;
					double fb = b.m_data.m_value.m_float;

					auto delta = std::abs(fa - fb);
					if (delta == 0 or std::isnan(delta))
						d = 0;
					else if (m_data.m_len and b.m_data.m_len)
					{
						auto epsilon = std::pow(10.0f, -1.0f * std::min(m_data.m_len, b.m_data.m_len));
						if (delta > epsilon)
							d = fa < fb ? -1 : 1;
						else
							d = 0;
					}
					else
					{
						auto dp = (m_data.m_value.m_float <=> b.m_data.m_value.m_float);
						if (dp == std::partial_ordering::less)
							d = -1;
						else if (dp == std::partial_ordering::greater)
							d = 1;
					}
				}
				else
				{
					auto dp = (m_data.m_value.m_float <=> b.m_data.m_value.m_float);
					if (dp == std::partial_ordering::less)
						d = -1;
					else if (dp == std::partial_ordering::greater)
						d = 1;
				}
				break;
			case TEXT:
				d = ignore_case
				        ? icompare(m_data.sv(), b.m_data.sv())
				        : m_data.sv().compare(b.m_data.sv());
				break;
			default:;
		}
	}
	else if (is_number() and b.is_number())
	{
		std::partial_ordering dp = std::partial_ordering::equivalent;

		if (is_number_float())
			dp = m_data.m_value.m_float <=> b.m_data.m_value.m_integer;
		else /*  if (is_number_int()) */
			dp = m_data.m_value.m_integer <=> b.m_data.m_value.m_float;

		if (dp == std::partial_ordering::less)
			d = -1;
		else if (dp == std::partial_ordering::greater)
			d = 1;
		else
			d = 0;
	}
	else if (is_number_int() and b.is_string())
		d = str().compare(b.m_data.sv());
	else if (is_string() and b.is_number_int())
		d = m_data.sv().compare(b.str());

	return d;
}

std::string item_value::str() const
{
	switch (m_data.m_type)
	{
		using enum item_value_type;

		case MISSING:
			return "?";

		case INAPPLICABLE:
			return ".";

		case TEXT:
			return std::string{ m_data.sv() };

		case INT:
		{
			char s[32];
			std::to_chars_result r = std::to_chars(s, s + sizeof(s), m_data.m_value.m_integer);
			return r.ec == std::errc{} ? std::string{ s, r.ptr } : "*****";
		}

		case FLOAT:
		{
			char s[32];

			std::to_chars_result r;

			if (m_data.m_len)
			{
				r = std::to_chars(s, s + sizeof(s), m_data.m_value.m_float, std::chars_format::fixed, m_data.m_len);
				if (r.ec != std::errc{})
					r = std::to_chars(s, s + sizeof(s), m_data.m_value.m_float);
			}
			else
				r = std::to_chars(s, s + sizeof(s), m_data.m_value.m_float);

			return r.ec == std::errc{} ? std::string{ s, r.ptr } : "*****";
		}
	}

	std::unreachable();
}

// void const_item_handle::assign_value(const item_value &value)
// {
// 	assert(not m_row_handle.empty());
// 	m_row_handle.assign(m_item_ix, value, true);
// }

std::ostream &operator<<(std::ostream &os, const item_value &v)
{
	switch (v.type())
	{
		using enum item_value_type;

		case INT: os << v.m_data.m_value.m_integer; break;
		case FLOAT: os << v.m_data.m_value.m_float; break;
		case TEXT: os << v.m_data.sv(); break;
		case MISSING: os << '?'; break;
		case INAPPLICABLE: os << '.'; break;
		default: os.setstate(std::ios::failbit);
	}

	return os;
}

void item_value::cast_to_int()
{
	switch (type())
	{
		using enum item_value_type;

		case INT:
			break;

		case FLOAT:
			*this = std::rint(m_data.m_value.m_float);
			break;

		case TEXT:
		{
			auto s = sv();
			int64_t v;

			auto sp = s.data();
			if (*sp == '+')
				++sp;

			auto [ptr, ec] = cif::from_chars(sp, s.data() + s.size(), v);
			if (ec != std::errc{})
				throw std::system_error(std::make_error_code(ec), "attempt to cast value to integer failed");
			if (ptr != s.data() + s.size())
				throw std::runtime_error("attempt to cast value to integer failed, trailing data");

			*this = v;
			break;
		}

		default:
			break;
	}
}

void item_value::cast_to_float()
{
	switch (type())
	{
		using enum item_value_type;

		case INT:
			*this = static_cast<double>(m_data.m_value.m_integer);
			break;

		case FLOAT:
			break;

		case TEXT:
		{
			auto s = sv();
			double v;
			auto [ptr, ec] = cif::from_chars(s.data(), s.data() + s.size(), v);
			if (ec != std::errc{})
				throw std::system_error(std::make_error_code(ec), "attempt to cast value to integer failed");
			if (ptr != s.data() + s.size())
				throw std::runtime_error("attempt to cast value to integer failed, trailing data");
			*this = v;
			break;
		}

		default:
			break;
	}
}

void item_value::cast_to_string()
{
	*this = str();
}

} // namespace cif
