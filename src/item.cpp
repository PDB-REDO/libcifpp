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

#include <cassert>
#include <compare>
#include <string_view>

namespace cif
{

bool item_handle::empty() const
{
	return m_item_ix >= m_row.size() or m_row[m_item_ix].empty();
}

item_value &item_handle::value()
{
	return m_row.operator[](m_item_ix);
}

const item_value &item_handle::value() const
{
	return m_row.operator[](m_item_ix);
}

void item_handle::set(item_value value, bool updateLinked)
{
	row_handle rh{ m_category, m_row };
	rh.assign(m_item_ix, std::move(value), updateLinked);
}

bool const_item_handle::empty() const
{
	return m_item_ix >= m_row.size() or m_row[m_item_ix].empty();
}

const item_value &const_item_handle::value() const
{
	return m_row.operator[](m_item_ix);
}

int item_value::compare(const item_value &b, bool ignore_case) const noexcept
{
	int d = static_cast<int>(m_data.m_type) - static_cast<int>(b.m_data.m_type);

	if (d == 0)
	{
		switch (m_data.m_type)
		{
			case cif::item_value_type::INT:
				d = m_data.m_value.m_integer - b.m_data.m_value.m_integer;
				break;
			case cif::item_value_type::FLOAT:
			{
				auto dp = (m_data.m_value.m_float <=> b.m_data.m_value.m_float);
				if (dp == std::partial_ordering::less)
					d = -1;
				else if (dp == std::partial_ordering::greater)
					d = 1;
				break;
			}
			case cif::item_value_type::TEXT:
				d = m_data.sv().compare(b.m_data.sv());
				break;
			default:;
		}
	}
	else if (is_number() and b.is_number())
	{
		std::partial_ordering dp = std::partial_ordering::equivalent;

		if (is_number_float())
			dp = m_data.m_value.m_float <=> b.m_data.m_value.m_integer;
		else/*  if (is_number_int()) */
			dp = m_data.m_value.m_integer <=> b.m_data.m_value.m_float;

		if (dp == std::partial_ordering::less)
			d = -1;
		else if (dp == std::partial_ordering::greater)
			d = 1;
		else
		 	d = 0;
	}

	return d;
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
		case cif::item_value_type::INT:
			os << v.m_data.m_value.m_integer;
			break;
		case cif::item_value_type::FLOAT:
			os << v.m_data.m_value.m_float;
			break;
		case cif::item_value_type::TEXT:
			os << v.m_data.sv();
			break;
		case cif::item_value_type::MISSING:
			os << '?';
			break;
		case cif::item_value_type::INAPPLICABLE:
			os << '.';
			break;
	}

	return os;
}

} // namespace cif
