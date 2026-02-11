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

#include "cif++/row.hpp"

#include "cif++/category.hpp"
#include "cif++/item.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace cif
{

// item_value &row_handle::operator[](uint16_t item_ix)
// {
// 	return empty() or item_ix >= m_row->size() ? s_null_item : m_row->operator[](item_ix);
// }

// const item_value &row_handle::operator[](uint16_t item_ix) const
// {
// 	return empty() or item_ix >= m_row->size() ? s_null_item : m_row->operator[](item_ix);
// }

// item_value &row_handle::operator[](std::string_view item_name)
// {
// 	return operator[](get_item_ix(item_name));
// }

// const item_value &row_handle::operator[](std::string_view item_name) const
// {
// 	return operator[](get_item_ix(item_name));
// }

uint16_t row_handle::get_item_ix(std::string_view name) const
{
	if (not m_category)
		throw std::runtime_error("uninitialized row");

	return m_category->get_item_ix(name);
}

std::string_view row_handle::get_item_name(uint16_t ix) const
{
	if (not m_category)
		throw std::runtime_error("uninitialized row");

	return m_category->get_item_name(ix);
}

uint16_t const_row_handle::get_item_ix(std::string_view name) const
{
	if (not m_category)
		throw std::runtime_error("uninitialized row");

	return m_category->get_item_ix(name);
}

std::string_view const_row_handle::get_item_name(uint16_t ix) const
{
	if (not m_category)
		throw std::runtime_error("uninitialized row");

	return m_category->get_item_name(ix);
}

// --------------------------------------------------------------------

void row_handle::assign(uint16_t item, item_value value, bool updateLinked, bool validate)
{
	if (not m_category)
		throw std::runtime_error("uninitialized row");

	m_category->update_value(m_row, item, std::move(value), updateLinked, validate);
}

uint16_t row_handle::add_item(std::string_view name)
{
	if (not m_category)
		throw std::runtime_error("uninitialized row");

	return m_category->add_item(name);
}

// --------------------------------------------------------------------

row_initializer::row_initializer(const_row_handle rh)
{
	if (not rh.m_category)
		throw std::runtime_error("uninitialized row");

	assert(rh.m_row);

	auto r = rh.get_row();
	auto &cat = *rh.m_category;

	for (uint16_t ix = 0; std::cmp_less(ix, r->size()); ++ix)
	{
		auto &i = r->operator[](ix);
		if (not i)
			continue;
		emplace_back(cat.get_item_name(ix), i);
	}
}

void row_initializer::set_value(std::string name, item_value value)
{
	for (auto &i : *this)
	{
		if (i.name() == name)
		{
			i.value(std::move(value));
			return;
		}
	}

	emplace_back(std::move(name), std::move(value));
}

void row_initializer::set_value_if_empty(std::string name, item_value value)
{
	if (std::ranges::find_if(*this, [name](auto &i)
			{ return i.name() == name; }) == end())
		emplace_back(std::move(name), std::move(value));
}

} // namespace cif
