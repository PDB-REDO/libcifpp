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

#include "cif++/category.hpp"

namespace cif
{

void row_handle::assign(uint16_t item, std::string_view value, bool updateLinked, bool validate)
{
	if (not m_category)
		throw std::runtime_error("uninitialized row");

	m_category->update_value(m_row, item, value, updateLinked, validate);
}

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

uint16_t row_handle::add_item(std::string_view name)
{
	if (not m_category)
		throw std::runtime_error("uninitialized row");

	return m_category->add_item(name);
}

void row_handle::swap(uint16_t item, row_handle &b)
{
	if (not m_category)
		throw std::runtime_error("uninitialized row");

	m_category->swap_item(item, *this, b);
}

// --------------------------------------------------------------------

row_initializer::row_initializer(row_handle rh)
{
	if (not rh.m_category)
		throw std::runtime_error("uninitialized row");

	assert(rh.m_row);

	row *r = rh.get_row();
	auto &cat = *rh.m_category;

	for (uint16_t ix = 0; ix < r->size(); ++ix)
	{
		auto &i = r->operator[](ix);
		if (not i)
			continue;
		emplace_back(cat.get_item_name(ix), i.text());
	}
}

void row_initializer::set_value(std::string_view name, std::string_view value)
{
	for (auto &i : *this)
	{
		if (i.name() == name)
		{
			i.value(value);
			return;
		}
	}

	emplace_back(name, value);
}

void row_initializer::set_value_if_empty(std::string_view name, std::string_view value)
{
	if (find_if(begin(), end(), [name](auto &i) { return i.name() == name; }) == end())
		emplace_back(name, value);
}

} // namespace cif