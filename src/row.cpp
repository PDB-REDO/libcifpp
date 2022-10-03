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

#include <cif++/category.hpp>

namespace cif
{

void row_handle::assign(size_t column, std::string_view value, bool updateLinked, bool validate)
{
	assert(m_category);
	m_category->update_value(m_row, column, value, updateLinked, validate);
}

uint16_t row_handle::get_column_ix(std::string_view name) const
{
	assert(m_category);
	return m_category->get_column_ix(name);
}

std::string_view row_handle::get_column_name(uint16_t ix) const
{
	assert(m_category);
	return m_category->get_column_name(ix);
}

uint16_t row_handle::add_column(std::string_view name)
{
	assert(m_category);
	return m_category->add_column(name);
}

void row_handle::swap(size_t column, row_handle &b)
{
	m_category->swap_item(column, *this, b);
}

// --------------------------------------------------------------------

row_initializer::row_initializer(row_handle rh)
{
	assert(rh.m_category);
	assert(rh.m_row);

	row *r = rh;
	auto &cat = *rh.m_category;

	for (auto i = r->m_head; i != nullptr; i = i->m_next)
		emplace_back(cat.get_column_name(i->m_column_ix), i->text());
}

} // namespace cif