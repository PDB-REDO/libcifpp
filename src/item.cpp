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

#include <cassert>

namespace cif
{

const item_handle item_handle::s_null_item;
row_handle s_null_row_handle;

item_handle::item_handle()
	: m_item_ix(std::numeric_limits<uint16_t>::max())
	, m_row_handle(s_null_row_handle)
{
}

std::string_view item_handle::text() const
{
	if (not m_row_handle.empty())
	{
		auto iv = m_row_handle.m_row->get(m_item_ix);
		if (iv != nullptr)
			return iv->text();
	}

	return {};
}

void item_handle::assign_value(std::string_view value)
{
	assert(not m_row_handle.empty());
	m_row_handle.assign(m_item_ix, value, true);
}

void item_handle::swap(item_handle &b)
{
	assert(m_item_ix == b.m_item_ix);
	// assert(&m_row_handle.m_category == &b.m_row_handle.m_category);
	m_row_handle.swap(m_item_ix, b.m_row_handle);
}

}
