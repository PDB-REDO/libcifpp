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

#include "item.hpp"

namespace cif::v2
{

// --------------------------------------------------------------------
/// \brief row_handle is the way to access data in rows

template<typename Category>
class row_handle
{
  public:

	using category_type = Category;

	row_handle(Category &cat)
		: m_cat(cat) {}

	// item_handle<row_t> operator[](uint32_t column_ix)
	// {
	// 	return item_handle<row_t>(column_ix, *this);
	// }

	// const item_handle<const row_t> operator[](uint32_t column_ix) const
	// {
	// 	return item_handle<const row_t>(column_ix, *this);
	// }

	// item_handle<row_t> operator[](std::string_view column_name)
	// {
	// 	return item_handle<row_t>(column_name, get_column_ix(column_name), *this);
	// }

	// const item_handle<const row_t> operator[](std::string_view column_name) const
	// {
	// 	return item_handle<const row_t>(column_name, get_column_ix(column_name), *this);
	// }



  private:

	uint32_t get_column_ix(std::string_view name) const
	{
		return m_cat.get_column_ix(name);
	}



	category_type &m_cat;
};


}