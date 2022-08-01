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

#include "row.hpp"

namespace cif::v2
{

// --------------------------------------------------------------------

template <typename CategoryType, typename... Ts>
class iterator_impl
{
  public:
	template <typename, typename...>
	friend class iterator_impl;

	static constexpr size_t N = sizeof...(Ts);

	using category_type = CategoryType;
	using row_type = typename category_type::row;

	// using row_type = std::conditional_t<std::is_const_v<category_type>, const category_type::row, category_type::row>;
	// using row_impl_type = std::conditional_t<std::is_const_v<category_type>, typename const category_type::row, typename category_type::row>;

	using iterator_category = std::forward_iterator_tag;
	using value_type = std::conditional_t<N == 0, row_handle<category_type>, std::tuple<Ts...>>;
	using difference_type = std::ptrdiff_t;
	using pointer = value_type *;
	using reference = std::conditional_t<N == 0, row_handle<category_type>, value_type>;

	friend class Category;

	// // default constructor, equal to end()
	// iterator_impl() {}

	iterator_impl(const iterator_impl &rhs)
		: m_category(rhs.m_category)
		, m_current(rhs.m_current)
		, m_value(rhs.m_value)
		, m_column_ix(rhs.m_column_ix)
	{
	}

	iterator_impl(category_type &cat, row_type *current)
		: m_category(&cat)
		, m_current(current)
		, m_value(cat, *current)
	{
		static_assert(N == 0, "Only valid if this is a row iterator, not a row<xxx> iterator");
	}

	// iterator_impl(ItemRow *data)
	// 	: m_current(data)
	// {
	// 	static_assert(N == 0, "Only valid if this is a row iterator, not a row<xxx> iterator");
	// }

	// iterator_impl(ItemRow *data, const std::array<size_t, N> &cix)
	// 	: m_current(data)
	// 	, m_column_ix(cix)
	// {
	// }

	template <typename IRowType>
	iterator_impl(iterator_impl<IRowType, Ts...> &rhs)
		: m_category(rhs.m_category)
		, m_current(rhs.m_current)
		, m_column_ix(rhs.m_column_ix)
	{
		if constexpr (N > 0)
			m_value = get(m_current, std::make_index_sequence<N>());
	}

	template <typename IRowType>
	iterator_impl(const iterator_impl<IRowType> &rhs, const std::array<size_t, N> &cix)
		: m_category(rhs.m_category)
		, m_current(rhs.m_current)
		, m_column_ix(cix)
	{
		if constexpr (N > 0)
			m_value = get(m_current, std::make_index_sequence<N>());
	}

	iterator_impl &operator=(const iterator_impl &i)
	{
		m_category = i.m_category;
		m_current = i.m_current;
		if constexpr (N != 0)
		{
			m_column_ix = i.m_column_ix;
			m_value = i.m_value;
		}
		return *this;
	}

	virtual ~iterator_impl() = default;

	reference operator*()
	{
		if constexpr (N == 0)
			return { *m_category, *m_current };
		else
			return m_value;
	}

	pointer operator->()
	{
		if constexpr (N == 0)
			return &m_current;
		else
			return &m_value;
	}

	row_type row() const
	{
		return m_current;
	}

	iterator_impl &operator++()
	{
		if (m_current != nullptr)
			m_current = m_current->m_next;

		if constexpr (N != 0)
			m_value = get(m_current, std::make_index_sequence<N>());

		return *this;
	}

	iterator_impl operator++(int)
	{
		iterator_impl result(*this);
		this->operator++();
		return result;
	}

	bool operator==(const iterator_impl &rhs) const { return m_current == rhs.m_current; }
	bool operator!=(const iterator_impl &rhs) const { return m_current != rhs.m_current; }

	template <typename IRowType, typename... ITs>
	bool operator==(const iterator_impl<IRowType, ITs...> &rhs) const
	{
		return m_current == rhs.m_current;
	}

	template <typename IRowType, typename... ITs>
	bool operator!=(const iterator_impl<IRowType, ITs...> &rhs) const
	{
		return m_current != rhs.m_current;
	}

  private:
	template <std::size_t... Is>
	std::tuple<Ts...> get(row_type row, std::index_sequence<Is...>) const
	{
		if (row)
			return std::tuple<Ts...>{row[m_column_ix[Is]].template as<Ts>()...};
		return {};
	}

	category_type *m_category;
	row_type *m_current;
	value_type m_value;
	std::array<size_t, N> m_column_ix;
};


} // namespace cif::v2