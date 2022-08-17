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

#include <cif++/item.hpp>

namespace cif
{

namespace detail
{

	// some helper classes to help create tuple result types
	template <typename... C>
	struct get_row_result
	{
		static constexpr size_t N = sizeof...(C);

		get_row_result(const row_handle &r, std::array<size_t, N> &&columns)
			: m_row(r)
			, m_columns(std::move(columns))
		{
		}

		const item_handle operator[](size_t ix) const
		{
			return m_row[m_columns[ix]];
		}

		template <typename... Ts, std::enable_if_t<N == sizeof...(Ts), int> = 0>
		operator std::tuple<Ts...>() const
		{
			return get<Ts...>(std::index_sequence_for<Ts...>{});
		}

		template <typename... Ts, std::size_t... Is>
		std::tuple<Ts...> get(std::index_sequence<Is...>) const
		{
			return std::tuple<Ts...>{ m_row[m_columns[Is]].template as<Ts>()... };
		}

		const row_handle &m_row;
		std::array<size_t, N> m_columns;
	};

	// we want to be able to tie some variables to a get_row_result, for this we use tiewraps
	template <typename... Ts>
	struct tie_wrap
	{
		tie_wrap(Ts... args)
			: m_value(args...)
		{
		}

		template <typename RR>
		void operator=(const RR &&rr)
		{
			// get_row_result will do the conversion, but only if the types
			// are compatible. That means the number of parameters to the get()
			// of the row should be equal to the number of items in the tuple
			// you are trying to tie.

			using RType = std::tuple<typename std::remove_reference<Ts>::type...>;

			m_value = static_cast<RType>(rr);
		}

		std::tuple<Ts...> m_value;
	};

} // namespace detail

template <typename... Ts>
auto tie(Ts &...v)
{
	return detail::tie_wrap<Ts &...>(std::forward<Ts &>(v)...);
}

// --------------------------------------------------------------------
/// \brief the row class, this one is not directly accessible from the outside

class row
{
  public:
	row() = default;

  private:
	friend class item_handle;
	friend class category_index;
	friend class row_initializer;

	template <typename, typename...>
	friend class iterator_impl;

	friend class category;

	void append(item_value *iv)
	{
		if (m_head == nullptr)
			m_head = m_tail = iv;
		else
			m_tail = m_tail->m_next = iv;
	}

	row *m_next = nullptr;
	item_value *m_head = nullptr, *m_tail = nullptr;
};

// --------------------------------------------------------------------
/// \brief row_handle is the way to access data stored in rows

class row_handle
{
  public:
	friend class item_handle;
	friend class category;
	friend class category_index;
	friend class row_initializer;

	row_handle() = default;

	row_handle(const row_handle &) = default;
	row_handle(row_handle &&) = default;

	row_handle &operator=(const row_handle &) = default;
	row_handle &operator=(row_handle &&) = default;

	row_handle(const category &cat, const row &r)
		: m_category(const_cast<category *>(&cat))
		, m_row(const_cast<row *>(&r))
	{
	}

	const category &get_category() const
	{
		return *m_category;
	}

	explicit operator bool() const
	{
		return m_category != nullptr and m_row != nullptr;
	}

	item_handle operator[](uint32_t column_ix)
	{
		return item_handle(column_ix, *this);
	}

	const item_handle operator[](uint32_t column_ix) const
	{
		return item_handle(column_ix, const_cast<row_handle &>(*this));
	}

	item_handle operator[](std::string_view column_name)
	{
		return item_handle(add_column(column_name), *this);
	}

	const item_handle operator[](std::string_view column_name) const
	{
		return item_handle(get_column_ix(column_name), const_cast<row_handle &>(*this));
	}

	// template <typename... Ts, size_t N>
	// std::tuple<Ts...> get(char const *const (&columns)[N]) const
	// {
	// 	static_assert(sizeof...(Ts) == N, "Number of columns should be equal to number of types to return");

	// 	std::array<size_t, N> cix;
	// 	for (size_t i = 0; i < N; ++i)
	// 		cix[i] = get_column_ix(columns[i]);
	// 	return detail::get_row_result<Ts...>(*this, std::move(cix));
	// }

	template <typename... C>
	auto get(C... columns) const
	{
		return detail::get_row_result<C...>(*this, { get_column_ix(columns)... });
	}

	template <typename... Ts, typename... C, std::enable_if_t<sizeof...(Ts) == sizeof...(C), int> = 0>
	std::tuple<Ts...> get(C... columns) const
	{
		return detail::get_row_result<Ts...>(*this, { get_column_ix(columns)... });
	}

	template<typename T>
	T get(const char *column)
	{
		return operator[](get_column_ix(column)).template as<T>();
	}

	void assign(const std::vector<item> &values)
	{
		for (auto &value : values)
			assign(value, true);
	}

	void assign(std::string_view name, std::string_view value, bool updateLinked, bool validate = true)
	{
		assign(add_column(name), value, updateLinked, validate);
	}

	void assign(size_t column, std::string_view value, bool updateLinked, bool validate = true);

	bool operator==(const row_handle &rhs) const { return m_category == rhs.m_category and m_row == rhs.m_row; }
	bool operator!=(const row_handle &rhs) const { return m_category != rhs.m_category or m_row != rhs.m_row; }

  private:
	uint16_t get_column_ix(std::string_view name) const;
	std::string_view get_column_name(uint16_t ix) const;

	uint16_t add_column(std::string_view name);

	operator row *()
	{
		return m_row;
	}

	void assign(const item &i, bool updateLinked)
	{
		assign(i.name(), i.value(), updateLinked);
	}

	category *m_category = nullptr;
	row *m_row = nullptr;
};

// --------------------------------------------------------------------

class row_initializer
{
  public:
	friend class category;

	row_initializer() = default;
	row_initializer(const row_initializer &) = default;
	row_initializer(row_initializer &&) = default;
	row_initializer &operator=(const row_initializer &) = default;
	row_initializer &operator=(row_initializer &&) = default;

	row_initializer(std::initializer_list<item> items)
		: m_items(items)
	{
	}

	template <typename ItemIter, std::enable_if_t<std::is_same_v<typename ItemIter::value_type, item>, int> = 0>
	row_initializer(ItemIter b, ItemIter e)
		: m_items(b, e)
	{
	}

	row_initializer(row_handle rh);

  private:
	std::vector<item> m_items;
};

} // namespace cif