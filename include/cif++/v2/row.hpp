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

template<typename> class row_handle;

namespace detail
{

	// some helper classes to help create tuple result types
	template <
		typename Category,
		typename... C
	>
	struct get_row_result
	{
		using category_type = Category;
		using row_type = category_type::row;
		using row_handle_type = row_handle<category_type>;
		using item_handle_type = item_handle<row_type>;

		static constexpr size_t N = sizeof...(C);

		get_row_result(const row_handle_type &r, std::array<size_t, N> &&columns)
			: m_row(r)
			, m_columns(std::move(columns))
		{
		}

		const item_handle_type operator[](size_t ix) const
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
			return std::tuple<Ts...>{m_row[m_columns[Is]].template as<Ts>()...};
		}

		const row_handle_type &m_row;
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

}

template <typename... Ts>
auto tie(Ts &...v)
{
	return detail::tie_wrap<Ts &...>(std::forward<Ts &>(v)...);
}

// --------------------------------------------------------------------
/// \brief row_handle is the way to access data in rows

template<typename Category>
class row_handle
{
  public:

	using category_type = Category;
	using row_type = std::conditional_t<std::is_const_v<category_type>, const typename category_type::row, typename category_type::row>;

	row_handle(category_type &cat, row_type &row)
		: m_cat(cat)
		, m_row(row) {}

	item_handle<row_type> operator[](uint32_t column_ix)
	{
		return item_handle<row_type>(column_ix, m_row);
	}

	const item_handle<const row_type> operator[](uint32_t column_ix) const
	{
		return item_handle<const row_type>(column_ix, m_row);
	}

	item_handle<row_type> operator[](std::string_view column_name)
	{
		return item_handle<row_type>(get_column_ix(column_name), m_row);
	}

	const item_handle<const row_type> operator[](std::string_view column_name) const
	{
		return item_handle<const row_type>(get_column_ix(column_name), m_row);
	}

	template <typename... Ts, size_t N>
	std::tuple<Ts...> get(char const *const (&columns)[N]) const
	{
		static_assert(sizeof...(Ts) == N, "Number of columns should be equal to number of types to return");

		std::array<size_t, N> cix;
		for (size_t i = 0; i < N; ++i)
			cix[i] = get_column_ix(columns[i]);
		return detail::get_row_result<category_type, Ts...>(*this, std::move(cix));
	}

	template <typename... C>
	auto get(C... columns) const
	{
		return detail::get_row_result<category_type, C...>(*this, {get_column_ix(columns)...});
	}

  private:

	uint32_t get_column_ix(std::string_view name) const
	{
		return m_cat.get_column_ix(name);
	}

	category_type &m_cat;
	row_type &m_row;
};


}