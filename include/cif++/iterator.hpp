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

#include "cif++/condition.hpp"
#include "cif++/row.hpp"

#include <array>
#include <cstdint>
#include <numeric>
#include <type_traits>

/**
 * @file iterator.hpp
 *
 * This file contains several implementations of generic iterators.
 *
 * Using partial specialization we can have implementation for
 * iterators that return row_handles, a single value or tuples of
 * multiple values.
 *
 */

namespace cif
{

class category;

// --------------------------------------------------------------------

/**
 * @brief Implementation of an iterator that can return
 * multiple values in a tuple. Of course, that tuple can
 * then be used in structured binding to receive the values
 * in a for loop e.g.
 *
 * @tparam Category The category for this iterator
 * @tparam Ts The types this iterator can be dereferenced to
 */
template <bool Const, typename... Ts>
class iterator_impl_base
{
  public:
	/** @cond */
	template <bool, typename...>
	friend class iterator_impl_base;

	friend class category;
	/** @endcond */

	/** variable that contains the number of elements in the tuple */
	static constexpr std::size_t N = sizeof...(Ts);

	/** @cond */
	using tuple_type = std::tuple<Ts...>;

	using row_handle_type = std::conditional_t<Const, const_row_handle, row_handle>;

	using iterator_category = std::forward_iterator_tag;
	using value_type = std::conditional_t<Const, const tuple_type, tuple_type>;
	using difference_type = std::ptrdiff_t;
	using pointer = value_type *;
	using reference = value_type &;

	iterator_impl_base() = default;

	iterator_impl_base(const iterator_impl_base &rhs) = default;
	iterator_impl_base(iterator_impl_base &&rhs) = default;

	template <bool C, typename... T2s>
	iterator_impl_base(const iterator_impl_base<C, T2s...> &rhs)
		: m_current(rhs.m_current)
		, m_value(rhs.m_value)
		, m_item_ix(rhs.m_item_ix)
	{
	}

	template <bool C>
	iterator_impl_base(iterator_impl_base<C, Ts...> &rhs)
		: m_current(rhs.m_current)
		, m_value(rhs.m_value)
		, m_item_ix(rhs.m_item_ix)
	{
		m_value = get(std::make_index_sequence<N>());
	}

	template <bool C>
	iterator_impl_base(const iterator_impl_base<C> &rhs, const std::array<uint16_t, N> &cix)
		: m_current(rhs.m_current)
		, m_item_ix(cix)
	{
		m_value = get(std::make_index_sequence<N>());
	}

	iterator_impl_base &operator=(iterator_impl_base i)
	{
		std::swap(m_current, i.m_current);
		std::swap(m_item_ix, i.m_item_ix);
		std::swap(m_value, i.m_value);
		return *this;
	}

	virtual ~iterator_impl_base() = default;

	auto operator*()
	{
		return m_value;
	}

	auto operator*() const
	{
		return m_value;
	}

	auto operator->()
	{
		return &m_value;
	}

	auto operator->() const
	{
		return &m_value;
	}

	operator const_row_handle() const
	{
		return m_current;
	}

	operator row_handle_type()
	{
		return m_current;
	}

	iterator_impl_base &operator++()
	{
		if (m_current)
			m_current.m_row = m_current.m_row->m_next;

		m_value = get(std::make_index_sequence<N>());

		return *this;
	}

	iterator_impl_base operator++(int)
	{
		iterator_impl_base result(*this);
		this->operator++();
		return result;
	}

	bool operator==(const iterator_impl_base &rhs) const { return m_current == rhs.m_current; }
	bool operator!=(const iterator_impl_base &rhs) const { return m_current != rhs.m_current; }

	template <bool C, typename... ITs>
	bool operator==(const iterator_impl_base<C, ITs...> &rhs) const
	{
		return m_current == rhs.m_current;
	}

	template <bool C, typename... ITs>
	bool operator!=(const iterator_impl_base<C, ITs...> &rhs) const
	{
		return m_current != rhs.m_current;
	}

	/** @endcond */

  private:
	template <std::size_t... Is>
	[[nodiscard]] tuple_type get(std::index_sequence<Is...>) const
	{
		return m_current ? tuple_type{ m_current[m_item_ix[Is]].template as<Ts>()... } : tuple_type{};
	}

	row_handle_type m_current;
	tuple_type m_value;
	std::array<uint16_t, N> m_item_ix;
};

/**
 * @brief Implementation of an iterator that returns
 * only row_handles
 *
 * @tparam Category The category for this iterator
 */
template <bool Const>
class iterator_impl_base<Const>
{
  public:
	/** @cond */

	template <bool, typename...>
	friend class iterator_impl_base;

	friend class category;

	using category_type = std::conditional_t<Const, const category, category>;
	using row_type = std::conditional_t<Const, const row, row>;
	using row_handle_type = std::conditional_t<Const, const_row_handle, row_handle>;

	using iterator_category = std::forward_iterator_tag;

	using value_type = std::conditional_t<Const, const_row_handle, row_handle>;
	using difference_type = std::ptrdiff_t;
	using pointer = value_type *;
	using reference = value_type &;

	iterator_impl_base() = default;

	iterator_impl_base(const iterator_impl_base &rhs) = default;
	iterator_impl_base(iterator_impl_base &&rhs) = default;

	template <bool C>
	iterator_impl_base(const iterator_impl_base<C> &rhs)
		: m_current(rhs.m_current)
	{
	}

	iterator_impl_base(category_type &cat, row_type *current)
		: m_current(cat, *current)
	{
	}

	template <bool C>
	iterator_impl_base(const iterator_impl_base<C> &rhs, const std::array<uint16_t, 0> &)
		: m_current(rhs.m_current)
	{
	}

	iterator_impl_base &operator=(iterator_impl_base i)
	{
		std::swap(m_current, i.m_current);
		return *this;
	}

	virtual ~iterator_impl_base() = default;

	auto operator*()
	{
		return m_current;
	}

	auto operator*() const
	{
		return m_current;
	}

	auto operator->()
	{
		return &m_current;
	}

	auto operator->() const
	{
		return &m_current;
	}

	operator const_row_handle() const
	{
		return m_current;
	}

	operator row_handle_type()
	{
		return m_current;
	}

	[[nodiscard]] int64_t row_id() const
	{
		return reinterpret_cast<int64_t>(m_current.m_row);
	}

	iterator_impl_base &operator++()
	{
		if (m_current)
			m_current.m_row = m_current.m_row->m_next;

		return *this;
	}

	iterator_impl_base operator++(int)
	{
		iterator_impl_base result(*this);
		this->operator++();
		return result;
	}

	bool operator==(const iterator_impl_base &rhs) const { return m_current == rhs.m_current; }
	bool operator!=(const iterator_impl_base &rhs) const { return m_current != rhs.m_current; }

	template <bool C, typename... ITs>
	bool operator==(const iterator_impl_base<C, ITs...> &rhs) const
	{
		return m_current == rhs.m_current;
	}

	template <bool C, typename... ITs>
	bool operator!=(const iterator_impl_base<C, ITs...> &rhs) const
	{
		return m_current != rhs.m_current;
	}

	/** @endcond */

  private:
	row_handle_type m_current;
};

/**
 * @brief Implementation of an iterator that can return
 * a single value.
 *
 * @tparam Category The category for this iterator
 * @tparam T The type this iterator can be dereferenced to
 */

template <bool Const, typename T>
class iterator_impl_base<Const, T>
{
  public:
	/** @cond */
	template <bool, typename...>
	friend class iterator_impl_base;

	friend class category;

	using category_type = std::conditional_t<Const, const category, category>;
	using row_handle_type = std::conditional_t<Const, const_row_handle, row_handle>;

	using iterator_category = std::forward_iterator_tag;
	using value_type = T;
	using difference_type = std::ptrdiff_t;
	using pointer = value_type *;
	using reference = value_type &;

	iterator_impl_base() = default;

	iterator_impl_base(const iterator_impl_base &rhs) = default;
	iterator_impl_base(iterator_impl_base &&rhs) = default;

	template <bool C, typename T2>
	iterator_impl_base(const iterator_impl_base<C, T2> &rhs)
		: m_current(rhs.m_current)
		, m_value(rhs.m_value)
		, m_item_ix(rhs.m_item_ix)
	{
	}

	template <bool C>
	iterator_impl_base(iterator_impl_base<C, T> &rhs)
		: m_current(rhs.m_current)
		, m_value(rhs.m_value)
		, m_item_ix(rhs.m_item_ix)
	{
		m_value = get();
	}

	template <bool C>
	iterator_impl_base(const iterator_impl_base<C> &rhs, const std::array<uint16_t, 1> &cix)
		: m_current(rhs.m_current)
		, m_item_ix(cix[0])
	{
		m_value = get();
	}

	iterator_impl_base &operator=(iterator_impl_base i)
	{
		std::swap(m_current, i.m_current);
		std::swap(m_item_ix, i.m_item_ix);
		std::swap(m_value, i.m_value);
		return *this;
	}

	virtual ~iterator_impl_base() = default;

	auto operator*()
	{
		return m_value;
	}

	auto operator*() const
	{
		return m_value;
	}

	auto operator->()
	{
		return &m_value;
	}

	auto operator->() const
	{
		return &m_value;
	}

	operator const_row_handle() const
	{
		return m_current;
	}

	operator row_handle_type()
	{
		return m_current;
	}

	iterator_impl_base &operator++()
	{
		if (m_current)
			m_current.m_row = m_current.m_row->m_next;

		m_value = get();

		return *this;
	}

	iterator_impl_base operator++(int)
	{
		iterator_impl_base result(*this);
		this->operator++();
		return result;
	}

	bool operator==(const iterator_impl_base &rhs) const { return m_current == rhs.m_current; }
	bool operator!=(const iterator_impl_base &rhs) const { return m_current != rhs.m_current; }

	template <bool C, typename... ITs>
	bool operator==(const iterator_impl_base<C, ITs...> &rhs) const
	{
		return m_current == rhs.m_current;
	}

	template <bool C, typename... ITs>
	bool operator!=(const iterator_impl_base<C, ITs...> &rhs) const
	{
		return m_current != rhs.m_current;
	}

	/** @endcond */

  private:
	[[nodiscard]] value_type get() const
	{
		return m_current ?  m_current[m_item_ix].template get<value_type>() : value_type{};
	}

	row_handle_type m_current;
	value_type m_value;
	uint16_t m_item_ix;
};

// --------------------------------------------------------------------

template<typename ... Ts>
using iterator_impl = iterator_impl_base<false, Ts...>;

template<typename ... Ts>
using const_iterator_impl = iterator_impl_base<true, Ts...>;


// --------------------------------------------------------------------
// iterator proxy

/**
 * @brief An iterator_proxy is used as a result type for methods that
 * return a range of values you want to iterate over.
 *
 * E.g. the class cif::category contains the method cif::category::rows()
 * that returns an iterator_proxy that allows you to iterate over
 * all the rows in the category.
 *
 * @tparam Category The category for the iterators
 * @tparam Ts The types the iterators return. See class: iterator
 */

template <bool Const, typename... Ts>
class iterator_proxy_base
{
  public:
	/** @cond */
	static constexpr const std::size_t N = sizeof...(Ts);

	using category_type = std::conditional_t<Const, const category, category>;

	using iterator = iterator_impl_base<Const, Ts...>;
	using row_iterator = iterator_impl_base<Const>;

	iterator_proxy_base(category_type &cat, row_iterator pos, char const *const items[N]);
	iterator_proxy_base(category_type &cat, row_iterator pos, std::initializer_list<char const *> items); // NOLINT(modernize-pass-by-value)

	iterator_proxy_base(iterator_proxy_base &&p);
	iterator_proxy_base &operator=(iterator_proxy_base &&p);

	iterator_proxy_base(const iterator_proxy_base &) = delete;
	iterator_proxy_base &operator=(const iterator_proxy_base &) = delete;
	/** @endcond */

	[[nodiscard]] iterator begin() const { return iterator(m_begin, m_item_ix); } ///< Return the iterator pointing to the first row
	[[nodiscard]] iterator end() const { return iterator(m_end, m_item_ix); }     ///< Return the iterator pointing past the last row

	[[nodiscard]] bool empty() const { return m_begin == m_end; }               ///< Return true if the range is empty
	explicit operator bool() const { return not empty(); }        ///< Easy way to detect if the range is empty
	[[nodiscard]] std::size_t size() const { return std::distance(begin(), end()); } ///< Return size of the range

	// row front() { return *begin(); }
	// row back() { return *(std::prev(end())); }

	[[nodiscard]] category_type &get_category() const { return *m_category; } ///< Return the category the iterator belong to

	/** swap */
	void swap(iterator_proxy_base &rhs)
	{
		std::swap(m_category, rhs.m_category);
		std::swap(m_begin, rhs.m_begin);
		std::swap(m_end, rhs.m_end);
		std::swap(m_item_ix, rhs.m_item_ix);
	}

  protected:
	iterator_proxy_base(category_type &cat);

  private:
	category_type *m_category;
	row_iterator m_begin, m_end;
	std::array<uint16_t, N> m_item_ix;
};

// --------------------------------------------------------------------

template <typename... Ts>
using iterator_proxy = iterator_proxy_base<false, Ts...>;

template <typename... Ts>
using const_iterator_proxy = iterator_proxy_base<true, Ts...>;

// --------------------------------------------------------------------
// conditional iterator proxy

/**
 * @brief A conditional iterator proxy is similar to an iterator_proxy
 * in that it can be used to return a range of rows you can iterate over.
 * In the case of an conditional_iterator_proxy a cif::condition is used
 * to filter out only those rows that match the condition.
 *
 * @tparam category_type The category the iterators belong to
 * @tparam Ts The types to which the iterators can be dereferenced
 */
template <bool Const, typename... Ts>
class conditional_iterator_proxy_base
{
  public:
	/** @cond */
	static constexpr const std::size_t N = sizeof...(Ts);

	using category_type = std::conditional_t<Const, const category, category>;
	using base_iterator = iterator_impl_base<Const, Ts...>;
	using value_type = typename base_iterator::value_type;
	using row_iterator = iterator_impl_base<Const>;

	class conditional_iterator_impl
	{
	  public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = conditional_iterator_proxy_base::value_type;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type *;
		using reference = value_type;

		conditional_iterator_impl() = default;
		conditional_iterator_impl(category_type &cat, row_iterator pos, const condition &cond, const std::array<uint16_t, N> &cix);
		conditional_iterator_impl(const conditional_iterator_impl &i) = default;
		conditional_iterator_impl &operator=(const conditional_iterator_impl &i) = default;

		virtual ~conditional_iterator_impl() = default;

		auto operator*()
		{
			return *m_begin;
		}

		auto operator*() const
		{
			return *m_begin;
		}

		auto operator->()
		{
			m_current = *m_begin;
			return &m_current;
		}

		auto operator->() const
		{
			m_current = *m_begin;
			return &m_current;
		}

		conditional_iterator_impl &operator++()
		{
			while (m_begin != m_end)
			{
				if (++m_begin == m_end)
					break;
				
				if (m_condition->operator()(m_begin))
					break;
			}

			return *this;
		}

		conditional_iterator_impl operator++(int)
		{
			conditional_iterator_impl result(*this);
			this->operator++();
			return result;
		}

		bool operator==(const conditional_iterator_impl &rhs) const { return m_begin == rhs.m_begin; }
		bool operator!=(const conditional_iterator_impl &rhs) const { return m_begin != rhs.m_begin; }

		bool operator==(const row_iterator &rhs) const { return m_begin == rhs; }
		bool operator!=(const row_iterator &rhs) const { return m_begin != rhs; }

		template <bool C, typename... ITs>
		bool operator==(const iterator_impl_base<C, ITs...> &rhs) const { return m_begin == rhs; }

		template <bool C, typename... ITs>
		bool operator!=(const iterator_impl_base<C, ITs...> &rhs) const { return m_begin != rhs; }

	  private:
		category_type *m_cat = nullptr;
		base_iterator m_begin, m_end;
		std::remove_cv_t<value_type> m_current;
		const condition *m_condition;
	};

	using iterator = conditional_iterator_impl;
	using reference = typename iterator::reference;

	template <typename... Ns>
	conditional_iterator_proxy_base(category_type &cat, row_iterator pos, condition &&cond, Ns... names); // NOLINT(modernize-pass-by-value)

	conditional_iterator_proxy_base(conditional_iterator_proxy_base &&p)
	{
		swap(*this, p);
	}

	conditional_iterator_proxy_base &operator=(conditional_iterator_proxy_base &&p)
	{
		swap(*this, p);
		return *this;
	}

	conditional_iterator_proxy_base(const conditional_iterator_proxy_base &) = delete;
	conditional_iterator_proxy_base &operator=(const conditional_iterator_proxy_base &) = delete;

	/** @endcond */

	[[nodiscard]] iterator begin() const; ///< Return the iterator pointing to the first row
	[[nodiscard]] iterator end() const;   ///< Return the iterator pointing past the last row

	[[nodiscard]] bool empty() const;                                           ///< Return true if the range is empty
	explicit operator bool() const { return not empty(); }        ///< Easy way to detect if the range is empty
	[[nodiscard]] std::size_t size() const { return std::distance(begin(), end()); } ///< Return size of the range

	auto front() { return *begin(); } ///< Return reference to the first row
	// row_handle back() { return *begin(); }

	[[nodiscard]] category_type &get_category() const { return *m_cat; } ///< Category the iterators belong to

	/** swap */
	template <bool C2, typename ... T2s>
	friend void swap(conditional_iterator_proxy_base<C2, T2s...> &lhs, conditional_iterator_proxy_base<C2, T2s...> &rhs);

  private:
	category_type *m_cat;
	condition m_condition;
	row_iterator mCBegin, mCEnd;
	std::array<uint16_t, N> mCix;
};

// --------------------------------------------------------------------

template <typename... Ts>
using conditional_iterator_proxy = conditional_iterator_proxy_base<false, Ts...>;

template <typename... Ts>
using const_conditional_iterator_proxy = conditional_iterator_proxy_base<true, Ts...>;

// --------------------------------------------------------------------

/** @cond */
template <bool Const, typename... Ts>
iterator_proxy_base<Const, Ts...>::iterator_proxy_base(category_type &cat, row_iterator pos, char const *const items[N])
	: m_category(&cat)
	, m_begin(pos)
	, m_end(cat.end())
{
	for (uint16_t i = 0; i < N; ++i)
		m_item_ix[i] = m_category->get_item_ix(items[i]);
}

template <bool Const, typename... Ts>
iterator_proxy_base<Const, Ts...>::iterator_proxy_base(category_type &cat, row_iterator pos, std::initializer_list<char const *> items)
	: m_category(&cat)
	, m_begin(pos)
	, m_end(cat.end())
{
	// static_assert(items.size() == N, "The list of item names should be exactly the same as the list of requested items");

	std::uint16_t i = 0;
	for (auto item : items)
		m_item_ix[i++] = m_category->get_item_ix(item);
}

template <bool Const, typename... Ts>
iterator_proxy_base<Const, Ts...>::iterator_proxy_base(category_type &cat)
	: m_category(&cat)
	, m_begin(cat.begin())
	, m_end(cat.end())
{
	std::iota(m_item_ix.begin(), m_item_ix.end(), 0);
}

// --------------------------------------------------------------------

template <bool Const, typename... Ts>
conditional_iterator_proxy_base<Const, Ts...>::conditional_iterator_impl::conditional_iterator_impl(
	category_type &cat, row_iterator pos, const condition &cond, const std::array<uint16_t, N> &cix)
	: m_cat(&cat)
	, m_begin(pos, cix)
	, m_end(cat.end(), cix)
	, m_condition(&cond)
{
	if (m_condition == nullptr or m_condition->empty())
		m_begin = m_end;
	else
		m_current = *m_begin;
}

template <bool Const, typename... Ts>
template <typename... Ns>
conditional_iterator_proxy_base<Const, Ts...>::conditional_iterator_proxy_base(category_type &cat, row_iterator pos, condition &&cond, Ns... names)
	: m_cat(&cat)
	, m_condition(std::move(cond))
	, mCBegin(pos)
	, mCEnd(cat.end())
{
	static_assert(sizeof...(Ts) == sizeof...(Ns), "Number of item names should be equal to number of requested value types");

	if (m_condition and m_condition.prepare(cat))
	{
		while (mCBegin != mCEnd and not m_condition(*mCBegin))
			++mCBegin;
	}
	else
		mCBegin = mCEnd;

	uint16_t i = 0;
	((mCix[i++] = m_cat->get_item_ix(names)), ...);
}

template <bool Const, typename... Ts>
auto conditional_iterator_proxy_base<Const, Ts...>::begin() const -> iterator
{
	return iterator{ *m_cat, mCBegin, m_condition, mCix };
}

template <bool Const, typename... Ts>
auto conditional_iterator_proxy_base<Const, Ts...>::end() const -> iterator
{
	return iterator{ *m_cat, mCEnd, m_condition, mCix };
}

template <bool Const, typename... Ts>
bool conditional_iterator_proxy_base<Const, Ts...>::empty() const
{
	return mCBegin == mCEnd;
}

template <bool Const, typename... Ts>
void swap(conditional_iterator_proxy_base<Const, Ts...> &lhs, conditional_iterator_proxy_base<Const, Ts...> &rhs)
{
	std::swap(lhs.m_cat, rhs.m_cat);
	std::swap(lhs.m_condition, rhs.m_condition);
	std::swap(lhs.mCBegin, rhs.mCBegin);
	std::swap(lhs.mCEnd, rhs.mCEnd);
	std::swap(lhs.mCix, rhs.mCix);
}

// --------------------------------------------------------------------

// template <bool Const, typename... Ts>


/** @endcond */

} // namespace cif
