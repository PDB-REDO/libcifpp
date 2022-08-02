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

#include <cif++/v2/row.hpp>

namespace cif::v2
{

// --------------------------------------------------------------------

template <typename Category, typename... Ts>
class iterator_impl
{
  public:
	template <typename, typename...>
	friend class iterator_impl;

	template <typename>
	friend class category_t;

	static constexpr size_t N = sizeof...(Ts);

	using category_type = Category;
	using row_type = std::conditional_t<std::is_const_v<category_type>, const typename category_type::row, typename category_type::row>;
	using row_handle_type = std::conditional_t<std::is_const_v<category_type>, const row_handle<category_type>, row_handle<category_type>>;

	using iterator_category = std::forward_iterator_tag;
	using value_type = std::conditional_t<N == 0, row_handle_type, std::tuple<Ts...>>;
	using difference_type = std::ptrdiff_t;
	using pointer = std::conditional_t<N == 0, row_handle_type, value_type *>;
	using reference = std::conditional_t<N == 0, row_handle_type, value_type &>;

	iterator_impl(const iterator_impl &rhs) = default;

	template<typename C2, typename... T2s>
	iterator_impl(const iterator_impl<C2, T2s...> &rhs)
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
		, m_value(rhs.m_value)
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
			m_value = get(std::make_index_sequence<N>());
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
			return {*m_category, *m_current};
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
			m_value = get(std::make_index_sequence<N>());

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
	std::tuple<Ts...> get(std::index_sequence<Is...>) const
	{
		if (m_current != nullptr)
		{
			row_handle_type rh{*m_category, *m_current};
			return std::tuple<Ts...>{rh[m_column_ix[Is]].template as<Ts>()...};
		}

		return {};
	}

	category_type *m_category;
	row_type *m_current;
	value_type m_value;
	std::array<size_t, N> m_column_ix;
};

// --------------------------------------------------------------------
// iterator proxy

template <typename Category, typename... Ts>
class iterator_proxy
{
  public:
	static constexpr const size_t N = sizeof...(Ts);

	using category_type = Category;
	using row_type = std::conditional_t<std::is_const_v<category_type>, const typename category_type::row, typename category_type::row>;
	using row_handle_type = std::conditional_t<std::is_const_v<category_type>, const row_handle<category_type>, row_handle<category_type>>;

	using iterator = iterator_impl<category_type, Ts...>;
	using row_iterator = iterator_impl<category_type>;

	iterator_proxy(category_type &cat, row_iterator pos, char const *const columns[N]);
	iterator_proxy(category_type &cat, row_iterator pos, std::initializer_list<char const *> columns);

	iterator_proxy(iterator_proxy &&p);
	iterator_proxy &operator=(iterator_proxy &&p);

	iterator_proxy(const iterator_proxy &) = delete;
	iterator_proxy &operator=(const iterator_proxy &) = delete;

	iterator begin() const { return iterator(m_begin, m_column_ix); }
	iterator end() const { return iterator(m_end, m_column_ix); }

	bool empty() const { return m_begin == m_end; }

	explicit operator bool() const { return not empty(); }

	size_t size() const { return std::distance(begin(), end()); }

	row_type front() { return *begin(); }
	row_type back() { return *(std::prev(end())); }

	category_type &category() const { return *m_category; }

	void swap(iterator_proxy &rhs)
	{
		std::swap(m_category, rhs.m_category);
		std::swap(m_begin, rhs.m_begin);
		std::swap(m_end, rhs.m_end);
		std::swap(m_column_ix, rhs.m_column_ix);
	}

  private:
	category_type *m_category;
	row_iterator m_begin, m_end;
	std::array<size_t, N> m_column_ix;
};

// // --------------------------------------------------------------------
// // conditional iterator proxy

// template <typename CategoryType, typename... Ts>
// class conditional_iterator_proxy
// {
//   public:
// 	static constexpr const size_t N = sizeof...(Ts);

// 	using base_iterator = iterator_impl<CategoryType, Ts...>;
// 	using value_type = typename base_iterator::value_type;
// 	using row_type = typename base_iterator::row_type;
// 	using row_iterator = iterator_impl<row_type>;

// 	class conditional_iterator_impl
// 	{
// 	  public:
// 		using iterator_category = std::forward_iterator_tag;
// 		using value_type = conditional_iterator_proxy::value_type;
// 		using difference_type = std::ptrdiff_t;
// 		using pointer = value_type *;
// 		using reference = value_type &;

// 		conditional_iterator_impl(CategoryType &cat, row_iterator pos, const Condition &cond, const std::array<size_t, N> &cix);
// 		conditional_iterator_impl(const conditional_iterator_impl &i) = default;
// 		conditional_iterator_impl &operator=(const conditional_iterator_impl &i) = default;

// 		virtual ~conditional_iterator_impl() = default;

// 		reference operator*()
// 		{
// 			return *mBegin;
// 		}

// 		pointer operator->()
// 		{
// 			return &*mBegin;
// 		}

// 		conditional_iterator_impl &operator++()
// 		{
// 			while (mBegin != mEnd)
// 			{
// 				if (++mBegin == mEnd)
// 					break;

// 				if ((*mCondition)(*mCat, mBegin.row()))
// 					break;
// 			}

// 			return *this;
// 		}

// 		conditional_iterator_impl operator++(int)
// 		{
// 			conditional_iterator_impl result(*this);
// 			this->operator++();
// 			return result;
// 		}

// 		bool operator==(const conditional_iterator_impl &rhs) const { return mBegin == rhs.mBegin; }
// 		bool operator!=(const conditional_iterator_impl &rhs) const { return mBegin != rhs.mBegin; }

// 		template <typename IRowType, typename... ITs>
// 		bool operator==(const iterator_impl<IRowType, ITs...> &rhs) const { return mBegin == rhs; }

// 		template <typename IRowType, typename... ITs>
// 		bool operator!=(const iterator_impl<IRowType, ITs...> &rhs) const { return mBegin != rhs; }

// 	  private:
// 		CategoryType *mCat;
// 		base_iterator mBegin, mEnd;
// 		const Condition *mCondition;
// 	};

// 	using iterator = conditional_iterator_impl;
// 	using reference = typename iterator::reference;

// 	template <typename... Ns>
// 	conditional_iterator_proxy(CategoryType &cat, row_iterator pos, Condition &&cond, Ns... names);

// 	conditional_iterator_proxy(conditional_iterator_proxy &&p);
// 	conditional_iterator_proxy &operator=(conditional_iterator_proxy &&p);

// 	conditional_iterator_proxy(const conditional_iterator_proxy &) = delete;
// 	conditional_iterator_proxy &operator=(const conditional_iterator_proxy &) = delete;

// 	iterator begin() const;
// 	iterator end() const;

// 	bool empty() const;

// 	explicit operator bool() const { return not empty(); }

// 	size_t size() const { return std::distance(begin(), end()); }

// 	row_type front() { return *begin(); }

// 	CategoryType &category() const { return *mCat; }

// 	void swap(conditional_iterator_proxy &rhs);

//   private:
// 	CategoryType *mCat;
// 	Condition mCondition;
// 	row_iterator mCBegin, mCEnd;
// 	std::array<size_t, N> mCix;
// };

// --------------------------------------------------------------------

template <typename Category, typename... Ts>
iterator_proxy<Category, Ts...>::iterator_proxy(Category &cat, row_iterator pos, char const *const columns[N])
	: m_category(&cat)
	, m_begin(pos)
	, m_end(cat.end())
{
	for (size_t i = 0; i < N; ++i)
		m_column_ix[i] = m_category->get_column_ix(columns[i]);
}

template <typename Category, typename... Ts>
iterator_proxy<Category, Ts...>::iterator_proxy(Category &cat, row_iterator pos, std::initializer_list<char const *> columns)
	: m_category(&cat)
	, m_begin(pos)
	, m_end(cat.end())
{
	// static_assert(columns.size() == N, "The list of column names should be exactly the same as the list of requested columns");

	std::size_t i = 0;
	for (auto column : columns)
		m_column_ix[i++] = m_category->get_column_ix(column);
}

// --------------------------------------------------------------------

// template <typename Category, typename... Ts>
// conditional_iterator_proxy<Category, Ts...>::conditional_iterator_impl::conditional_iterator_impl(
// 	Category &cat, row_iterator pos, const Condition &cond, const std::array<size_t, N> &cix)
// 	: mCat(&cat)
// 	, mBegin(pos, cix)
// 	, mEnd(cat.end(), cix)
// 	, mCondition(&cond)
// {
// }

// template <typename Category, typename... Ts>
// conditional_iterator_proxy<Category, Ts...>::conditional_iterator_proxy(conditional_iterator_proxy &&p)
// 	: mCat(nullptr)
// 	, mCBegin(p.mCBegin)
// 	, mCEnd(p.mCEnd)
// 	, mCix(p.mCix)
// {
// 	std::swap(mCat, p.mCat);
// 	std::swap(mCix, p.mCix);
// 	mCondition.swap(p.mCondition);
// }

// template <typename Category, typename... Ts>
// template <typename... Ns>
// conditional_iterator_proxy<Category, Ts...>::conditional_iterator_proxy(Category &cat, row_iterator pos, Condition &&cond, Ns... names)
// 	: mCat(&cat)
// 	, mCondition(std::move(cond))
// 	, mCBegin(pos)
// 	, mCEnd(cat.end())
// {
// 	static_assert(sizeof...(Ts) == sizeof...(Ns), "Number of column names should be equal to number of requested value types");

// 	mCondition.prepare(cat);

// 	while (mCBegin != mCEnd and not mCondition(*mCat, mCBegin.row()))
// 		++mCBegin;

// 	size_t i = 0;
// 	((mCix[i++] = mCat->getColumnIndex(names)), ...);
// }

// template <typename Category, typename... Ts>
// conditional_iterator_proxy<Category, Ts...> &conditional_iterator_proxy<Category, Ts...>::operator=(conditional_iterator_proxy &&p)
// {
// 	swap(p);
// 	return *this;
// }

// template <typename Category, typename... Ts>
// typename conditional_iterator_proxy<Category, Ts...>::iterator conditional_iterator_proxy<Category, Ts...>::begin() const
// {
// 	return iterator(*mCat, mCBegin, mCondition, mCix);
// }

// template <typename Category, typename... Ts>
// typename conditional_iterator_proxy<Category, Ts...>::iterator conditional_iterator_proxy<Category, Ts...>::end() const
// {
// 	return iterator(*mCat, mCEnd, mCondition, mCix);
// }

// template <typename Category, typename... Ts>
// bool conditional_iterator_proxy<Category, Ts...>::empty() const
// {
// 	return mCBegin == mCEnd;
// }

// template <typename Category, typename... Ts>
// void conditional_iterator_proxy<Category, Ts...>::swap(conditional_iterator_proxy &rhs)
// {
// 	std::swap(mCat, rhs.mCat);
// 	mCondition.swap(rhs.mCondition);
// 	std::swap(mCBegin, rhs.mCBegin);
// 	std::swap(mCEnd, rhs.mCEnd);
// 	std::swap(mCix, rhs.mCix);
// }

} // namespace cif::v2