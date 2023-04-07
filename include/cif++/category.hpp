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

#include "cif++/forward_decl.hpp"

#include "cif++/condition.hpp"
#include "cif++/iterator.hpp"
#include "cif++/row.hpp"
#include "cif++/validate.hpp"

#include <array>

// TODO: implement all of:
// https://en.cppreference.com/w/cpp/named_req/Container
// https://en.cppreference.com/w/cpp/named_req/SequenceContainer
// and more?

namespace cif
{

// --------------------------------------------------------------------
// special exception
class duplicate_key_error : public std::runtime_error
{
  public:
	duplicate_key_error(const std::string &msg)
		: std::runtime_error(msg)
	{
	}
};

class multiple_results_error : public std::runtime_error
{
  public:
	multiple_results_error()
		: std::runtime_error("query should have returned exactly one row")
	{
	}
};

// --------------------------------------------------------------------
// These should be moved elsewhere, one day.

template<typename _Tp> inline constexpr bool is_optional_v = false;
template<typename _Tp> inline constexpr bool is_optional_v<std::optional<_Tp>> = true;

// --------------------------------------------------------------------

class category
{
  public:
	friend class row_handle;

	template <typename, typename...>
	friend class iterator_impl;

	using value_type = row_handle;
	using reference = value_type;
	using const_reference = const value_type;
	using iterator = iterator_impl<category>;
	using const_iterator = iterator_impl<const category>;

	category() = default;

	category(std::string_view name);

	category(const category &rhs);

	category(category &&rhs);

	category &operator=(const category &rhs);

	category &operator=(category &&rhs);

	~category();

	// --------------------------------------------------------------------

	const std::string &name() const { return m_name; }

	iset key_fields() const;

	std::set<uint16_t> key_field_indices() const;

	void set_validator(const validator *v, datablock &db);
	void update_links(datablock &db);

	const validator *get_validator() const { return m_validator; }
	const category_validator *get_cat_validator() const { return m_cat_validator; }

	bool is_valid() const;
	bool validate_links() const;

	bool operator==(const category &rhs) const;
	bool operator!=(const category &rhs) const
	{
		return not operator==(rhs);
	}

	// --------------------------------------------------------------------

	reference front()
	{
		return { *this, *m_head };
	}

	const_reference front() const
	{
		return { const_cast<category &>(*this), const_cast<row &>(*m_head) };
	}

	reference back()
	{
		return { *this, *m_tail };
	}

	const_reference back() const
	{
		return { const_cast<category &>(*this), const_cast<row &>(*m_tail) };
	}

	iterator begin()
	{
		return { *this, m_head };
	}

	iterator end()
	{
		return { *this, nullptr };
	}

	const_iterator begin() const
	{
		return { *this, m_head };
	}

	const_iterator end() const
	{
		return { *this, nullptr };
	}

	const_iterator cbegin() const
	{
		return { *this, m_head };
	}

	const_iterator cend() const
	{
		return { *this, nullptr };
	}

	size_t size() const
	{
		return std::distance(cbegin(), cend());
	}

	bool empty() const
	{
		return m_head == nullptr;
	}

	// --------------------------------------------------------------------
	// A category can have a key, as defined by the validator/dictionary

	/// @brief The key type
	using key_type = row_initializer;

	/// @brief Return a row_handle for the row specified by \a key
	/// @param key The value for the key, fields specified in the dictionary should have a value
	/// @return The row found in the index, or an undefined row_handle
	row_handle operator[](const key_type &key);

	const row_handle operator[](const key_type &key) const
	{
		return const_cast<category *>(this)->operator[](key);
	}

	// --------------------------------------------------------------------

	template <typename... Ts, typename... Ns>
	iterator_proxy<const category, Ts...> rows(Ns... names) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of column titles should be equal to the number of types to return");
		return iterator_proxy<const category, Ts...>(*this, begin(), { names... });
	}

	template <typename... Ts, typename... Ns>
	iterator_proxy<category, Ts...> rows(Ns... names)
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of column titles should be equal to the number of types to return");
		return iterator_proxy<category, Ts...>(*this, begin(), { names... });
	}

	// --------------------------------------------------------------------

	conditional_iterator_proxy<category> find(condition &&cond)
	{
		return find(begin(), std::move(cond));
	}

	conditional_iterator_proxy<category> find(iterator pos, condition &&cond)
	{
		return { *this, pos, std::move(cond) };
	}

	conditional_iterator_proxy<const category> find(condition &&cond) const
	{
		return find(cbegin(), std::move(cond));
	}

	conditional_iterator_proxy<const category> find(const_iterator pos, condition &&cond) const
	{
		return conditional_iterator_proxy<const category>{ *this, pos, std::move(cond) };
	}

	template <typename... Ts, typename... Ns>
	conditional_iterator_proxy<category, Ts...> find(condition &&cond, Ns... names)
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of column titles should be equal to the number of types to return");
		return find<Ts...>(cbegin(), std::move(cond), std::forward<Ns>(names)...);
	}

	template <typename... Ts, typename... Ns>
	conditional_iterator_proxy<const category, Ts...> find(condition &&cond, Ns... names) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of column titles should be equal to the number of types to return");
		return find<Ts...>(cbegin(), std::move(cond), std::forward<Ns>(names)...);
	}

	template <typename... Ts, typename... Ns>
	conditional_iterator_proxy<category, Ts...> find(const_iterator pos, condition &&cond, Ns... names)
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of column titles should be equal to the number of types to return");
		return { *this, pos, std::move(cond), std::forward<Ns>(names)... };
	}

	template <typename... Ts, typename... Ns>
	conditional_iterator_proxy<const category, Ts...> find(const_iterator pos, condition &&cond, Ns... names) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of column titles should be equal to the number of types to return");
		return { *this, pos, std::move(cond), std::forward<Ns>(names)... };
	}

	// --------------------------------------------------------------------
	// if you only expect a single row

	row_handle find1(condition &&cond)
	{
		return find1(begin(), std::move(cond));
	}

	row_handle find1(iterator pos, condition &&cond)
	{
		auto h = find(pos, std::move(cond));

		if (h.size() != 1)
			throw multiple_results_error();

		return *h.begin();
	}

	const row_handle find1(condition &&cond) const
	{
		return find1(cbegin(), std::move(cond));
	}

	const row_handle find1(const_iterator pos, condition &&cond) const
	{
		auto h = find(pos, std::move(cond));

		if (h.size() != 1)
			throw multiple_results_error();

		return *h.begin();
	}

	template <typename T>
	T find1(condition &&cond, const char *column) const
	{
		return find1<T>(cbegin(), std::move(cond), column);
	}

	template <typename T, std::enable_if_t<not is_optional_v<T>, int> = 0>
	T find1(const_iterator pos, condition &&cond, const char *column) const
	{
		auto h = find<T>(pos, std::move(cond), column);

		if (h.size() != 1)
			throw multiple_results_error();

		return *h.begin();
	}

	template <typename T, std::enable_if_t<is_optional_v<T>, int> = 0>
	T find1(const_iterator pos, condition &&cond, const char *column) const
	{
		auto h = find<typename T::value_type>(pos, std::move(cond), column);

		if (h.size() > 1)
			throw multiple_results_error();

		if (h.empty())
			return {};

		return *h.begin();
	}

	template <typename... Ts, typename... Cs, typename U = std::enable_if_t<sizeof...(Ts) != 1>>
	std::tuple<Ts...> find1(condition &&cond, Cs... columns) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Cs), "The number of column titles should be equal to the number of types to return");
		// static_assert(std::is_same_v<Cs, const char*>..., "The column names should be const char");
		return find1<Ts...>(cbegin(), std::move(cond), std::forward<Cs>(columns)...);
	}

	template <typename... Ts, typename... Cs, typename U = std::enable_if_t<sizeof...(Ts) != 1>>
	std::tuple<Ts...> find1(const_iterator pos, condition &&cond, Cs... columns) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Cs), "The number of column titles should be equal to the number of types to return");
		auto h = find<Ts...>(pos, std::move(cond), std::forward<Cs>(columns)...);

		if (h.size() != 1)
			throw multiple_results_error();

		return *h.begin();
	}

	// --------------------------------------------------------------------
	// if you want only a first hit

	row_handle find_first(condition &&cond)
	{
		return find_first(begin(), std::move(cond));
	}

	row_handle find_first(iterator pos, condition &&cond)
	{
		auto h = find(pos, std::move(cond));

		return h.empty() ? row_handle{} : *h.begin();
	}

	const row_handle find_first(condition &&cond) const
	{
		return find_first(cbegin(), std::move(cond));
	}

	const row_handle find_first(const_iterator pos, condition &&cond) const
	{
		auto h = find(pos, std::move(cond));

		return h.empty() ? row_handle{} : *h.begin();
	}

	template <typename T>
	T find_first(condition &&cond, const char *column) const
	{
		return find_first<T>(cbegin(), std::move(cond), column);
	}

	template <typename T>
	T find_first(const_iterator pos, condition &&cond, const char *column) const
	{
		auto h = find<T>(pos, std::move(cond), column);

		return h.empty() ? T{} : *h.begin();
	}

	template <typename... Ts, typename... Cs, typename U = std::enable_if_t<sizeof...(Ts) != 1>>
	std::tuple<Ts...> find_first(condition &&cond, Cs... columns) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Cs), "The number of column titles should be equal to the number of types to return");
		// static_assert(std::is_same_v<Cs, const char*>..., "The column names should be const char");
		return find_first<Ts...>(cbegin(), std::move(cond), std::forward<Cs>(columns)...);
	}

	template <typename... Ts, typename... Cs, typename U = std::enable_if_t<sizeof...(Ts) != 1>>
	std::tuple<Ts...> find_first(const_iterator pos, condition &&cond, Cs... columns) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Cs), "The number of column titles should be equal to the number of types to return");
		auto h = find<Ts...>(pos, std::move(cond), std::forward<Cs>(columns)...);

		return h.empty() ? std::tuple<Ts...>{} : *h.begin();
	}

	// --------------------------------------------------------------------

	template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	T find_max(const char *column, condition &&cond) const
	{
		T result = std::numeric_limits<T>::min();

		for (auto v : find<T>(std::move(cond), column))
		{
			if (result < v)
				result = v;
		}

		return result;
	}

	template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	T find_max(const char *column) const
	{
		return find_max<T>(column, all());
	}

	template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	T find_min(const char *column, condition &&cond) const
	{
		T result = std::numeric_limits<T>::max();

		for (auto v : find<T>(std::move(cond), column))
		{
			if (result > v)
				result = v;
		}

		return result;
	}

	template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	T find_min(const char *column) const
	{
		return find_min<T>(column, all());
	}

	bool exists(condition &&cond) const
	{
		bool result = false;

		if (cond)
		{
			cond.prepare(*this);

			auto sh = cond.single();

			if (sh.has_value() and *sh)
				result = true;
			else
			{
				for (auto r : *this)
				{
					if (cond(r))
					{
						result = true;
						break;
					}
				}
			}
		}

		return result;
	}

	size_t count(condition &&cond) const
	{
		size_t result = 0;

		if (cond)
		{
			cond.prepare(*this);

			auto sh = cond.single();

			if (sh.has_value() and *sh)
				result = 1;
			else
			{
				for (auto r : *this)
				{
					if (cond(r))
						++result;
				}
			}
		}

		return result;
	}

	// --------------------------------------------------------------------

	bool has_children(row_handle r) const;
	bool has_parents(row_handle r) const;

	std::vector<row_handle> get_children(row_handle r, const category &childCat) const;
	std::vector<row_handle> get_parents(row_handle r, const category &parentCat) const;
	std::vector<row_handle> get_linked(row_handle r, const category &cat) const;

	// --------------------------------------------------------------------

	// void insert(const_iterator pos, const row_initializer &row)
	// {
	// 	insert_impl(pos, row);
	// }

	// void insert(const_iterator pos, row_initializer &&row)
	// {
	// 	insert_impl(pos, std::move(row));
	// }

	iterator erase(iterator pos);
	void erase(row_handle rh)
	{
		erase(iterator(*this, rh.m_row));
	}

	size_t erase(condition &&cond);
	size_t erase(condition &&cond, std::function<void(row_handle)> &&visit);

	iterator emplace(row_initializer &&ri)
	{
		return this->emplace(ri.begin(), ri.end());
	}

	template <typename ItemIter>
	iterator emplace(ItemIter b, ItemIter e)
	{
		row *r = this->create_row();

		try
		{
			for (auto i = b; i != e; ++i)
			{
				// item_value *new_item = this->create_item(*i);
				r->append(add_column(i->name()), { i->value() });
			}
		}
		catch (...)
		{
			if (r != nullptr)
				this->delete_row(r);
			throw;
		}

		return insert_impl(cend(), r);
	}

	void clear();

	// --------------------------------------------------------------------
	/// \brief generate a new, unique ID. Pass it an ID generating function
	/// based on a sequence number. This function will be called until the
	/// result is unique in the context of this category
	std::string get_unique_id(std::function<std::string(int)> generator = cif::cif_id_for_number);
	std::string get_unique_id(const std::string &prefix)
	{
		return get_unique_id([prefix](int nr)
			{ return prefix + std::to_string(nr + 1); });
	}

	// --------------------------------------------------------------------

	/// \brief Rename a single column in the rows that match \a cond to value \a value
	/// making sure the linked categories are updated according to the link.
	/// That means, child categories are updated if the links are absolute
	/// and unique. If they are not, the child category rows are split.

	void update_value(condition &&cond, std::string_view tag, std::string_view value)
	{
		auto rs = find(std::move(cond));
		std::vector<row_handle> rows;
		std::copy(rs.begin(), rs.end(), std::back_inserter(rows));
		update_value(rows, tag, value);
	}

	void update_value(const std::vector<row_handle> &rows, std::string_view tag, std::string_view value);

	// --------------------------------------------------------------------
	/// \brief Return the index number for \a column_name

	uint16_t get_column_ix(std::string_view column_name) const
	{
		uint16_t result;

		for (result = 0; result < m_columns.size(); ++result)
		{
			if (iequals(column_name, m_columns[result].m_name))
				break;
		}

		if (VERBOSE > 0 and result == m_columns.size() and m_cat_validator != nullptr) // validate the name, if it is known at all (since it was not found)
		{
			auto iv = m_cat_validator->get_validator_for_item(column_name);
			if (iv == nullptr)
				std::cerr << "Invalid name used '" << column_name << "' is not a known column in " + m_name << std::endl;
		}

		return result;
	}

	std::string_view get_column_name(uint16_t ix) const
	{
		if (ix >= m_columns.size())
			throw std::out_of_range("column index is out of range");

		return m_columns[ix].m_name;
	}

	uint16_t add_column(std::string_view column_name)
	{
		using namespace std::literals;

		uint16_t result = get_column_ix(column_name);

		if (result == m_columns.size())
		{
			const item_validator *item_validator = nullptr;

			if (m_cat_validator != nullptr)
			{
				item_validator = m_cat_validator->get_validator_for_item(column_name);
				if (item_validator == nullptr)
					m_validator->report_error("tag " + std::string(column_name) + " not allowed in category " + m_name, false);
			}

			m_columns.emplace_back(column_name, item_validator);
		}

		return result;
	}

	bool has_column(std::string_view name) const
	{
		return get_column_ix(name) < m_columns.size();
	}

	iset get_columns() const;

	// --------------------------------------------------------------------

	void sort(std::function<int(row_handle, row_handle)> f);
	void reorder_by_index();

	// --------------------------------------------------------------------

	std::vector<std::string> get_tag_order() const;

	void write(std::ostream &os) const;
	void write(std::ostream &os, const std::vector<std::string> &order, bool addMissingColumns = true);

  private:
	void write(std::ostream &os, const std::vector<uint16_t> &order, bool includeEmptyColumns) const;

  public:
	friend std::ostream &operator<<(std::ostream &os, const category &cat)
	{
		cat.write(os);
		return os;
	}

  private:
	void update_value(row *row, uint16_t column, std::string_view value, bool updateLinked, bool validate = true);

  private:
	void erase_orphans(condition &&cond, category &parent);

	using allocator_type = std::allocator<void>;

	constexpr allocator_type get_allocator() const
	{
		return {};
	}

	using char_allocator_type = typename std::allocator_traits<allocator_type>::template rebind_alloc<char>;
	using char_allocator_traits = std::allocator_traits<char_allocator_type>;

	using row_allocator_type = typename std::allocator_traits<allocator_type>::template rebind_alloc<row>;
	using row_allocator_traits = std::allocator_traits<row_allocator_type>;

	row_allocator_traits::pointer get_row()
	{
		row_allocator_type ra(get_allocator());
		return row_allocator_traits::allocate(ra, 1);
	}

	row *create_row()
	{
		auto p = this->get_row();
		row_allocator_type ra(get_allocator());
		row_allocator_traits::construct(ra, p);
		return p;
	}

	row *clone_row(const row &r);

	void delete_row(row *r);

	row_handle create_copy(row_handle r);

	struct item_column
	{
		std::string m_name;
		const item_validator *m_validator;

		item_column(std::string_view name, const item_validator *validator)
			: m_name(name)
			, m_validator(validator)
		{
		}
	};

	struct link
	{
		link(category *linked, const link_validator *v)
			: linked(linked)
			, v(v)
		{
		}

		category *linked;
		const link_validator *v;
	};

	// proxy methods for every insertion
	iterator insert_impl(const_iterator pos, row *n);
	iterator erase_impl(const_iterator pos);

	// --------------------------------------------------------------------

	condition get_parents_condition(row_handle rh, const category &parentCat) const;
	condition get_children_condition(row_handle rh, const category &childCat) const;

	// --------------------------------------------------------------------

	void swap_item(uint16_t column_ix, row_handle &a, row_handle &b);

	// --------------------------------------------------------------------

	std::string m_name;
	std::vector<item_column> m_columns;
	const validator *m_validator = nullptr;
	const category_validator *m_cat_validator = nullptr;
	std::vector<link> m_parent_links, m_child_links;
	bool m_cascade = true;
	uint32_t m_last_unique_num = 0;
	class category_index *m_index = nullptr;
	row *m_head = nullptr, *m_tail = nullptr;
};

} // namespace cif
