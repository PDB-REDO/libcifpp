/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2025  NKI/AVL, Netherlands Cancer Institute
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

/// @file cql.hpp
///
/// This file contains code to access stored data as if it were
/// a relation database. The underlying code uses SQLite as engine.
/// categories are exposed as virtual tables.
///

#include "cif++/category.hpp"
#include "cif++/item.hpp"
#include "cif++/iterator.hpp"
#include "cif++/row.hpp"

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

// --------------------------------------------------------------------

namespace cif::cql
{

class connection;
struct result_impl;

// --------------------------------------------------------------------

/// Reference to a field in the result set
class field_ref final
{
  public:
	/// The name of the field
	[[nodiscard]] std::string_view name() const &
	{
		return m_row.get_category().get_item_name(m_index);
	}

	/// The index number of the field
	[[nodiscard]] constexpr size_t num() const noexcept
	{
		return m_index;
	}

	/// Return the contents of this item as type @tparam T
	template <typename T = std::string>
	[[nodiscard]] auto get() const -> T
	{
		return m_row[m_index].get<T>();
	}

	/// Returns true if the field contains NULL
	[[nodiscard]] bool is_null() const
	{
		return m_row[m_index].is_null();
	}

	/// Return the contents of this item as type @tparam T or, if not
	/// set, use @a dv as the default value.
	///
	template <typename T>
	auto value_or(const T &dv) const
	{
		return m_row[m_index].value_or(dv);
	}

	/// Constructor
	field_ref(const_row_handle rh, uint16_t col, std::shared_ptr<result_impl> result_impl)
		: m_row(std::move(rh))
		, m_index(col)
		, m_result_impl(std::move(result_impl))
	{
	}

	/// Copy constructor
	field_ref(const field_ref &) = default;
	/// Move constructor
	field_ref(field_ref &&) = default;

	/// Copy assignment
	field_ref &operator=(const field_ref &) = default;
	/// Move assignment
	field_ref &operator=(field_ref &&) = default;

  private:
	const_row_handle m_row;
	uint16_t m_index;

	std::shared_ptr<result_impl> m_result_impl;
};

// --------------------------------------------------------------------

/// A reference to a row in the result set
class row_ref final
{
  public:
	/// Iterator for the items in this row
	class const_field_iterator
	{
	  public:
		friend class result;

		/// @cond
		using iterator_category = std::forward_iterator_tag;
		using value_type = const field_ref;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type *;
		using reference = value_type &;

		const_field_iterator(const const_field_iterator &) = default;
		const_field_iterator(const_field_iterator &&) = default;

		const_field_iterator &operator=(const const_field_iterator &) = default;
		const_field_iterator &operator=(const_field_iterator &&) = default;

		reference operator*()
		{
			return m_current;
		}

		pointer operator->()
		{
			return &m_current;
		}

		const_field_iterator &operator++()
		{
			if (m_row)
			{
				++m_col;
				m_current = field_ref(m_row, m_col, m_result_impl);
			}

			return *this;
		}

		const_field_iterator operator++(int)
		{
			const_field_iterator result(*this);
			this->operator++();
			return result;
		}

		bool operator==(const const_field_iterator &rhs) const
		{
			return m_row == rhs.m_row and m_col == rhs.m_col;
		}

		bool operator!=(const const_field_iterator &rhs) const
		{
			return m_row != rhs.m_row or m_col != rhs.m_col;
		}

	  private:
		friend class row_ref;

		const_field_iterator(const_row_handle row, uint16_t column, std::shared_ptr<result_impl> result_impl)
			: m_row(std::move(row))
			, m_col(column)
			, m_current(m_row, m_col, result_impl)
			, m_result_impl(result_impl)
		{
		}

		const_row_handle m_row;
		uint16_t m_col;
		field_ref m_current;

		std::shared_ptr<result_impl> m_result_impl;

		/// @endcond
	};

	// --------------------------------------------------------------------

	row_ref() = default;

	/// Constructor
	row_ref(const_row_handle rh, std::shared_ptr<result_impl> result_impl)
		: m_row(std::move(rh))
		, m_result_impl(std::move(result_impl))
	{
	}

	/// @cond
	row_ref(const row_ref &) = default;
	row_ref &operator=(const row_ref &) = default;
	/// @endcond

	// --------------------------------------------------------------------

	[[nodiscard]] const_field_iterator begin() const noexcept { return { m_row, 0, m_result_impl }; }                            ///< Return begin field iterator
	[[nodiscard]] const_field_iterator cbegin() const noexcept { return { m_row, 0, m_result_impl }; }                           ///< Return cbegin field iterator
	[[nodiscard]] const_field_iterator end() const noexcept { return { m_row, static_cast<uint16_t>(size()), m_result_impl }; }  ///< Return end field iterator
	[[nodiscard]] const_field_iterator cend() const noexcept { return { m_row, static_cast<uint16_t>(size()), m_result_impl }; } ///< Return cend field iterator

	[[nodiscard]] field_ref front() const noexcept { return { m_row, 0, m_result_impl }; }                                ///< return reference to front field
	[[nodiscard]] field_ref back() const noexcept { return { m_row, static_cast<uint16_t>(size() - 1), m_result_impl }; } ///< return reference to back field

	[[nodiscard]] size_t size() const noexcept;                       ///< return number of items in the row
	[[nodiscard]] bool empty() const noexcept { return size() == 0; } ///< return if the row contains no items at all

	[[nodiscard]] field_ref operator[](uint16_t index) const noexcept { return { m_row, index, m_result_impl }; } ///< access field by index
	[[nodiscard]] field_ref operator[](std::string_view name) const;                                              ///< access field by name

	// --------------------------------------------------------------------

	/// @cond
	bool operator==(const row_ref &rhs) const { return m_row == rhs.m_row; }
	bool operator!=(const row_ref &rhs) const { return m_row != rhs.m_row; }

  private:
	const_row_handle m_row;
	std::shared_ptr<result_impl> m_result_impl;

	/// @endcond
};

// --------------------------------------------------------------------

/// The result set, containing the result of a query
class result
{
  public:
	// --------------------------------------------------------------------

	/// iterator to the rows in the result set
	class iterator
	{
	  public:
		friend class view;

		/// @cond

		using iterator_category = std::forward_iterator_tag;
		using value_type = const row_ref;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type *;
		using reference = value_type &;

		// const_row_iterator() = default;

		iterator(std::shared_ptr<result_impl> result_impl, category::const_iterator cat_iter)
			: m_iter(std::move(cat_iter))
			, m_current(*m_iter, result_impl)
			, m_result_impl(result_impl)
		{
		}

		iterator(const iterator &) = default;
		iterator(iterator &&) = default;

		// const_row_iterator &operator=(const const_row_iterator &) = default;
		// const_row_iterator &operator=(const_row_iterator &&) = default;

		reference operator*()
		{
			return m_current;
		}

		pointer operator->()
		{
			return &m_current;
		}

		iterator &operator++()
		{
			++m_iter;
			m_current = { *m_iter, m_result_impl };
			return *this;
		}

		iterator operator++(int)
		{
			iterator result(*this);
			this->operator++();
			return result;
		}

		bool operator==(const iterator &rhs) const
		{
			return m_result_impl == rhs.m_result_impl and m_iter == rhs.m_iter;
		}

		bool operator!=(const iterator &rhs) const
		{
			return m_result_impl != rhs.m_result_impl or m_iter != rhs.m_iter;
		}

	  private:
		category::const_iterator m_iter;
		row_ref m_current;
		std::shared_ptr<result_impl> m_result_impl;

		/// @endcond
	};

	// --------------------------------------------------------------------

	/// @cond
	result() = delete;
	result(result const &rhs) noexcept = default;
	result(result &&rhs) noexcept = default;
	result &operator=(result const &rhs) noexcept = default;
	result &operator=(result &&rhs) noexcept = default;

	result(category &&data, const std::string &query = "");

	~result() = default;
	/// @endcond

	/// Return the row if and only if the result set contains exactly one row, throws otherwise
	[[nodiscard]] row_ref one_row() const
	{
		if (size() != 1)
			throw std::runtime_error("Expected one row");
		return front();
	}

	/// Return the row if and only if the result set contains exactly one row,
	/// and this row also contains only one field, throws otherwise
	[[nodiscard]] field_ref one_field() const
	{
		expect_columns(1);

		if (size() != 1)
			throw std::runtime_error("Expected one row");

		return one_row().front();
	}

	// --------------------------------------------------------------------

	/// @cond
	[[nodiscard]] iterator begin() const noexcept;
	[[nodiscard]] iterator cbegin() const noexcept;

	[[nodiscard]] iterator end() const noexcept;
	[[nodiscard]] iterator cend() const noexcept;

	[[nodiscard]] row_ref front() const;
	[[nodiscard]] row_ref back() const;

	[[nodiscard]] size_t size() const noexcept;
	[[nodiscard]] bool empty() const noexcept { return size() == 0; }
	/// @endcond

	/// Return the number of colums/fields in each row
	[[nodiscard]] size_t column_count() const;

	/// Return the result set as a cif::category
	[[nodiscard]] category &get_category() const;

	/// Test to see if the result set contains at least the number of fields/columns
	/// but only when not empty
	void expect_columns(size_t cols) const
	{
		if (auto actual = column_count(); size() > 0 and cols != actual)
			throw std::runtime_error("Unexpected number of columns");
	}

	// --------------------------------------------------------------------

	/// Print out the result set, for debugging mostly
	friend std::ostream &operator<<(std::ostream &os, const result &r)
	{
		os << r.get_category();
		return os;
	}

  private:
	friend class transaction;
	friend class SelectStatement;

	std::shared_ptr<result_impl> m_impl;
};

// --------------------------------------------------------------------

/// Helper class to allow access to the data as a stream
template <typename... Ts>
class cql_iterator_proxy : public cif::iterator_proxy<Ts...>
{
  public:
	/// Constructor
	cql_iterator_proxy(result &&res)
		: cif::iterator_proxy<Ts...>(res.get_category())
		, m_result(std::forward<result>(res))
	{
		m_result.expect_columns(cif::iterator_proxy<Ts...>::N);
	}

  private:
	result m_result;
};

// --------------------------------------------------------------------

/// Transaction class.
/// At construction, this class starts a transaction on the connection
/// and at exit an automatic rollback is done, unless commit was called.
class transaction final
{
  public:
	/// @brief Constructor
	transaction(connection &conn);

	/// @cond
	~transaction();

	transaction(const transaction &) = delete;
	transaction &operator=(const transaction &) = delete;
	/// @endcond

	/// @brief Execute the sql in @a query returning an iterable result
	result exec(std::string query);

	/// @brief Execute the sql in @a query returning an iterable result.
	/// Updates @a tail with what remains after the first statement in @a query
	result exec(std::string query, std::string &tail);

	/// Execute the sql in @a sql and return the result as a stream
	template <typename... Ts>
	cql_iterator_proxy<Ts...> stream(const std::string &sql)
	{
		return cql_iterator_proxy<Ts...>{ exec(sql) };
	}

	/// Commit the result of the operations
	void commit();

	/// Rollback the result of the operations, the underlying data is
	/// restored to the state before the construction of this transaction.
	void rollback();

  private:
	connection &m_conn;
	bool m_transaction_active = false;
};

// --------------------------------------------------------------------

/// This connection class creates a SQLite environment with the data in
/// the provided datablock as tables.
class connection final
{
  public:
	/// @brief Constructor
	connection(datablock &db);

	/// @brief Destructor
	~connection();

	friend class transaction;

	/// @brief Return true if the string @a sql contains a complete statement.
	[[nodiscard]] bool is_complete_statement(const std::string &sql) const;

	/// @brief Execute the sql in @a query returning an iterable result
	result exec(std::string query);

	/// @brief Execute the sql in @a query returning an iterable result.
	/// Updates @a tail with what remains after the first statement in @a query
	result exec(std::string query, std::string &tail);

	/// @brief Return true if the underlying data was modified by any query.
	[[nodiscard]] bool is_modified() const;

  private:
	struct connection_impl *m_impl;
};

} // namespace cif::cql