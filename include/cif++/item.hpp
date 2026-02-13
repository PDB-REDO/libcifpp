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

#include "cif++/text.hpp"

#include <algorithm>
#include <cassert>
#include <catch2/catch_tostring.hpp>
#include <cctype>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>

/** \file item.hpp
 *
 * This file contains the declaration of item but also the item_value and item_handle
 * These handle the storage of and access to the data for a single data item.
 */

namespace cif
{

class category;
class row;

// --------------------------------------------------------------------
/** @brief item is a transient class that is used to pass data into rows
 * but it also takes care of formatting data.
 *
 *
 *
 * The class cif::item is often used implicitly when creating a row in a category
 * using the emplace function.
 *
 * @code{.cpp}
 * cif::category cat("my-cat");
 * cat.emplace({
 *   { "item-1", 1 },                             // <- stores an item with value 1
 *   { "item-2", 1.0, 2 },                        // <- stores an item with value 1.00
 *   { "item-3", std::optional<int>() },          // <- stores an item with value ?
 *   { "item-4", std::make_optional<int>(42) },   // <- stores an item with value 42
 *   { "item-5" }                                 // <- stores an item with value .
 * });
 *
 * std::cout << cat << '\n';
 * @endcode
 *
 * Will result in:
 *
 * @code{.txt}
 * _my-cat.item-1 1
 * _my-cat.item-2 1.00
 * _my-cat.item-3 ?
 * _my-cat.item-4 42
 * _my-cat.item-5 .
 * @endcode
 */

enum class item_value_type
{
	INT,
	FLOAT,
	TEXT,

	INAPPLICABLE,
	MISSING
};

template <typename T>
concept IntegralType = (std::is_integral_v<std::remove_cvref_t<T>>);

template <typename T>
concept FloatType = std::is_floating_point_v<std::remove_cvref_t<T>>;

template <typename T>
concept StringType = (std::is_assignable_v<std::string, T> and not std::is_integral_v<T> and not std::is_floating_point_v<T>);

// --------------------------------------------------------------------

/// \cond
template <typename T>
inline constexpr bool is_optional_v = false;
template <typename T>
inline constexpr bool is_optional_v<std::optional<T>> = true;
/// \endcond

class item_value
{
  public:
	item_value() noexcept
	{
		m_data.m_type = item_value_type::MISSING;
	}

	item_value(item_value_type type) noexcept
		: m_data(type)
	{
	}

	item_value(const item_value &rhs)
	{
		m_data.m_type = rhs.m_data.m_type;
		switch (m_data.m_type)
		{
			case item_value_type::INT:
				m_data.m_value = rhs.m_data.m_value.m_integer;
				break;
			case item_value_type::FLOAT:
				m_data.m_len = rhs.m_data.m_len;
				m_data.m_value = rhs.m_data.m_value.m_float;
				break;
			case item_value_type::TEXT:
				m_data.m_len = rhs.m_data.m_len;
				m_data.m_value = rhs.m_data.sv();
				break;
			default: break;
		}
	}

	item_value(std::nullptr_t)
	{
		m_data.m_type = item_value_type::MISSING;
	}

	item_value(std::string_view s)
	{
		if (s == ".")
			m_data.m_type = item_value_type::INAPPLICABLE;
		else if (s == "?")
			m_data.m_type = item_value_type::MISSING;
		else
		{
			m_data.m_type = item_value_type::TEXT;
			m_data.m_len = s.length();
			m_data.m_value = s;
		}
	}

	template <size_t N>
	item_value(const char(s)[N])
		: item_value(std::string_view{ s, N })
	{
	}

	item_value(const char *s)
		: item_value(std::string_view{ s })
	{
	}

	item_value(const std::string &s)
		: item_value(std::string_view{ s })
	{
	}

	template <IntegralType T>
	item_value(T v)
	{
		m_data.m_type = item_value_type::INT;
		m_data.m_value = static_cast<int64_t>(v);
	}

	template <FloatType T>
	item_value(T v, int precision = 0)
	{
		m_data.m_type = item_value_type::FLOAT;
		m_data.m_value = static_cast<double>(v);
		m_data.m_len = precision;
	}

	template <typename T>
	item_value(std::optional<T> v)
	{
		if (v.has_value())
		{
			item_value iv{ *v };
			swap(*this, iv);
		}
		else
			m_data.m_type = item_value_type::MISSING;
	}

	item_value(item_value &&rhs) noexcept
	{
		swap(*this, rhs);
	}

	item_value &operator=(item_value rhs) noexcept
	{
		swap(*this, rhs);
		return *this;
	}

	// --------------------------------------------------------------------

	[[nodiscard]] constexpr bool is_inapplicable() const noexcept { return m_data.m_type == item_value_type::INAPPLICABLE; }
	[[nodiscard]] constexpr bool is_missing() const noexcept { return m_data.m_type == item_value_type::MISSING; }
	[[nodiscard]] constexpr bool is_null() const noexcept { return is_inapplicable() or is_missing(); }

	[[nodiscard]] constexpr bool is_string() const noexcept { return m_data.m_type == item_value_type::TEXT; }

	[[nodiscard]] constexpr bool is_number_int() const noexcept { return m_data.m_type == item_value_type::INT; }
	[[nodiscard]] constexpr bool is_number_float() const noexcept { return m_data.m_type == item_value_type::FLOAT; }
	[[nodiscard]] constexpr bool is_number() const noexcept { return is_number_int() or is_number_float(); }

	[[nodiscard]] constexpr item_value_type type() const { return m_data.m_type; }

	explicit operator bool() const noexcept
	{
		bool result;
		switch (m_data.m_type)
		{
			case item_value_type::INT: result = m_data.m_value.m_integer != 0; break;
			case item_value_type::FLOAT: result = m_data.m_value.m_float != 0; break;
			case item_value_type::TEXT: result = m_data.m_len != 0; break;
			case item_value_type::INAPPLICABLE:
			case item_value_type::MISSING: result = false; break;
		}
		return result;
	}

	[[nodiscard]] bool empty() const noexcept
	{
		switch (m_data.m_type)
		{
			case item_value_type::INAPPLICABLE:
			case item_value_type::MISSING:
				return true;

			case item_value_type::TEXT:
				return m_data.sv().empty();

			default:
				return false;
		}
	}

	// --------------------------------------------------------------------

	template <StringType T>
	[[nodiscard]] inline std::string get() const
	{
		return str();
	}

	template <IntegralType T>
	[[nodiscard]] std::remove_cvref_t<T> get() const
	{
		static_assert(not std::is_same_v<std::remove_cvref_t<T>, bool>, "bool is no longer supported");

		switch (m_data.m_type)
		{
			case item_value_type::INT:
				return m_data.m_value.m_integer;
			case item_value_type::FLOAT:
				return m_data.m_value.m_float;
			case item_value_type::TEXT:
			{
				auto sv = m_data.sv();
				int64_t v;
				auto &&[ptr, ec] = std::from_chars(sv.data(), sv.data() + sv.length(), v);
				if (ec != std::errc{})
					throw std::system_error(std::make_error_code(ec));
				if (ptr != sv.data() + sv.length())
					throw std::invalid_argument("String value does not contain only an integer");

				return v;
			}
			default:
				return not empty();
		}
	}

	template <FloatType T>
	[[nodiscard]] std::remove_cvref_t<T> get() const
	{
		switch (m_data.m_type)
		{
			case item_value_type::INT:
				return m_data.m_value.m_integer;
			case item_value_type::FLOAT:
				return m_data.m_value.m_float;
			case item_value_type::TEXT:
			{
				auto sv = m_data.sv();
				double v;
				auto &&[ptr, ec] = std::from_chars(sv.data(), sv.data() + sv.length(), v);
				if (ec != std::errc{})
					throw std::system_error(std::make_error_code(ec));
				if (ptr != sv.data() + sv.length())
					throw std::invalid_argument("String value does not contain only a floating point number");
				return v;
			}
			default:
				return not empty();
		}
	}

	template <typename T>
		requires is_optional_v<T>
	[[nodiscard]] auto get() const
	{
		switch (m_data.m_type)
		{
			case item_value_type::INAPPLICABLE:
			case item_value_type::MISSING:
				return T{};

			default:
			{
				auto v = get<typename T::value_type>();
				return T{ v };
			}
		}
	}

	[[nodiscard]] std::string str() const;

	// --------------------------------------------------------------------

	friend void swap(item_value &a, item_value &b) noexcept
	{
		std::swap(a.m_data.m_type, b.m_data.m_type);
		std::swap(a.m_data.m_len, b.m_data.m_len);
		std::swap(a.m_data.m_value, b.m_data.m_value);
	}

	// --------------------------------------------------------------------
	// std::partial_ordering operator<=>(const item_value &rhs) const
	// {
	// 	if (m_data.m_type == rhs.m_data.m_type)
	// 	{
	// 		switch (m_data.m_type)
	// 		{
	// 			case item_value_type::INT: return m_data.m_value.m_integer <=> rhs.m_data.m_value.m_integer;
	// 			case item_value_type::FLOAT: return m_data.m_value.m_float <=> rhs.m_data.m_value.m_float;
	// 			case item_value_type::TEXT: return m_data.sv() <=> rhs.m_data.sv();
	// 			case item_value_type::MISSING:
	// 			case item_value_type::EMPTY: return std::strong_ordering::equivalent;
	// 		}
	// 	}
	// 	else
	// 		return m_data.m_type <=> rhs.m_data.m_type;
	// }

	bool operator==(const item_value &rhs) const
	{
		if (m_data.m_type == rhs.m_data.m_type)
		{
			switch (m_data.m_type)
			{
				case item_value_type::INT: return m_data.m_value.m_integer == rhs.m_data.m_value.m_integer;
				case item_value_type::FLOAT: return m_data.m_value.m_float == rhs.m_data.m_value.m_float;
				case item_value_type::TEXT: return m_data.sv() == rhs.m_data.sv();
				case item_value_type::INAPPLICABLE:
				case item_value_type::MISSING: return true;
			}
		}

		return false;
	}

	[[nodiscard]] int compare(const item_value &b, bool ignore_case = false) const noexcept;

	friend std::ostream &operator<<(std::ostream &os, const item_value &v);

  private:
	union value
	{
		int64_t m_integer{};
		double m_float;
		char m_local_str[8];
		char *m_str;

		value() = default;

		value(int64_t v)
			: m_integer(v)
		{
		}

		value(double v)
			: m_float(v)
		{
		}

		value(std::string_view s)
		{
			if (s.length() >= sizeof(m_local_str))
			{
				m_str = new char[s.length() + 1];
				std::copy(s.data(), s.data() + s.length(), m_str);
				m_str[s.length()] = 0;
			}
			else
				std::memcpy(m_local_str, s.data(), s.length() + 1);
		}

		value(item_value_type t)
		{
			m_integer = 0;
		}

		void destroy(item_value_type t, size_t len)
		{
			if (t == item_value_type::TEXT and len >= sizeof(m_local_str))
				delete[] m_str;
		}
	};

	struct data
	{
		item_value_type m_type = item_value_type::MISSING;
		uint32_t m_len{};
		value m_value{};

		data(item_value_type t)
			: m_type(t)
			, m_value(t)
		{
		}

		data() noexcept = default;
		data(data &&rhs) noexcept
		{
			std::swap(m_type, rhs.m_type);
			std::swap(m_len, rhs.m_len);
			std::swap(m_value, rhs.m_value);
		}

		data(const data &) noexcept = delete;
		data &operator=(data &&) noexcept = delete;
		data &operator=(const data &) noexcept = delete;

		~data()
		{
			m_value.destroy(m_type, m_len);
		}

		[[nodiscard]] std::string_view sv() const noexcept
		{
			assert(m_type == item_value_type::TEXT);
			return m_type == item_value_type::TEXT ? std::string_view(m_len >= sizeof(m_value.m_local_str) ? m_value.m_str : m_value.m_local_str, m_len) : std::string_view{};
		}

		[[nodiscard]] const char *c_str() const noexcept
		{
			assert(m_type == item_value_type::TEXT);
			return m_type == item_value_type::TEXT ? (m_len >= sizeof(m_value.m_local_str) ? m_value.m_str : m_value.m_local_str) : nullptr;
		}
	} m_data{};
};

static_assert(sizeof(item_value) == 16, "item_value should be 16 bytes");

class item
{
  public:
	/// \brief Default constructor, empty item
	item() = default;

	/// \brief constructor for an item with name \a name and as
	/// content the character '.', i.e. an inapplicable value.
	item(std::string name)
		: m_name(std::move(name))
		, m_value(item_value_type::MISSING)
	{
	}

	item(std::string name, item_value value)
		: m_name(std::move(name))
		, m_value(std::move(value))
	{
	}

	/** @cond */
	item(const item &rhs) = default;

	item(item &&rhs)
	{
		swap(*this, rhs);
	}

	item &operator=(item rhs) noexcept
	{
		swap(*this, rhs);
		return *this;
	}
	/** @endcond */

	friend void swap(item &a, item &b) noexcept
	{
		std::swap(a.m_name, b.m_name);
		std::swap(a.m_value, b.m_value);
	}

	[[nodiscard]] const std::string &name() const { return m_name; }    ///< Return the name of the item
	[[nodiscard]] const item_value &value() const & { return m_value; } ///< Return the value of the item
	item_value &value() & { return m_value; }                           ///< Return the value of the item

	/// \brief replace the content of the stored value with \a v
	void value(item_value v) { m_value = std::move(v); }

	/// \brief empty means either null or unknown
	[[nodiscard]] bool empty() const { return m_value.empty(); }

	/// \brief returns true if the item contains '.' or '?'
	[[nodiscard]] bool is_null() const { return m_value.is_null(); }

	/// \brief returns true if the item contains '?'
	[[nodiscard]] bool is_unknown() const { return m_value.is_missing(); }

	// /// \brief the length of the value string
	// std::size_t length() const { return m_value.length(); }

	/// \brief support for structured binding
	template <std::size_t N>
	decltype(auto) get() const
	{
		if constexpr (N == 0)
			return name();
		else if constexpr (N == 1)
			return value();
	}

	// auto operator<=>(const item &rhs) const = default;

  private:
	std::string m_name;
	item_value m_value;
};

// --------------------------------------------------------------------
/// \brief This is item_handle, it is used to access the data stored in
/// item_value's in rows

struct item_handle
{
  public:
	item_handle() = delete;

	/**
	 * @brief Assign value @a value to the item referenced
	 *
	 * @tparam T Type of the value
	 * @param value The value
	 * @return reference to this item_handle
	 */
	item_handle &operator=(item_value value)
	{
		set(std::move(value), true);
		return *this;
	}

	[[nodiscard]] item_value &value();
	[[nodiscard]] const item_value &value() const;

	[[nodiscard]] constexpr bool is_inapplicable() const noexcept { return value().type() == item_value_type::INAPPLICABLE; }
	[[nodiscard]] constexpr bool is_missing() const noexcept { return value().type() == item_value_type::MISSING; }
	[[nodiscard]] constexpr bool is_null() const noexcept { return is_inapplicable() or is_missing(); }

	[[nodiscard]] constexpr bool is_string() const noexcept { return value().type() == item_value_type::TEXT; }

	[[nodiscard]] constexpr bool is_number_int() const noexcept { return value().type() == item_value_type::INT; }
	[[nodiscard]] constexpr bool is_number_float() const noexcept { return value().type() == item_value_type::FLOAT; }
	[[nodiscard]] constexpr bool is_number() const noexcept { return is_number_int() or is_number_float(); }

	[[nodiscard]] auto type() const { return value().type(); }

	template <typename T>
	[[nodiscard]] auto get() const
	{
		if (empty())
			return T{};
		else
			return value().template get<T>();
	}

	template <typename T>
	[[nodiscard]] auto as() const
	{
		if (empty())
			return T{};
		else
			return value().template get<T>();
	}

	[[nodiscard]] auto str() const
	{
		return value().str();
	}

	/** Swap contents of @a a and @a b */
	friend void swap(item_handle &a, item_handle &b) noexcept;

	/** Return the contents of this item as type @tparam T or, if not
	 * set, use @a dv as the default value.
	 */
	template <typename T>
	[[nodiscard]] auto value_or(const T &dv) const
	{
		return empty() ? dv : this->get<T>();
	}

	/**
	 * @brief Compare the contents of this item with value @a value
	 * optionally ignoring character case, if @a icase is true.
	 * Returns 0 if both are equal, -1 if this sorts before @a value
	 * and 1 if this sorts after @a value
	 *
	 * @tparam T Type of the value @a value
	 * @param value The value to compare with
	 * @param icase Flag indicating if we should compare character case sensitive
	 * @return -1, 0 or 1
	 */

	[[nodiscard]] int compare(const item_value &value, bool icase = true) const noexcept
	{
		return this->value().compare(value, icase);
	}

	[[nodiscard]] int compare(const item_handle &value, bool icase = true) const noexcept
	{
		return compare(value.value(), icase);
	}

	/**
	 * @brief Compare the value contained with the value @a value and
	 * return true if both are equal.
	 */
	[[nodiscard]] bool operator==(const item_value &value) const noexcept
	{
		// TODO: icase or not icase?
		return this->value().compare(value) != 0;
	}

	// We may not have C++20 yet...

	/**
	 * @brief Compare the value contained with the value @a value and
	 * return true if both are not equal.
	 */
	template <typename T>
	[[nodiscard]] bool operator!=(const T &value) const noexcept
	{
		return not operator==(value);
	}

	/**
	 * @brief Returns true if the content string is empty or
	 * only contains '.' meaning null or '?' meaning unknown
	 * in a mmCIF context
	 */
	[[nodiscard]] bool empty() const;

	/** Easy way to test for an empty item */
	explicit operator bool() const { return not empty(); }

	/** Return a std::string_view for the contents */
	[[nodiscard]] std::string_view text_() const;

	/**
	 * @brief Construct a new item handle object
	 *
	 * @param item Item index
	 * @param row Reference to the row
	 */
	item_handle(category &cat, row &row, uint16_t item_ix)
		: m_category(cat)
		, m_row(row)
		, m_item_ix(item_ix)
	{
	}

  private:
	category &m_category;
	row &m_row;
	uint16_t m_item_ix;

	friend class parser;

	void set(item_value value, bool updateLinked);
};

struct const_item_handle
{
  public:
	const_item_handle() = delete;

	[[nodiscard]] const item_value &value() const;

	[[nodiscard]] constexpr bool is_inapplicable() const noexcept { return value().type() == item_value_type::INAPPLICABLE; }
	[[nodiscard]] constexpr bool is_missing() const noexcept { return value().type() == item_value_type::MISSING; }
	[[nodiscard]] constexpr bool is_null() const noexcept { return is_inapplicable() or is_missing(); }

	[[nodiscard]] constexpr bool is_string() const noexcept { return value().type() == item_value_type::TEXT; }

	[[nodiscard]] constexpr bool is_number_int() const noexcept { return value().type() == item_value_type::INT; }
	[[nodiscard]] constexpr bool is_number_float() const noexcept { return value().type() == item_value_type::FLOAT; }
	[[nodiscard]] constexpr bool is_number() const noexcept { return is_number_int() or is_number_float(); }

	[[nodiscard]] auto type() const { return value().type(); }

	template <typename T>
	[[nodiscard]] auto get() const
	{
		if (empty())
			return T{};
		else
			return value().template get<T>();
	}

	template <typename T>
	[[nodiscard]] auto as() const
	{
		if (empty())
			return T{};
		else
			return value().template get<T>();
	}

	[[nodiscard]] auto str() const
	{
		return value().str();
	}

	/** Return the contents of this item as type @tparam T or, if not
	 * set, use @a dv as the default value.
	 */
	template <typename T>
	[[nodiscard]] auto value_or(const T &dv) const
	{
		return empty() ? dv : this->get<T>();
	}

	/**
	 * @brief Compare the contents of this item with value @a value
	 * optionally ignoring character case, if @a icase is true.
	 * Returns 0 if both are equal, -1 if this sorts before @a value
	 * and 1 if this sorts after @a value
	 *
	 * @tparam T Type of the value @a value
	 * @param value The value to compare with
	 * @param icase Flag indicating if we should compare character case sensitive
	 * @return -1, 0 or 1
	 */

	[[nodiscard]] int compare(const item_value &value, bool icase = true) const noexcept
	{
		return this->value().compare(value, icase);
	}

	[[nodiscard]] int compare(const const_item_handle &value, bool icase = true) const noexcept
	{
		if (empty() and value.empty())
			return 0;
		else if (empty())
			return -1;
		else if (value.empty())
			return 1;
		else
			return compare(value.value(), icase);
	}

	/**
	 * @brief Compare the value contained with the value @a value and
	 * return true if both are equal.
	 */
	[[nodiscard]] bool operator==(const item_value &value) const noexcept
	{
		// TODO: icase or not icase?
		return this->value().compare(value) != 0;
	}

	// We may not have C++20 yet...

	/**
	 * @brief Compare the value contained with the value @a value and
	 * return true if both are not equal.
	 */
	template <typename T>
	[[nodiscard]] bool operator!=(const T &value) const noexcept
	{
		return not operator==(value);
	}

	/**
	 * @brief Returns true if the content string is empty or
	 * only contains '.' meaning null or '?' meaning unknown
	 * in a mmCIF context
	 */
	[[nodiscard]] bool empty() const;

	/** Easy way to test for an empty item */
	explicit operator bool() const { return not empty(); }

	/**
	 * @brief Construct a new item handle object
	 *
	 * @param item Item index
	 * @param row Reference to the row
	 */
	const_item_handle(const category &cat, const row &row, uint16_t item_ix)
		: m_category(cat)
		, m_row(row)
		, m_item_ix(item_ix)
	{
	}

  private:
	const category &m_category;
	const row &m_row;
	uint16_t m_item_ix;
};

} // namespace cif

namespace std
{

/** @cond */

template <>
struct tuple_size<::cif::item>
	: public std::integral_constant<std::size_t, 2>
{
};

template <>
struct tuple_element<0, ::cif::item>
{
	using type = decltype(std::declval<::cif::item>().name());
};

template <>
struct tuple_element<1, ::cif::item>
{
	using type = decltype(std::declval<::cif::item>().value());
};

} // namespace std
