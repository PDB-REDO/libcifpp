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

#include "cif++/exports.hpp"
#include "cif++/forward_decl.hpp"
#include "cif++/text.hpp"
#include "cif++/utilities.hpp"

#include <cassert>
#include <charconv>
#include <compare>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <system_error>
#include <utility>

/** \file item.hpp
 *
 * This file contains the declaration of item but also the item_value and item_handle
 * These handle the storage of and access to the data for a single data item.
 */

namespace cif
{

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
	BOOLEAN,
	INT,
	FLOAT,
	TEXT,
	MISSING,
	EMPTY // This is the real NULL in SQL terms
};

template <typename T>
concept BooleanType = std::is_same_v<std::remove_cvref_t<T>, bool>;

template <typename T>
concept IntegralType = (std::is_integral_v<std::remove_cvref_t<T>> and not std::is_same_v<std::remove_cvref_t<T>, bool>);

template <typename T>
concept FloatType = std::is_floating_point_v<std::remove_cvref_t<T>>;

template <typename T>
concept StringType = (std::is_assignable_v<std::string, T> and not std::is_integral_v<T> and not std::is_floating_point_v<T>);

// --------------------------------------------------------------------

/// \cond
template <typename _Tp>
inline constexpr bool is_optional_v = false;
template <typename _Tp>
inline constexpr bool is_optional_v<std::optional<_Tp>> = true;
/// \endcond



class item_value
{
  public:
	item_value() noexcept
	{
		m_data.m_type = item_value_type::EMPTY;
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
			case item_value_type::BOOLEAN: m_data.m_value = rhs.m_data.m_value.m_boolean; break;
			case item_value_type::INT: m_data.m_value = rhs.m_data.m_value.m_integer; break;
			case item_value_type::FLOAT: m_data.m_value = rhs.m_data.m_value.m_float; break;
			case item_value_type::TEXT:
				m_data.m_len = rhs.m_data.m_len;
				m_data.m_value = rhs.m_data.sv();
				break;
			default: break;
		}
	}

	item_value(std::nullptr_t)
	{
		m_data.m_type = item_value_type::EMPTY;
	}

	template <BooleanType T>
	item_value(T v)
	{
		m_data.m_type = item_value_type::BOOLEAN;
		m_data.m_value = v;
	}

	item_value(std::string_view s)
	{
		m_data.m_type = item_value_type::TEXT;
		m_data.m_len = s.length();
		m_data.m_value = s;
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
			item_value iv{ *v  };
			swap(*this, iv);
		}
		else
			m_data.m_type = item_value_type::EMPTY;
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

	constexpr bool is_null() const noexcept { return m_data.m_type == item_value_type::MISSING; }
	constexpr bool is_empty() const noexcept { return m_data.m_type == item_value_type::EMPTY; }
	constexpr bool is_string() const noexcept { return m_data.m_type == item_value_type::TEXT; }
	constexpr bool is_number() const noexcept { return is_number_int() or is_number_float(); }
	constexpr bool is_number_int() const noexcept { return m_data.m_type == item_value_type::INT; }
	constexpr bool is_number_float() const noexcept { return m_data.m_type == item_value_type::FLOAT; }
	constexpr bool is_boolean() const noexcept { return m_data.m_type == item_value_type::BOOLEAN; }

	constexpr item_value_type type() const { return m_data.m_type; }

	explicit operator bool() const noexcept
	{
		bool result;
		switch (m_data.m_type)
		{
			case item_value_type::BOOLEAN: result = m_data.m_value.m_boolean; break;
			case item_value_type::INT: result = m_data.m_value.m_integer != 0; break;
			case item_value_type::FLOAT: result = m_data.m_value.m_float != 0; break;
			case item_value_type::TEXT: result = m_data.m_len != 0; break;
			case item_value_type::MISSING:
			case item_value_type::EMPTY: result = false; break;
		}
		return result;
	}

	bool empty() const noexcept
	{
		switch (m_data.m_type)
		{
			case item_value_type::MISSING:
			case item_value_type::EMPTY:
				return true;

			case item_value_type::TEXT:
				return m_data.sv().empty();

			default:
				return false;
		}
	}

	// --------------------------------------------------------------------

	template <StringType T>
	inline std::string get() const
	{
		switch (m_data.m_type)
		{
			case item_value_type::EMPTY:
			case item_value_type::MISSING:
				return "";

			case item_value_type::TEXT:
				return std::string{ m_data.sv() };

			case cif::item_value_type::BOOLEAN:
				return m_data.m_value.m_boolean ? "y" : "n";

			default:
			{
				char b[32];

				const auto &[ptr, ec] =
					m_data.m_type == item_value_type::INT ? std::to_chars(b, b + sizeof(b), m_data.m_value.m_integer)
					: m_data.m_len                        ? std::to_chars(b, b + sizeof(b), m_data.m_value.m_float, std::chars_format::fixed, m_data.m_len)
														  : std::to_chars(b, b + sizeof(b), m_data.m_value.m_float, std::chars_format::general);

				if (ec != std::errc{})
					throw std::system_error(std::make_error_code(ec));

				return std::string{ b, ptr };
			}
		}
	}

	template <IntegralType T>
	std::remove_cvref_t<T> get() const
	{
		switch (m_data.m_type)
		{
			case cif::item_value_type::BOOLEAN:
				return m_data.m_value.m_boolean;
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
	std::remove_cvref_t<T> get() const
	{
		switch (m_data.m_type)
		{
			case cif::item_value_type::BOOLEAN:
				return m_data.m_value.m_boolean;
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

	template <BooleanType T>
	std::remove_cvref_t<T> get() const
	{
		switch (m_data.m_type)
		{
			case cif::item_value_type::BOOLEAN:
				return m_data.m_value.m_boolean;
			case item_value_type::INT:
				return m_data.m_value.m_integer != 0;
			case item_value_type::FLOAT:
				return m_data.m_value.m_float != 0.;
			case item_value_type::TEXT:
				return iequals(m_data.sv(), "y") or iequals(m_data.sv(), "yes") or iequals(m_data.sv(), "true");
			default:
				return not empty();
		}
	}

	template <typename T>
		requires is_optional_v<T>
	auto get() const
	{
		using value_type = T::value_type;

		switch (m_data.m_type)
		{
			case item_value_type::MISSING:
			case item_value_type::EMPTY:
				return T{};

			default:
				value_type v = get<value_type>();
				return T{ v };
		}
	}

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
	// 			case item_value_type::BOOLEAN: return m_data.m_value.m_boolean <=> rhs.m_data.m_value.m_boolean;
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
				case item_value_type::BOOLEAN: return m_data.m_value.m_boolean == rhs.m_data.m_value.m_boolean;
				case item_value_type::INT: return m_data.m_value.m_integer == rhs.m_data.m_value.m_integer;
				case item_value_type::FLOAT: return m_data.m_value.m_float == rhs.m_data.m_value.m_float;
				case item_value_type::TEXT: return m_data.sv() == rhs.m_data.sv();
				case item_value_type::MISSING:
				case item_value_type::EMPTY: return true;
			}
		}

		return false;
	}

	int compare(const item_value &b, bool ignore_case = false) const noexcept;

	friend std::ostream &operator<<(std::ostream &os, const item_value &v);

  private:
	union value
	{
		bool m_boolean;
		int64_t m_integer;
		double m_float;
		char m_local_str[8];
		char *m_str;

		value()
			: m_integer(0)
		{
		}

		value(bool v)
			: m_boolean(v)
		{
		}

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
				memcpy(m_local_str, s.data(), s.length() + 1);
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
		item_value_type m_type = item_value_type::EMPTY;
		uint32_t m_len{};
		value m_value{};

		data(item_value_type t)
			: m_type(t)
			, m_value(t)
		{
		}

		data() noexcept
		{
		}
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

		std::string_view sv() const noexcept
		{
			return m_type == item_value_type::TEXT ? std::string_view(m_len >= sizeof(m_value.m_local_str) ? m_value.m_str : m_value.m_local_str, m_len) : std::string_view{};
		}

		const char *c_str() const noexcept
		{
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
		, m_value(item_value_type::EMPTY)
	{
	}

	item(std::string name, item_value value)
		: m_name(std::move(name))
		, m_value(std::move(value))
	{
	}

	// /// \brief constructor for an item with name \a name and as
	// /// content the character '.', i.e. an inapplicable value.
	// item(std::string_view name, std::nullptr_t)
	// 	: m_name(name)
	// 	, m_value(item_value_type::EMPTY)
	// {
	// }

	// /// \brief constructor for an item with name \a name and as
	// /// content a single character string with content \a value
	// item(std::string_view name, char value)
	// 	: m_name(name)
	// 	, m_value(std::string_view{ &value, 1 })
	// {
	// }

	// /// \brief constructor for an item with name \a name and as
	// /// content the formatted floating point value \a value
	// template <FloatType T>
	// item(std::string_view name, T value)
	// 	: m_name(name)
	// 	, m_value(value)
	// {
	// }

	// /// \brief constructor for an item with name \a name and as
	// /// content the formatted floating point value \a value with
	// /// precision \a precision
	// template <FloatType T>
	// item(std::string_view name, T value, int precision)
	// 	: m_name(name)
	// 	, m_value(value, precision)
	// {
	// }

	// /// \brief constructor for an item with name \a name and as
	// /// content the formatted integral value \a value
	// template <IntegralType T>
	// item(const std::string_view name, T value)
	// 	: m_name(name)
	// 	, m_value(value)
	// {
	// }

	// // TODO: Perhaps introduce a real boolean type?
	// /// \brief constructor for an item with name \a name and as
	// /// content the formatted boolean value \a value
	// template <BooleanType T>
	// item(const std::string_view name, T value)
	// 	: m_name(name)
	// 	, m_value(value)
	// {
	// }

	// /// \brief constructor for an item with name \a name and as
	// /// content value \a value
	// item(const std::string_view name, std::string_view value)
	// 	: m_name(name)
	// 	, m_value(value)
	// {
	// }

	// /// \brief constructor for an item with name \a name and as
	// /// content the optional value \a value
	// template <typename T>
	// item(const std::string_view name, const std::optional<T> &value)
	// 	: m_name(name)
	// 	, m_value(item_value_type::MISSING)
	// {
	// 	if (value.has_value())
	// 		m_value = *value;
	// }

	// /// \brief constructor for an item with name \a name and as
	// /// content the formatted floating point value \a value with
	// /// precision \a precision
	// template <typename T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
	// item(std::string_view name, const std::optional<T> &value, int precision)
	// 	: m_name(name)
	// 	, m_value(item_value_type::MISSING)
	// {
	// 	if (value.has_value())
	// 		m_value = item_value(*value, precision);
	// }

	/** @cond */
	item(const item &rhs)
		: m_name(rhs.m_name)
		, m_value(rhs.m_value)
	{
	}

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

	const std::string &name() const { return m_name; }    ///< Return the name of the item
	const item_value &value() const & { return m_value; } ///< Return the value of the item
	item_value &value() & { return m_value; }             ///< Return the value of the item

	/// \brief replace the content of the stored value with \a v
	void value(item_value v) { m_value = std::move(v); }

	/// \brief empty means either null or unknown
	bool empty() const { return m_value.empty(); }

	/// \brief returns true if the item contains '.'
	bool is_null() const { return m_value.is_null(); }

	/// \brief returns true if the item contains '?'
	bool is_unknown() const { return m_value.is_empty(); }

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

// // --------------------------------------------------------------------
// /// \brief the internal storage for items in a category
// ///
// /// Internal storage, strictly forward linked list with minimal space
// /// requirements. Strings of size 7 or shorter are stored internally.
// /// Typically, more than 99% of the strings in an mmCIF file are less
// /// than 8 bytes in length.

// struct item_value
// {
// 	/** @cond */
// 	item_value() = default;
// 	/** @endcond */

// 	/// \brief constructor
// 	item_value(std::string_view text)
// 		: m_length(text.length())
// 		, m_storage(0)
// 	{
// 		if (m_length >= kBufferSize)
// 		{
// 			m_data = new char[m_length + 1];
// 			std::copy(text.begin(), text.end(), m_data);
// 			m_data[m_length] = 0;
// 		}
// 		else
// 		{
// 			std::copy(text.begin(), text.end(), m_local_data);
// 			m_local_data[m_length] = 0;
// 		}
// 	}

// 	/** @cond */
// 	item_value(item_value &&rhs) noexcept
// 		: m_length(std::exchange(rhs.m_length, 0))
// 		, m_storage(std::exchange(rhs.m_storage, 0))
// 	{
// 	}

// 	item_value &operator=(item_value &&rhs) noexcept
// 	{
// 		std::swap(m_length, rhs.m_length);
// 		std::swap(m_storage, rhs.m_storage);
// 		return *this;
// 	}

// 	~item_value()
// 	{
// 		if (m_length >= kBufferSize)
// 			delete[] m_data;
// 		m_storage = 0;
// 		m_length = 0;
// 	}

// 	item_value(const item_value &) = delete;
// 	item_value &operator=(const item_value &) = delete;
// 	/** @endcond */

// 	/** operator bool, allows easy checking for empty items */
// 	explicit operator bool() const
// 	{
// 		return m_length != 0;
// 	}

// 	std::size_t m_length = 0; ///< Length of the data
// 	union
// 	{
// 		char m_local_data[8]; ///< Storage area for small strings (strings smaller than kBufferSize)
// 		char *m_data;         ///< Pointer to a string stored in the heap
// 		uint64_t m_storage;   ///< Alternative storage of the data, used in move operations
// 	};

// 	/** The maximum length of locally stored strings */
// 	static constexpr std::size_t kBufferSize = sizeof(m_local_data);

// 	// By using std::string_view instead of c_str we obain a
// 	// nice performance gain since we avoid many calls to strlen.

// 	/** Return the content of the item as a std::string_view */
// 	constexpr inline std::string_view text() const
// 	{
// 		return { m_length >= kBufferSize ? m_data : m_local_data, m_length };
// 	}
// };

// // --------------------------------------------------------------------
// // Transient object to access stored data

// /// \brief This is item_handle, it is used to access the data stored in item_value.

// struct item_handle
// {
//   public:
// // 	/** @cond */
// // 	// conversion helper class
// // 	template <typename T, typename = void>
// // 	struct item_value_as;
// // 	/** @endcond */

// // 	/**
// // 	 * @brief Assign value @a value to the item referenced
// // 	 *
// // 	 * @tparam T Type of the value
// // 	 * @param value The value
// // 	 * @return reference to this item_handle
// // 	 */
// // 	template <typename T>
// // 	item_handle &operator=(const T &value)
// // 	{
// // 		assign_value(item{ "", value }.value());
// // 		return *this;
// // 	}

// // 	/**
// // 	 * @brief Assign value @a value to the item referenced
// // 	 *
// // 	 * @tparam T Type of the value
// // 	 * @param value The value
// // 	 * @return reference to this item_handle
// // 	 */
// // 	template <typename T>
// // 	item_handle &operator=(T &&value)
// // 	{
// // 		assign_value(item{ "", std::forward<T>(value) }.value());
// // 		return *this;
// // 	}

// // 	/**
// // 	 * @brief Assign value @a value to the item referenced
// // 	 *
// // 	 * @tparam T Type of the value
// // 	 * @param value The value
// // 	 * @return reference to this item_handle
// // 	 */
// // 	template <std::size_t N>
// // 	item_handle &operator=(const char (&value)[N])
// // 	{
// // 		assign_value(item{ "", std::move(value) }.value());
// // 		return *this;
// // 	}

// // 	/**
// // 	 * @brief A method with a variable number of arguments that will be concatenated and
// // 	 * assigned as a string. Use it like this:
// // 	 *
// // 	 * @code{.cpp}
// // 	 * cif::item_handle ih;
// // 	 * is.os("The result of ", 1, " * ", 42, " is of course ", 42);
// // 	 * @endcode
// // 	 *
// // 	 * And the content will then be `The result of 1 * 42 is of course 42`.
// // 	 *
// // 	 * @tparam Ts Types of the parameters
// // 	 * @param v The parameters to concatenate
// // 	 */
// // 	template <typename... Ts>
// // 	void os(const Ts &...v)
// // 	{
// // 		std::ostringstream ss;
// // 		((ss << v), ...);
// // 		this->operator=(ss.str());
// // 	}

// // 	/** Swap contents of this and @a b */
// // 	void swap(item_handle &b);

// // 	/** Return the contents of this item as type @tparam T */
// // 	template <typename T = std::string>
// // 	auto as() const -> T
// // 	{
// // 		using value_type = std::remove_cv_t<std::remove_reference_t<T>>;
// // 		return item_value_as<value_type>::convert(*this);
// // 	}

// // 	/** Return the contents of this item as type @tparam T or, if not
// // 	 * set, use @a dv as the default value.
// // 	 */
// // 	template <typename T>
// // 	auto value_or(const T &dv) const
// // 	{
// // 		return empty() ? dv : this->as<T>();
// // 	}

// // 	/**
// // 	 * @brief Compare the contents of this item with value @a value
// // 	 * optionally ignoring character case, if @a icase is true.
// // 	 * Returns 0 if both are equal, -1 if this sorts before @a value
// // 	 * and 1 if this sorts after @a value
// // 	 *
// // 	 * @tparam T Type of the value @a value
// // 	 * @param value The value to compare with
// // 	 * @param icase Flag indicating if we should compare character case sensitive
// // 	 * @return -1, 0 or 1
// // 	 */
// // 	template <typename T>
// // 	int compare(const T &value, bool icase = true) const
// // 	{
// // 		return item_value_as<T>::compare(*this, value, icase);
// // 	}

// // 	/**
// // 	 * @brief Compare the value contained with the value @a value and
// // 	 * return true if both are equal.
// // 	 */
// // 	template <typename T>
// // 	bool operator==(const T &value) const
// // 	{
// // 		// TODO: icase or not icase?
// // 		return item_value_as<T>::compare(*this, value, true) == 0;
// // 	}

// // 	// We may not have C++20 yet...

// // 	/**
// // 	 * @brief Compare the value contained with the value @a value and
// // 	 * return true if both are not equal.
// // 	 */
// // 	template <typename T>
// // 	bool operator!=(const T &value) const
// // 	{
// // 		return not operator==(value);
// // 	}

// // 	/**
// // 	 * @brief Returns true if the content string is empty or
// // 	 * only contains '.' meaning null or '?' meaning unknown
// // 	 * in a mmCIF context
// // 	 */
// // 	bool empty() const
// // 	{
// // 		auto txt = text();
// // 		return txt.empty() or (txt.length() == 1 and (txt.front() == '.' or txt.front() == '?'));
// // 	}

// // 	/** Easy way to test for an empty item */
// // 	explicit operator bool() const { return not empty(); }

// // 	/// is_null return true if the item contains '.'
// // 	bool is_null() const
// // 	{
// // 		auto txt = text();
// // 		return txt.length() == 1 and txt.front() == '.';
// // 	}

// // 	/// is_unknown returns true if the item contains '?'
// // 	bool is_unknown() const
// // 	{
// // 		auto txt = text();
// // 		return txt.length() == 1 and txt.front() == '?';
// // 	}

// // 	/** Return a std::string_view for the contents */
// // 	std::string_view text() const;

// 	/**
// 	 * @brief Construct a new item handle object
// 	 *
// 	 * @param item Item index
// 	 * @param row Reference to the row
// 	 */
// 	item_handle(uint16_t item, row_handle &row)
// 		: m_item_ix(item)
// 		, m_row_handle(row)
// 	{
// 	}

// 	/** A variable holding an empty item */
// 	CIFPP_EXPORT static const item_handle s_null_item;

// // 	/** friend to swap two item handles */
// // 	friend void swap(item_handle a, item_handle b)
// // 	{
// // 		a.swap(b);
// // 	}

//   private:
// 	item_handle();

// 	uint16_t m_item_ix;
// 	row_handle &m_row_handle;

// 	// void assign_value(std::string_view value);
// };

// // So sad that older gcc implementations of from_chars did not support floats yet...

// /** @cond */
// template <typename T>
// struct item_handle::item_value_as<T, std::enable_if_t<std::is_arithmetic_v<T> and not std::is_same_v<T, bool>>>
// {
// 	using value_type = std::remove_reference_t<std::remove_cv_t<T>>;

// 	static value_type convert(const item_handle &ref)
// 	{
// 		value_type result = {};

// 		if (not ref.empty())
// 		{
// 			auto txt = ref.text();

// 			auto b = txt.data();
// 			auto e = txt.data() + txt.size();

// 			std::from_chars_result r = (b + 1 < e and *b == '+' and std::isdigit(b[1])) //
// 			                               ? from_chars(b + 1, e, result)
// 			                               : from_chars(b, e, result);

// 			if ((bool)r.ec or r.ptr != e)
// 			{
// 				result = {};
// 				if (cif::VERBOSE)
// 				{
// 					if (r.ec == std::errc::invalid_argument)
// 						std::cerr << "Attempt to convert " << std::quoted(txt) << " into a number\n";
// 					else if (r.ec == std::errc::result_out_of_range)
// 						std::cerr << "Conversion of " << std::quoted(txt) << " into a type that is too small\n";
// 					else
// 						std::cerr << "Not a valid number " << std::quoted(txt) << '\n';
// 				}
// 			}
// 		}

// 		return result;
// 	}

// 	static int compare(const item_handle &ref, const T &value, bool icase)
// 	{
// 		int result = 0;

// 		auto txt = ref.text();

// 		if (ref.empty())
// 			result = 1;
// 		else
// 		{
// 			value_type v = {};

// 			auto b = txt.data();
// 			auto e = txt.data() + txt.size();

// 			std::from_chars_result r = (b + 1 < e and *b == '+' and std::isdigit(b[1]))
// 			                               ? from_chars(b + 1, e, v)
// 			                               : from_chars(b, e, v);

// 			if ((bool)r.ec or r.ptr != e)
// 			{
// 				if (cif::VERBOSE)
// 				{
// 					if (r.ec == std::errc::invalid_argument)
// 						std::cerr << "Attempt to convert " << std::quoted(txt) << " into a number\n";
// 					else if (r.ec == std::errc::result_out_of_range)
// 						std::cerr << "Conversion of " << std::quoted(txt) << " into a type that is too small\n";
// 					else
// 						std::cerr << "Not a valid number " << std::quoted(txt) << '\n';
// 				}
// 				result = 1;
// 			}
// 			else if (std::abs(v - value) <= std::numeric_limits<value_type>::epsilon())
// 				result = 0;
// 			else if (v < value)
// 				result = -1;
// 			else if (v > value)
// 				result = 1;
// 		}

// 		return result;
// 	}
// };

// template <typename T>
// struct item_handle::item_value_as<std::optional<T>>
// {
// 	static std::optional<T> convert(const item_handle &ref)
// 	{
// 		std::optional<T> result;
// 		if (ref)
// 			result = ref.as<T>();
// 		return result;
// 	}

// 	static int compare(const item_handle &ref, std::optional<T> value, bool icase)
// 	{
// 		if (ref.empty() and not value)
// 			return 0;

// 		if (ref.empty())
// 			return -1;
// 		else if (not value)
// 			return 1;
// 		else
// 			return ref.compare(*value, icase);
// 	}
// };

// template <typename T>
// struct item_handle::item_value_as<T, std::enable_if_t<std::is_same_v<T, bool>>>
// {
// 	static bool convert(const item_handle &ref)
// 	{
// 		bool result = false;
// 		if (not ref.empty())
// 			result = iequals(ref.text(), "y");
// 		return result;
// 	}

// 	static int compare(const item_handle &ref, bool value, bool icase)
// 	{
// 		bool rv = convert(ref);
// 		return value && rv ? 0
// 		                   : (rv < value ? -1 : 1);
// 	}
// };

// template <std::size_t N>
// struct item_handle::item_value_as<char[N]>
// {
// 	static std::string convert(const item_handle &ref)
// 	{
// 		if (ref.empty())
// 			return {};
// 		return { ref.text().data(), ref.text().size() };
// 	}

// 	static int compare(const item_handle &ref, const char (&value)[N], bool icase)
// 	{
// 		return icase ? cif::icompare(ref.text(), value) : ref.text().compare(value);
// 	}
// };

// template <typename T>
// struct item_handle::item_value_as<T, std::enable_if_t<std::is_same_v<T, const char *>>>
// {
// 	static std::string convert(const item_handle &ref)
// 	{
// 		if (ref.empty())
// 			return {};
// 		return { ref.text().data(), ref.text().size() };
// 	}

// 	static int compare(const item_handle &ref, const char *value, bool icase)
// 	{
// 		return icase ? cif::icompare(ref.text(), value) : ref.text().compare(value);
// 	}
// };

// template <typename T>
// struct item_handle::item_value_as<T, std::enable_if_t<std::is_same_v<T, std::string_view>>>
// {
// 	static std::string convert(const item_handle &ref)
// 	{
// 		if (ref.empty())
// 			return {};
// 		return { ref.text().data(), ref.text().size() };
// 	}

// 	static int compare(const item_handle &ref, const std::string_view &value, bool icase)
// 	{
// 		return icase ? cif::icompare(ref.text(), value) : ref.text().compare(value);
// 	}
// };

// template <typename T>
// struct item_handle::item_value_as<T, std::enable_if_t<std::is_same_v<T, std::string>>>
// {
// 	static std::string convert(const item_handle &ref)
// 	{
// 		if (ref.empty())
// 			return {};
// 		return { ref.text().data(), ref.text().size() };
// 	}

// 	static int compare(const item_handle &ref, const std::string &value, bool icase)
// 	{
// 		return icase ? cif::icompare(ref.text(), value) : ref.text().compare(value);
// 	}
// };

/** @endcond */

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

/** @endcond */

} // namespace std