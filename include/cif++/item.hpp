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
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
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
class item
{
  public:
	/// \brief Default constructor, empty item
	item() = default;

	/// \brief constructor for an item with name \a name and as
	/// content the character '.', i.e. an inapplicable value.
	item(std::string_view name)
		: m_name(name)
		, m_value({ '.' })
	{
	}

	/// \brief constructor for an item with name \a name and as
	/// content a single character string with content \a value
	item(std::string_view name, char value)
		: m_name(name)
		, m_value({ value })
	{
	}

	/// \brief constructor for an item with name \a name and as
	/// content the formatted floating point value \a value with
	/// precision \a precision
	template <typename T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
	item(std::string_view name, const T &value, int precision)
		: m_name(name)
	{
		using namespace std;
		using namespace cif;

		char buffer[32];

		auto r = to_chars(buffer, buffer + sizeof(buffer) - 1, value, chars_format::fixed, precision);
		if ((bool)r.ec)
			throw std::runtime_error("Could not format number");

		m_value.assign(buffer, r.ptr - buffer);
	}

	/// \brief constructor for an item with name \a name and as
	/// content a formatted floating point value \a value with
	/// so-called general formatting
	template <typename T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
	item(const std::string_view name, const T &value)
		: m_name(name)
	{
		using namespace std;
		using namespace cif;

		char buffer[32];

		auto r = to_chars(buffer, buffer + sizeof(buffer) - 1, value, chars_format::general);
		if ((bool)r.ec)
			throw std::runtime_error("Could not format number");

		m_value.assign(buffer, r.ptr - buffer);
	}

	/// \brief constructor for an item with name \a name and as
	/// content the formatted integral value \a value
	template <typename T, std::enable_if_t<std::is_integral_v<T> and not std::is_same_v<T, bool>, int> = 0>
	item(const std::string_view name, const T &value)
		: m_name(name)
	{
		char buffer[32];

		auto r = std::to_chars(buffer, buffer + sizeof(buffer) - 1, value);
		if ((bool)r.ec)
			throw std::runtime_error("Could not format number");

		m_value.assign(buffer, r.ptr - buffer);
	}

	/// \brief constructor for an item with name \a name and as
	/// content the formatted boolean value \a value
	template <typename T, std::enable_if_t<std::is_same_v<T, bool>, int> = 0>
	item(const std::string_view name, const T &value)
		: m_name(name)
	{
		m_value.assign(value ? "y" : "n");
	}

	/// \brief constructor for an item with name \a name and as
	/// content value \a value
	item(const std::string_view name, std::string_view value)
		: m_name(name)
		, m_value(value)
	{
	}

	/// \brief constructor for an item with name \a name and as
	/// content value \a value
	template<typename T, std::enable_if_t<std::is_same_v<T, std::string>, int> = 0>
	item(const std::string_view name, T &&value)
		: m_name(name)
		, m_value(std::move(value))
	{
	}

	/// \brief constructor for an item with name \a name and as
	/// content the optional value \a value
	template <typename T>
	item(const std::string_view name, const std::optional<T> &value)
		: m_name(name)
	{
		if (value.has_value())
		{
			item tmp(name, *value);
			std::swap(tmp.m_value, m_value);
		}
		else
			m_value.assign("?");
	}

	/// \brief constructor for an item with name \a name and as
	/// content the formatted floating point value \a value with
	/// precision \a precision
	template <typename T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
	item(std::string_view name, const std::optional<T> &value, int precision)
		: m_name(name)
	{
		if (value.has_value())
		{
			item tmp(name, *value, precision);
			std::swap(tmp.m_value, m_value);
		}
		else
			m_value.assign("?");
	}

	/** @cond */
	item(const item &rhs) = default;
	item(item &&rhs) noexcept = default;
	item &operator=(const item &rhs) = default;
	item &operator=(item &&rhs) noexcept = default;
	/** @endcond */

	std::string_view name() const { return m_name; }   ///< Return the name of the item
	std::string_view value() const & { return m_value; } ///< Return the value of the item
	std::string value() const && { return std::move(m_value); } ///< Return the value of the item

	/// \brief replace the content of the stored value with \a v
	void value(std::string_view v) { m_value = v; }

	/// \brief empty means either null or unknown
	bool empty() const { return m_value.empty(); }

	/// \brief returns true if the item contains '.'
	bool is_null() const { return m_value == "."; }

	/// \brief returns true if the item contains '?'
	bool is_unknown() const { return m_value == "?"; }

	/// \brief the length of the value string
	std::size_t length() const { return m_value.length(); }

	/// \brief support for structured binding
	template <std::size_t N>
	decltype(auto) get() const
	{
		if constexpr (N == 0)
			return name();
		else if constexpr (N == 1)
			return value();
	}

  private:
	std::string_view m_name;
	std::string m_value;
};

// --------------------------------------------------------------------
/// \brief the internal storage for items in a category
///
/// Internal storage, strictly forward linked list with minimal space
/// requirements. Strings of size 7 or shorter are stored internally.
/// Typically, more than 99% of the strings in an mmCIF file are less
/// than 8 bytes in length.

struct item_value
{
	/** @cond */
	item_value() = default;
	/** @endcond */

	/// \brief constructor
	item_value(std::string_view text)
		: m_length(text.length())
		, m_storage(0)
	{
		if (m_length >= kBufferSize)
		{
			m_data = new char[m_length + 1];
			std::copy(text.begin(), text.end(), m_data);
			m_data[m_length] = 0;
		}
		else
		{
			std::copy(text.begin(), text.end(), m_local_data);
			m_local_data[m_length] = 0;
		}
	}

	/** @cond */
	item_value(item_value &&rhs) noexcept
		: m_length(std::exchange(rhs.m_length, 0))
		, m_storage(std::exchange(rhs.m_storage, 0))
	{
	}

	item_value &operator=(item_value &&rhs) noexcept
	{
		std::swap(m_length, rhs.m_length);
		std::swap(m_storage, rhs.m_storage);
		return *this;
	}

	~item_value()
	{
		if (m_length >= kBufferSize)
			delete[] m_data;
		m_storage = 0;
		m_length = 0;
	}

	item_value(const item_value &) = delete;
	item_value &operator=(const item_value &) = delete;
	/** @endcond */

	/** operator bool, allows easy checking for empty items */
	explicit operator bool() const
	{
		return m_length != 0;
	}

	std::size_t m_length = 0; ///< Length of the data
	union
	{
		char m_local_data[8]; ///< Storage area for small strings (strings smaller than kBufferSize)
		char *m_data;         ///< Pointer to a string stored in the heap
		uint64_t m_storage;   ///< Alternative storage of the data, used in move operations
	};

	/** The maximum length of locally stored strings */
	static constexpr std::size_t kBufferSize = sizeof(m_local_data);

	// By using std::string_view instead of c_str we obain a
	// nice performance gain since we avoid many calls to strlen.

	/** Return the content of the item as a std::string_view */
	constexpr inline std::string_view text() const
	{
		return { m_length >= kBufferSize ? m_data : m_local_data, m_length };
	}
};

// --------------------------------------------------------------------
// Transient object to access stored data

/// \brief This is item_handle, it is used to access the data stored in item_value.

struct item_handle
{
  public:
	/** @cond */
	// conversion helper class
	template <typename T, typename = void>
	struct item_value_as;
	/** @endcond */

	/**
	 * @brief Assign value @a value to the item referenced
	 *
	 * @tparam T Type of the value
	 * @param value The value
	 * @return reference to this item_handle
	 */
	template <typename T>
	item_handle &operator=(const T &value)
	{
		assign_value(item{ "", value }.value());
		return *this;
	}

	/**
	 * @brief Assign value @a value to the item referenced
	 *
	 * @tparam T Type of the value
	 * @param value The value
	 * @return reference to this item_handle
	 */
	template <typename T>
	item_handle &operator=(T &&value)
	{
		assign_value(item{ "", std::forward<T>(value) }.value());
		return *this;
	}

	/**
	 * @brief Assign value @a value to the item referenced
	 *
	 * @tparam T Type of the value
	 * @param value The value
	 * @return reference to this item_handle
	 */
	template <std::size_t N>
	item_handle &operator=(const char (&value)[N])
	{
		assign_value(item{ "", std::move(value) }.value());
		return *this;
	}

	/**
	 * @brief A method with a variable number of arguments that will be concatenated and
	 * assigned as a string. Use it like this:
	 *
	 * @code{.cpp}
	 * cif::item_handle ih;
	 * is.os("The result of ", 1, " * ", 42, " is of course ", 42);
	 * @endcode
	 *
	 * And the content will then be `The result of 1 * 42 is of course 42`.
	 *
	 * @tparam Ts Types of the parameters
	 * @param v The parameters to concatenate
	 */
	template <typename... Ts>
	void os(const Ts &...v)
	{
		std::ostringstream ss;
		((ss << v), ...);
		this->operator=(ss.str());
	}

	/** Swap contents of this and @a b */
	void swap(item_handle &b);

	/** Return the contents of this item as type @tparam T */
	template <typename T = std::string>
	auto as() const -> T
	{
		using value_type = std::remove_cv_t<std::remove_reference_t<T>>;
		return item_value_as<value_type>::convert(*this);
	}

	/** Return the contents of this item as type @tparam T or, if not
	 * set, use @a dv as the default value.
	 */
	template <typename T>
	auto value_or(const T &dv) const
	{
		return empty() ? dv : this->as<T>();
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
	template <typename T>
	int compare(const T &value, bool icase = true) const
	{
		return item_value_as<T>::compare(*this, value, icase);
	}

	/**
	 * @brief Compare the value contained with the value @a value and
	 * return true if both are equal.
	 */
	template <typename T>
	bool operator==(const T &value) const
	{
		// TODO: icase or not icase?
		return item_value_as<T>::compare(*this, value, true) == 0;
	}

	// We may not have C++20 yet...

	/**
	 * @brief Compare the value contained with the value @a value and
	 * return true if both are not equal.
	 */
	template <typename T>
	bool operator!=(const T &value) const
	{
		return not operator==(value);
	}

	/**
	 * @brief Returns true if the content string is empty or
	 * only contains '.' meaning null or '?' meaning unknown
	 * in a mmCIF context
	 */
	bool empty() const
	{
		auto txt = text();
		return txt.empty() or (txt.length() == 1 and (txt.front() == '.' or txt.front() == '?'));
	}

	/** Easy way to test for an empty item */
	explicit operator bool() const { return not empty(); }

	/// is_null return true if the item contains '.'
	bool is_null() const
	{
		auto txt = text();
		return txt.length() == 1 and txt.front() == '.';
	}

	/// is_unknown returns true if the item contains '?'
	bool is_unknown() const
	{
		auto txt = text();
		return txt.length() == 1 and txt.front() == '?';
	}

	/** Return a std::string_view for the contents */
	std::string_view text() const;

	/**
	 * @brief Construct a new item handle object
	 *
	 * @param item Item index
	 * @param row Reference to the row
	 */
	item_handle(uint16_t item, row_handle &row)
		: m_item_ix(item)
		, m_row_handle(row)
	{
	}

	/** A variable holding an empty item */
	CIFPP_EXPORT static const item_handle s_null_item;

	/** friend to swap two item handles */
	friend void swap(item_handle a, item_handle b)
	{
		a.swap(b);
	}

  private:
	item_handle();

	uint16_t m_item_ix;
	row_handle &m_row_handle;

	void assign_value(std::string_view value);
};

// So sad that older gcc implementations of from_chars did not support floats yet...

/** @cond */
template <typename T>
struct item_handle::item_value_as<T, std::enable_if_t<std::is_arithmetic_v<T> and not std::is_same_v<T, bool>>>
{
	using value_type = std::remove_reference_t<std::remove_cv_t<T>>;

	static value_type convert(const item_handle &ref)
	{
		value_type result = {};

		if (not ref.empty())
		{
			auto txt = ref.text();

			auto b = txt.data();
			auto e = txt.data() + txt.size();

			std::from_chars_result r = (b + 1 < e and *b == '+' and std::isdigit(b[1])) ? selected_charconv<value_type>::from_chars(b + 1, e, result) : selected_charconv<value_type>::from_chars(b, e, result);

			if ((bool)r.ec or r.ptr != e)
			{
				result = {};
				if (cif::VERBOSE)
				{
					if (r.ec == std::errc::invalid_argument)
						std::cerr << "Attempt to convert " << std::quoted(txt) << " into a number\n";
					else if (r.ec == std::errc::result_out_of_range)
						std::cerr << "Conversion of " << std::quoted(txt) << " into a type that is too small\n";
					else
						std::cerr << "Not a valid number " << std::quoted(txt) << '\n';
				}
			}
		}

		return result;
	}

	static int compare(const item_handle &ref, const T &value, bool icase)
	{
		int result = 0;

		auto txt = ref.text();

		if (ref.empty())
			result = 1;
		else
		{
			value_type v = {};

			auto b = txt.data();
			auto e = txt.data() + txt.size();

			std::from_chars_result r = (b + 1 < e and *b == '+' and std::isdigit(b[1])) ? selected_charconv<value_type>::from_chars(b + 1, e, v) : selected_charconv<value_type>::from_chars(b, e, v);

			if ((bool)r.ec or r.ptr != e)
			{
				if (cif::VERBOSE)
				{
					if (r.ec == std::errc::invalid_argument)
						std::cerr << "Attempt to convert " << std::quoted(txt) << " into a number\n";
					else if (r.ec == std::errc::result_out_of_range)
						std::cerr << "Conversion of " << std::quoted(txt) << " into a type that is too small\n";
					else
						std::cerr << "Not a valid number " << std::quoted(txt) << '\n';
				}
				result = 1;
			}
			else if (std::abs(v - value) <= std::numeric_limits<value_type>::epsilon())
				result = 0;
			else if (v < value)
				result = -1;
			else if (v > value)
				result = 1;
		}

		return result;
	}
};

template <typename T>
struct item_handle::item_value_as<std::optional<T>>
{
	static std::optional<T> convert(const item_handle &ref)
	{
		std::optional<T> result;
		if (ref)
			result = ref.as<T>();
		return result;
	}

	static int compare(const item_handle &ref, std::optional<T> value, bool icase)
	{
		if (ref.empty() and not value)
			return 0;

		if (ref.empty())
			return -1;
		else if (not value)
			return 1;
		else
			return ref.compare(*value, icase);
	}
};

template <typename T>
struct item_handle::item_value_as<T, std::enable_if_t<std::is_same_v<T, bool>>>
{
	static bool convert(const item_handle &ref)
	{
		bool result = false;
		if (not ref.empty())
			result = iequals(ref.text(), "y");
		return result;
	}

	static int compare(const item_handle &ref, bool value, bool icase)
	{
		bool rv = convert(ref);
		return value && rv ? 0
		                   : (rv < value ? -1 : 1);
	}
};

template <std::size_t N>
struct item_handle::item_value_as<char[N]>
{
	static std::string convert(const item_handle &ref)
	{
		if (ref.empty())
			return {};
		return { ref.text().data(), ref.text().size() };
	}

	static int compare(const item_handle &ref, const char (&value)[N], bool icase)
	{
		return icase ? cif::icompare(ref.text(), value) : ref.text().compare(value);
	}
};

template <typename T>
struct item_handle::item_value_as<T, std::enable_if_t<std::is_same_v<T, const char *>>>
{
	static std::string convert(const item_handle &ref)
	{
		if (ref.empty())
			return {};
		return { ref.text().data(), ref.text().size() };
	}

	static int compare(const item_handle &ref, const char *value, bool icase)
	{
		return icase ? cif::icompare(ref.text(), value) : ref.text().compare(value);
	}
};

template <typename T>
struct item_handle::item_value_as<T, std::enable_if_t<std::is_same_v<T, std::string_view>>>
{
	static std::string convert(const item_handle &ref)
	{
		if (ref.empty())
			return {};
		return { ref.text().data(), ref.text().size() };
	}

	static int compare(const item_handle &ref, const std::string_view &value, bool icase)
	{
		return icase ? cif::icompare(ref.text(), value) : ref.text().compare(value);
	}
};

template <typename T>
struct item_handle::item_value_as<T, std::enable_if_t<std::is_same_v<T, std::string>>>
{
	static std::string convert(const item_handle &ref)
	{
		if (ref.empty())
			return {};
		return { ref.text().data(), ref.text().size() };
	}

	static int compare(const item_handle &ref, const std::string &value, bool icase)
	{
		return icase ? cif::icompare(ref.text(), value) : ref.text().compare(value);
	}
};

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