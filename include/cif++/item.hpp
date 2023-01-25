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

#include <cif++/exports.hpp>
#include <cif++/forward_decl.hpp>
#include <cif++/text.hpp>
#include <cif++/utilities.hpp>

#include <cassert>
#include <charconv>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <utility>

/// \file item.hpp
/// This file contains the declaration of item but also the item_value and item_handle
/// These handle the storage of and access to the data for a single data field. 

namespace cif
{

// --------------------------------------------------------------------
/// \brief item is a transient class that is used to pass data into rows
///        but it also takes care of formatting data. 
class item
{
  public:
	/// \brief Default constructor, empty item
	item() = default;

	/// \brief constructor for an item with name \a name and as
	/// content a single character string with content \a value
	item(std::string_view name, char value)
		: m_name(name)
		, m_value({ value })
	{
	}

	/// \brief constructor for an item with name \a name and as
	/// content a the formatted floating point value \a value with
	/// precision \a precision
	template <typename T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
	item(std::string_view name, const T &value, int precision)
		: m_name(name)
	{
		using namespace std;
		using namespace cif;

		char buffer[32];

		auto r = to_chars(buffer, buffer + sizeof(buffer) - 1, value, chars_format::fixed, precision);
		if (r.ec != std::errc())
			throw std::runtime_error("Could not format number");

		assert(r.ptr >= buffer and r.ptr < buffer + sizeof(buffer));
		*r.ptr = 0;
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
		if (r.ec != std::errc())
			throw std::runtime_error("Could not format number");

		assert(r.ptr >= buffer and r.ptr < buffer + sizeof(buffer));
		*r.ptr = 0;
		m_value.assign(buffer, r.ptr - buffer);
	}

	/// \brief constructor for an item with name \a name and as
	/// content a the formatted integral value \a value
	template <typename T, std::enable_if_t<std::is_integral_v<T> and not std::is_same_v<T,bool>, int> = 0>
	item(const std::string_view name, const T &value)
		: m_name(name)
	{
		char buffer[32];

		auto r = std::to_chars(buffer, buffer + sizeof(buffer) - 1, value);
		if (r.ec != std::errc())
			throw std::runtime_error("Could not format number");

		assert(r.ptr >= buffer and r.ptr < buffer + sizeof(buffer));
		*r.ptr = 0;
		m_value.assign(buffer, r.ptr - buffer);
	}

	/// \brief constructor for an item with name \a name and as
	/// content a the formatted boolean value \a value
	template <typename T, std::enable_if_t<std::is_same_v<T,bool>, int> = 0>
	item(const std::string_view name, const T &value)
		: m_name(name)
	{
		m_value.assign(value ? "y" : "n");
	}

	/// \brief constructor for an item with name \a name and as
	/// content value \a value
	item(const std::string_view name, const std::string_view value)
		: m_name(name)
		, m_value(value)
	{
	}

	item(const item &rhs) = default;

	item(item &&rhs) noexcept = default;

	item &operator=(const item &rhs) = default;

	item &operator=(item &&rhs) noexcept = default;

	std::string_view name() const { return m_name; }
	std::string_view value() const { return m_value; }

	/// \brief replace the content of the stored value with \a v
	void value(std::string_view v) { m_value = v; }

	/// \brief empty means either null or unknown
	bool empty() const { return m_value.empty(); }

	/// \brief returns true if the field contains '.'
	bool is_null() const { return m_value == "."; }

	/// \brief returns true if the field contains '?'
	bool is_unknown() const { return m_value == "?"; }

	/// \brief the length of the value string
	size_t length() const { return m_value.length(); }

	/// \brief support for structured binding
	template<size_t N>
	decltype(auto) get() const
	{
		     if constexpr (N == 0) return name();
		else if constexpr (N == 1) return value();
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
	item_value() = default;

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

	item_value(item_value &&rhs)
		: m_length(std::exchange(rhs.m_length, 0))
		, m_storage(std::exchange(rhs.m_storage, 0))
	{
	}

	item_value &operator=(item_value &&rhs)
	{
		if (this != &rhs)
		{
			m_length = std::exchange(rhs.m_length, m_length);
			m_storage = std::exchange(rhs.m_storage, m_storage);
		}
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

	explicit operator bool() const
	{
		return m_length != 0;
	}

	size_t m_length = 0;
	union
	{
		char m_local_data[8];
		char *m_data;
		uint64_t m_storage;
	};

	static constexpr size_t kBufferSize = sizeof(m_local_data);

	// By using std::string_view instead of c_str we obain a
	// nice performance gain since we avoid many calls to strlen.
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
	// conversion helper class
	template <typename T, typename = void>
	struct item_value_as;

	template <typename T>
	item_handle &operator=(const T &value)
	{
		item v{ "", value };
		assign_value(v);
		return *this;
	}

	template <typename... Ts>
	void os(const Ts &...v)
	{
		std::ostringstream ss;
		((ss << v), ...);
		this->operator=(ss.str());
	}

	void swap(item_handle &b);

	template <typename T = std::string>
	auto as() const -> T
	{
		using value_type = std::remove_cv_t<std::remove_reference_t<T>>;
		return item_value_as<value_type>::convert(*this);
	}

	template <typename T>
	auto value_or(const T &dv) const
	{
		return empty() ? dv : this->as<T>();
	}

	template <typename T>
	int compare(const T &value, bool icase = true) const
	{
		return item_value_as<T>::compare(*this, value, icase);
	}

	template <typename T>
	bool operator==(const T &value) const
	{
		// TODO: icase or not icase?
		return item_value_as<T>::compare(*this, value, true) == 0;
	}

	// We may not have C++20 yet...
	template <typename T>
	bool operator!=(const T &value) const
	{
		return not operator==(value);
	}

	// empty means either null or unknown
	bool empty() const
	{
		auto txt = text();
		return txt.empty() or (txt.length() == 1 and (txt.front() == '.' or txt.front() == '?'));
	}

	explicit operator bool() const { return not empty(); }

	// is_null means the field contains '.'
	bool is_null() const
	{
		auto txt = text();
		return txt.length() == 1 and txt.front() == '.';
	}

	// is_unknown means the field contains '?'
	bool is_unknown() const
	{
		auto txt = text();
		return txt.length() == 1 and txt.front() == '?';
	}

	std::string_view text() const;

	item_handle(uint16_t column, row_handle &row)
		: m_column(column)
		, m_row_handle(row)
	{
	}

	static CIFPP_EXPORT const item_handle s_null_item;

	friend void swap(item_handle a, item_handle b)
	{
		a.swap(b);
	}

  private:
	item_handle();

	uint16_t m_column;
	row_handle &m_row_handle;

	void assign_value(const item &value);
};

// So sad that older gcc implementations of from_chars did not support floats yet...

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

			std::from_chars_result r = selected_charconv<value_type>::from_chars(txt.data(), txt.data() + txt.size(), result);

			if (r.ec != std::errc())
			{
				result = {};
				if (cif::VERBOSE)
				{
					if (r.ec == std::errc::invalid_argument)
						std::cerr << "Attempt to convert " << std::quoted(txt) << " into a number" << std::endl;
					else if (r.ec == std::errc::result_out_of_range)
						std::cerr << "Conversion of " << std::quoted(txt) << " into a type that is too small" << std::endl;
				}
			}
		}

		return result;
	}

	static int compare(const item_handle &ref, const T &value, bool icase)
	{
		int result = 0;

		auto txt = ref.text();

		if (txt.empty())
			result = 1;
		else
		{
			value_type v = {};

			std::from_chars_result r = selected_charconv<value_type>::from_chars(txt.data(), txt.data() + txt.size(), v);

			if (r.ec != std::errc())
			{
				if (cif::VERBOSE)
				{
					if (r.ec == std::errc::invalid_argument)
						std::cerr << "Attempt to convert " << std::quoted(txt) << " into a number" << std::endl;
					else if (r.ec == std::errc::result_out_of_range)
						std::cerr << "Conversion of " << std::quoted(txt) << " into a type that is too small" << std::endl;
				}
				result = 1;
			}
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

template <size_t N>
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

} // namespace cif

namespace std
{

template<> struct tuple_size<::cif::item>
            : public std::integral_constant<std::size_t, 2> {};

template<> struct tuple_element<0, ::cif::item>
{
	using type = decltype(std::declval<::cif::item>().name());
};

template<> struct tuple_element<1, ::cif::item>
{
	using type = decltype(std::declval<::cif::item>().value());
};

}