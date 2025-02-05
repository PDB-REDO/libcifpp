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

#include "cif++/item.hpp"

#include <array>

/**
 * @file row.hpp
 * 
 * The class cif::row should be an opaque type. It is used to store the
 * internal data per row in a category. You should use cif::row_handle
 * to get access to the contents in a row.
 * 
 * One could think of rows as vectors of cif::item. But internally
 * that's not the case.
 * 
 * You can access the values of stored items by name or index.
 * The return value of operator[] is an cif::item_handle object.
 * 
 * @code {.cpp}
 * cif::category &atom_site = my_db["atom_site"];
 * cif::row_handle rh = atom_site.front();
 * 
 * // by name:
 * std::string name = rh["label_atom_id"].as<std::string>();
 * 
 * // by index:
 * uint16_t ix = atom_site.get_item_ix("label_atom_id");
 * assert(rh[ix].as<std::string() == name);
 * @endcode
 * 
 * There some template magic here to allow easy extracting of data
 * from rows. This can be done using cif::tie e.g.:
 * 
 * @code {.cpp}
 * std::string name;
 * float x, y, z;
 * 
 * cif::tie(name, x, y, z) = rh.get("label_atom_id", "cartn_x", "cartn_y", "cartn_z");
 * @endcode
 * 
 * However, a more modern way uses structured binding:
 * 
 * @code {.cpp}
 * const auto &[name, x, y, z] = rh.get<std::string,float,float,float>("label_atom_id", "cartn_x", "cartn_y", "cartn_z");
 * @endcode
 * 
 * 
 * 
 */

namespace cif
{

namespace detail
{

	// some helper classes to help create tuple result types
	template <typename... C>
	struct get_row_result
	{
		static constexpr std::size_t N = sizeof...(C);

		get_row_result(const row_handle &r, std::array<uint16_t, N> &&items)
			: m_row(r)
			, m_items(std::move(items))
		{
		}

		const item_handle operator[](uint16_t ix) const
		{
			return m_row[m_items[ix]];
		}

		template <typename... Ts, std::enable_if_t<N == sizeof...(Ts), int> = 0>
		operator std::tuple<Ts...>() const
		{
			return get<Ts...>(std::index_sequence_for<Ts...>{});
		}

		template <typename... Ts, std::size_t... Is>
		std::tuple<Ts...> get(std::index_sequence<Is...>) const
		{
			return std::tuple<Ts...>{ m_row[m_items[Is]].template as<Ts>()... };
		}

		const row_handle &m_row;
		std::array<uint16_t, N> m_items;
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

/// \brief similar to std::tie, assign values to each element in @a v from the 
/// result of a get on a row_handle.
template <typename... Ts>
auto tie(Ts &...v)
{
	return detail::tie_wrap<Ts &...>(std::forward<Ts &>(v)...);
}

// --------------------------------------------------------------------
/// \brief the row class, this one is not directly accessible from the outside

class row : public std::vector<item_value>
{
  public:
	row() = default;

	/**
	 * @brief Return the item_value pointer for item at index @a ix
	 */
	item_value* get(uint16_t ix)
	{
		return ix < size() ? &data()[ix] : nullptr;
	}

	/**
	 * @brief Return the const item_value pointer for item at index @a ix
	 */
	const item_value* get(uint16_t ix) const
	{
		return ix < size() ? &data()[ix] : nullptr;
	}

  private:
	friend class category;
	friend class category_index;

	template <typename, typename...>
	friend class iterator_impl;

	void append(uint16_t ix, item_value &&iv)
	{
		if (ix >= size())
			resize(ix + 1);
		
		at(ix) = std::move(iv);
	}

	void remove(uint16_t ix)
	{
		if (ix < size())
			at(ix) = item_value{};
	}

	row *m_next = nullptr;
};

// --------------------------------------------------------------------
/// \brief row_handle is the way to access data stored in rows

class row_handle
{
  public:
	/** @cond */
	friend struct item_handle;
	friend class category;
	friend class category_index;
	friend class row_initializer;
	template <typename, typename...> friend class iterator_impl;

	row_handle() = default;

	row_handle(const row_handle &) = default;
	row_handle(row_handle &&) = default;

	row_handle &operator=(const row_handle &) = default;
	row_handle &operator=(row_handle &&) = default;

	/** @endcond */

	/// \brief constructor taking a category @a cat and a row @a r
	row_handle(const category &cat, const row &r)
		: m_category(const_cast<category *>(&cat))
		, m_row(const_cast<row *>(&r))
	{
	}

	/// \brief return the category this row belongs to
	const category &get_category() const
	{
		return *m_category;
	}

	/// \brief Return true if the row is empty or uninitialised
	bool empty() const
	{
		return m_category == nullptr or m_row == nullptr;
	}

	/// \brief convenience method to test for empty()
	explicit operator bool() const
	{
		return not empty();
	}

	/// \brief return a cif::item_handle to the item in item @a item_ix
	item_handle operator[](uint16_t item_ix)
	{
		return empty() ? item_handle::s_null_item : item_handle(item_ix, *this);
	}

	/// \brief return a const cif::item_handle to the item in item @a item_ix
	const item_handle operator[](uint16_t item_ix) const
	{
		return empty() ? item_handle::s_null_item : item_handle(item_ix, const_cast<row_handle &>(*this));
	}

	/// \brief return a cif::item_handle to the item in the item named @a item_name
	item_handle operator[](std::string_view item_name)
	{
		return empty() ? item_handle::s_null_item : item_handle(add_item(item_name), *this);
	}

	/// \brief return a const cif::item_handle to the item in the item named @a item_name
	const item_handle operator[](std::string_view item_name) const
	{
		return empty() ? item_handle::s_null_item : item_handle(get_item_ix(item_name), const_cast<row_handle &>(*this));
	}

	/// \brief Return an object that can be used in combination with cif::tie
	/// to assign the values for the items @a items
	template <typename... C>
	auto get(C... items) const
	{
		return detail::get_row_result<C...>(*this, { get_item_ix(items)... });
	}

	/// \brief Return a tuple of values of types @a Ts for the items @a items
	template <typename... Ts, typename... C, std::enable_if_t<sizeof...(Ts) == sizeof...(C) and sizeof...(C) != 1, int> = 0>
	std::tuple<Ts...> get(C... items) const
	{
		return detail::get_row_result<Ts...>(*this, { get_item_ix(items)... });
	}

	/// \brief Get the value of item @a item cast to type @a T
	template <typename T>
	T get(const char *item) const
	{
		return operator[](get_item_ix(item)).template as<T>();
	}

	/// \brief Get the value of item @a item cast to type @a T
	template <typename T>
	T get(std::string_view item) const
	{
		return operator[](get_item_ix(item)).template as<T>();
	}

	/// \brief assign each of the items named in @a values to their respective value
	void assign(const std::vector<item> &values)
	{
		for (auto &value : values)
			assign(value, true);
	}

	/** \brief assign the value @a value to the item named @a name 
	 * 
	 * If updateLinked it true, linked records are updated as well.
	 * That means that if item @a name is part of the link definition
	 * and the link results in a linked record in another category
	 * this record in the linked category is updated as well.
	 * 
	 * If validate is true, which is default, the assigned value is
	 * checked to see if it conforms to the rules defined in the dictionary
	 */

	void assign(std::string_view name, std::string_view value, bool updateLinked, bool validate = true)
	{
		assign(add_item(name), value, updateLinked, validate);
	}

	/** \brief assign the value @a value to item at index @a item
	 * 
	 * If updateLinked it true, linked records are updated as well.
	 * That means that if item @a item is part of the link definition
	 * and the link results in a linked record in another category
	 * this record in the linked category is updated as well.
	 * 
	 * If validate is true, which is default, the assigned value is
	 * checked to see if it conforms to the rules defined in the dictionary
	 */

	void assign(uint16_t item, std::string_view value, bool updateLinked, bool validate = true);

	/// \brief compare two rows
	bool operator==(const row_handle &rhs) const { return m_category == rhs.m_category and m_row == rhs.m_row; }

	/// \brief compare two rows
	bool operator!=(const row_handle &rhs) const { return m_category != rhs.m_category or m_row != rhs.m_row; }

  private:
	uint16_t get_item_ix(std::string_view name) const;
	std::string_view get_item_name(uint16_t ix) const;

	uint16_t add_item(std::string_view name);

	row *get_row()
	{
		return m_row;
	}

	const row *get_row() const
	{
		return m_row;
	}

	void assign(const item &i, bool updateLinked)
	{
		assign(i.name(), i.value(), updateLinked);
	}

	void swap(uint16_t item, row_handle &r);

	category *m_category = nullptr;
	row *m_row = nullptr;
};

// --------------------------------------------------------------------

/**
 * @brief The class row_initializer is a list of cif::item's.
 * 
 * This class is used to construct new rows, it allows to
 * group a list of item name and value pairs and pass it
 * in one go to the constructing function.
 */
class row_initializer : public std::vector<item>
{
  public:
	/** @cond */
	friend class category;

	row_initializer() = default;
	row_initializer(const row_initializer &) = default;
	row_initializer(row_initializer &&) = default;
	row_initializer &operator=(const row_initializer &) = default;
	row_initializer &operator=(row_initializer &&) = default;

	/** @endcond */

	/// \brief constructor taking a std::initializer_list of items
	row_initializer(std::initializer_list<item> items)
		: std::vector<item>(items)
	{
	}

	/// \brief constructor taking a range of items
	template <typename ItemIter, std::enable_if_t<std::is_same_v<typename ItemIter::value_type, item>, int> = 0>
	row_initializer(ItemIter b, ItemIter e)
		: std::vector<item>(b, e)
	{
	}

	/// \brief constructor taking the values of an existing row
	row_initializer(row_handle rh);


	/// \brief set the value for item name @a name to @a value
	void set_value(std::string_view name, std::string_view value);

	/// \brief set the value for item based on @a i
	void set_value(const item &i)
	{
		set_value(i.name(), i.value());
	}

	/// \brief set the value for item name @a name to @a value, but only if the item did not have a value already
	void set_value_if_empty(std::string_view name, std::string_view value);

	/// \brief set the value for item @a i, but only if the item did not have a value already
	void set_value_if_empty(const item &i)
	{
		set_value_if_empty(i.name(), i.value());
	}
};

} // namespace cif