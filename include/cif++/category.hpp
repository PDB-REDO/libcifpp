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
#include "cif++/text.hpp"
#include "cif++/validate.hpp"

#include <array>

/** \file category.hpp
 * Documentation for the cif::category class
 *
 * The category class should meet the requirements of Container and
 * SequenceContainer.
 *
 * TODO: implement all of:
 * https://en.cppreference.com/w/cpp/named_req/Container
 * https://en.cppreference.com/w/cpp/named_req/SequenceContainer
 * and more?
 */

namespace cif
{

// --------------------------------------------------------------------
// special exceptions

/// @brief A duplicate_key_error is thrown when an attempt is made
/// to insert a row with values that would introduce a duplicate key
/// in the index. Of course, this can only happen if a @ref category_validator
/// has been defined for this category.
class duplicate_key_error : public std::runtime_error
{
  public:
	/**
	 * @brief Construct a new duplicate key error object
	 */
	duplicate_key_error(const std::string &msg)
		: std::runtime_error(msg)
	{
	}
};

/// @brief A missing_key_error is thrown when an attempt is made
/// to create an index when one of the key items is missing.
class missing_key_error : public std::runtime_error
{
  public:
	/**
	 * @brief Construct a new duplicate key error object
	 */
	missing_key_error(const std::string &msg, const std::string &key)
		: std::runtime_error(msg)
		, m_key(key)
	{
	}

	const std::string &get_key() const noexcept { return m_key; }

  private:
	std::string m_key;
};

/// @brief A multiple_results_error is throw when you request a single
/// row using a query but the query contains more than exactly one row.
class multiple_results_error : public std::runtime_error
{
  public:
	/**
	 * @brief Construct a new multiple results error object
	 */
	multiple_results_error()
		: std::runtime_error("query should have returned exactly one row")
	{
	}
};

// --------------------------------------------------------------------
// These should be moved elsewhere, one day.

/// \cond
template <typename _Tp>
inline constexpr bool is_optional_v = false;
template <typename _Tp>
inline constexpr bool is_optional_v<std::optional<_Tp>> = true;
/// \endcond

// --------------------------------------------------------------------

/// The class category is a sequence container for rows of data values.
/// You could think of it as a std::vector<cif::row_handle> like class.
///
/// A @ref category_validator can be assigned to an object of category
/// after which this class can validate contained data and use an
/// index to keep key values unique.

class category
{
  public:
	/// \cond

	friend class row_handle;

	template <typename, typename...>
	friend class iterator_impl;

	using value_type = row_handle;
	using reference = value_type;
	using const_reference = const value_type;
	using iterator = iterator_impl<category>;
	using const_iterator = iterator_impl<const category>;

	/// \endcond

	category() = default;            ///< Default constructor
	category(std::string_view name); ///< Constructor taking a \a name
	category(const category &rhs);   ///< Copy constructor

	category(category &&rhs) noexcept ///< Move constructor
	{
		swap(*this, rhs);
	}

	category &operator=(category rhs) ///< assignement operator
	{
		swap(*this, rhs);
		return *this;
	}

	/// @brief Destructor
	/// @note Please note that the destructor is not virtual. It is assumed that
	/// you will not derive from this class.
	~category();

	friend void swap(category &a, category &b) noexcept;

	// --------------------------------------------------------------------

	const std::string &name() const { return m_name; } ///< Returns the name of the category

	[[deprecated("use key_items instead")]] iset key_fields() const; ///< Returns the cif::iset of key item names. Retrieved from the @ref category_validator for this category

	iset key_items() const; ///< Returns the cif::iset of key item names. Retrieved from the @ref category_validator for this category

	[[deprecated("use key_item_indices instead")]] std::set<uint16_t> key_field_indices() const; ///< Returns a set of indices for the key items.

	std::set<uint16_t> key_item_indices() const; ///< Returns a set of indices for the key items.

	/// @brief Set the validator for this category to @a v
	/// @param v The category_validator to assign. A nullptr value is allowed.
	/// @param db The enclosing @ref datablock
	void set_validator(const validator *v, datablock &db);

	/// @brief Update the links in this category
	/// @param db The enclosing @ref datablock
	void update_links(const datablock &db);

	/// @brief Return the global @ref validator for the data
	/// @return The @ref validator or nullptr if not assigned
	const validator *get_validator() const { return m_validator; }

	/// @brief Return the category validator for this category
	/// @return The @ref category_validator or nullptr if not assigned
	const category_validator *get_cat_validator() const { return m_cat_validator; }

	/// @brief Validate the data stored using the assigned @ref category_validator
	/// @return Returns true is all validations pass
	bool is_valid() const;

	/// @brief Validate links, that means, values in this category should have an
	/// accompanying value in parent categories.
	///
	/// @note
	/// The code makes one exception when validating missing links and that's between
	/// *atom_site* and a parent *pdbx_poly_seq_scheme* or *entity_poly_seq*.
	/// This particular case should be skipped because it is wrong:
	/// there are atoms that are not part of a polymer, and thus will have no
	/// parent in those categories.
	///
	/// @return Returns true is all validations pass
	bool validate_links() const;

	/// @brief Equality operator, returns true if @a rhs is equal to this
	/// @param rhs The object to compare with
	/// @return True if the data contained is equal
	bool operator==(const category &rhs) const;

	/// @brief Unequality operator, returns true if @a rhs is not equal to this
	/// @param rhs The object to compare with
	/// @return True if the data contained is not equal
	bool operator!=(const category &rhs) const
	{
		return not operator==(rhs);
	}

	// --------------------------------------------------------------------

	/// @brief Return a reference to the first row in this category.
	/// @return Reference to the first row in this category. The result is undefined if
	/// the category is empty.
	reference front()
	{
		return { *this, *m_head };
	}

	/// @brief Return a const reference to the first row in this category.
	/// @return const reference to the first row in this category. The result is undefined if
	/// the category is empty.
	const_reference front() const
	{
		return { const_cast<category &>(*this), const_cast<row &>(*m_head) };
	}

	/// @brief Return a reference to the last row in this category.
	/// @return Reference to the last row in this category. The result is undefined if
	/// the category is empty.
	reference back()
	{
		return { *this, *m_tail };
	}

	/// @brief Return a const reference to the last row in this category.
	/// @return const reference to the last row in this category. The result is undefined if
	/// the category is empty.
	const_reference back() const
	{
		return { const_cast<category &>(*this), const_cast<row &>(*m_tail) };
	}

	/// Return an iterator to the first row
	iterator begin()
	{
		return { *this, m_head };
	}

	/// Return an iterator pointing past the last row
	iterator end()
	{
		return { *this, nullptr };
	}

	/// Return a const iterator to the first row
	const_iterator begin() const
	{
		return { *this, m_head };
	}

	/// Return a const iterator pointing past the last row
	const_iterator end() const
	{
		return { *this, nullptr };
	}

	/// Return a const iterator to the first row
	const_iterator cbegin() const
	{
		return { *this, m_head };
	}

	/// Return an iterator pointing past the last row
	const_iterator cend() const
	{
		return { *this, nullptr };
	}

	/// Return a count of the rows in this container
	std::size_t size() const
	{
		return std::distance(cbegin(), cend());
	}

	/// Return the theoretical maximum number or rows that can be stored
	std::size_t max_size() const
	{
		return std::numeric_limits<std::size_t>::max(); // this is a bit optimistic, I guess
	}

	/// Return true if the category is empty
	bool empty() const
	{
		return m_head == nullptr;
	}

	// --------------------------------------------------------------------
	// A category can have a key, as defined by the validator/dictionary

	/// @brief The key type
	using key_type = row_initializer;

	/// @brief Return a row_handle for the row specified by \a key
	/// @param key The value for the key, items specified in the dictionary should have a value
	/// @return The row found in the index, or an undefined row_handle
	row_handle operator[](const key_type &key);

	/// @brief Return a const row_handle for the row specified by \a key
	/// @param key The value for the key, items specified in the dictionary should have a value
	/// @return The row found in the index, or an undefined row_handle
	const row_handle operator[](const key_type &key) const
	{
		return const_cast<category *>(this)->operator[](key);
	}

	// --------------------------------------------------------------------

	/// @brief Return a special const iterator for all rows in this category.
	/// This iterator can be used in a structured binding context. E.g.:
	///
	/// @code{.cpp}
	/// for (const auto &[name, value] : cat.rows<std::string,int>("item_name", "item_value"))
	///   std::cout << name << ": " << value << '\n';
	/// @endcode
	///
	/// @tparam Ts The types for the items requested
	/// @param names The names for the items requested

	template <typename... Ts, typename... Ns>
	iterator_proxy<const category, Ts...> rows(Ns... names) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of item names should be equal to the number of types to return");
		return iterator_proxy<const category, Ts...>(*this, begin(), { names... });
	}

	/// @brief Return a special iterator for all rows in this category.
	/// This iterator can be used in a structured binding context. E.g.:
	///
	/// @code{.cpp}
	/// for (const auto &[name, value] : cat.rows<std::string,int>("item_name", "item_value"))
	///   std::cout << name << ": " << value << '\n';
	///
	/// // or in case we only need one item:
	///
	/// for (int id : cat.rows<int>("id"))
	///   std::cout << id << '\n';
	/// @endcode
	///
	/// @tparam Ts The types for the items requested
	/// @param names The names for the items requested

	template <typename... Ts, typename... Ns>
	iterator_proxy<category, Ts...> rows(Ns... names)
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of item names should be equal to the number of types to return");
		return iterator_proxy<category, Ts...>(*this, begin(), { names... });
	}

	// --------------------------------------------------------------------

	/// @brief Return a special iterator to loop over all rows that conform to @a cond
	///
	/// @code{.cpp}
	/// for (row_handle rh : cat.find(cif::key("first_name") == "John" and cif::key("last_name") == "Doe"))
	///    .. // do something with rh
	/// @endcode
	///
	/// @param cond The condition for the query
	/// @return A special iterator that loops over all elements that match. The iterator can be dereferenced
	/// to a @ref row_handle

	conditional_iterator_proxy<category> find(condition &&cond)
	{
		return find(begin(), std::move(cond));
	}

	/// @brief Return a special iterator to loop over all rows that conform to @a cond
	/// starting at @a pos
	///
	/// @param pos Where to start searching
	/// @param cond The condition for the query
	/// @return A special iterator that loops over all elements that match. The iterator can be dereferenced
	/// to a @ref row_handle

	conditional_iterator_proxy<category> find(iterator pos, condition &&cond)
	{
		return { *this, pos, std::move(cond) };
	}

	/// @brief Return a special const iterator to loop over all rows that conform to @a cond
	///
	/// @param cond The condition for the query
	/// @return A special iterator that loops over all elements that match. The iterator can be dereferenced
	/// to a const @ref row_handle

	conditional_iterator_proxy<const category> find(condition &&cond) const
	{
		return find(cbegin(), std::move(cond));
	}

	/// @brief Return a special const iterator to loop over all rows that conform to @a cond
	/// starting at @a pos
	///
	/// @param pos Where to start searching
	/// @param cond The condition for the query
	/// @return A special iterator that loops over all elements that match. The iterator can be dereferenced
	/// to a const @ref row_handle

	conditional_iterator_proxy<const category> find(const_iterator pos, condition &&cond) const
	{
		return conditional_iterator_proxy<const category>{ *this, pos, std::move(cond) };
	}

	/// @brief Return a special iterator to loop over all rows that conform to @a cond. The resulting
	/// iterator can be used in a structured binding context.
	///
	/// @code{.cpp}
	/// for (const auto &[name, value] : cat.find<std::string,int>(cif::key("item_value") > 10, "item_name", "item_value"))
	///    std::cout << name << ": " << value << '\n';
	/// @endcode
	///
	/// @param cond The condition for the query
	/// @tparam Ts The types for the items requested
	/// @param names The names for the items requested
	/// @return A special iterator that loops over all elements that match.

	template <typename... Ts, typename... Ns>
	conditional_iterator_proxy<category, Ts...> find(condition &&cond, Ns... names)
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of item names should be equal to the number of types to return");
		return find<Ts...>(cbegin(), std::move(cond), std::forward<Ns>(names)...);
	}

	/// @brief Return a special const iterator to loop over all rows that conform to @a cond. The resulting
	/// iterator can be used in a structured binding context.
	///
	/// @param cond The condition for the query
	/// @tparam Ts The types for the items requested
	/// @param names The names for the items requested
	/// @return A special iterator that loops over all elements that match.

	template <typename... Ts, typename... Ns>
	conditional_iterator_proxy<const category, Ts...> find(condition &&cond, Ns... names) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of item names should be equal to the number of types to return");
		return find<Ts...>(cbegin(), std::move(cond), std::forward<Ns>(names)...);
	}

	/// @brief Return a special iterator to loop over all rows that conform to @a cond starting at @a pos.
	/// The resulting iterator can be used in a structured binding context.
	///
	/// @param pos Iterator pointing to the location where to start
	/// @param cond The condition for the query
	/// @tparam Ts The types for the items requested
	/// @param names The names for the items requested
	/// @return A special iterator that loops over all elements that match.

	template <typename... Ts, typename... Ns>
	conditional_iterator_proxy<category, Ts...> find(const_iterator pos, condition &&cond, Ns... names)
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of item names should be equal to the number of types to return");
		return { *this, pos, std::move(cond), std::forward<Ns>(names)... };
	}

	/// @brief Return a special const iterator to loop over all rows that conform to @a cond starting at @a pos.
	/// The resulting iterator can be used in a structured binding context.
	///
	/// @param pos Iterator pointing to the location where to start
	/// @param cond The condition for the query
	/// @tparam Ts The types for the items requested
	/// @param names The names for the items requested
	/// @return A special iterator that loops over all elements that match.

	template <typename... Ts, typename... Ns>
	conditional_iterator_proxy<const category, Ts...> find(const_iterator pos, condition &&cond, Ns... names) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Ns), "The number of item names should be equal to the number of types to return");
		return { *this, pos, std::move(cond), std::forward<Ns>(names)... };
	}

	// --------------------------------------------------------------------
	// if you only expect a single row

	/// @brief Return the row handle for the row that matches @a cond Throws @a multiple_results_error if
	/// there are is not exactly one row matching @a cond
	/// @param cond The condition to search for
	/// @return Row handle to the row found
	row_handle find1(condition &&cond)
	{
		return find1(begin(), std::move(cond));
	}

	/// @brief Return the row handle for the row that matches @a cond starting at @a pos
	/// Throws @a multiple_results_error if there are is not exactly one row matching @a cond
	/// @param pos The position to start the search
	/// @param cond The condition to search for
	/// @return Row handle to the row found
	row_handle find1(iterator pos, condition &&cond)
	{
		auto h = find(pos, std::move(cond));

		if (h.size() != 1)
			throw multiple_results_error();

		return *h.begin();
	}

	/// @brief Return the const row handle for the row that matches @a cond Throws @a multiple_results_error if
	/// there are is not exactly one row matching @a cond
	/// @param cond The condition to search for
	/// @return Row handle to the row found
	const row_handle find1(condition &&cond) const
	{
		return find1(cbegin(), std::move(cond));
	}

	/// @brief Return const the row handle for the row that matches @a cond starting at @a pos
	/// Throws @a multiple_results_error if there are is not exactly one row matching @a cond
	/// @param pos The position to start the search
	/// @param cond The condition to search for
	/// @return Row handle to the row found
	const row_handle find1(const_iterator pos, condition &&cond) const
	{
		auto h = find(pos, std::move(cond));

		if (h.size() != 1)
			throw multiple_results_error();

		return *h.begin();
	}

	/// @brief Return value for the item named @a item for the single row that
	/// matches @a cond. Throws @a multiple_results_error if there are is not exactly one row
	/// @tparam The type to use for the result
	/// @param cond The condition to search for
	/// @param item The name of the item to return the value for
	/// @return The value found
	template <typename T>
	T find1(condition &&cond, std::string_view item) const
	{
		return find1<T>(cbegin(), std::move(cond), item);
	}

	/// @brief Return value for the item named @a item for the single row that
	/// matches @a cond when starting to search at @a pos.
	/// Throws @a multiple_results_error if there are is not exactly one row
	/// @tparam The type to use for the result
	/// @param pos The location to start the search
	/// @param cond The condition to search for
	/// @param item The name of the item to return the value for
	/// @return The value found
	template <typename T, std::enable_if_t<not is_optional_v<T>, int> = 0>
	T find1(const_iterator pos, condition &&cond, std::string_view item) const
	{
		auto h = find<T>(pos, std::move(cond), item);

		if (h.size() != 1)
			throw multiple_results_error();

		return *h.begin();
	}

	/// @brief Return a value of type std::optional<T> for the item named @a item for the single row that
	/// matches @a cond when starting to search at @a pos.
	/// If the row was not found, an empty value is returned.
	/// @tparam The type to use for the result
	/// @param pos The location to start the search
	/// @param cond The condition to search for
	/// @param item The name of the item to return the value for
	/// @return The value found, can be empty if no row matches the condition
	template <typename T, std::enable_if_t<is_optional_v<T>, int> = 0>
	T find1(const_iterator pos, condition &&cond, std::string_view item) const
	{
		auto h = find<typename T::value_type>(pos, std::move(cond), item);

		if (h.size() > 1)
			throw multiple_results_error();

		if (h.empty())
			return {};

		return *h.begin();
	}

	/// @brief Return a std::tuple for the values for the items named in @a items
	/// for the single row that matches @a cond
	/// Throws @a multiple_results_error if there are is not exactly one row
	/// @tparam The types to use for the resulting tuple
	/// @param cond The condition to search for
	/// @param items The names of the items to return the value for
	/// @return The values found as a single tuple of type std::tuple<Ts...>
	template <typename... Ts, typename... Cs, typename U = std::enable_if_t<sizeof...(Ts) != 1>>
	std::tuple<Ts...> find1(condition &&cond, Cs... items) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Cs), "The number of item names should be equal to the number of types to return");
		// static_assert(std::is_same_v<Cs, const char*>..., "The item names should be const char");
		return find1<Ts...>(cbegin(), std::move(cond), std::forward<Cs>(items)...);
	}

	/// @brief Return a std::tuple for the values for the items named in @a items
	/// for the single row that matches @a cond when starting to search at @a pos
	/// Throws @a multiple_results_error if there are is not exactly one row
	/// @tparam The types to use for the resulting tuple
	/// @param pos The location to start the search
	/// @param cond The condition to search for
	/// @param items The names of the items to return the value for
	/// @return The values found as a single tuple of type std::tuple<Ts...>
	template <typename... Ts, typename... Cs, typename U = std::enable_if_t<sizeof...(Ts) != 1>>
	std::tuple<Ts...> find1(const_iterator pos, condition &&cond, Cs... items) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Cs), "The number of item names should be equal to the number of types to return");
		auto h = find<Ts...>(pos, std::move(cond), std::forward<Cs>(items)...);

		if (h.size() != 1)
			throw multiple_results_error();

		return *h.begin();
	}

	// --------------------------------------------------------------------
	// if you want only a first hit

	/// @brief Return a row handle to the first row that matches @a cond
	/// @param cond The condition to search for
	/// @return The handle to the row that matches or an empty row_handle
	row_handle find_first(condition &&cond)
	{
		return find_first(begin(), std::move(cond));
	}

	/// @brief Return a row handle to the first row that matches @a cond starting at @a pos
	/// @param pos The location to start searching
	/// @param cond The condition to search for
	/// @return The handle to the row that matches or an empty row_handle
	row_handle find_first(iterator pos, condition &&cond)
	{
		auto h = find(pos, std::move(cond));

		return h.empty() ? row_handle{} : *h.begin();
	}

	/// @brief Return a const row handle to the first row that matches @a cond
	/// @param cond The condition to search for
	/// @return The const handle to the row that matches or an empty row_handle
	const row_handle find_first(condition &&cond) const
	{
		return find_first(cbegin(), std::move(cond));
	}

	/// @brief Return a const row handle to the first row that matches @a cond starting at @a pos
	/// @param pos The location to start searching
	/// @param cond The condition to search for
	/// @return The const handle to the row that matches or an empty row_handle
	const row_handle find_first(const_iterator pos, condition &&cond) const
	{
		auto h = find(pos, std::move(cond));

		return h.empty() ? row_handle{} : *h.begin();
	}

	/// @brief Return the value for item @a item for the first row that matches condition @a cond
	/// @tparam The type of the value to return
	/// @param cond The condition to search for
	/// @param item The item for which the value should be returned
	/// @return The value found or a default constructed value if not found
	template <typename T>
	T find_first(condition &&cond, std::string_view item) const
	{
		return find_first<T>(cbegin(), std::move(cond), item);
	}

	/// @brief Return the value for item @a item for the first row that matches condition @a cond
	/// when starting the search at @a pos
	/// @tparam The type of the value to return
	/// @param pos The location to start searching
	/// @param cond The condition to search for
	/// @param item The item for which the value should be returned
	/// @return The value found or a default constructed value if not found
	template <typename T>
	T find_first(const_iterator pos, condition &&cond, std::string_view item) const
	{
		auto h = find<T>(pos, std::move(cond), item);

		return h.empty() ? T{} : *h.begin();
	}

	/// @brief Return a tuple containing the values for the items @a items for the first row that matches condition @a cond
	/// @tparam The types of the values to return
	/// @param cond The condition to search for
	/// @param items The items for which the values should be returned
	/// @return The values found or default constructed values if not found
	template <typename... Ts, typename... Cs, typename U = std::enable_if_t<sizeof...(Ts) != 1>>
	std::tuple<Ts...> find_first(condition &&cond, Cs... items) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Cs), "The number of item names should be equal to the number of types to return");
		// static_assert(std::is_same_v<Cs, const char*>..., "The item names should be const char");
		return find_first<Ts...>(cbegin(), std::move(cond), std::forward<Cs>(items)...);
	}

	/// @brief Return a tuple containing the values for the items @a items for the first row that matches condition @a cond
	/// when starting the search at @a pos
	/// @tparam The types of the values to return
	/// @param pos The location to start searching
	/// @param cond The condition to search for
	/// @param items The items for which the values should be returned
	/// @return The values found or default constructed values if not found
	template <typename... Ts, typename... Cs, typename U = std::enable_if_t<sizeof...(Ts) != 1>>
	std::tuple<Ts...> find_first(const_iterator pos, condition &&cond, Cs... items) const
	{
		static_assert(sizeof...(Ts) == sizeof...(Cs), "The number of item names should be equal to the number of types to return");
		auto h = find<Ts...>(pos, std::move(cond), std::forward<Cs>(items)...);

		return h.empty() ? std::tuple<Ts...>{} : *h.begin();
	}

	// --------------------------------------------------------------------

	/// @brief Return the maximum value for item @a item for all rows that match condition @a cond
	/// @tparam The type of the value to return
	/// @param item The item to use for the value
	/// @param cond The condition to search for
	/// @return The value found or the minimal value for the type
	template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	T find_max(std::string_view item, condition &&cond) const
	{
		T result = std::numeric_limits<T>::min();

		for (auto v : find<T>(std::move(cond), item))
		{
			if (result < v)
				result = v;
		}

		return result;
	}

	/// @brief Return the maximum value for item @a item for all rows
	/// @tparam The type of the value to return
	/// @param item The item to use for the value
	/// @return The value found or the minimal value for the type
	template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	T find_max(std::string_view item) const
	{
		return find_max<T>(item, all());
	}

	/// @brief Return the minimum value for item @a item for all rows that match condition @a cond
	/// @tparam The type of the value to return
	/// @param item The item to use for the value
	/// @param cond The condition to search for
	/// @return The value found or the maximum value for the type
	template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	T find_min(std::string_view item, condition &&cond) const
	{
		T result = std::numeric_limits<T>::max();

		for (auto v : find<T>(std::move(cond), item))
		{
			if (result > v)
				result = v;
		}

		return result;
	}

	/// @brief Return the maximum value for item @a item for all rows
	/// @tparam The type of the value to return
	/// @param item The item to use for the value
	/// @return The value found or the maximum value for the type
	template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	T find_min(std::string_view item) const
	{
		return find_min<T>(item, all());
	}

	/// @brief Return whether a row exists that matches condition @a cond
	/// @param cond The condition to match
	/// @return True if a row exists
	[[deprecated("Use contains instead")]] bool exists(condition &&cond) const
	{
		return contains(std::move(cond));
	}

	/// @brief Return whether a row exists that matches condition @a cond
	/// @param cond The condition to match
	/// @return True if a row exists
	bool contains(condition &&cond) const
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

	/// @brief Return the total number of rows that match condition @a cond
	/// @param cond The condition to match
	/// @return The count
	std::size_t count(condition &&cond) const
	{
		std::size_t result = 0;

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

	/// Using the relations defined in the validator, return whether the row
	/// in @a r has any children in other categories
	bool has_children(row_handle r) const;

	/// Using the relations defined in the validator, return whether the row
	/// in @a r has any parents in other categories
	bool has_parents(row_handle r) const;

	/// Using the relations defined in the validator, return the row handles
	/// for all rows in @a childCat that are linked to row @a r
	std::vector<row_handle> get_children(row_handle r, const category &childCat) const;

	/// Using the relations defined in the validator, return the row handles
	/// for all rows in @a parentCat that are linked to row @a r
	std::vector<row_handle> get_parents(row_handle r, const category &parentCat) const;

	/// Using the relations defined in the validator, return the row handles
	/// for all rows in @a cat that are in any way linked to row @a r
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

	/// Erase the row pointed to by @a pos and return the iterator to the
	/// row following pos.
	iterator erase(iterator pos);

	/// Erase row @a rh
	void erase(row_handle rh)
	{
		erase(iterator(*this, rh.m_row));
	}

	/// @brief Erase all rows that match condition @a cond
	/// @param cond The condition
	/// @return The number of rows that have been erased
	std::size_t erase(condition &&cond);

	/// @brief Erase all rows that match condition @a cond calling
	/// the visitor function @a visit for each before actually erasing it.
	/// @param cond The condition
	/// @param visit The visitor function
	/// @return The number of rows that have been erased
	std::size_t erase(condition &&cond, std::function<void(row_handle)> &&visit);

	/// @brief Emplace the values in @a ri in a new row
	/// @param ri An object containing the values to insert
	/// @return iterator to the newly created row
	iterator emplace(row_initializer &&ri)
	{
		return this->emplace(ri.begin(), ri.end());
	}

	/// @brief Create a new row and emplace the values in the range @a b to @a e in it
	/// @param b Iterator to the beginning of the range of @ref item_value
	/// @param e Iterator to the end of the range of @ref item_value
	/// @return iterator to the newly created row
	template <typename ItemIter>
	iterator emplace(ItemIter b, ItemIter e)
	{
		row *r = this->create_row();

		try
		{
			for (auto i = b; i != e; ++i)
			{
				// item_value *new_item = this->create_item(*i);
				r->append(add_item(i->name()), { i->value() });
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

	/// @brief Completely erase all rows contained in this category
	void clear();

	// --------------------------------------------------------------------
	/// \brief generate a new, unique ID. Pass it an ID generating function
	/// based on a sequence number. This function will be called until the
	/// result is unique in the context of this category
	std::string get_unique_id(std::function<std::string(int)> generator = cif::cif_id_for_number);

	/// @brief Generate a new, unique ID based on a string prefix followed by a number
	/// @param prefix The string prefix
	/// @return a new unique ID
	std::string get_unique_id(const std::string &prefix)
	{
		return get_unique_id([prefix](int nr)
			{ return prefix + std::to_string(nr + 1); });
	}

	/// @brief Generate a new, unique value for a item named @a item_name
	/// @param item_name The name of the item
	/// @return a new unique value
	std::string get_unique_value(std::string_view item_name);

	// --------------------------------------------------------------------

	using value_provider_type = std::function<std::string_view(std::string_view)>;

	/// \brief Update a single item named @a item_name in the rows that match
	/// \a cond to values provided by a callback function \a value_provider
	/// making sure the linked categories are updated according to the link.
	/// That means, child categories are updated if the links are absolute
	/// and unique. If they are not, the child category rows are split.

	void update_value(condition &&cond, std::string_view item_name,
		value_provider_type &&value_provider)
	{
		auto rs = find(std::move(cond));
		std::vector<row_handle> rows;
		std::copy(rs.begin(), rs.end(), std::back_inserter(rows));
		update_value(rows, item_name, std::move(value_provider));
	}

	/// \brief Update a single item named @a item_name in the rows \a rows
	/// to values provided by a callback function \a value_provider
	/// making sure the linked categories are updated according to the link.
	/// That means, child categories are updated if the links are absolute
	/// and unique. If they are not, the child category rows are split.

	void update_value(const std::vector<row_handle> &rows, std::string_view item_name,
		value_provider_type &&value_provider);

	/// \brief Update a single item named @a item_name in the rows that match \a cond to value \a value
	/// making sure the linked categories are updated according to the link.
	/// That means, child categories are updated if the links are absolute
	/// and unique. If they are not, the child category rows are split.

	void update_value(condition &&cond, std::string_view item_name, std::string_view value)
	{
		auto rs = find(std::move(cond));
		std::vector<row_handle> rows;
		std::copy(rs.begin(), rs.end(), std::back_inserter(rows));
		update_value(rows, item_name, value);
	}

	/// \brief Update a single item named @a item_name in @a rows to value \a value
	/// making sure the linked categories are updated according to the link.
	/// That means, child categories are updated if the links are absolute
	/// and unique. If they are not, the child category rows are split.

	void update_value(const std::vector<row_handle> &rows, std::string_view item_name, std::string_view value)
	{
		update_value(rows, item_name, [value](std::string_view)
			{ return value; });
	}

	// --------------------------------------------------------------------
	// Naming used to be very inconsistent. For backward compatibility,
	// the old function names are here as deprecated variants.

	/// \brief Return the index number for \a column_name
	[[deprecated("Use get_item_ix instead")]] uint16_t get_column_ix(std::string_view column_name) const
	{
		return get_item_ix(column_name);
	}

	/// @brief Return the name for column with index @a ix
	/// @param ix The index number
	/// @return The name of the column
	[[deprecated("use get_item_name instead")]] std::string_view get_column_name(uint16_t ix) const
	{
		return get_item_name(ix);
	}

	/// @brief Make sure a item with name @a item_name is known and return its index number
	/// @param item_name The name of the item
	/// @return The index number of the item
	[[deprecated("use add_item instead")]] uint16_t add_column(std::string_view item_name)
	{
		return add_item(item_name);
	}

	/** @brief Remove column name @a colum_name
	 * @param column_name The column to be removed
	 */
	[[deprecated("use remove_item instead")]] void remove_column(std::string_view column_name)
	{
		remove_item(column_name);
	}

	/** @brief Rename column @a from_name to @a to_name */
	[[deprecated("use rename_item instead")]] void rename_column(std::string_view from_name, std::string_view to_name)
	{
		rename_item(from_name, to_name);
	}

	/// @brief Return whether a column with name @a name exists in this category
	/// @param name The name of the column
	/// @return True if the column exists
	[[deprecated("use has_item instead")]] bool has_column(std::string_view name) const
	{
		return has_item(name);
	}

	/// @brief Return the cif::iset of columns in this category
	[[deprecated("use get_items instead")]] iset get_columns() const
	{
		return get_items();
	}

	// --------------------------------------------------------------------
	/// \brief Return the index number for \a item_name

	uint16_t get_item_ix(std::string_view item_name) const
	{
		uint16_t result;

		for (result = 0; result < m_items.size(); ++result)
		{
			if (iequals(item_name, m_items[result].m_name))
				break;
		}

		if (VERBOSE > 0 and result == m_items.size() and m_cat_validator != nullptr) // validate the name, if it is known at all (since it was not found)
		{
			auto iv = m_cat_validator->get_validator_for_item(item_name);
			if (iv == nullptr)
				std::cerr << "Invalid name used '" << item_name << "' is not a known item in " + m_name << '\n';
		}

		return result;
	}

	/// @brief Return the name for item with index @a ix
	/// @param ix The index number
	/// @return The name of the item
	std::string_view get_item_name(uint16_t ix) const
	{
		if (ix >= m_items.size())
			throw std::out_of_range("item index is out of range");

		return m_items[ix].m_name;
	}

	/// @brief Make sure a item with name @a item_name is known and return its index number
	/// @param item_name The name of the item
	/// @return The index number of the item
	uint16_t add_item(std::string_view item_name)
	{
		using namespace std::literals;

		uint16_t result = get_item_ix(item_name);

		if (result == m_items.size())
		{
			const item_validator *item_validator = nullptr;

			if (m_cat_validator != nullptr)
			{
				item_validator = m_cat_validator->get_validator_for_item(item_name);
				if (item_validator == nullptr)
					m_validator->report_error(validation_error::item_not_allowed_in_category, m_name, item_name, false);
			}

			m_items.emplace_back(item_name, item_validator);
		}

		return result;
	}

	/** @brief Remove item name @a colum_name
	 * @param item_name The item to be removed
	 */
	void remove_item(std::string_view item_name);

	/** @brief Rename item @a from_name to @a to_name */
	void rename_item(std::string_view from_name, std::string_view to_name);

	/// @brief Return whether a item with name @a name exists in this category
	/// @param name The name of the item
	/// @return True if the item exists
	bool has_item(std::string_view name) const
	{
		return get_item_ix(name) < m_items.size();
	}

	/// @brief Return the cif::iset of items in this category
	iset get_items() const;

	// --------------------------------------------------------------------

	/// @brief Sort the rows using comparator function @a f
	/// @param f The comparator function taking two row_handles and returning
	/// an int indicating whether the first is smaller, equal or larger than
	/// the second. ( respectively a value <0, 0, or >0 )
	void sort(std::function<int(row_handle, row_handle)> f);

	/// @brief Reorder the rows in the category using the index defined by
	/// the @ref category_validator
	void reorder_by_index();

	// --------------------------------------------------------------------

	/// This function returns effectively the list of fully qualified item
	/// names, that is category_name + '.' + item_name for each item
	[[deprecated("use get_item_order instead")]] std::vector<std::string> get_tag_order() const
	{
		return get_item_order();
	}

	/// This function returns effectively the list of fully qualified item
	/// names, that is category_name + '.' + item_name for each item
	std::vector<std::string> get_item_order() const;

	/// Write the contents of the category to the std::ostream @a os
	void write(std::ostream &os) const;

	/// @brief Write the contents of the category to the std::ostream @a os and
	/// use @a order as the order of the items. If @a addMissingItems is
	/// false, items that do not contain any value will be suppressed
	/// @param os The std::ostream to write to
	/// @param order The order in which the items should appear
	/// @param addMissingItems When false, empty items are suppressed from the output
	void write(std::ostream &os, const std::vector<std::string> &order, bool addMissingItems = true);

  private:
	void write(std::ostream &os, const std::vector<uint16_t> &order, bool includeEmptyItems) const;

  public:
	/// friend function to make it possible to do:
	/// @code {.cpp}
	/// std::cout << my_category;
	/// @endcode
	friend std::ostream &operator<<(std::ostream &os, const category &cat)
	{
		cat.write(os);
		return os;
	}

  private:
	void update_value(row *row, uint16_t item, std::string_view value, bool updateLinked, bool validate = true);

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

	struct item_entry
	{
		std::string m_name;
		const item_validator *m_validator;

		item_entry(std::string_view name, const item_validator *validator)
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

// TODO: NEED TO FIX THIS!
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

	void swap_item(uint16_t item_ix, row_handle &a, row_handle &b);

	// --------------------------------------------------------------------

	std::string m_name;
	std::vector<item_entry> m_items;
	const validator *m_validator = nullptr;
	const category_validator *m_cat_validator = nullptr;
	std::vector<link> m_parent_links, m_child_links;
	bool m_cascade = true;
	uint32_t m_last_unique_num = 0;
	class category_index *m_index = nullptr;
	row *m_head = nullptr, *m_tail = nullptr;
};

} // namespace cif
