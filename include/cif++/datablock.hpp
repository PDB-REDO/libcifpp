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

#include "cif++/category.hpp"
#include "cif++/forward_decl.hpp"

/** \file datablock.hpp
 * Each valid mmCIF file contains at least one @ref cif::datablock.
 * A datablock has a name and can contain one or more @ref cif::category "categories"
 */

namespace cif
{

// --------------------------------------------------------------------

/**
 * @brief A datablock is a list of category objects with some additional features
 * 
 */

class datablock : public std::list<category>
{
  public:
	datablock() = default;

	/**
	 * @brief Construct a new datablock object with name @a name
	 * 
	 * @param name The name for the new datablock
	 */
	datablock(std::string_view name)
		: m_name(name)
	{
	}

	/** @cond */
	datablock(const datablock &);

	datablock(datablock &&db) noexcept
	{
		swap_(*this, db);
	}

	datablock &operator=(datablock db)
	{
		swap_(*this, db);
		return *this;
	}
	/** @endcond */

	friend void swap_(datablock &a, datablock &b) noexcept
	{
		std::swap(a.m_name, b.m_name);
		std::swap(a.m_validator, b.m_validator);
		std::swap(static_cast<std::list<category>&>(a), static_cast<std::list<category>&>(b));
	}

	// --------------------------------------------------------------------

	/**
	 * @brief Return the name of this datablock
	 */
	const std::string &name() const { return m_name; }

	/**
	 * @brief Set the name of this datablock to @a name
	 * 
	 * @param name The new name
	 */
	void set_name(std::string_view name)
	{
		m_name = name;
	}

	/**
	 * @brief Set the validator object to @a v
	 * 
	 * @param v The new validator object, may be null
	 */
	void set_validator(const validator *v);

	/**
	 * @brief Get the validator object
	 * 
	 * @return const validator* The validator or nullptr if there is none
	 */
	const validator *get_validator() const;

	/**
	 * @brief Validates the content of this datablock and all its content
	 * 
	 * @return true If the content is valid
	 * @return false If the content is not valid
	 */
	bool is_valid() const;

	/**
	 * @brief Validates the content of this datablock and all its content
	 * and updates or removes the audit_conform category to match the result.
	 * 
	 * @return true If the content is valid
	 * @return false If the content is not valid
	 */
	bool is_valid();

	/**
	 * @brief Validates all contained data for valid links between parents and children
	 * as defined in the validator
	 * 
	 * @return true If all links are valid
	 * @return false If all links are not valid
	 */
	bool validate_links() const;

	// --------------------------------------------------------------------

	/**
	 * @brief Return the category named @a name, will create a new and empty
	 * category named @a name if it does not exist.
	 * 
	 * @param name The name of the category to return
	 * @return category& Reference to the named category
	 */
	category &operator[](std::string_view name);

	/**
	 * @brief Return the const category named @a name, will return a reference
	 * to a static empty category if it was not found.
	 * 
	 * @param name The name of the category to return
	 * @return category& Reference to the named category
	 */
	const category &operator[](std::string_view name) const;

	/**
	 * @brief Return a pointer to the category named @a name or nullptr if
	 * it does not exist.
	 * 
	 * @param name The name of the category
	 * @return category* Pointer to the category found or nullptr
	 */
	category *get(std::string_view name);

	/**
	 * @brief Return a pointer to the category named @a name or nullptr if
	 * it does not exist.
	 * 
	 * @param name The name of the category
	 * @return category* Pointer to the category found or nullptr
	 */
	const category *get(std::string_view name) const;

	/**
	 * @brief Tries to find a category with name @a name and will create a
	 * new one if it is not found. The result is a tuple of an iterator
	 * pointing to the category and a boolean indicating whether the category
	 * was created or not.
	 * 
	 * @param name The name for the category
	 * @return std::tuple<iterator, bool> A tuple containing an iterator pointing
	 * at the category and a boolean indicating whether the category was newly
	 * created.
	 */
	std::tuple<iterator, bool> emplace(std::string_view name);

	/**
	 * @brief Get the preferred order of the categories when writing them
	 */
	[[deprecated("use get_item_order instead")]]
	std::vector<std::string> get_tag_order() const
	{
		return get_item_order();
	}

	/**
	 * @brief Get the preferred order of the categories when writing them
	 */
	std::vector<std::string> get_item_order() const;

	/**
	 * @brief Write out the contents to @a os
	 */
	void write(std::ostream &os) const;

	/**
	 * @brief Write out the contents to @a os using the order defined in @a item_name_order
	 */
	void write(std::ostream &os, const std::vector<std::string> &item_name_order);

	/**
	 * @brief Friend operator<< to write datablock @a db to std::ostream @a os
	 */
	friend std::ostream &operator<<(std::ostream &os, const datablock &db)
	{
		db.write(os);
		return os;
	}

	// --------------------------------------------------------------------

	/**
	 * @brief Comparison operator to compare two datablock for equal content
	 */
	bool operator==(const datablock &rhs) const;

  private:
	std::string m_name;
	const validator *m_validator = nullptr;
};

} // namespace cif