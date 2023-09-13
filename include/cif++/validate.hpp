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

#include <filesystem>
#include <list>
#include <mutex>
#include <utility>

/**
 * @file validate.hpp
 *
 * Support for validating mmCIF files based on a dictionary. These dictionaries
 * contain information about the categories and items therein, what they may
 * contain and how this should be formatted. There's also information on links
 * between parent and child categories.
 *
 */

namespace cif
{

struct category_validator;

// --------------------------------------------------------------------

/**
 * @brief The exception thrown when a validation error occurs
 *
 */
class validation_error : public std::exception
{
  public:
	/// @brief Constructor
	validation_error(const std::string &msg);

	/// @brief Constructor
	validation_error(const std::string &cat, const std::string &item,
		const std::string &msg);

	/// @brief The description of the error
	const char *what() const noexcept { return m_msg.c_str(); }

	/// @cond
	std::string m_msg;
	/// @endcond
};

// --------------------------------------------------------------------

/** @brief the primitive types known */
enum class DDL_PrimitiveType
{
	Char,  ///< Text
	UChar, ///< Text that is compared ignoring the character case
	Numb   ///< Nummeric values
};

/// @brief Return the DDL_PrimitiveType encoded in @a s
DDL_PrimitiveType map_to_primitive_type(std::string_view s);

struct regex_impl;

/**
 * @brief For each defined type in a dictionary a type_validator is created
 *
 * A type validator can check if the contents of an item are conforming the
 * specification. The check is done using regular expressions.
 *
 * A type_validator can also be used to compare two values that conform to
 * this type. Comparison is of course based on the primitive type.
 *
 */
struct type_validator
{
	std::string m_name;                 ///< The name of the type
	DDL_PrimitiveType m_primitive_type; ///< The primitive_type of the type
	regex_impl *m_rx;                   ///< The regular expression for the type

	type_validator() = delete;

	/// @brief Constructor
	type_validator(std::string_view name, DDL_PrimitiveType type, std::string_view rx);

	type_validator(const type_validator &) = delete;

	/// @brief Copy constructor
	type_validator(type_validator &&rhs)
		: m_name(std::move(rhs.m_name))
		, m_primitive_type(rhs.m_primitive_type)
	{
		m_rx = std::exchange(rhs.m_rx, nullptr);
	}

	type_validator &operator=(const type_validator &) = delete;

	/// @brief Move constructor
	type_validator &operator=(type_validator &&rhs)
	{
		m_name = std::move(rhs.m_name);
		m_primitive_type = rhs.m_primitive_type;
		m_rx = std::exchange(rhs.m_rx, nullptr);

		return *this;
	}

	/// @brief Destructor
	~type_validator();

	/// @brief Return the sorting order
	bool operator<(const type_validator &rhs) const
	{
		return icompare(m_name, rhs.m_name) < 0;
	}

	/// @brief Compare the contents of @a a and @a b based on the
	/// primitive type of this type. A value of zero indicates the
	/// values are equal. Less than zero means @a a sorts before @a b
	/// and a value larger than zero likewise means the opposite
	int compare(std::string_view a, std::string_view b) const;
};

/**
 * @brief An item_validator binds a type_validator to an item in
 * a category along with other information found in the dictionary.
 *
 * mmCIF dictionaries may indicate an item is e.g. mandatory or
 * consists of a certain list of allowed values. Even default
 * values can be provided.
 *
 */
struct item_validator
{
	std::string m_tag;                        ///< The item name
	bool m_mandatory;                         ///< Flag indicating this item is mandatory
	const type_validator *m_type;             ///< The type for this item
	cif::iset m_enums;                        ///< If filled, the set of allowed values
	std::string m_default;                    ///< If filled, a default value for this item
	category_validator *m_category = nullptr; ///< The category_validator this item_validator belongs to

	/// @brief Compare based on the name
	bool operator<(const item_validator &rhs) const
	{
		return icompare(m_tag, rhs.m_tag) < 0;
	}

	/// @brief Compare based on the name
	bool operator==(const item_validator &rhs) const
	{
		return iequals(m_tag, rhs.m_tag);
	}

	/// @brief Validate the value in @a value for this item
	/// Will throw a validation_error exception if it fails
	void operator()(std::string_view value) const;
};

/**
 * @brief A validator for categories
 *
 * Categories can have a key, a set of items that in combination
 * should be unique.
 */
struct category_validator
{
	std::string m_name;                         ///< The name of the category
	std::vector<std::string> m_keys;            ///< The list of items that make up the key
	cif::iset m_groups;							///< The category groups this category belongs to
	cif::iset m_mandatory_fields;               ///< The mandatory fields for this category
	std::set<item_validator> m_item_validators; ///< The item validators for the items in this category

	/// @brief return true if this category sorts before @a rhs
	bool operator<(const category_validator &rhs) const
	{
		return icompare(m_name, rhs.m_name) < 0;
	}

	/// @brief Add item_validator @a v to the list of item validators
	void addItemValidator(item_validator &&v);

	/// @brief Return the item_validator for item @a tag, may return nullptr
	const item_validator *get_validator_for_item(std::string_view tag) const;
};

/**
 * @brief A validator for links between categories
 *
 * Links are defined as a set of pairs of item names in a
 * parent category and a corresponding item in a child
 * category. This means that the size of m_parent_keys
 * is always equal to the size of m_child_keys.
 *
 * Multiple links may be defined between two categories.
 *
 */
struct link_validator
{
	int m_link_group_id;                    ///< The link group ID
	std::string m_parent_category;          ///< The name of the parent category
	std::vector<std::string> m_parent_keys; ///< The items in the parent category making up the set of linked items
	std::string m_child_category;           ///< The name of the child category
	std::vector<std::string> m_child_keys;  ///< The items in the child category making up the set of linked items
	std::string m_link_group_label;         ///< The group label assigned to this link
};

// --------------------------------------------------------------------

/**
 * @brief The validator class combines all the link, category and item validator classes
 *
 */
class validator
{
  public:
	/**
	 * @brief Construct a new validator object
	 *
	 * @param name The name of the underlying dictionary
	 */
	validator(std::string_view name)
		: m_name(name)
	{
	}

	/// @brief destructor
	~validator() = default;

	validator(const validator &rhs) = delete;
	validator &operator=(const validator &rhs) = delete;

	/// @brief move constructor
	validator(validator &&rhs) = default;

	/// @brief move assignment operator
	validator &operator=(validator &&rhs) = default;

	friend class dictionary_parser;

	/// @brief Add type_validator @a v to the list of type validators
	void add_type_validator(type_validator &&v);

	/// @brief Return the type validator for @a type_code, may return nullptr
	const type_validator *get_validator_for_type(std::string_view type_code) const;

	/// @brief Add category_validator @a v to the list of category validators
	void add_category_validator(category_validator &&v);

	/// @brief Return the category validator for @a category, may return nullptr
	const category_validator *get_validator_for_category(std::string_view category) const;

	/// @brief Add link_validator @a v to the list of link validators
	void add_link_validator(link_validator &&v);

	/// @brief Return the list of link validators for which the parent is @a category
	std::vector<const link_validator *> get_links_for_parent(std::string_view category) const;

	/// @brief Return the list of link validators for which the child is @a category
	std::vector<const link_validator *> get_links_for_child(std::string_view category) const;

	/// @brief Bottleneck function to report an error in validation
	void report_error(const std::string &msg, bool fatal) const;

	const std::string &name() const { return m_name; }        ///< Get the name of this validator
	void set_name(const std::string &name) { m_name = name; } ///< Set the name of this validator

	const std::string &version() const { return m_version; }              ///< Get the version of this validator
	void set_version(const std::string &version) { m_version = version; } ///< Set the version of this validator

  private:
	// name is fully qualified here:
	item_validator *get_validator_for_item(std::string_view name) const;

	std::string m_name;
	std::string m_version;
	bool m_strict = false;
	std::set<type_validator> m_type_validators;
	std::set<category_validator> m_category_validators;
	std::vector<link_validator> m_link_validators;
};

// --------------------------------------------------------------------

/**
 * @brief Validators are globally unique objects, use the validator_factory
 * class to construct them. This class is a singleton.
 */

class validator_factory
{
  public:
	/// @brief Return the singleton instance
	static validator_factory &instance()
	{
		static validator_factory s_instance;
		return s_instance;
	}

	/// @brief Return the validator with name @a dictionary_name
	const validator &operator[](std::string_view dictionary_name);

	/// @brief Construct a new validator with name @a name from the data in @a is
	const validator &construct_validator(std::string_view name, std::istream &is);

  private:
	// --------------------------------------------------------------------

	validator_factory() = default;

	std::mutex m_mutex;
	std::list<validator> m_validators;
};

} // namespace cif
