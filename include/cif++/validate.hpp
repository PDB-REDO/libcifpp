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

#include <cassert>
#include <filesystem>
#include <list>
#include <mutex>
#include <system_error>
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
// New: error_code

/**
 * @enum validation_error
 *
 * @brief A stronly typed class containing the error codes reported by @ref cif::validator and friends
 */
enum class validation_error
{
	value_does_not_match_rx = 1,      /**< The value of an item does not conform to the regular expression specified for it */
	value_is_not_in_enumeration_list, /**< The value of an item is not in the list of values allowed */
	not_a_known_primitive_type,       /**< The type is not a known primitive type */
	undefined_category,               /**< Category has no definition in the dictionary */
	unknown_item,                     /**< The item is not defined to be part of the category */
	incorrect_item_validator,         /**< Incorrectly specified validator for item */
	missing_mandatory_items,          /**< Missing mandatory items */
	missing_key_items,                /**< An index could not be constructed due to missing key items */
	item_not_allowed_in_category,     /**< Requested item allowed in category according to dictionary */
	empty_file,                       /**< The file contains no datablocks */
	empty_datablock,                  /**< The datablock contains no categories */
	empty_category,                   /**< The category is empty */
	not_valid_pdbx,                   /**< The file is not a valid PDBx file */
};
/**
 * @brief The implementation for @ref validation_category error messages
 *
 */
class validation_category_impl : public std::error_category
{
  public:
	/**
	 * @brief User friendly name
	 *
	 * @return const char*
	 */

	const char *name() const noexcept override
	{
		return "cif::validation";
	}

	/**
	 * @brief Provide the error message as a string for the error code @a ev
	 *
	 * @param ev The error code
	 * @return std::string
	 */

	std::string message(int ev) const override
	{
		switch (static_cast<validation_error>(ev))
		{
			case validation_error::value_does_not_match_rx:
				return "Value in item does not match regular expression";
			case validation_error::value_is_not_in_enumeration_list:
				return "Value is not in the enumerated list of valid values";
			case validation_error::not_a_known_primitive_type:
				return "The type is not a known primitive type";
			case validation_error::undefined_category:
				return "Category has no definition in the dictionary";
			case validation_error::unknown_item:
				return "Item is not defined to be part of the category";
			case validation_error::incorrect_item_validator:
				return "Incorrectly specified validator for item";
			case validation_error::missing_mandatory_items:
				return "Missing mandatory items";
			case validation_error::missing_key_items:
				return "An index could not be constructed due to missing key items";
			case validation_error::item_not_allowed_in_category:
				return "Requested item not allowed in category according to dictionary";
			case validation_error::empty_file:
				return "The file contains no datablocks";
			case validation_error::empty_datablock:
				return "The datablock contains no categories";
			case validation_error::empty_category:
				return "The category is empty";
			case validation_error::not_valid_pdbx:
				return "The file is not a valid PDBx file";

			default:
				assert(false);
				return "unknown error code";
		}
	}

	/**
	 * @brief Return whether two error codes are equivalent, always false in this case
	 *
	 */

	bool equivalent(const std::error_code & /*code*/, int /*condition*/) const noexcept override
	{
		return false;
	}
};

/**
 * @brief Return the implementation for the validation_category
 *
 * @return std::error_category&
 */
inline std::error_category &validation_category()
{
	static validation_category_impl instance;
	return instance;
}

inline std::error_code make_error_code(validation_error e)
{
	return std::error_code(static_cast<int>(e), validation_category());
}

inline std::error_condition make_error_condition(validation_error e)
{
	return std::error_condition(static_cast<int>(e), validation_category());
}

// --------------------------------------------------------------------

class validation_exception : public std::runtime_error
{
  public:
	validation_exception(validation_error err)
		: validation_exception(make_error_code(err))
	{
	}

	validation_exception(validation_error err, std::string_view category)
		: validation_exception(make_error_code(err), category)
	{
	}

	validation_exception(validation_error err, std::string_view category, std::string_view item)
		: validation_exception(make_error_code(err), category, item)
	{
	}

	validation_exception(std::error_code ec);

	validation_exception(std::error_code ec, std::string_view category);

	validation_exception(std::error_code ec, std::string_view category, std::string_view item);
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

/// @brief Return the DDL_PrimitiveType encoded in @a s, error reporting variant
DDL_PrimitiveType map_to_primitive_type(std::string_view s, std::error_code &ec) noexcept;

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

/** @brief Item alias, items can be renamed over time
 */

struct item_alias
{
	item_alias(const std::string &alias_name, const std::string &dictionary, const std::string &version)
		: m_name(alias_name)
		, m_dict(dictionary)
		, m_vers(version)
	{
	}

	item_alias(const item_alias &) = default;
	item_alias &operator=(const item_alias &) = default;

	std::string m_name; ///< The alias_name
	std::string m_dict; ///< The dictionary in which it was known
	std::string m_vers; ///< The version of the dictionary
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
	std::string m_item_name;                  ///< The item name
	bool m_mandatory;                         ///< Flag indicating this item is mandatory
	const type_validator *m_type;             ///< The type for this item
	cif::iset m_enums;                        ///< If filled, the set of allowed values
	std::string m_default;                    ///< If filled, a default value for this item
	category_validator *m_category = nullptr; ///< The category_validator this item_validator belongs to
	std::vector<item_alias> m_aliases;        ///< The aliases for this item

	/// @brief Compare based on the name
	bool operator<(const item_validator &rhs) const
	{
		return icompare(m_item_name, rhs.m_item_name) < 0;
	}

	/// @brief Compare based on the name
	bool operator==(const item_validator &rhs) const
	{
		return iequals(m_item_name, rhs.m_item_name);
	}

	/// @brief Validate the value in @a value for this item
	/// Will throw a std::system_error exception if it fails
	void operator()(std::string_view value) const;

	/// @brief A more gentle version of value validation
	bool validate_value(std::string_view value, std::error_code &ec) const noexcept;
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
	cif::iset m_groups;                         ///< The category groups this category belongs to
	cif::iset m_mandatory_items;                ///< The mandatory items for this category
	std::set<item_validator> m_item_validators; ///< The item validators for the items in this category

	/// @brief return true if this category sorts before @a rhs
	bool operator<(const category_validator &rhs) const
	{
		return icompare(m_name, rhs.m_name) < 0;
	}

	/// @brief Add item_validator @a v to the list of item validators
	void add_item_validator(item_validator &&v);

	/// @brief Return the item_validator for item @a item_name, may return nullptr
	const item_validator *get_validator_for_item(std::string_view item_name) const;

	/// @brief Return the item_validator for an item that has as alias name @a item_name, may return nullptr
	const item_validator *get_validator_for_aliased_item(std::string_view item_name) const;
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
	void report_error(validation_error err, bool fatal = true) const
	{
		report_error(make_error_code(err), fatal);
	}

	/// @brief Bottleneck function to report an error in validation
	void report_error(std::error_code ec, bool fatal = true) const;

	/// @brief Bottleneck function to report an error in validation
	void report_error(validation_error err, std::string_view category,
		std::string_view item, bool fatal = true) const
	{
		report_error(make_error_code(err), category, item, fatal);
	}

	/// @brief Bottleneck function to report an error in validation
	void report_error(std::error_code ec, std::string_view category,
		std::string_view item, bool fatal = true) const;

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
	static validator_factory &instance();

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
