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

namespace cif
{

struct category_validator;

// --------------------------------------------------------------------

class validation_error : public std::exception
{
  public:
	validation_error(const std::string &msg);
	validation_error(const std::string &cat, const std::string &item,
		const std::string &msg);
	const char *what() const noexcept { return m_msg.c_str(); }
	std::string m_msg;
};

// --------------------------------------------------------------------

enum class DDL_PrimitiveType
{
	Char,
	UChar,
	Numb
};

DDL_PrimitiveType map_to_primitive_type(std::string_view s);

struct regex_impl;

struct type_validator
{
	std::string m_name;
	DDL_PrimitiveType m_primitive_type;
	regex_impl *m_rx;

	type_validator() = delete;
	type_validator(std::string_view name, DDL_PrimitiveType type, std::string_view rx);

	type_validator(const type_validator &) = delete;
	type_validator(type_validator &&rhs)
		: m_name(std::move(rhs.m_name))
		, m_primitive_type(rhs.m_primitive_type)
	{
		m_rx = std::exchange(rhs.m_rx, nullptr);
	}

	type_validator &operator=(const type_validator &) = delete;
	type_validator &operator=(type_validator &&rhs)
	{
		m_name = std::move(rhs.m_name);
		m_primitive_type = rhs.m_primitive_type;
		m_rx = std::exchange(rhs.m_rx, nullptr);

		return *this;
	}

	~type_validator();

	bool operator<(const type_validator &rhs) const
	{
		return icompare(m_name, rhs.m_name) < 0;
	}

	int compare(std::string_view a, std::string_view b) const;
};

struct item_validator
{
	std::string m_tag;
	bool m_mandatory;
	const type_validator *m_type;
	cif::iset m_enums;
	std::string m_default;
	bool m_default_is_null;
	category_validator *m_category = nullptr;

	// ItemLinked is used for non-key links
	struct item_link
	{
		item_validator *m_parent;
		std::string m_parent_item;
		std::string m_child_item;
	};

	std::vector<item_link> mLinked;

	bool operator<(const item_validator &rhs) const
	{
		return icompare(m_tag, rhs.m_tag) < 0;
	}

	bool operator==(const item_validator &rhs) const
	{
		return iequals(m_tag, rhs.m_tag);
	}

	void operator()(std::string_view value) const;
};

struct category_validator
{
	std::string m_name;
	std::vector<std::string> m_keys;
	cif::iset m_groups;
	cif::iset m_mandatory_fields;
	std::set<item_validator> m_item_validators;

	bool operator<(const category_validator &rhs) const
	{
		return icompare(m_name, rhs.m_name) < 0;
	}

	void addItemValidator(item_validator &&v);

	const item_validator *get_validator_for_item(std::string_view tag) const;

	const std::set<item_validator> &item_validators() const
	{
		return m_item_validators;
	}
};

struct link_validator
{
	int m_link_group_id;
	std::string m_parent_category;
	std::vector<std::string> m_parent_keys;
	std::string m_child_category;
	std::vector<std::string> m_child_keys;
	std::string m_link_group_label;
};

// --------------------------------------------------------------------

class validator
{
  public:
	validator(std::string_view name)
		: m_name(name)
	{
	}

	~validator() = default;

	validator(const validator &rhs) = delete;
	validator &operator=(const validator &rhs) = delete;

	validator(validator &&rhs) = default;
	validator &operator=(validator &&rhs) = default;

	friend class dictionary_parser;

	void add_type_validator(type_validator &&v);
	const type_validator *get_validator_for_type(std::string_view type_code) const;

	void add_category_validator(category_validator &&v);
	const category_validator *get_validator_for_category(std::string_view category) const;

	void add_link_validator(link_validator &&v);
	std::vector<const link_validator *> get_links_for_parent(std::string_view category) const;
	std::vector<const link_validator *> get_links_for_child(std::string_view category) const;

	void report_error(const std::string &msg, bool fatal) const;

	const std::string &name() const { return m_name; }
	void set_name(const std::string &name) { m_name = name; }

	const std::string &version() const { return m_version; }
	void version(const std::string &version) { m_version = version; }

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
class validator_factory
{
  public:
	static validator_factory &instance()
	{
		static validator_factory s_instance;
		return s_instance;
	}

	const validator &operator[](std::string_view dictionary_name);

	const validator &construct_validator(std::string_view name, std::istream &is);

  private:

	// --------------------------------------------------------------------

	validator_factory() = default;

	std::mutex m_mutex;
	std::list<validator> m_validators;
};

} // namespace cif
