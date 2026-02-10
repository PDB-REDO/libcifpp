/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
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

#include "cif++/validate.hpp"

#include "cif++/category.hpp"
#include "cif++/dictionary_parser.hpp"
#include "cif++/utilities.hpp"

#include <algorithm>
#include <cassert>
#include <format>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <stdexcept>

// The validator depends on regular expressions. Unfortunately,
// the implementation of std::regex in g++ is buggy and crashes
// on reading the pdbx dictionary. We used to use boost regex
// instead but using pcre2 is even easier and faster.

#define PCRE2_CODE_UNIT_WIDTH 8
#include <pcre2.h>

namespace cif
{

validation_exception::validation_exception(std::error_code ec)
	: runtime_error(ec.message())
{
}

validation_exception::validation_exception(std::error_code ec, std::string_view category)
	: runtime_error((ec.message() + "; category: ").append(category))
{
}

validation_exception::validation_exception(std::error_code ec, std::string_view category, std::string_view item)
	: runtime_error((ec.message() + "; category: ").append(category).append("; item: ").append(item))
{
}

// --------------------------------------------------------------------

struct regex_impl
{
	regex_impl(std::string_view rx);
	~regex_impl();

	regex_impl(const regex_impl &) = delete;
	regex_impl &operator=(const regex_impl &) = delete;

	bool match(std::string_view v) const;

  private:
	pcre2_code_8 *m_rx = nullptr;
	pcre2_match_data *m_data = nullptr;
	mutable std::mutex m_mutex;
};

regex_impl::regex_impl(std::string_view rx)
{
	int err_code;
	size_t err_offset;
	m_rx = pcre2_compile(reinterpret_cast<PCRE2_SPTR>(rx.data()), rx.length(), 0, &err_code, &err_offset, nullptr);
	if (m_rx == nullptr)
	{
		char buffer[256];
		int n = pcre2_get_error_message(err_code, reinterpret_cast<unsigned char *>(&buffer[0]), sizeof(buffer));

		throw std::runtime_error(std::string("PCRE2 compilation failed: ") + std::string{ buffer, buffer + n });
	}

	m_data = pcre2_match_data_create_from_pattern(m_rx, nullptr);
}

regex_impl::~regex_impl()
{
	std::unique_lock lock(m_mutex);

	if (m_data)
		pcre2_match_data_free(m_data);

	if (m_rx)
		pcre2_code_free(m_rx);
}

bool regex_impl::match(std::string_view v) const
{
	std::unique_lock lock(m_mutex);

	bool result = false;

	if (int rc = pcre2_match(m_rx, reinterpret_cast<PCRE2_SPTR>(v.data()), v.length(), 0, 0, m_data, nullptr); rc >= 0)
		result = true;
	else if (rc != PCRE2_ERROR_NOMATCH)
		std::cerr << "Error matching with pcre\n";

	return result;
}

// --------------------------------------------------------------------

DDL_PrimitiveType map_to_primitive_type(std::string_view s, std::error_code &ec) noexcept
{
	ec = {};
	DDL_PrimitiveType result = DDL_PrimitiveType::Char;
	if (iequals(s, "char"))
		result = DDL_PrimitiveType::Char;
	else if (iequals(s, "uchar"))
		result = DDL_PrimitiveType::UChar;
	else if (iequals(s, "numb"))
		result = DDL_PrimitiveType::Numb;
	else
		ec = make_error_code(validation_error::not_a_known_primitive_type);
	return result;
}

DDL_PrimitiveType map_to_primitive_type(std::string_view s)
{
	std::error_code ec;
	auto result = map_to_primitive_type(s, ec);
	if (ec)
		throw std::system_error(ec, std::string{ s });
	return result;
}

// --------------------------------------------------------------------

type_validator::type_validator(std::string_view name, DDL_PrimitiveType type, std::string_view rx)
	: m_name(name)
	, m_primitive_type(type)
	, m_rx(new regex_impl(rx.empty() ? ".+" : rx)) /// Empty regular expressions are not allowed, in libcpp's std::regex (POSIX?)
{
}

int type_validator::compare(std::string_view a, std::string_view b) const
{
	int result = 0;

	if (a.empty())
		result = b.empty() ? 0 : -1;
	else if (b.empty())
		result = a.empty() ? 0 : +1;
	else
	{
		switch (m_primitive_type)
		{
			case DDL_PrimitiveType::Numb:
			{
				double da, db;

				using namespace cif;
				using namespace std;

				std::from_chars_result ra, rb;

				ra = from_chars(a.data(), a.data() + a.length(), da);
				rb = from_chars(b.data(), b.data() + b.length(), db);

				if (ra.ec == std::errc{} and rb.ec == std::errc{})
				{
					auto d = da - db;
					if (std::abs(d) > std::numeric_limits<double>::epsilon())
					{
						if (d > 0)
							result = 1;
						else if (d < 0)
							result = -1;
					}
				}
				else if (ra.ec != std::errc{})
					result = 1;
				else
					result = -1;
				break;
			}

			case DDL_PrimitiveType::UChar:
			case DDL_PrimitiveType::Char:
			{
				// CIF is guaranteed to have ascii only, therefore this primitive code will do
				// also, we're collapsing spaces

				auto ai = a.begin(), bi = b.begin();
				for (;;)
				{
					if (ai == a.end())
					{
						if (bi != b.end())
							result = -1;
						break;
					}
					else if (bi == b.end())
					{
						result = 1;
						break;
					}

					char ca = *ai;
					char cb = *bi;

					if (m_primitive_type == DDL_PrimitiveType::UChar)
					{
						ca = tolower(ca);
						cb = tolower(cb);
					}

					result = ca - cb;

					if (result != 0)
						break;

					if (ca == ' ')
					{
						while (ai[1] == ' ')
							++ai;
						while (bi[1] == ' ')
							++bi;
					}

					++ai;
					++bi;
				}

				break;
			}
		}
	}

	return result;
}

// --------------------------------------------------------------------

void item_validator::operator()(std::string_view value) const
{
	std::error_code ec;
	if (not validate_value(value, ec))
		throw std::system_error(ec, std::format("'{}' is not a valid value for {}", value, m_item_name));
}

bool item_validator::validate_value(std::string_view value, std::error_code &ec) const noexcept
{
	ec.clear();

	if (not value.empty() and value != "?" and value != ".")
	{
		if (m_type != nullptr and not m_type->m_rx->match(value))
			ec = make_error_code(validation_error::value_does_not_match_rx);
		else if (not m_enums.empty() and m_enums.count(std::string{ value }) == 0)
			ec = make_error_code(validation_error::value_is_not_in_enumeration_list);
	}

	return ec == std::errc{};
}

// --------------------------------------------------------------------

void category_validator::add_item_validator(item_validator &&v)
{
	if (v.m_mandatory)
		m_mandatory_items.insert(v.m_item_name);

	v.m_category = m_name;

	auto i = std::ranges::find(m_item_validators, v);
	if (i != m_item_validators.end())
	{
		if (VERBOSE >= 4)
			std::cout << "Could not add validator for item " << v.m_item_name << " to category " << m_name << '\n';
	}
	else
		m_item_validators.emplace_back(std::move(v));
}

const item_validator *category_validator::get_validator_for_item(std::string_view item_name) const
{
	const item_validator *result = nullptr;
	auto i = std::ranges::find(m_item_validators, item_validator{ std::string(item_name) });
	if (i != m_item_validators.end())
		result = &*i;
	else if (VERBOSE > 4)
		std::cout << "No validator for item " << item_name << '\n';
	return result;
}

const item_validator *category_validator::get_validator_for_aliased_item(std::string_view item_name) const
{
	const item_validator *result = nullptr;

	for (auto &iv : m_item_validators)
	{
		for (auto &ai : iv.m_aliases)
		{
			const auto &[cat, name] = split_item_name(ai.m_name);
			if (iequals(name, item_name) and iequals(cat, m_name))
			{
				result = &iv;
				break;
			}
		}
		if (result)
			break;
	}

	return result;
}

// --------------------------------------------------------------------

void swap(validator &a, validator &b) noexcept
{
	std::swap(a.m_audit_conform, b.m_audit_conform);
	std::swap(a.m_strict, b.m_strict);
	std::swap(a.m_type_validators, b.m_type_validators);
	std::swap(a.m_category_validators, b.m_category_validators);
	std::swap(a.m_link_validators, b.m_link_validators);
}

void validator::parse(std::istream &is)
{
	parse_dictionary(*this, is);
}

void validator::add_type_validator(type_validator &&v)
{
	m_type_validators.emplace(v);
}

const type_validator *validator::get_validator_for_type(std::string_view typeCode) const
{
	const type_validator *result = nullptr;

	auto i = m_type_validators.find(type_validator{ std::string(typeCode), DDL_PrimitiveType::Char, {} });
	if (i != m_type_validators.end())
		result = &*i;
	return result;
}

void validator::add_category_validator(category_validator &&v)
{
	m_category_validators.emplace(v);
}

const category_validator *validator::get_validator_for_category(std::string_view category) const
{
	const category_validator *result = nullptr;
	auto i = m_category_validators.find(category_validator{ std::string(category) });
	if (i != m_category_validators.end())
		result = &*i;
	return result;
}

item_validator *validator::get_validator_for_item(std::string_view item_name) const
{
	item_validator *result = nullptr;

	std::string cat, item;
	std::tie(cat, item) = split_item_name(item_name);

	auto *cv = get_validator_for_category(cat);
	if (cv != nullptr)
		result = const_cast<item_validator *>(cv->get_validator_for_item(item));

	if (result == nullptr and VERBOSE > 4)
		std::cout << "No validator for item " << item_name << '\n';

	return result;
}

void validator::add_link_validator(link_validator &&v)
{
	assert(v.m_parent_keys.size() == v.m_child_keys.size());
	if (v.m_parent_keys.size() != v.m_child_keys.size())
		throw std::runtime_error("unequal number of keys for parent and child in link");

	auto pcv = get_validator_for_category(v.m_parent_category);
	auto ccv = get_validator_for_category(v.m_child_category);

	if (pcv == nullptr)
		throw std::runtime_error("unknown parent category " + v.m_parent_category);

	if (ccv == nullptr)
		throw std::runtime_error("unknown child category " + v.m_child_category);

	for (std::size_t i = 0; i < v.m_parent_keys.size(); ++i)
	{
		auto piv = pcv->get_validator_for_item(v.m_parent_keys[i]);

		if (piv == nullptr)
			throw std::runtime_error("unknown parent item _" + v.m_parent_category + '.' + v.m_parent_keys[i]);

		auto civ = ccv->get_validator_for_item(v.m_child_keys[i]);
		if (civ == nullptr)
			throw std::runtime_error("unknown child item _" + v.m_child_category + '.' + v.m_child_keys[i]);

		if (civ->m_type == nullptr and piv->m_type != nullptr)
			const_cast<item_validator *>(civ)->m_type = piv->m_type;
	}

	m_link_validators.emplace_back(std::move(v));
}

std::vector<const link_validator *> validator::get_links_for_parent(std::string_view category) const
{
	std::vector<const link_validator *> result;

	for (auto &l : m_link_validators)
	{
		if (iequals(l.m_parent_category, category))
			result.push_back(&l);
	}

	return result;
}

std::vector<const link_validator *> validator::get_links_for_child(std::string_view category) const
{
	std::vector<const link_validator *> result;

	for (auto &l : m_link_validators)
	{
		if (iequals(l.m_child_category, category))
			result.push_back(&l);
	}

	return result;
}

void validator::report_error(std::error_code ec, bool fatal) const
{
	if (m_strict or fatal)
		throw validation_exception(ec);
	else if (VERBOSE > 0)
		std::cerr << ec.message() << '\n';
}

void validator::report_error(std::error_code ec, std::string_view category,
	std::string_view item, bool fatal) const
{
	if (m_strict or fatal)
	{
		if (item.empty())
			throw validation_exception(ec, category);
		else
		 	throw validation_exception(ec, category, item);
	}

	if (VERBOSE > 0)
		std::cerr << ec.message() << " category: " << std::quoted(category) << " item: " << std::quoted(item) << '\n';
}

void validator::fill_audit_conform(category &audit_conform) const
{
	audit_conform.clear();
	audit_conform.emplace(m_audit_conform.begin(), m_audit_conform.end());
}

bool validator::matches_audit_conform(const category &audit_conform) const
{
	if (audit_conform.empty())
		return false;

	auto ai = m_audit_conform.begin();
	auto bi = audit_conform.begin();

	while (ai != m_audit_conform.end() and bi != audit_conform.end())
	{
		const auto &[name_a, version_a] = ai->get<std::string, std::optional<std::string>>("dict_name", "dict_version");
		const auto &[name_b, version_b] = bi->get<std::string, std::optional<std::string>>("dict_name", "dict_version");

		++ai;
		++bi;

		if (name_a != name_b)
			return false;

		if (not version_b.has_value() or not version_a.has_value())
			continue;

		if (validator_factory::check_version(name_a, *version_b, *version_a) == false)
			return false;
	}

	return ai == m_audit_conform.end() and bi == audit_conform.end();
}

void validator::append_audit_conform(const std::string &name, const std::optional<std::string> &version)
{
	m_audit_conform.emplace({ //
		{ "dict_name", name },
		{ "dict_version", version } });
}

// --------------------------------------------------------------------

validator_factory &validator_factory::instance()
{
	static validator_factory s_instance;
	return s_instance;
}

const validator *validator_factory::get(std::string_view dictionary_name)
{
	category audit_conform("audit_conform");
	for (auto part : cif::split(dictionary_name, ";,", true))
		audit_conform.emplace({ { "dict_name", part } });

	return get(audit_conform);
}

const validator *validator_factory::get(const category &audit_conform)
{
	const validator *result = nullptr;

	std::scoped_lock lock(m_mutex);

	// Check existing first
	for (auto &v : m_validators)
	{
		if (v.matches_audit_conform(audit_conform))
			result = &v;
	}

	// If the audit conform contains only one record, this is easy
	if (result == nullptr and audit_conform.size() == 1)
	{
		const auto &[name, version] =
			audit_conform.front().get<std::string, std::optional<std::string>>("dict_name", "dict_version");
		if (not name.empty())
			result = &m_validators.emplace_back(construct_validator(name, version));
	}

	if (result == nullptr)
	{
		// A new, merged dictionary
		std::optional<validator> v;
		for (const auto &[name, version] : audit_conform.rows<std::string, std::optional<std::string>>("dict_name", "dict_version"))
		{
			if (name.empty())
				continue;

			if (not v) // first dict
				v = construct_validator(name, version);
			else // additional/extending dict
			{
				auto data = load_resource(name);
				if (not data)
					throw std::runtime_error("Could not load dictionary " + std::string{ name });

				v->parse(*data);
			}
		}

		if (v)
			result = &m_validators.emplace_back(std::move(*v));
	}

	return result;
}

const validator &validator_factory::operator[](const category &audit_conform)
{
	auto v = get(audit_conform);
	if (v == nullptr)
		throw std::runtime_error("Could not load dictionary for audit_conform");
	return *v;
}

const validator &validator_factory::operator[](std::string_view dictionary_name)
{
	auto v = get(dictionary_name);
	if (v == nullptr)
		throw std::runtime_error("Could not load dictionary for " + std::string{ dictionary_name });
	return *v;
}

validator validator_factory::construct_validator(std::string_view name, std::optional<std::string> version)
{
	auto data = load_resource(name);
	if (not data and name == "mmcif_pdbx_v50")
		data = load_resource("mmcif_pdbx.dic");

	if (not data)
		throw std::runtime_error("Could not load dictionary " + std::string{ name });

	validator v;
	v.parse(*data);

	if (version.has_value() and VERBOSE >= 0 and
		not v.matches_audit_conform(category{ "audit_conform", //
			{ { "dict_name", name }, { "dict_version", version } } }))
	{
		std::clog << "Loaded dictionary does not match name=" << name << " and version=" << version.value_or("''") << "\n";
	}

	return v;
}

bool validator_factory::check_version(std::string_view name, std::string_view expected, std::string_view found)
{
	bool result = true;
	auto el = cif::split(expected, ".");
	auto fl = cif::split(found, ".");

	auto eli = el.begin();
	auto fli = fl.begin();

	while (eli != el.end() and fli != fl.end())
	{
		int e_int, f_int;
		if (auto [ptr, ec] = std::from_chars(eli->data(), eli->data() + eli->length(), e_int); ec != std::errc{})
		{
			std::clog << "Could not parse requested version string for dictionary " << std::quoted(expected) << "\n";
			result = false;
			break;
		}

		if (auto [ptr, ec] = std::from_chars(fli->data(), fli->data() + fli->length(), f_int); ec != std::errc{})
		{
			std::clog << "Could not parse version string in dictionary " << name << " " << std::quoted(found) << "\n";
			result = false;
			break;
		}

		if (f_int > e_int) // newer version, assume this is ok
			break;

		if (f_int < e_int)
		{
			std::clog << "The version in dictionary " << name << " is lower than requested, this may cause validation errors\n";
			result = false;
			break;
		}

		++eli;
		++fli;
	}

	return result;
}

} // namespace cif
