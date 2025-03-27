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
#include "cif++/gzio.hpp"
#include "cif++/utilities.hpp"

#include <cassert>
#include <fstream>
#include <iostream>

// The validator depends on regular expressions. Unfortunately,
// the implementation of std::regex in g++ is buggy and crashes
// on reading the pdbx dictionary. Therefore, in case g++ is used
// the code will use boost::regex instead.

#if USE_BOOST_REGEX
# include <boost/regex.hpp>
using boost::regex;
#else
# include <regex>
using std::regex;
#endif

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

struct regex_impl : public regex
{
	regex_impl(std::string_view rx)
		: regex(rx.begin(), rx.end(), regex::extended | regex::optimize)
	{
	}
};

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

type_validator::~type_validator()
{
	delete m_rx;
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

				ra = selected_charconv<double>::from_chars(a.data(), a.data() + a.length(), da);
				rb = selected_charconv<double>::from_chars(b.data(), b.data() + b.length(), db);

				if (not(bool) ra.ec and not(bool) rb.ec)
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
				else if ((bool)ra.ec)
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
		throw std::system_error(ec, std::string{ value } + " does not match rx for " + m_item_name);
}

bool item_validator::validate_value(std::string_view value, std::error_code &ec) const noexcept
{
	ec.clear();

	if (not value.empty() and value != "?" and value != ".")
	{
		if (m_type != nullptr and not regex_match(value.begin(), value.end(), *m_type->m_rx))
			ec = make_error_code(validation_error::value_does_not_match_rx);
		else if (not m_enums.empty() and m_enums.count(std::string{ value }) == 0)
			ec = make_error_code(validation_error::value_is_not_in_enumeration_list);
	}

	return not(bool) ec;
}

// --------------------------------------------------------------------

void category_validator::add_item_validator(item_validator &&v)
{
	if (v.m_mandatory)
		m_mandatory_items.insert(v.m_item_name);

	v.m_category = this;

	auto r = m_item_validators.insert(std::move(v));
	if (not r.second and VERBOSE >= 4)
		std::cout << "Could not add validator for item " << v.m_item_name << " to category " << m_name << '\n';
}

const item_validator *category_validator::get_validator_for_item(std::string_view item_name) const
{
	const item_validator *result = nullptr;
	auto i = m_item_validators.find(item_validator{ std::string(item_name) });
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

void validator_base::report_error(std::error_code ec, bool fatal) const
{
	if (m_strict or fatal)
		throw validation_exception(ec);
	else if (VERBOSE > 0)
		std::cerr << ec.message() << '\n';
}

void validator_base::report_error(std::error_code ec, std::string_view category,
	std::string_view item, bool fatal) const
{
	auto ex = item.empty() ? validation_exception(ec, category) : validation_exception(ec, category, item);

	if (m_strict or fatal)
		throw ex;
	else if (VERBOSE > 0)
		std::cerr << ex.what() << '\n';
}

// --------------------------------------------------------------------

void validator::add_type_validator(type_validator &&v)
{
	auto r = m_type_validators.insert(std::move(v));
	if (not r.second and VERBOSE > 4)
		std::cout << "Could not add validator for type " << v.m_name << '\n';
}

const type_validator *validator::get_validator_for_type(std::string_view typeCode) const
{
	const type_validator *result = nullptr;

	auto i = m_type_validators.find(type_validator{ std::string(typeCode), DDL_PrimitiveType::Char, {} });
	if (i != m_type_validators.end())
		result = &*i;
	else if (VERBOSE > 4)
		std::cout << "No validator for type " << typeCode << '\n';
	return result;
}

void validator::add_category_validator(category_validator &&v)
{
	auto r = m_category_validators.insert(std::move(v));
	if (not r.second and VERBOSE > 4)
		std::cout << "Could not add validator for category " << v.m_name << '\n';
}

const category_validator *validator::get_validator_for_category(std::string_view category) const
{
	const category_validator *result = nullptr;
	auto i = m_category_validators.find(category_validator{ std::string(category) });
	if (i != m_category_validators.end())
		result = &*i;
	else if (VERBOSE > 4)
		std::cout << "No validator for category " << category << '\n';
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

// --------------------------------------------------------------------

extended_validator::extended_validator(std::vector<const validator *> validators)
	: m_validators(validators)
{
	std::vector<std::string> names, versions;

	for (auto v : m_validators)
	{
		names.emplace_back(v->name());
		versions.emplace_back(v->version());
		m_strict = m_strict or v->is_strict();
	}

	m_name = cif::join(names, "; ");
	m_version = cif::join(versions, "; ");
}

const type_validator *extended_validator::get_validator_for_type(std::string_view type_code) const
{
	const type_validator *result = nullptr;

	for (auto v : m_validators)
	{
		result = v->get_validator_for_type(type_code);
		if (result)
			break;
	}

	return result;
}

const category_validator *extended_validator::get_validator_for_category(std::string_view category) const
{
	const category_validator *result = nullptr;

	for (auto v : m_validators)
	{
		result = v->get_validator_for_category(category);
		if (result)
			break;
	}

	return result;
}

std::vector<const link_validator *> extended_validator::get_links_for_parent(std::string_view category) const
{
	std::vector<const link_validator *> result;

	for (auto v : m_validators)
	{
		auto links = v->get_links_for_parent(category);
		result.insert(result.end(), links.begin(), links.end());
	}

	return result;
}

std::vector<const link_validator *> extended_validator::get_links_for_child(std::string_view category) const
{
	std::vector<const link_validator *> result;

	for (auto v : m_validators)
	{
		auto links = v->get_links_for_child(category);
		result.insert(result.end(), links.begin(), links.end());
	}

	return result;
}

// --------------------------------------------------------------------

validator_factory &validator_factory::instance()
{
	static validator_factory s_instance;
	return s_instance;
}

const validator_base &validator_factory::operator[](std::string_view dictionary_name)
{
	category audit_conform("audit_conform");
	for (auto part : cif::split(dictionary_name, ";", true))
		audit_conform.emplace({ { "dict_name", part } });

	return construct_validator(audit_conform);
}

const validator_base &validator_factory::construct_validator(const category &audit_conform)
{
	if (audit_conform.empty())
		throw std::runtime_error("Empty audit_conform category, cannot create a validator");

	std::lock_guard lock(m_mutex);

	std::vector<const validator *> validators;

	for (const auto &[name, version] : audit_conform.rows<std::string, std::optional<std::string>>("dict_name", "dict_version"))
	{
		auto &v = construct_validator(name, version);
		validators.emplace_back(&v);
	}

	if (validators.size() == 1)
		return *validators.front();

	// override mode, last dictionary is most important
	std::reverse(validators.begin(), validators.end());
		
	for (auto &ev : m_extended_validators)
	{
		if (ev.m_validators == validators)
			return ev;
	}

	return m_extended_validators.emplace_back(validators);
}

const validator &validator_factory::construct_validator(std::string_view name,
	std::optional<std::string> version)
{
	for (auto &v : m_validators)
	{
		if (version.has_value())
			check_version(name, *version, v.version());

		if (v.name() == name)
			return v;
	}

	std::filesystem::path dictionary(name);

	auto data = load_resource(name);

	if (not data and dictionary.extension().string() != ".dic")
		data = load_resource(dictionary.parent_path() / (dictionary.filename().string() + ".dic"));

	if (not data)
	{
		std::error_code ec;

		// might be a compressed dictionary on disk
		std::filesystem::path p = dictionary;
		if (p.extension() == ".dic")
			p = p.parent_path() / (p.filename().string() + ".gz");
		else
			p = p.parent_path() / (p.filename().string() + ".dic.gz");

#if defined(CACHE_DIR) or defined(DATA_DIR)
		if (not std::filesystem::exists(p, ec) or ec)
		{
			for (const char *dir : {
# if defined(CACHE_DIR)
					 CACHE_DIR,
# endif
# if defined(DATA_DIR)
					 DATA_DIR
# endif
				 })
			{
				auto p2 = std::filesystem::path(dir) / p;
				if (std::filesystem::exists(p2, ec) and not ec)
				{
					swap(p, p2);
					break;
				}
			}
		}
#endif

		if (std::filesystem::exists(p, ec) and not ec)
		{
			auto in = std::make_unique<gzio::ifstream>(p);

			if (not in->is_open())
				throw std::runtime_error("Could not open dictionary (" + p.string() + ")");

			data.reset(in.release());
		}
		else
			throw std::runtime_error("Dictionary not found or defined (" + dictionary.string() + ")");
	}

	return construct_validator(name, version, *data);
}

const validator &validator_factory::construct_validator(std::string_view name,
	std::optional<std::string> version, std::istream &is)
{
	auto v = parse_dictionary(name, is);

	if (version.has_value() and VERBOSE >= 0)
	{
		auto vv = v.version();

		if (vv.empty())
			std::clog << "Could not check version of dictionary " << name << " since this info is missing\n";
		else
			check_version(name, *version, vv);
	}

	return m_validators.emplace_back(std::move(v));
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
		if (auto [ptr, ec] = std::from_chars(eli->begin(), eli->end(), e_int); ec != std::errc{})
		{
			std::clog << "Could not parse requested version string for dictionary " << std::quoted(expected) << "\n";
			result = false;
			break;
		}

		if (auto [ptr, ec] = std::from_chars(fli->begin(), fli->end(), f_int); ec != std::errc{})
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
