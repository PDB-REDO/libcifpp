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

#include <cassert>
#include <fstream>
#include <iostream>

#include <cif++/v2/dictionary_parser.hpp>
#include <cif++/v2/validate.hpp>

#include <cif++/CifUtils.hpp>

namespace cif
{
extern int VERBOSE;
}

namespace cif::v2
{

using cif::VERBOSE;

validation_error::validation_error(const std::string &msg)
	: m_msg(msg)
{
}

validation_error::validation_error(const std::string &cat, const std::string &item, const std::string &msg)
	: m_msg("When validating _" + cat + '.' + item + ": " + msg)
{
}

// --------------------------------------------------------------------

DDL_PrimitiveType map_to_primitive_type(std::string_view s)
{
	DDL_PrimitiveType result;
	if (iequals(s, "char"))
		result = DDL_PrimitiveType::Char;
	else if (iequals(s, "uchar"))
		result = DDL_PrimitiveType::UChar;
	else if (iequals(s, "numb"))
		result = DDL_PrimitiveType::Numb;
	else
		throw validation_error("Not a known primitive type");
	return result;
}

// --------------------------------------------------------------------

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

				auto ra = cif::from_chars(a.begin(), a.end(), da);
				auto rb = cif::from_chars(b.begin(), b.end(), db);

				if (ra.ec == std::errc() and rb.ec == std::errc())
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
				else if (ra.ec == std::errc())
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

// void ValidateItem::addLinked(ValidateItem* parent, const std::string& parentItem, const std::string& childItem)
//{
////	if (mParent != nullptr and VERBOSE)
////		cerr << "replacing parent in " << mCategory->m_name << " from " << mParent->mCategory->m_name << " to " << parent->mCategory->m_name << endl;
////	mParent = parent;
//
//	if (mType == nullptr and parent != nullptr)
//		mType = parent->mType;
//
//	if (parent != nullptr)
//	{
//		mLinked.push_back({parent, parentItem, childItem});
//
//		parent->mChildren.insert(this);
////
////		if (mCategory->mKeys == std::vector<std::string>{mTag})
////			parent->mForeignKeys.insert(this);
//	}
//}

void item_validator::operator()(std::string_view value) const
{
	if (not value.empty() and value != "?" and value != ".")
	{
		if (m_type != nullptr and not regex_match(value.begin(), value.end(), m_type->m_rx))
			throw validation_error(m_category->m_name, m_tag, "Value '" + std::string{value} + "' does not match type expression for type " + m_type->m_name);

		if (not m_enums.empty())
		{
			if (m_enums.count(std::string{value}) == 0)
				throw validation_error(m_category->m_name, m_tag, "Value '" + std::string{value} + "' is not in the list of allowed values");
		}
	}
}

// --------------------------------------------------------------------

void category_validator::addItemValidator(item_validator &&v)
{
	if (v.m_mandatory)
		m_mandatory_fields.insert(v.m_tag);

	v.m_category = this;

	auto r = m_item_validators.insert(std::move(v));
	if (not r.second and VERBOSE >= 4)
		std::cout << "Could not add validator for item " << v.m_tag << " to category " << m_name << std::endl;
}

const item_validator *category_validator::get_validator_for_item(std::string_view tag) const
{
	const item_validator *result = nullptr;
	auto i = m_item_validators.find(item_validator{std::string(tag)});
	if (i != m_item_validators.end())
		result = &*i;
	else if (VERBOSE > 4)
		std::cout << "No validator for tag " << tag << std::endl;
	return result;
}

// --------------------------------------------------------------------

void validator::add_type_validator(type_validator &&v)
{
	auto r = m_type_validators.insert(std::move(v));
	if (not r.second and VERBOSE > 4)
		std::cout << "Could not add validator for type " << v.m_name << std::endl;
}

const type_validator *validator::get_validator_for_type(std::string_view typeCode) const
{
	const type_validator *result = nullptr;

	auto i = m_type_validators.find(type_validator{std::string(typeCode), DDL_PrimitiveType::Char, boost::regex()});
	if (i != m_type_validators.end())
		result = &*i;
	else if (VERBOSE > 4)
		std::cout << "No validator for type " << typeCode << std::endl;
	return result;
}

void validator::add_category_validator(category_validator &&v)
{
	auto r = m_category_validators.insert(std::move(v));
	if (not r.second and VERBOSE > 4)
		std::cout << "Could not add validator for category " << v.m_name << std::endl;
}

const category_validator *validator::get_validator_for_category(std::string_view category) const
{
	const category_validator *result = nullptr;
	auto i = m_category_validators.find(category_validator{std::string(category)});
	if (i != m_category_validators.end())
		result = &*i;
	else if (VERBOSE > 4)
		std::cout << "No validator for category " << category << std::endl;
	return result;
}

item_validator *validator::get_validator_for_item(std::string_view tag) const
{
	item_validator *result = nullptr;

	std::string cat, item;
	std::tie(cat, item) = splitTagName(tag);

	auto *cv = get_validator_for_category(cat);
	if (cv != nullptr)
		result = const_cast<item_validator *>(cv->get_validator_for_item(item));

	if (result == nullptr and VERBOSE > 4)
		std::cout << "No validator for item " << tag << std::endl;

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

	for (size_t i = 0; i < v.m_parent_keys.size(); ++i)
	{
		auto piv = pcv->get_validator_for_item(v.m_parent_keys[i]);

		if (piv == nullptr)
			throw std::runtime_error("unknown parent tag _" + v.m_parent_category + '.' + v.m_parent_keys[i]);

		auto civ = ccv->get_validator_for_item(v.m_child_keys[i]);
		if (civ == nullptr)
			throw std::runtime_error("unknown child tag _" + v.m_child_category + '.' + v.m_child_keys[i]);

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
		if (l.m_parent_category == category)
			result.push_back(&l);
	}

	return result;
}

std::vector<const link_validator *> validator::get_links_for_child(std::string_view category) const
{
	std::vector<const link_validator *> result;

	for (auto &l : m_link_validators)
	{
		if (l.m_child_category == category)
			result.push_back(&l);
	}

	return result;
}

void validator::report_error(const std::string &msg, bool fatal) const
{
	if (m_strict or fatal)
		throw validation_error(msg);
	else if (VERBOSE > 0)
		std::cerr << msg << std::endl;
}

// --------------------------------------------------------------------

const validator &validator_factory::operator[](std::string_view dictionary_name)
{
	std::lock_guard lock(m_mutex);

	for (auto &validator : m_validators)
	{
		if (iequals(validator.name(), dictionary_name))
			return validator;
	}

	// not found, add it

	// too bad clang version 10 did not have a constructor for std::filesystem::path that accepts a std::string_view
	std::filesystem::path dictionary(dictionary_name.data(), dictionary_name.data() + dictionary_name.length());

	auto data = loadResource(dictionary_name);

	if (not data and dictionary.extension().string() != ".dic")
		data = loadResource(dictionary.parent_path() / (dictionary.filename().string() + ".dic"));

	if (data)
		construct_validator(dictionary_name, *data);
	else
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
#if defined(CACHE_DIR)
					 CACHE_DIR,
#endif
#if defined(DATA_DIR)
						 DATA_DIR
#endif
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
			std::ifstream file(p, std::ios::binary);
			if (not file.is_open())
				throw std::runtime_error("Could not open dictionary (" + p.string() + ")");

			io::filtering_stream<io::input> in;
			in.push(io::gzip_decompressor());
			in.push(file);

			construct_validator(dictionary_name, in);
		}
		else
			throw std::runtime_error("Dictionary not found or defined (" + dictionary.string() + ")");
	}

	assert(iequals(m_validators.back().name(), dictionary_name));

	return m_validators.back();
}

void validator_factory::construct_validator(std::string_view name, std::istream &is)
{
	parse_dictionary(name, is);
}

} // namespace cif::v2
