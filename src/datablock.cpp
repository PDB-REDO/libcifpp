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

#include "cif++/datablock.hpp"

#include "cif++/validate.hpp"
#include "cif++/utilities.hpp"

#include <algorithm>
#include <exception>
#include <ranges>

namespace cif
{

datablock::datablock(const datablock &db)
	: std::list<category>(db)
	, m_name(db.m_name)
	, m_validator(db.m_validator)
{
	for (auto &cat : *this)
		cat.update_links(*this);
}

void datablock::load_dictionary()
{
	if (auto *audit_conform = get("audit_conform"); audit_conform and not audit_conform->empty())
	{
		try
		{
			set_validator(validator_factory::instance().get(*audit_conform));
		}
		catch (const std::exception &ex)
		{
			std::clog << ex.what() << '\n';
		}
	}
}

void datablock::load_dictionary(std::string_view dict)
{
	try
	{
		set_validator(validator_factory::instance().get(dict));
	}
	catch (const std::exception &ex)
	{
		std::clog << ex.what() << '\n';
	}
}

void datablock::set_validator(const validator *v)
{
	m_validator = v;

	try
	{
		for (auto &cat : *this)
			cat.set_validator(v, *this);
	}
	catch (const std::exception &)
	{
		m_validator = nullptr;
		throw_with_nested(std::runtime_error("Error while setting validator in datablock " + m_name));
	}
}

const validator *datablock::get_validator() const
{
	return m_validator;
}

bool datablock::is_valid() const
{
	if (m_validator == nullptr)
		throw std::runtime_error("Validator not specified for datablock data_" + name());

	bool result = true;
	for (auto &cat : *this)
		result = cat.is_valid() and result;

	return result;
}

bool datablock::validate_links() const
{
	bool result = true;

	for (auto &cat : *this)
		const_cast<category &>(cat).update_links(*this);

	for (auto &cat : *this)
		result = cat.validate_links() and result;

	return result;
}

bool datablock::strip()
{
	bool result = true;

	// remove all categories that have no validator
	std::erase_if(*this, [](category &c)
		{
		bool result = false;
		if (c.get_cat_validator() == nullptr)
		{
			if (VERBOSE > 0)
				std::clog << "Dropping category " << c.name() << '\n';
			result = true;
		}
		return result; });

	// then strip the remaining categories
	for (auto &cat : *this)
		cat.strip();

	// Add or remove the audit_conform block here.
	if (is_valid())
	{
		// If the dictionary declares an audit_conform category, put it in,
		// but only if it does not exist already!

		if (auto audit_conform = get("audit_conform");
			audit_conform != nullptr and m_validator->get_validator_for_category("audit_conform") != nullptr)
		{
			audit_conform->clear();
			m_validator->fill_audit_conform(*audit_conform);
		}
	}
	else
		result = false;

	return result;
}

// --------------------------------------------------------------------

category &datablock::operator[](std::string_view name)
{
	auto i = std::ranges::find_if(*this, [name](const category &c)
		{ return iequals(c.name(), name); });

	if (i != end())
		return *i;

	auto &cat = emplace_back(name);

	if (m_validator)
		cat.set_validator(m_validator, *this);

	return back();
}

const category &datablock::operator[](std::string_view name) const
{
	static const category s_empty;
	auto i = std::ranges::find_if(*this, [name](const category &c)
		{ return iequals(c.name(), name); });
	return i == end() ? s_empty : *i;
}

category *datablock::get(std::string_view name)
{
	auto i = std::ranges::find_if(*this, [name](const category &c)
		{ return iequals(c.name(), name); });
	return i == end() ? nullptr : &*i;
}

const category *datablock::get(std::string_view name) const
{
	return const_cast<datablock *>(this)->get(name);
}

std::tuple<datablock::iterator, bool> datablock::emplace(std::string_view name)
{
	bool is_new = true;

	auto i = begin();
	while (i != end())
	{
		if (iequals(name, i->name()))
		{
			is_new = false;
			break;
		}

		++i;
	}

	if (is_new)
	{
		i = insert(end(), { name });
		i->set_validator(m_validator, *this);
	}

	assert(i != end());

	// links may have changed...
	for (auto &cat : *this)
		cat.update_links(*this);

	return std::make_tuple(i, is_new);
}

std::vector<std::string> datablock::get_item_order() const
{
	std::vector<std::string> result;

	// for entry and audit_conform on top
	auto ci = std::ranges::find_if(*this, [](const category &cat)
		{ return cat.name() == "entry"; });
	if (ci != end())
	{
		auto cto = ci->get_item_order();
		result.insert(result.end(), cto.begin(), cto.end());
	}

	ci = std::ranges::find_if(*this, [](const category &cat)
		{ return cat.name() == "audit_conform"; });
	if (ci != end())
	{
		auto cto = ci->get_item_order();
		result.insert(result.end(), cto.begin(), cto.end());
	}

	for (auto &cat : *this)
	{
		if (cat.name() == "entry" or cat.name() == "audit_conform")
			continue;
		auto cto = cat.get_item_order();
		result.insert(result.end(), cto.begin(), cto.end());
	}

	return result;
}

namespace
{
	using elem_t = std::tuple<std::string, int, bool>;
	using cat_order_t = std::vector<elem_t>;
	using iter_t = cat_order_t::iterator;

	inline int get_count(iter_t i)
	{
		return std::get<1>(*i);
	}

	inline bool is_on_stack(iter_t i)
	{
		return std::get<2>(*i);
	}

	void calculate_cat_order(cat_order_t &cat_order, iter_t i, const validator &validator)
	{
		if (i == cat_order.end() or get_count(i) >= 0)
			return;

		auto &&[cat, count, on_stack] = *i;

		on_stack = true;

		int parent_count = 0;

		for (auto link : validator.get_links_for_child(cat))
		{
			auto ei = std::ranges::find_if(cat_order, [parent = link->m_parent_category](elem_t &a)
				{ return std::get<0>(a) == parent; });

			if (ei == cat_order.end())
				continue;

			if (not is_on_stack(ei))
				calculate_cat_order(cat_order, ei, validator);

			parent_count += get_count(ei);
		}

		count = parent_count + 1;
	}
} // namespace

void datablock::write(std::ostream &os) const
{
	os << "data_" << m_name << '\n'
	   << "# \n";

	if (m_validator and size() > 0)
	{
		// base order on parent child relationships, parents first

		cat_order_t cat_order;

		for (auto &cat : *this)
		{
			if (cat.name() == "entry" or cat.name() == "audit_conform")
				continue;
			cat_order.emplace_back(cat.name(), -1, false);
		}

		for (auto i = cat_order.begin(); i != cat_order.end(); ++i)
			calculate_cat_order(cat_order, i, *m_validator);

		std::ranges::sort(cat_order, [](const elem_t &a, const elem_t &b)
			{
			const auto &[cat_a, count_a, on_stack_a] = a;
			const auto &[cat_b, count_b, on_stack_b] = b;

			int d = std::get<1>(a) - std::get<1>(b);
			if (d == 0)
				d = cat_b.compare(cat_a);

			return d < 0; });

		if (auto entry = get("entry"); entry != nullptr)
			entry->write(os);

		if (auto audit_conform = get("audit_conform"); audit_conform != nullptr)
			audit_conform->write(os);

		for (auto &&[cat, count, on_stack] : cat_order)
			get(cat)->write(os);
	}
	else
	{
		// mmcif support, sort of. First write the 'entry' Category
		// and if it exists, _AND_ we have a Validator, write out the
		// audit_conform record.

		if (auto entry = get("entry"); entry != nullptr)
			entry->write(os);

		// If the dictionary declares an audit_conform category, put it in,
		// but only if it does not exist already!
		if (auto audit_conform = get("audit_conform"); audit_conform != nullptr)
			audit_conform->write(os);

		for (auto &cat : *this)
		{
			if (cat.name() != "entry" and cat.name() != "audit_conform")
				cat.write(os);
		}
	}
}

void datablock::write(std::ostream &os, const std::vector<std::string> &item_name_order)
{
	os << "data_" << m_name << '\n'
	   << "# \n";

	std::vector<std::string> cat_order{ "entry", "audit_conform" };
	for (auto &o : item_name_order)
	{
		std::string cat_name, item_name;
		std::tie(cat_name, item_name) = split_item_name(o);
		if (std::ranges::find_if(std::ranges::reverse_view(cat_order), [cat_name](const std::string &s) -> bool
				{ return iequals(cat_name, s); }) == cat_order.rend())
			cat_order.push_back(cat_name);
	}

	for (auto &c : cat_order)
	{
		auto cat = get(c);
		if (cat == nullptr)
			continue;

		std::vector<std::string> items;
		for (auto &o : item_name_order)
		{
			std::string cat_name, item_name;
			std::tie(cat_name, item_name) = split_item_name(o);

			if (cat_name == c)
				items.push_back(item_name);
		}

		cat->write(os, items);
	}

	// for any Category we missed in the catOrder
	for (auto &cat : *this)
	{
		if (std::ranges::find_if(cat_order, [&](const std::string &s) -> bool
				{ return iequals(cat.name(), s); }) != cat_order.end())
			continue;

		cat.write(os);
	}
}

bool datablock::operator==(const datablock &rhs) const
{
	// shortcut
	if (this == &rhs)
		return true;

	auto &dbA = *this;
	auto &dbB = rhs;

	std::vector<std::string> catA, catB;

	for (auto &cat : dbA)
	{
		if (not cat.empty())
			catA.push_back(cat.name());
	}
	std::ranges::sort(catA);

	for (auto &cat : dbB)
	{
		if (not cat.empty())
			catB.push_back(cat.name());
	}
	std::ranges::sort(catB);

	// loop over categories twice, to group output
	// First iteration is to list missing categories.

	std::vector<std::string> missingA, missingB;

	auto catA_i = catA.begin(), catB_i = catB.begin();

	while (catA_i != catA.end() and catB_i != catB.end())
	{
		if (not iequals(*catA_i, *catB_i))
			return false;

		++catA_i, ++catB_i;
	}

	if (catA_i != catA.end() or catB_i != catB.end())
		return false;

	// Second loop, now compare category values
	catA_i = catA.begin(), catB_i = catB.begin();

	while (catA_i != catA.end() and catB_i != catB.end())
	{
		std::string nA = *catA_i;
		to_lower(nA);

		std::string nB = *catB_i;
		to_lower(nB);

		int d = nA.compare(nB);
		if (d > 0)
			++catB_i;
		else if (d < 0)
			++catA_i;
		else
		{
			if (not(*dbA.get(*catA_i) == *dbB.get(*catB_i)))
				return false;
			++catA_i;
			++catB_i;
		}
	}

	return true;
}

} // namespace cif
