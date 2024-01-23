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

#include "cif++/category.hpp"
#include "cif++/condition.hpp"

namespace cif
{

iset get_category_items(const category &cat)
{
	return cat.key_items();
}

uint16_t get_item_ix(const category &cat, std::string_view col)
{
	return cat.get_item_ix(col);
}

bool is_item_type_uchar(const category &cat, std::string_view col)
{
	bool result = false;

	auto cv = cat.get_cat_validator();
	if (cv)
	{
		auto iv = cv->get_validator_for_item(col);
		if (iv != nullptr and iv->m_type != nullptr)
		{
			auto type = iv->m_type;
			result = type->m_primitive_type == DDL_PrimitiveType::UChar;
		}
	}

	return result;
}

namespace detail
{

	condition_impl *key_equals_condition_impl::prepare(const category &c)
	{
		m_item_ix = c.get_item_ix(m_item_name);
		m_icase = is_item_type_uchar(c, m_item_name);

		if (c.get_cat_validator() != nullptr and
			c.key_item_indices().contains(m_item_ix) and
			c.key_item_indices().size() == 1)
		{
			m_single_hit = c[{ { m_item_name, m_value } }];
		}

		return this;
	}

	bool found_in_range(condition_impl *c, std::vector<and_condition_impl *>::iterator b, std::vector<and_condition_impl *>::iterator e)
	{
		bool result = true;

		for (auto s = b; s != e; ++s)
		{
			auto &cs = (*s)->m_sub;

			if (find_if(cs.begin(), cs.end(), [c](const condition_impl *i) { return i->equals(c); }) == cs.end())
			{
				result = false;
				break;
			}
		}

		return result;
	}

	condition_impl *and_condition_impl::combine_equal(std::vector<and_condition_impl *> &subs, or_condition_impl *oc)
	{
		and_condition_impl *and_result = nullptr;

		auto first = subs.front();
		auto &fc = first->m_sub;

		for (auto c : fc)
		{
			if (not found_in_range(c, subs.begin() + 1, subs.end()))
				continue;

			if (and_result == nullptr)
				and_result = new and_condition_impl();

			and_result->m_sub.push_back(c);
			fc.erase(remove(fc.begin(), fc.end(), c), fc.end());

			for (auto sub : subs)
			{
				auto &ssub = sub->m_sub;

				for (auto sc : ssub)
				{
					if (not sc->equals(c))
						continue;
					
					ssub.erase(remove(ssub.begin(), ssub.end(), sc), ssub.end());
					delete sc;
					break;
				}
			}
		}

		if (and_result != nullptr)
		{
			and_result->m_sub.push_back(oc);
			return and_result;
		}

		return oc;
	}

	condition_impl *or_condition_impl::prepare(const category &c)
	{
		std::vector<and_condition_impl *> and_conditions;

		for (auto &sub : m_sub)
		{
			sub = sub->prepare(c);
			if (typeid(*sub) == typeid(and_condition_impl))
				and_conditions.push_back(static_cast<and_condition_impl *>(sub));
		}

		if (and_conditions.size() == m_sub.size())
			return and_condition_impl::combine_equal(and_conditions, this);

		return this;
	}

} // namespace detail

void condition::prepare(const category &c)
{
	if (m_impl)
		m_impl = m_impl->prepare(c);
	
	m_prepared = true;
}

} // namespace cif
