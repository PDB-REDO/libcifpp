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

#include "cif++/condition.hpp"

#include "cif++/category.hpp"
#include "cif++/validate.hpp"

#include <algorithm>

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
	// 	// index lookup
	// 	struct index_lookup_condition_impl : public condition_impl
	// 	{
	// 		index_lookup_condition_impl(row_initializer &&key_values)
	// 			: m_key_values(std::move(key_values))
	// 		{
	// 		}
	//
	// 		condition_impl *prepare(const category &c) override
	// 		{
	// 			m_single_hit = c[m_key_values];
	// 			return this;
	// 		}
	//
	// 		bool test(row_handle r) const override
	// 		{
	// 			return m_single_hit == r;
	// 		}
	//
	// 		void str(std::ostream &os) const override
	// 		{
	// 			os << "index scan";
	// 		}
	//
	// 		virtual std::optional<row_handle> single() const override
	// 		{
	// 			return m_single_hit;
	// 		}
	//
	// 		virtual bool equals(const condition_impl *rhs) const override
	// 		{
	// 			if (typeid(*rhs) == typeid(index_lookup_condition_impl))
	// 			{
	// 				auto ri = static_cast<const index_lookup_condition_impl *>(rhs);
	// 				if (m_single_hit or ri->m_single_hit)
	// 					return m_single_hit == ri->m_single_hit;
	// 				else
	// 					// watch out, both m_item_ix might be the same while item_names might be diffent (in case they both do not exist in the category)
	// 					return m_key_values == ri->m_key_values;
	// 			}
	// 			return this == rhs;
	// 		}
	//
	// 		row_initializer m_key_values;
	// 		row_handle m_single_hit;
	// 	};

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

	condition_impl *key_equals_number_condition_impl::prepare(const category &c)
	{
		m_item_ix = c.get_item_ix(m_item_name);

		if (c.get_cat_validator() != nullptr and
			c.key_item_indices().contains(m_item_ix) and
			c.key_item_indices().size() == 1)
		{
			item v(m_item_name, m_value);
			m_single_hit = c[{ { m_item_name, std::string{ v.value() }, false } }];
		}

		return this;
	}

	bool found_in_range(condition_impl *c, std::vector<and_condition_impl *>::iterator b, std::vector<and_condition_impl *>::iterator e)
	{
		bool result = true;

		for (auto s = b; s != e; ++s)
		{
			auto &cs = (*s)->m_sub;

			if (std::ranges::find_if(cs, [c](const condition_impl *i)
					{ return i->equals(c); }) == cs.end())
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

		for (size_t fc_i = 0; fc_i < fc.size();)
		{
			auto c = fc[fc_i];
			if (not found_in_range(c, subs.begin() + 1, subs.end()))
			{
				++fc_i;
				continue;
			}

			if (and_result == nullptr)
				and_result = new and_condition_impl();

			and_result->m_sub.push_back(c);
			fc.erase(fc.begin() + static_cast<std::string::difference_type>(fc_i));

			for (auto sub : subs)
			{
				auto &ssub = sub->m_sub;

				for (size_t ssub_i = 0; ssub_i < ssub.size();)
				{
					auto sc = ssub[ssub_i];
					if (not sc->equals(c))
					{
						++ssub_i;
						continue;
					}

					ssub.erase(ssub.begin() + static_cast<std::string::difference_type>(ssub_i));
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

	condition_impl *and_condition_impl::prepare(const category &c)
	{
		for (auto &sub : m_sub)
			sub = sub->prepare(c);

		if (auto cv = c.get_cat_validator(); cv != nullptr)
		{
			// See if we can collapse a search part of this and_condition into a single index lookup

			cif::iset keys{ cv->m_keys.begin(), cv->m_keys.end() };
			category::key_type lookup;
			std::vector<condition_impl *> subs;
			std::vector<std::string> may_be_empty;

			for (auto &sub : m_sub)
			{
				if (auto s = dynamic_cast<const key_equals_condition_impl *>(sub); s != nullptr)
				{
					if (keys.contains(s->m_item_name))
					{
						lookup.emplace_back(s->m_item_name, s->m_value);
						subs.emplace_back(sub);
					}
					continue;
				}

				if (auto s = dynamic_cast<const key_equals_number_condition_impl *>(sub); s != nullptr)
				{
					if (keys.contains(s->m_item_name))
					{
						item v{ s->m_item_name, s->m_value };
						lookup.emplace_back(s->m_item_name, std::string{ v.value() });
						subs.emplace_back(sub);
					}
					continue;
				}

				if (auto s = dynamic_cast<const key_equals_or_empty_condition_impl *>(sub); s != nullptr)
				{
					if (keys.contains(s->m_item_name))
					{
						lookup.emplace_back(s->m_item_name, s->m_value, true);
						subs.emplace_back(sub);
						may_be_empty.emplace_back(s->m_item_name);
					}
					continue;
				}

				if (auto s = dynamic_cast<const key_equals_number_or_empty_condition_impl *>(sub); s != nullptr)
				{
					if (keys.contains(s->m_item_name))
					{
						item v{ s->m_item_name, s->m_value };
						lookup.emplace_back(s->m_item_name, std::string{ v.value() }, true);
						subs.emplace_back(sub);
					}
					continue;
				}
			}

			if (lookup.size() == keys.size())
			{
				m_single = c[lookup];

				for (auto s : subs)
					std::erase(m_sub, s);
			}
		}

		return this;
	}

	bool and_condition_impl::test(row_handle r) const
	{
		bool result = true;

		if (m_single.has_value() and *m_single != r)
			result = false;
		else
		{
			for (auto sub : m_sub)
			{
				if (sub->test(r))
					continue;

				result = false;
				break;
			}
		}

		return result;
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
