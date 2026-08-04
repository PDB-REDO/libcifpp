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
#include "cif++/cif++.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <typeinfo>
#include <vector>

namespace cif
{

iset get_category_items(const category &cat)
{
	return cat.key_items();
}

std::optional<uint16_t> get_item_ix(const category &cat, std::string_view col)
{
	auto ix = cat.get_item_ix(col);
	std::optional<uint16_t> result;
	if (ix < cat.get_item_count())
		result = ix;
	return result;
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
	// 		bool test(const_row_handle r) const override
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

	bool key_equals_condition_impl::prepare(const category &c)
	{
		bool result = false;

		if (auto ix = get_item_ix(c, m_item_name); ix.has_value())
		{
			m_item_ix = *ix;
			m_icase = is_item_type_uchar(c, m_item_name);

			if (auto cv = c.get_cat_validator();
				cv != nullptr and cv->m_keys.size() == 1 and
				cv->m_keys.front() == m_item_name)
			{
				m_single_hit = c[{ { m_item_name, m_value } }];
			}

			result = true;
		}

		return result;
	}

	bool found_in_range(condition_impl *c, std::vector<and_condition_impl *>::iterator b, std::vector<and_condition_impl *>::iterator e)
	{
		bool result = true;

		for (auto s : std::span(b, e))
		{
			auto &cs = s->m_sub;

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
		auto and_result = std::make_unique<and_condition_impl>();

		auto first = subs.front();

		for (auto &fc : first->m_sub)
		{
			if (not found_in_range(fc, subs.begin() + 1, subs.end()))
				continue;

			and_result->m_sub.push_back(fc);

			for (auto sub : std::span(subs.begin() + 1, subs.end()))
			{
				auto &ssub = sub->m_sub;

				for (auto &sc : ssub)
				{
					if (not sc->equals(fc))
						continue;

					delete sc;
					sc = nullptr;
					break;
				}

				std::erase(ssub, nullptr);
			}

			fc = nullptr;
		}

		std::erase(first->m_sub, nullptr);

		auto new_or = std::make_unique<or_condition_impl>();
		new_or->m_sub = std::move(oc->m_sub);
		and_result->m_sub.push_back(new_or.release());
		return and_result.release();
	}

	bool and_condition_impl::prepare(const category &c)
	{
		for (auto &sub : m_sub)
		{
			if (not sub->prepare(c))
				return false;
		}

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
			}

			if (lookup.size() == keys.size())
			{
				m_single = c[lookup];

				for (auto s : subs)
					std::erase(m_sub, s);
			}
		}

		return true;
	}

	bool and_condition_impl::test(const_row_handle r) const
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

	bool or_condition_impl::prepare(const category &c)
	{
		std::vector<and_condition_impl *> and_conditions;

		for (auto &sub : m_sub)
		{
			if (not sub->prepare(c))
			{
				delete sub;
				sub = nullptr;
				continue;
			}

			if (typeid(*sub) == typeid(and_condition_impl))
				and_conditions.push_back(static_cast<and_condition_impl *>(sub));
		}

		std::erase(m_sub, nullptr);

		if (not m_sub.empty() and and_conditions.size() == m_sub.size())
			m_sub = { and_condition_impl::combine_equal(and_conditions, this) };

		return not m_sub.empty();
	}

} // namespace detail

bool condition::prepare(const category &c)
{
	return m_impl != nullptr and m_impl->prepare(c);
}

} // namespace cif
