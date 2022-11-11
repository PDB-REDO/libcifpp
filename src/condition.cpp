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

#include <cif++/category.hpp>
#include <cif++/condition.hpp>

namespace cif
{

iset get_category_fields(const category &cat)
{
	return cat.key_fields();
}

uint16_t get_column_ix(const category &cat, std::string_view col)
{
	return cat.get_column_ix(col);
}

bool is_column_type_uchar(const category &cat, std::string_view col)
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
		m_item_ix = get_column_ix(c, m_item_tag);
		m_icase = is_column_type_uchar(c, m_item_tag);

		if (c.get_cat_validator() != nullptr and
			c.key_field_indices().contains(m_item_ix) and
			c.key_field_indices().size() == 1)
		{
			m_single_hit = c[{ { m_item_tag, m_value } }];
		}

		return this;
	}

	condition_impl *and_condition_impl::prepare(const category &c)
	{
		for (auto &sub : mSub)
			sub = sub->prepare(c);

		for (;;)
		{
			auto si = find_if(mSub.begin(), mSub.end(), [](condition_impl *sub) { return dynamic_cast<and_condition_impl *>(sub) != nullptr; });
			if (si == mSub.end())
				break;
			
			and_condition_impl *sub_and = static_cast<and_condition_impl *>(*si);

			mSub.erase(si);

			mSub.insert(mSub.end(), sub_and->mSub.begin(), sub_and->mSub.end());
			sub_and->mSub.clear();
			delete sub_and;
		}

		return this;
	}

	condition_impl *or_condition_impl::prepare(const category &c)
	{
		condition_impl *result = this;

		mA = mA->prepare(c);
		mB = mB->prepare(c);

		key_equals_condition_impl *equals = dynamic_cast<key_equals_condition_impl*>(mA);
		key_is_empty_condition_impl *empty = dynamic_cast<key_is_empty_condition_impl*>(mB);

		if (equals == nullptr and empty == nullptr)
		{
			equals = dynamic_cast<key_equals_condition_impl*>(mB);
			empty = dynamic_cast<key_is_empty_condition_impl*>(mA);			
		}

		if (equals != nullptr and empty != nullptr and equals->m_item_tag == empty->m_item_tag)
		{
			result = new detail::key_equals_or_empty_condition_impl(equals);
			result = result->prepare(c);
			delete this;
		}

		return result;
	}

} // namespace detail

void condition::prepare(const category &c)
{
	if (m_impl)
		m_impl = m_impl->prepare(c);
	
	m_prepared = true;
}

} // namespace cif
