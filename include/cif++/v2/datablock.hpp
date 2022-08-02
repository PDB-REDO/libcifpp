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

#include <cif++/v2/category.hpp>

namespace cif::v2
{

// --------------------------------------------------------------------

template <
	typename Alloc = std::allocator<void>,
	typename Category = category_t<Alloc>>
class datablock_t
{
  public:
	using category_type = Category;
	using allocator_type = Alloc;

	using category_allocator_type = typename std::allocator_traits<Alloc>::template rebind_alloc<category_type>;
	using category_type_list = std::list<category_type, category_allocator_type>;

	using iterator = category_type_list::iterator;
	using const_iterator = category_type_list::const_iterator;

	datablock_t(std::string_view name, const allocator_type &alloc = allocator_type())
		: m_categories(alloc)
		, m_name(name)
	{
	}

	datablock_t(const datablock_t &) = default;

	datablock_t(datablock_t &&) = default;

	// template <typename Alloc2>
	// datablock_t(const datablock_t &db, const Alloc2 &a)
	// 	: m_categories(db, a)
	// 	, m_name(db.m_name)
	// {
	// }

	// template <typename Alloc2>
	// datablock_t(datablock_t &&db, const Alloc2 &a)
	// 	: base_type(std::move(db), a)
	// 	, m_name(db.m_name)
	// {
	// }

	datablock_t &operator=(const datablock_t &) = default;
	datablock_t &operator=(datablock_t &&) = default;

	// --------------------------------------------------------------------

	const std::string &name() const { return m_name; }

	// --------------------------------------------------------------------

	category_type &operator[](std::string_view name)
	{
		auto i = std::find_if(m_categories.begin(), m_categories.end(), [name](const category_type &c)
			{ return iequals(c.name(), name); });
		
		if (i != m_categories.end())
			return *i;

		m_categories.emplace_back(name);
		return m_categories.back();
	}

	const category_type &operator[](std::string_view name) const
	{
		static const category_type s_empty;
		auto i = std::find_if(m_categories.begin(), m_categories.end(), [name](const category_type &c)
			{ return iequals(c.name(), name); });
		return i == m_categories.end() ? s_empty : *i;
	}

	std::tuple<iterator, bool> emplace(std::string_view name)
	{
		bool is_new = true;

		auto i = m_categories.begin();
		while (i != m_categories.end())
		{
			if (iequals(name, i->name()))
			{
				is_new = false;

				if (i != m_categories.begin())
				{
					auto n = std::next(i);
					m_categories.splice(m_categories.begin(), m_categories, i, n);
				}

				break;
			}

			++i;
		}

		if (is_new)
		{
			m_categories.emplace(m_categories.begin(), name);
			// m_categories.emplace(begin(), *this, std::string(name), mValidator);

			// for (auto &cat : mCategories)
			// 	cat.updateLinks();
		}

		return std::make_tuple(m_categories.begin(), is_new);		
	}

	void write(std::ostream &os) const
	{
		// std::shared_lock lock(mLock);

		os << "data_" << m_name << std::endl
		   << "# " << std::endl;

		// mmcif support, sort of. First write the 'entry' Category
		// and if it exists, _AND_ we have a Validator, write out the
		// audit_conform record.

		for (auto &cat : m_categories)
		{
			if (cat.name() != "entry")
				continue;

			cat.write(os);

			// if (mValidator != nullptr)
			// {
			// 	Category auditConform(*this, "audit_conform", nullptr);
			// 	auditConform.emplace({{"dict_name", mValidator->dictName()},
			// 		{"dict_version", mValidator->dictVersion()}});
			// 	auditConform.write(os);
			// }

			break;
		}

		for (auto &cat : m_categories)
		{
			if (cat.name() != "entry" and cat.name() != "audit_conform")
				cat.write(os);
		}
	}

	friend std::ostream &operator<<(std::ostream &os, const datablock_t &db)
	{
		db.write(os);
		return os;
	}

  private:
	category_type_list m_categories;
	std::string m_name;
};

using datablock = datablock_t<>;

}