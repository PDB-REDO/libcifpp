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

#include <cif++/v2/forward_decl.hpp>

#include <cif++/v2/category.hpp>

namespace cif::v2
{

// --------------------------------------------------------------------

class datablock
{
  public:
	using category_list = std::list<category>;

	using iterator = category_list::iterator;
	using const_iterator = category_list::const_iterator;

	using reference = typename category_list::reference;

	datablock() = default;

	datablock(std::string_view name)
		: m_name(name)
	{
	}

	datablock(const datablock &) = default;

	datablock(datablock &&) = default;

	datablock &operator=(const datablock &) = default;
	datablock &operator=(datablock &&) = default;

	// --------------------------------------------------------------------

	const std::string &name() const { return m_name; }

	void set_validator(const validator *v)
	{
		m_validator = v;

		for (auto &cat : *this)
			cat.set_validator(v, *this);
	}

	const validator *get_validator() const
	{
		return m_validator;
	}

	bool is_valid() const
	{
		if (m_validator == nullptr)
			throw std::runtime_error("Validator not specified");

		bool result = true;
		for (auto &cat : *this)
			result = cat.is_valid() and result;
		
		return result;
	}

	// --------------------------------------------------------------------

	bool empty() const { return m_categories.empty(); }
	size_t size() const { return m_categories.size(); }

	reference front() { return m_categories.front(); }
	reference back() { return m_categories.back(); }

	iterator begin() { return m_categories.begin(); }
	iterator end() { return m_categories.end(); }

	const_iterator cbegin() { return m_categories.cbegin(); }
	const_iterator cend() { return m_categories.cend(); }

	const_iterator begin() const { return m_categories.begin(); }
	const_iterator end() const { return m_categories.end(); }

	// --------------------------------------------------------------------

	category &operator[](std::string_view name)
	{
		auto i = std::find_if(m_categories.begin(), m_categories.end(), [name](const category &c)
			{ return iequals(c.name(), name); });

		if (i != m_categories.end())
			return *i;

		m_categories.emplace_back(name);
		return m_categories.back();
	}

	const category &operator[](std::string_view name) const
	{
		static const category s_empty;
		auto i = std::find_if(m_categories.begin(), m_categories.end(), [name](const category &c)
			{ return iequals(c.name(), name); });
		return i == m_categories.end() ? s_empty : *i;
	}

	category *get(std::string_view name)
	{
		auto i = std::find_if(m_categories.begin(), m_categories.end(), [name](const category &c)
			{ return iequals(c.name(), name); });
		return i == m_categories.end() ? nullptr : &*i;
	}

	const category *get(std::string_view name) const
	{
		return const_cast<datablock *>(this)->get(name);
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

	// void write(std::ostream &os) const
	// {
	// 	// std::shared_lock lock(mLock);

	// 	os << "data_" << m_name << std::endl
	// 	   << "# " << std::endl;

	// 	// mmcif support, sort of. First write the 'entry' Category
	// 	// and if it exists, _AND_ we have a Validator, write out the
	// 	// audit_conform record.

	// 	for (auto &cat : m_categories)
	// 	{
	// 		if (cat.name() != "entry")
	// 			continue;

	// 		cat.write(os);

	// 		// if (mValidator != nullptr)
	// 		// {
	// 		// 	Category auditConform(*this, "audit_conform", nullptr);
	// 		// 	auditConform.emplace({{"dict_name", mValidator->dictName()},
	// 		// 		{"dict_version", mValidator->dictVersion()}});
	// 		// 	auditConform.write(os);
	// 		// }

	// 		break;
	// 	}

	// 	for (auto &cat : m_categories)
	// 	{
	// 		if (cat.name() != "entry" and cat.name() != "audit_conform")
	// 			cat.write(os);
	// 	}
	// }

	// friend std::ostream &operator<<(std::ostream &os, const datablock &db)
	// {
	// 	db.write(os);
	// 	return os;
	// }

  private:
	category_list m_categories;
	std::string m_name;
	const validator *m_validator = nullptr;
};

} // namespace cif::v2