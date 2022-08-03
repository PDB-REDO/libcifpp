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

#include <cif++/v2/datablock.hpp>
#include <cif++/v2/parser.hpp>

namespace cif::v2
{

// --------------------------------------------------------------------

class file
{
  public:
	using datablock_list = std::list<datablock>;

	using value_type = datablock_list::value_type;
	using reference = datablock_list::reference;
	using pointer = datablock_list::pointer;

	using iterator = datablock_list::iterator;
	using const_iterator = datablock_list::const_iterator;

	file() = default;

	file(std::istream &is)
	{
		load(is);
	}

	file(const file &) = default;
	file(file &&) = default;
	file &operator=(const file &) = default;
	file &operator=(file &&) = default;

	datablock &operator[](std::string_view name)
	{
		auto i = std::find_if(m_datablocks.begin(), m_datablocks.end(), [name](const datablock &c)
			{ return iequals(c.name(), name); });
		
		if (i != m_datablocks.end())
			return *i;

		m_datablocks.emplace_back(name);
		return m_datablocks.back();
	}

	const datablock &operator[](std::string_view name) const
	{
		static const datablock s_empty;
		auto i = std::find_if(m_datablocks.begin(), m_datablocks.end(), [name](const datablock &c)
			{ return iequals(c.name(), name); });
		return i == m_datablocks.end() ? s_empty : *i;
	}

	std::tuple<iterator, bool> emplace(std::string_view name)
	{
		bool is_new = true;

		auto i = m_datablocks.begin();
		while (i != m_datablocks.end())
		{
			if (iequals(name, i->name()))
			{
				is_new = false;

				if (i != m_datablocks.begin())
				{
					auto n = std::next(i);
					m_datablocks.splice(m_datablocks.begin(), m_datablocks, i, n);
				}

				break;
			}

			++i;
		}

		if (is_new)
			m_datablocks.emplace(m_datablocks.begin(), name);

		return std::make_tuple(m_datablocks.begin(), is_new);		
	}

	bool empty() const { return m_datablocks.empty(); }
	size_t size() const { return m_datablocks.size(); }

	reference front() { return m_datablocks.front(); }
	reference back() { return m_datablocks.back(); }

	void load(std::istream &is)
	{
		// auto saved = mValidator;
		// setValidator(nullptr);

		parser p(is, *this);
		p.parse_file();

		// if (saved != nullptr)
		// {
		// 	setValidator(saved);
		// 	(void)isValid();
		// }
	}

  private:
	datablock_list m_datablocks;
	std::unique_ptr<Validator> m_validator;
};

}