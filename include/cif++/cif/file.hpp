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

#include <cif++/cif/forward_decl.hpp>

#include <cif++/cif/datablock.hpp>
#include <cif++/cif/parser.hpp>

namespace cif::v2
{

// --------------------------------------------------------------------

class file : public std::list<datablock>
{
  public:

	file() = default;

	file(std::istream &is)
	{
		load(is);
	}

	file(const file &) = default;
	file(file &&) = default;
	file &operator=(const file &) = default;
	file &operator=(file &&) = default;

	void set_validator(const validator *v)
	{
		m_validator = v;
		for (auto &db : *this)
			db.set_validator(v);
	}

	const validator *get_validator() const
	{
		return m_validator;
	}

	bool is_valid() const
	{
		if (m_validator == nullptr)
			std::runtime_error("No validator loaded explicitly, cannot continue");

		bool result = true;
		for (auto &d : *this)
			result = d.is_valid() and result;

		return result;
	}
	
	bool is_valid()
	{
		if (m_validator == nullptr)
		{
			if (VERBOSE > 0)
				std::cerr << "No dictionary loaded explicitly, loading default" << std::endl;

			load_dictionary();
		}

		bool result = true;
		for (auto &d : *this)
			result = d.is_valid() and result;

		return result;
	}
	
	void load_dictionary()
	{
		load_dictionary("mmcif_ddl");
	}

	void load_dictionary(std::string_view name)
	{
		set_validator(&validator_factory::instance()[name]);
	}

	datablock &operator[](std::string_view name)
	{
		auto i = std::find_if(begin(), end(), [name](const datablock &c)
			{ return iequals(c.name(), name); });
		
		if (i != end())
			return *i;

		emplace_back(name);
		return back();
	}

	const datablock &operator[](std::string_view name) const
	{
		static const datablock s_empty;
		auto i = std::find_if(begin(), end(), [name](const datablock &c)
			{ return iequals(c.name(), name); });
		return i == end() ? s_empty : *i;
	}

	std::tuple<iterator, bool> emplace(std::string_view name)
	{
		bool is_new = true;

		auto i = begin();
		while (i != end())
		{
			if (iequals(name, i->name()))
			{
				is_new = false;

				if (i != begin())
				{
					auto n = std::next(i);
					splice(begin(), *this, i, n);
				}

				break;
			}

			++i;
		}

		if (is_new)
		{
			auto &db = emplace_front(name);
			db.set_validator(m_validator);
		}

		return std::make_tuple(begin(), is_new);		
	}

	void load(std::istream &is)
	{
		auto saved = m_validator;
		set_validator(nullptr);

		parser p(is, *this);
		p.parse_file();

		if (saved != nullptr)
		{
			set_validator(saved);
			(void)is_valid();
		}
	}

  private:
	const validator* m_validator = nullptr;
};

}