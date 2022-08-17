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

#include <gxrio.hpp>

#include <cif++/file.hpp>

namespace cif
{

// --------------------------------------------------------------------
void file::set_validator(const validator *v)
{
	m_validator = v;
	for (auto &db : *this)
		db.set_validator(v);
}

bool file::is_valid() const
{
	if (m_validator == nullptr)
		std::runtime_error("No validator loaded explicitly, cannot continue");

	bool result = true;
	for (auto &d : *this)
		result = d.is_valid() and result;

	return result;
}

bool file::is_valid()
{
	if (m_validator == nullptr)
	{
		if (VERBOSE > 0)
			std::cerr << "No dictionary loaded explicitly, loading default" << std::endl;

		load_dictionary();
	}

	bool result = not empty();

	for (auto &d : *this)
		result = d.is_valid() and result;

	return result;
}

void file::load_dictionary()
{
	if (not empty())
	{
		auto *audit_conform = front().get("audit_conform");
		if (audit_conform and not audit_conform->empty())
		{
			std::string name;
			cif::tie(name) = audit_conform->front().get("dict_name");

			if (not name.empty())
				load_dictionary(name);
		}
	}

	if (not m_validator)
		load_dictionary("mmcif_ddl");
}

void file::load_dictionary(std::string_view name)
{
	set_validator(&validator_factory::instance()[name]);
}

datablock &file::operator[](std::string_view name)
{
	auto i = std::find_if(begin(), end(), [name](const datablock &c)
		{ return iequals(c.name(), name); });

	if (i != end())
		return *i;

	emplace_back(name);
	return back();
}

const datablock &file::operator[](std::string_view name) const
{
	static const datablock s_empty;
	auto i = std::find_if(begin(), end(), [name](const datablock &c)
		{ return iequals(c.name(), name); });
	return i == end() ? s_empty : *i;
}

std::tuple<file::iterator, bool> file::emplace(std::string_view name)
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

void file::load(const std::filesystem::path &p)
{
	gxrio::ifstream in(p);
	load(in);
}

void file::load(std::istream &is)
{
	auto saved = m_validator;
	set_validator(nullptr);

	parser p(is, *this);
	p.parse_file();

	if (saved != nullptr)
		set_validator(saved);
	else
		load_dictionary();
}

void file::save(const std::filesystem::path &p) const
{
	gxrio::ofstream outFile(p);
	save(outFile);
}

void file::save(std::ostream &os) const
{
	for (auto &db : *this)
		db.write(os);
}

} // namespace cif