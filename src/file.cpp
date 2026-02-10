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

#include "cif++/file.hpp"

#include "cif++/gzio.hpp"
#include "cif++/parser.hpp"
#include "cif++/text.hpp"

#include <exception>
#include <ranges>
#include <stdexcept>
#include <string>

namespace cif
{

class validator;

bool file::is_valid() const
{
	bool result = true;
	for (auto &d : *this)
		result = d.is_valid() and result;

	if (result)
		result = validate_links();

	return result;
}

bool file::is_valid()
{
	bool result = not empty();

	for (bool first = true; auto &d : *this)
	{
		if (first)
		{
			result = d.is_valid() and result;
			first = false;
		}
		else if (d.get_validator() != nullptr)
			result = d.is_valid() and result;
	}

	if (result)
		result = validate_links();

	return result;
}

bool file::validate_links() const
{
	bool result = true;

	for (auto &db : *this)
		result = db.validate_links() and result;

	return result;
}

bool file::contains(std::string_view name) const
{
	return std::ranges::find_if(*this, [name](const datablock &db)
			   { return iequals(db.name(), name); }) != end();
}

datablock &file::operator[](std::string_view name)
{
	auto i = std::ranges::find_if(*this, [name](const datablock &c)
		{ return iequals(c.name(), name); });

	if (i != end())
		return *i;

	emplace_back(name);
	return back();
}

const datablock &file::operator[](std::string_view name) const
{
	static const datablock s_empty;
	auto i = std::ranges::find_if(*this, [name](const datablock &c)
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
			break;
		}

		++i;
	}

	if (is_new)
		i = insert(end(), { name });

	assert(i != end());
	return std::make_tuple(i, is_new);
}

void file::load(const std::filesystem::path &p)
{
	gzio::ifstream in(p);
	if (not in.is_open())
		throw std::runtime_error("Could not open file '" + p.string() + '\'');

	try
	{
		load(in);
	}
	catch (const std::exception &)
	{
		throw_with_nested(std::runtime_error("Error reading file '" + p.string() + '\''));
	}
}

void file::load(const std::filesystem::path &p, const validator &v)
{
	gzio::ifstream in(p);
	if (not in.is_open())
		throw std::runtime_error("Could not open file '" + p.string() + '\'');

	try
	{
		load(in, v);
	}
	catch (const std::exception &)
	{
		throw_with_nested(std::runtime_error("Error reading file '" + p.string() + '\''));
	}
}

void file::load(std::istream &is, const validator &v)
{
	parser p(is, *this);
	p.parse_file();
	for (auto &db : *this)
		db.set_validator(&v);
}

void file::load(std::istream &is)
{
	parser p(is, *this);
	p.parse_file();
}

void file::save(const std::filesystem::path &p) const
{
	gzio::ofstream outFile(p);
	save(outFile);
}

void file::save(std::ostream &os) const
{
	// if (not is_valid())
	// 	std::cout << "File is not valid!\n";

	for (auto &db : *this)
		db.write(os);
}

} // namespace cif
