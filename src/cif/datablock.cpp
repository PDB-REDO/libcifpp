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

#include <cif++/cif/datablock.hpp>

namespace cif::v2
{

void datablock::set_validator(const validator *v)
{
	m_validator = v;

	for (auto &cat : *this)
		cat.set_validator(v, *this);
}

const validator *datablock::get_validator() const
{
	return m_validator;
}

bool datablock::is_valid() const
{
	if (m_validator == nullptr)
		throw std::runtime_error("Validator not specified");

	bool result = true;
	for (auto &cat : *this)
		result = cat.is_valid() and result;

	return result;
}

// --------------------------------------------------------------------

category &datablock::operator[](std::string_view name)
{
	auto i = std::find_if(begin(), end(), [name](const category &c)
		{ return iequals(c.name(), name); });

	if (i != end())
		return *i;

	emplace_back(name);
	return back();
}

const category &datablock::operator[](std::string_view name) const
{
	static const category s_empty;
	auto i = std::find_if(begin(), end(), [name](const category &c)
		{ return iequals(c.name(), name); });
	return i == end() ? s_empty : *i;
}

category *datablock::get(std::string_view name)
{
	auto i = std::find_if(begin(), end(), [name](const category &c)
		{ return iequals(c.name(), name); });
	return i == end() ? nullptr : &*i;
}

const category *datablock::get(std::string_view name) const
{
	return const_cast<datablock *>(this)->get(name);
}

std::tuple<datablock::iterator, bool> datablock::emplace(std::string_view name)
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
		auto &c = emplace_front(name);
		c.set_validator(m_validator, *this);
	}

	return std::make_tuple(begin(), is_new);
}

std::vector<std::string> datablock::get_tag_order() const
{
	std::vector<std::string> result;

	for (auto &cat : *this)
	{
		auto cto = cat.get_tag_order();
		result.insert(result.end(), cto.begin(), cto.end());
	}

	return result;
}

void datablock::write(std::ostream &os) const
{
	os << "data_" << m_name << std::endl
	   << "# " << std::endl;

	// mmcif support, sort of. First write the 'entry' Category
	// and if it exists, _AND_ we have a Validator, write out the
	// audit_conform record.

	for (auto &cat : *this)
	{
		if (cat.name() != "entry")
			continue;

		cat.write(os);

		if (m_validator != nullptr)
		{
			category auditConform("audit_conform");
			auditConform.emplace({
				{"dict_name", m_validator->name()},
				{"dict_version", m_validator->version()}});
			auditConform.write(os);
		}

		break;
	}

	for (auto &cat : *this)
	{
		if (cat.name() != "entry" and cat.name() != "audit_conform")
			cat.write(os);
	}
}

} // namespace cif::cif