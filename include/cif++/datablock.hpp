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

#include <cif++/forward_decl.hpp>

#include <cif++/category.hpp>

namespace cif
{

// --------------------------------------------------------------------

class datablock : public std::list<category>
{
  public:
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

	void set_name(std::string_view name)
	{
		m_name = name;
	}

	void set_validator(const validator *v);

	const validator *get_validator() const;

	bool is_valid() const;
	void validate_links() const;

	// --------------------------------------------------------------------

	category &operator[](std::string_view name);
	const category &operator[](std::string_view name) const;

	category *get(std::string_view name);
	const category *get(std::string_view name) const;

	std::tuple<iterator, bool> emplace(std::string_view name);

	std::vector<std::string> get_tag_order() const;
	void write(std::ostream &os) const;
	void write(std::ostream &os, const std::vector<std::string> &tag_order);

	friend std::ostream &operator<<(std::ostream &os, const datablock &db)
	{
		db.write(os);
		return os;
	}

	// --------------------------------------------------------------------
	
	bool operator==(const datablock &rhs) const;

  private:
	std::string m_name;
	const validator *m_validator = nullptr;
};

} // namespace cif