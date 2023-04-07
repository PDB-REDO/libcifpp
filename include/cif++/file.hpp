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

#include <list>

#include "exports.hpp"
#include "datablock.hpp"
#include "parser.hpp"

namespace cif
{

// --------------------------------------------------------------------

class file : public std::list<datablock>
{
  public:
	file() = default;

	explicit file(const std::filesystem::path &p)
	{
		load(p);
	}

	explicit file(std::istream &is)
	{
		load(is);
	}

	explicit file(const char *data, size_t length)
	{
		struct membuf : public std::streambuf
		{
			membuf(char *text, size_t length)
			{
				this->setg(text, text, text + length);
			}
		} buffer(const_cast<char *>(data), length);

		std::istream is(&buffer);
		load(is);
	}

	file(const file &) = default;
	file(file &&) = default;
	file &operator=(const file &) = default;
	file &operator=(file &&) = default;

	void set_validator(const validator *v);

	const validator *get_validator() const
	{
		return m_validator;
	}

	bool is_valid() const;
	bool is_valid();
	bool validate_links() const;

	void load_dictionary();
	void load_dictionary(std::string_view name);

	bool contains(std::string_view name) const;

	datablock &front()
	{
		assert(not empty());
		return std::list<datablock>::front();
	}

	const datablock &front() const
	{
		assert(not empty());
		return std::list<datablock>::front();
	}

	datablock &operator[](std::string_view name);
	const datablock &operator[](std::string_view name) const;

	std::tuple<iterator, bool> emplace(std::string_view name);

	void load(const std::filesystem::path &p);
	void load(std::istream &is);

	void save(const std::filesystem::path &p) const;
	void save(std::ostream &os) const;

	friend std::ostream &operator<<(std::ostream &os, const file &f)
	{
		f.save(os);
		return os;
	}

  private:
	const validator *m_validator = nullptr;
};

} // namespace cif