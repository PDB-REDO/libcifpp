/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
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

#include <filesystem>
#include <stdexcept>
#include <unordered_map>

#include <cif++/structure/Structure.hpp>

namespace mmcif
{

class BondMapException : public std::runtime_error
{
  public:
	BondMapException(const std::string &msg)
		: runtime_error(msg)
	{
	}
};

class BondMap
{
  public:
	BondMap(const Structure &p);

	BondMap(const BondMap &) = delete;
	BondMap &operator=(const BondMap &) = delete;

	bool operator()(const Atom &a, const Atom &b) const
	{
		return isBonded(index.at(a.id()), index.at(b.id()));
	}

	bool is1_4(const Atom &a, const Atom &b) const
	{
		uint32_t ixa = index.at(a.id());
		uint32_t ixb = index.at(b.id());

		return bond_1_4.count(key(ixa, ixb));
	}

	// links coming from the struct_conn records:
	std::vector<std::string> linked(const Atom &a) const;

	// This list of atomID's is comming from either CCD or the CCP4 dictionaries loaded
	static std::vector<std::string> atomIDsForCompound(const std::string &compoundID);

  private:
	bool isBonded(uint32_t ai, uint32_t bi) const
	{
		return bond.count(key(ai, bi)) != 0;
	}

	uint64_t key(uint32_t a, uint32_t b) const
	{
		if (a > b)
			std::swap(a, b);
		return static_cast<uint64_t>(a) | (static_cast<uint64_t>(b) << 32);
	}

	std::tuple<uint32_t, uint32_t> dekey(uint64_t k) const
	{
		return std::make_tuple(
			static_cast<uint32_t>(k >> 32),
			static_cast<uint32_t>(k));
	}

	uint32_t dim;
	std::unordered_map<std::string, uint32_t> index;
	std::set<uint64_t> bond, bond_1_4;

	std::map<std::string, std::set<std::string>> link;
};

} // namespace mmcif
