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

#include <cstdint>
#include <array>

namespace mmcif
{

// --------------------------------------------------------------------

struct Spacegroup
{
	const char* name;
	const char* xHM;
	const char* Hall;
	int nr;
};

extern const Spacegroup kSpaceGroups[];
extern const std::size_t kNrOfSpaceGroups;

// --------------------------------------------------------------------

struct SymopData
{
	constexpr SymopData(const std::array<int,15>& data)
		: m_packed((static_cast<uint64_t>(data[ 0]) & 0x03) << 34 bitor
				   (static_cast<uint64_t>(data[ 1]) & 0x03) << 32 bitor
				   (static_cast<uint64_t>(data[ 2]) & 0x03) << 30 bitor
				   (static_cast<uint64_t>(data[ 3]) & 0x03) << 28 bitor
				   (static_cast<uint64_t>(data[ 4]) & 0x03) << 26 bitor
				   (static_cast<uint64_t>(data[ 5]) & 0x03) << 24 bitor
				   (static_cast<uint64_t>(data[ 6]) & 0x03) << 22 bitor
				   (static_cast<uint64_t>(data[ 7]) & 0x03) << 20 bitor
				   (static_cast<uint64_t>(data[ 8]) & 0x03) << 18 bitor
				   (static_cast<uint64_t>(data[ 9]) & 0x07) << 15 bitor
				   (static_cast<uint64_t>(data[10]) & 0x07) << 12 bitor
				   (static_cast<uint64_t>(data[11]) & 0x07) <<  9 bitor
				   (static_cast<uint64_t>(data[12]) & 0x07) <<  6 bitor
				   (static_cast<uint64_t>(data[13]) & 0x07) <<  3 bitor
				   (static_cast<uint64_t>(data[14]) & 0x07) <<  0)
	{
	}

	bool operator==(const SymopData& rhs) const
	{
		return m_packed == rhs.m_packed;
	}

  private:

	friend struct SymopDataBlock;

	const uint64_t kPackMask = (~0UL >> (64-36));

	SymopData(uint64_t v)
		: m_packed(v bitand kPackMask) {}

	uint64_t m_packed;
};

struct SymopDataBlock
{
	constexpr SymopDataBlock(int spacegroup, int rotational_number, const std::array<int,15>& rt_data)
		: m_v((spacegroup & 0xffffULL) << 48 bitor
			  (rotational_number & 0xffULL) << 40 bitor
			  SymopData(rt_data).m_packed)
	{
	}

	uint16_t spacegroup() const			{ return m_v >> 48; }
	SymopData symop() const				{ return SymopData(m_v); }
	uint8_t rotational_number() const	{ return (m_v >> 40) bitand 0xff; }

  private:
	uint64_t m_v;
};

static_assert(sizeof(SymopDataBlock) == sizeof(uint64_t), "Size of SymopData is wrong");

extern const SymopDataBlock kSymopNrTable[];
extern const std::size_t kSymopNrTableSize;

// --------------------------------------------------------------------

int GetSpacegroupNumber(std::string spacegroup);	// alternative for clipper's parsing code

}
