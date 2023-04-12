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

#include "cif++/exports.hpp"
#include "cif++/matrix.hpp"
#include "cif++/point.hpp"

#include <array>
#include <cstdint>
#include <string>

namespace cif
{

// --------------------------------------------------------------------

enum class space_group_name
{
	full,
	xHM,
	Hall
};

struct space_group
{
	const char *name;
	const char *xHM;
	const char *Hall;
	int nr;
};

extern CIFPP_EXPORT const space_group kSpaceGroups[];
extern CIFPP_EXPORT const std::size_t kNrOfSpaceGroups;

// --------------------------------------------------------------------

struct symop_data
{
	constexpr symop_data(const std::array<int, 15> &data)
		: m_packed((data[0] bitand 0x03ULL) << 34 bitor
				   (data[1] bitand 0x03ULL) << 32 bitor
				   (data[2] bitand 0x03ULL) << 30 bitor
				   (data[3] bitand 0x03ULL) << 28 bitor
				   (data[4] bitand 0x03ULL) << 26 bitor
				   (data[5] bitand 0x03ULL) << 24 bitor
				   (data[6] bitand 0x03ULL) << 22 bitor
				   (data[7] bitand 0x03ULL) << 20 bitor
				   (data[8] bitand 0x03ULL) << 18 bitor
				   (data[9] bitand 0x07ULL) << 15 bitor
				   (data[10] bitand 0x07ULL) << 12 bitor
				   (data[11] bitand 0x07ULL) << 9 bitor
				   (data[12] bitand 0x07ULL) << 6 bitor
				   (data[13] bitand 0x07ULL) << 3 bitor
				   (data[14] bitand 0x07ULL) << 0)
	{
	}

	bool operator==(const symop_data &rhs) const
	{
		return m_packed == rhs.m_packed;
	}

	bool operator<(const symop_data &rhs) const
	{
		return m_packed < rhs.m_packed;
	}

	inline constexpr int unpack3(int offset) const
	{
		int result = (m_packed >> offset) bitand 0x03;
		return result == 3 ? -1 : result;
	}

	inline constexpr int unpack7(int offset) const
	{
		return (m_packed >> offset) bitand 0x07;
	}

	constexpr std::array<int, 15> data() const
	{
		return {
			unpack3(34),
			unpack3(32),
			unpack3(30),
			unpack3(28),
			unpack3(26),
			unpack3(24),
			unpack3(22),
			unpack3(20),
			unpack3(18),
			unpack7(15),
			unpack7(12),
			unpack7(9),
			unpack7(6),
			unpack7(3),
			unpack7(0)
		};
	}

  private:
	friend struct symop_datablock;

	const uint64_t kPackMask = (~0ULL >> (64 - 36));

	symop_data(uint64_t v)
		: m_packed(v bitand kPackMask)
	{
	}

	uint64_t m_packed;
};

struct symop_datablock
{
	constexpr symop_datablock(int spacegroup, int rotational_number, const std::array<int, 15> &rt_data)
		: m_v((spacegroup bitand 0xffffULL) << 48 bitor
			  (rotational_number bitand 0xffULL) << 40 bitor
			  symop_data(rt_data).m_packed)
	{
	}

	uint16_t spacegroup() const { return m_v >> 48; }
	symop_data symop() const { return symop_data(m_v); }
	uint8_t rotational_number() const { return (m_v >> 40) bitand 0xff; }

  private:
	uint64_t m_v;
};

static_assert(sizeof(symop_datablock) == sizeof(uint64_t), "Size of symop_data is wrong");

extern CIFPP_EXPORT const symop_datablock kSymopNrTable[];
extern CIFPP_EXPORT const std::size_t kSymopNrTableSize;

// --------------------------------------------------------------------
// Some more symmetry related stuff here.

class datablock;

class cell;
class spacegroup;
class rtop;
class sym_op;

// --------------------------------------------------------------------
// class cell

class cell
{
  public:
	cell(float a, float b, float c, float alpha = 90.f, float beta = 90.f, float gamma = 90.f);
	cell(const datablock &db);

	float get_a() const { return m_a; }
	float get_b() const { return m_b; }
	float get_c() const { return m_c; }

	float get_alpha() const { return m_alpha; }
	float get_beta() const { return m_beta; }
	float get_gamma() const { return m_gamma; }

	matrix3x3<float> get_orthogonal_matrix() const { return m_orthogonal; }
	matrix3x3<float> get_fractional_matrix() const { return m_fractional; }

  private:
	void init();

	float m_a, m_b, m_c, m_alpha, m_beta, m_gamma;
	matrix3x3<float> m_orthogonal, m_fractional;
};

/// @brief A class that encapsulates the symmetry operations as used in PDB files, i.e. a rotational number and a translation vector
class sym_op
{
  public:
	sym_op(uint8_t nr = 1, uint8_t ta = 5, uint8_t tb = 5, uint8_t tc = 5)
		: m_nr(nr)
		, m_ta(ta)
		, m_tb(tb)
		, m_tc(tc)
	{
	}

	explicit sym_op(std::string_view s);

	sym_op(const sym_op &) = default;
	sym_op(sym_op &&) = default;
	sym_op &operator=(const sym_op &) = default;
	sym_op &operator=(sym_op &&) = default;

	constexpr bool is_identity() const
	{
		return m_nr == 1 and m_ta == 5 and m_tb == 5 and m_tc == 5;
	}

	explicit constexpr operator bool() const
	{
		return not is_identity();
	}

	std::string string() const;

	friend class spacegroup;

  private:
	uint8_t m_nr;
	uint8_t m_ta, m_tb, m_tc;
};

namespace literals
{
	inline sym_op operator""_symop(const char *text, size_t length)
	{
		return sym_op({ text, length });
	}
} // namespace literals

class transformation
{
  public:
	transformation(const symop_data &data);
	transformation(const matrix3x3<float> &r, const cif::point &t)
		: m_rotation(r)
		, m_translation(t)
	{
	}

	transformation(const transformation &) = default;
	transformation(transformation &&) = default;
	transformation &operator=(const transformation &) = default;
	transformation &operator=(transformation &&) = default;

	point operator()(const cell &c, const point &pt) const;

  private:
	matrix3x3<float> m_rotation;
	point m_translation;
};

// --------------------------------------------------------------------

int get_space_group_number(std::string_view spacegroup);                        // alternative for clipper's parsing code, using space_group_name::full
int get_space_group_number(std::string_view spacegroup, space_group_name type); // alternative for clipper's parsing code

class spacegroup : public std::vector<transformation>
{
  public:
	spacegroup(std::string_view name)
		: spacegroup(get_space_group_number(name))
	{
	}

	spacegroup(std::string_view name, space_group_name type)
		: spacegroup(get_space_group_number(name, type))
	{
	}

	spacegroup(int nr);

	int get_nr() const { return m_nr; }
	std::string get_name() const;

	point operator()(const point &pt, const cell &c, sym_op symop);

  private:
	int m_nr;
	size_t m_index;
};

// --------------------------------------------------------------------

int get_space_group_number(std::string spacegroup);                        // alternative for clipper's parsing code, using space_group_name::full
int get_space_group_number(std::string spacegroup, space_group_name type); // alternative for clipper's parsing code

// class rtop
// {
//   public:
// 	rtop(const spacegroup &sg, const cell &c, int nr);

// 	friend rtop operator+(rtop rt, cif::point t);
// 	friend cif::point operator*(cif::point p, rtop rt);

//   private:
// 	cell m_c;
// 	cif::quaternion m_q;
// 	cif::point m_t;
// };

// class spacegroup
// {
//   public:
// 	spacegroup(const cif::datablock &db);

//   private:
// 	std::vector<
// };

static_assert(sizeof(sym_op) == 4, "Sym_op should be four bytes");

// --------------------------------------------------------------------
// Symmetry operations on points

template <typename T>
inline point_type<T> operator*(const matrix3x3<T> &m, const point_type<T> &pt)
{
	return point_type(m(0, 0) * pt.m_x + m(0, 1) * pt.m_y + m(0, 2) * pt.m_z,
		m(1, 0) * pt.m_x + m(1, 1) * pt.m_y + m(1, 2) * pt.m_z,
		m(2, 0) * pt.m_x + m(2, 1) * pt.m_y + m(2, 2) * pt.m_z);
}

inline point orthogonal(const point &pt, const cell &c)
{
	return c.get_orthogonal_matrix() * pt;
}

inline point fractional(const point &pt, const cell &c)
{
	return c.get_fractional_matrix() * pt;
}

} // namespace cif
