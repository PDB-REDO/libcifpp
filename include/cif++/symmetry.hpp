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

#if defined(__cpp_impl_three_way_comparison)
#include <compare>
#endif

/** \file cif++/symmetry.hpp
 *
 * This file contains code to do symmetry operations based on the
 * operations as specified in the International Tables.
 */

namespace cif
{

// --------------------------------------------------------------------

/// \brief Apply matrix transformation @a m on point @a pt and return the result
inline point operator*(const matrix3x3<float> &m, const point &pt)
{
	return {
		m(0, 0) * pt.m_x + m(0, 1) * pt.m_y + m(0, 2) * pt.m_z,
		m(1, 0) * pt.m_x + m(1, 1) * pt.m_y + m(1, 2) * pt.m_z,
		m(2, 0) * pt.m_x + m(2, 1) * pt.m_y + m(2, 2) * pt.m_z
	};
}

// --------------------------------------------------------------------

/// \brief the space groups we know
enum class space_group_name
{
	full, ///< The *full* spacegroup
	xHM,  ///< The *xHM* spacegroup
	Hall  ///< The *Hall* spacegroup
};

/// \brief For each known spacegroup we define a structure like this
struct space_group
{
	const char *name; ///< The name according to *full*
	const char *xHM;  ///< The name according to *xHM*
	const char *Hall; ///< The name according to *Hall*
	int nr;           ///< The number for this spacegroup
};

/// \brief Global list of spacegroups
extern CIFPP_EXPORT const space_group kSpaceGroups[];

/// \brief Global for the size of the list of spacegroups
extern CIFPP_EXPORT const std::size_t kNrOfSpaceGroups;

// --------------------------------------------------------------------

/**
 * @brief Helper class to efficiently pack the data that
 * makes up a symmetry operation
 *
 */

struct symop_data
{
	/// \brief constructor
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

	/// \brief compare
	bool operator==(const symop_data &rhs) const
	{
		return m_packed == rhs.m_packed;
	}

	/// \brief sorting order
	bool operator<(const symop_data &rhs) const
	{
		return m_packed < rhs.m_packed;
	}

	/// \brief return an int representing the value stored in the two bits at offset @a offset
	inline constexpr int unpack3(int offset) const
	{
		int result = (m_packed >> offset) bitand 0x03;
		return result == 3 ? -1 : result;
	}

	/// \brief return an int representing the value stored in the three bits at offset @a offset
	inline constexpr int unpack7(int offset) const
	{
		return (m_packed >> offset) bitand 0x07;
	}

	/// \brief return an array of 15 ints representing the values stored
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

/**
 * @brief For each symmetry operator defined in the international tables
 * we have an entry in this struct type. It contains the spacegroup
 * number, the symmetry operations and the rotational number.
 */
struct symop_datablock
{
	/// \brief constructor
	constexpr symop_datablock(int spacegroup, int rotational_number, const std::array<int, 15> &rt_data)
		: m_v((spacegroup bitand 0xffffULL) << 48 bitor
			  (rotational_number bitand 0xffULL) << 40 bitor
			  symop_data(rt_data).m_packed)
	{
	}

	uint16_t spacegroup() const { return m_v >> 48; }                     ///< Return the spacegroup
	symop_data symop() const { return symop_data(m_v); }                  ///< Return the symmetry operation
	uint8_t rotational_number() const { return (m_v >> 40) bitand 0xff; } ///< Return the rotational_number

  private:
	uint64_t m_v;
};

static_assert(sizeof(symop_datablock) == sizeof(uint64_t), "Size of symop_data is wrong");

/// \brief Global containing the list of known symmetry operations
extern CIFPP_EXPORT const symop_datablock kSymopNrTable[];

/// \brief Size of the list of known symmetry operations
extern CIFPP_EXPORT const std::size_t kSymopNrTableSize;

// --------------------------------------------------------------------
// Some more symmetry related stuff here.

class datablock;

class cell;
class spacegroup;
class rtop;
struct sym_op;

/** @brief A class that encapsulates the symmetry operations as used in PDB files,
 * i.e. a rotational number and a translation vector.
 *
 * The syntax in string format follows the syntax as used in mmCIF files, i.e.
 * rotational number followed by underscore and the three translations where 5 is
 * no movement.
 *
 * So the string 1_555 means no symmetry movement at all since the rotational number
 * 1 always corresponds to the symmetry operation [x, y, z].
 */

struct sym_op
{
  public:
	/// \brief constructor
	sym_op(uint8_t nr = 1, uint8_t ta = 5, uint8_t tb = 5, uint8_t tc = 5)
		: m_nr(nr)
		, m_ta(ta)
		, m_tb(tb)
		, m_tc(tc)
	{
	}

	/// \brief construct a sym_op based on the contents encoded in string @a s
	explicit sym_op(std::string_view s);

	/** @cond */
	sym_op(const sym_op &) = default;
	sym_op(sym_op &&) = default;
	sym_op &operator=(const sym_op &) = default;
	sym_op &operator=(sym_op &&) = default;
	/** @endcond */

	/// \brief return true if this sym_op is the identity operator
	constexpr bool is_identity() const
	{
		return m_nr == 1 and m_ta == 5 and m_tb == 5 and m_tc == 5;
	}

	/// \brief quick test for unequal to identity
	explicit constexpr operator bool() const
	{
		return not is_identity();
	}

	/// \brief return the content encoded in a string
	std::string string() const;

#if defined(__cpp_impl_three_way_comparison)
	/// \brief a default spaceship operator
	constexpr auto operator<=>(const sym_op &rhs) const = default;
#else
	/// \brief a default equals operator
	constexpr bool operator==(const sym_op &rhs) const
	{
		return m_nr == rhs.m_nr and m_ta == rhs.m_ta and m_tb == rhs.m_tb and m_tc == rhs.m_tc;
	}

	/// \brief a default not-equals operator
	constexpr bool operator!=(const sym_op &rhs) const
	{
		return not operator==(rhs);
	}
#endif

	/// @cond
	uint8_t m_nr;
	uint8_t m_ta, m_tb, m_tc;
	/// @endcond
};

static_assert(sizeof(sym_op) == 4, "Sym_op should be four bytes");

namespace literals
{
	/**
	 * @brief This operator allows you to write code like this:
	 *
	 * @code {.cpp}
	 * using namespace cif::literals;
	 *
	 * cif::sym_op so = "1_555"_symop;
	 * @endcode
	 *
	 */
	inline sym_op operator""_symop(const char *text, std::size_t length)
	{
		return sym_op({ text, length });
	}
} // namespace literals

// --------------------------------------------------------------------
// The transformation class

/**
 * @brief A class you can use to apply symmetry transformations on points
 *
 * Transformations consist of two operations, a matrix transformation which
 * is often a rotation followed by a translation.
 *
 * In case the matrix transformation is a pure rotation a quaternion
 * is created to do the actual calculations. That's faster and more
 * precise.
 */
class transformation
{
  public:
	/// \brief constructor taking a symop_data object @a data
	transformation(const symop_data &data);

	/// \brief constructor taking a rotation matrix @a r and a translation vector @a t
	transformation(const matrix3x3<float> &r, const cif::point &t);

	/** @cond */
	transformation(const transformation &) = default;
	transformation(transformation &&) = default;
	transformation &operator=(const transformation &) = default;
	transformation &operator=(transformation &&) = default;
	/** @endcond */

	/// \brief operator() to perform the transformation on point @a pt and return the result
	point operator()(point pt) const
	{
		if (m_q)
			pt.rotate(m_q);
		else
			pt = m_rotation * pt;

		return pt + m_translation;
	}

	/// \brief return a transformation object that is the result of applying @a rhs after @a lhs
	friend transformation operator*(const transformation &lhs, const transformation &rhs);

	/// \brief return the inverse transformation for @a t
	friend transformation inverse(const transformation &t);

	/// \brief return the inverse tranformation for this
	transformation operator-() const
	{
		return inverse(*this);
	}

	friend class spacegroup;

  private:
	// Most rotation matrices provided by the International Tables
	// are really rotation matrices, in those cases we can construct
	// a quaternion. Unfortunately, that doesn't work for all of them

	void try_create_quaternion();

	matrix3x3<float> m_rotation;
	quaternion m_q;
	point m_translation;
};

// --------------------------------------------------------------------
// class cell

/**
 * @brief The cell class describes the dimensions and angles of a unit cell
 * in a crystal
 */

class cell
{
  public:
	/// \brief constructor
	cell(float a, float b, float c, float alpha = 90.f, float beta = 90.f, float gamma = 90.f);

	/// \brief constructor that takes the appropriate values from the *cell* category in datablock @a db
	cell(const datablock &db);

	float get_a() const { return m_a; } ///< return dimension a
	float get_b() const { return m_b; } ///< return dimension b
	float get_c() const { return m_c; } ///< return dimension c

	float get_alpha() const { return m_alpha; } ///< return angle alpha
	float get_beta() const { return m_beta; }   ///< return angle beta
	float get_gamma() const { return m_gamma; } ///< return angle gamma

	float get_volume() const; ///< return the calculated volume for this cell

	matrix3x3<float> get_orthogonal_matrix() const { return m_orthogonal; } ///< return the matrix to use to transform coordinates from fractional to orthogonal
	matrix3x3<float> get_fractional_matrix() const { return m_fractional; } ///< return the matrix to use to transform coordinates from orthogonal to fractional

  private:
	void init();

	float m_a, m_b, m_c, m_alpha, m_beta, m_gamma;
	matrix3x3<float> m_orthogonal, m_fractional;
};

// --------------------------------------------------------------------

/// \brief Return the spacegroup number from the *symmetry* category in datablock @a db
int get_space_group_number(const datablock &db);

/// \brief Return the spacegroup number for spacegroup named @a spacegroup
int get_space_group_number(std::string_view spacegroup);

/// \brief Return the spacegroup number for spacegroup named @a spacegroup assuming space_group_name @a type
int get_space_group_number(std::string_view spacegroup, space_group_name type);

/**
 * @brief class to encapsulate the list of transformations making up a spacegroup
 *
 */
class spacegroup : public std::vector<transformation>
{
  public:
	/// \brief constructor using the information in the *symmetry* category in datablock @a db
	spacegroup(const datablock &db)
		: spacegroup(get_space_group_number(db))
	{
	}

	/// \brief constructor using the spacegroup named @a name
	spacegroup(std::string_view name)
		: spacegroup(get_space_group_number(name))
	{
	}

	/// \brief constructor using the spacegroup named @a name assuming space_group_name @a type
	spacegroup(std::string_view name, space_group_name type)
		: spacegroup(get_space_group_number(name, type))
	{
	}

	/// \brief constructor using the spacegroup number @a nr
	spacegroup(int nr);

	int get_nr() const { return m_nr; } ///< Return the nr
	std::string get_name() const;       ///< Return the name

	/** \brief perform a spacegroup operation on point @a pt using
	 * cell @a c and sym_op @a symop
	 */

	point operator()(const point &pt, const cell &c, sym_op symop) const;

	/** \brief perform an inverse spacegroup operation on point @a pt using
	 * cell @a c and sym_op @a symop
	 */
	point inverse(const point &pt, const cell &c, sym_op symop) const;

  private:
	int m_nr;
	std::size_t m_index;
};

// --------------------------------------------------------------------
/**
 * @brief A crystal combines a cell and a spacegroup.
 *
 * The information in cell and spacegroup together make up all
 * information you need to do symmetry calculations in a crystal
 */

class crystal
{
  public:
	/// \brief constructor using the information found in datablock @a db
	crystal(const datablock &db)
		: m_cell(db)
		, m_spacegroup(db)
	{
	}

	/// \brief constructor using cell @a c and spacegroup @a sg
	crystal(const cell &c, const spacegroup &sg)
		: m_cell(c)
		, m_spacegroup(sg)
	{
	}

	/** @cond */
	crystal(const crystal &) = default;
	crystal(crystal &&) = default;
	crystal &operator=(const crystal &) = default;
	crystal &operator=(crystal &&) = default;
	/** @endcond */

	const cell &get_cell() const { return m_cell; }                   ///< Return the cell
	const spacegroup &get_spacegroup() const { return m_spacegroup; } ///< Return the spacegroup

	/// \brief Return the symmetry copy of point @a pt using symmetry operation @a symop
	point symmetry_copy(const point &pt, sym_op symop) const
	{
		return m_spacegroup(pt, m_cell, symop);
	}

	/// \brief Return the symmetry copy of point @a pt using the inverse of symmetry operation @a symop
	point inverse_symmetry_copy(const point &pt, sym_op symop) const
	{
		return m_spacegroup.inverse(pt, m_cell, symop);
	}

	/// \brief Return a tuple consisting of distance, new location and symmetry operation
	/// for the point @a b with respect to point @a a.
	std::tuple<float, point, sym_op> closest_symmetry_copy(point a, point b) const;

  private:
	cell m_cell;
	spacegroup m_spacegroup;
};

// --------------------------------------------------------------------
// Symmetry operations on points

/// \brief convenience function returning the fractional point @a pt in orthogonal coordinates for cell @a c
inline point orthogonal(const point &pt, const cell &c)
{
	return c.get_orthogonal_matrix() * pt;
}

/// \brief convenience function returning the orthogonal point @a pt in fractional coordinates for cell @a c
inline point fractional(const point &pt, const cell &c)
{
	return c.get_fractional_matrix() * pt;
}

// --------------------------------------------------------------------

} // namespace cif
