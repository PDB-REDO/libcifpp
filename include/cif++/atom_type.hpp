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

/** \file atom_type.hpp
 * 
 * This file contains information about all known elements
 */

#pragma once

#include "cif++/exports.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

namespace cif
{

/** Atom type as an integer. All known elements are available as a constant. */

enum atom_type : uint8_t
{
	Nn = 0, ///< Unknown

	H = 1,  ///< Hydro­gen
	He = 2, ///< He­lium

	Li = 3,  ///< Lith­ium
	Be = 4,  ///< Beryl­lium
	B = 5,   ///< Boron
	C = 6,   ///< Carbon
	N = 7,   ///< Nitro­gen
	O = 8,   ///< Oxy­gen
	F = 9,   ///< Fluor­ine
	Ne = 10, ///< Neon

	Na = 11, ///< So­dium
	Mg = 12, ///< Magne­sium
	Al = 13, ///< Alumin­ium
	Si = 14, ///< Sili­con
	P = 15,  ///< Phos­phorus
	S = 16,  ///< Sulfur
	Cl = 17, ///< Chlor­ine
	Ar = 18, ///< Argon

	K = 19,  ///< Potas­sium
	Ca = 20, ///< Cal­cium
	Sc = 21, ///< Scan­dium
	Ti = 22, ///< Tita­nium
	V = 23,  ///< Vana­dium
	Cr = 24, ///< Chrom­ium
	Mn = 25, ///< Manga­nese
	Fe = 26, ///< Iron
	Co = 27, ///< Cobalt
	Ni = 28, ///< Nickel
	Cu = 29, ///< Copper
	Zn = 30, ///< Zinc
	Ga = 31, ///< Gallium
	Ge = 32, ///< Germa­nium
	As = 33, ///< Arsenic
	Se = 34, ///< Sele­nium
	Br = 35, ///< Bromine
	Kr = 36, ///< Kryp­ton

	Rb = 37, ///< Rubid­ium
	Sr = 38, ///< Stront­ium
	Y = 39,  ///< Yttrium
	Zr = 40, ///< Zirco­nium
	Nb = 41, ///< Nio­bium
	Mo = 42, ///< Molyb­denum
	Tc = 43, ///< Tech­netium
	Ru = 44, ///< Ruthe­nium
	Rh = 45, ///< Rho­dium
	Pd = 46, ///< Pallad­ium
	Ag = 47, ///< Silver
	Cd = 48, ///< Cad­mium
	In = 49, ///< Indium
	Sn = 50, ///< Tin
	Sb = 51, ///< Anti­mony
	Te = 52, ///< Tellurium
	I = 53,  ///< Iodine
	Xe = 54, ///< Xenon
	Cs = 55, ///< Cae­sium
	Ba = 56, ///< Ba­rium
	La = 57, ///< Lan­thanum

	Hf = 72, ///< Haf­nium
	Ta = 73, ///< Tanta­lum
	W = 74,  ///< Tung­sten
	Re = 75, ///< Rhe­nium
	Os = 76, ///< Os­mium
	Ir = 77, ///< Iridium
	Pt = 78, ///< Plat­inum
	Au = 79, ///< Gold
	Hg = 80, ///< Mer­cury
	Tl = 81, ///< Thallium
	Pb = 82, ///< Lead
	Bi = 83, ///< Bis­muth
	Po = 84, ///< Polo­nium
	At = 85, ///< Asta­tine
	Rn = 86, ///< Radon
	Fr = 87, ///< Fran­cium
	Ra = 88, ///< Ra­dium
	Ac = 89, ///< Actin­ium

	Rf = 104, ///< Ruther­fordium
	Db = 105, ///< Dub­nium
	Sg = 106, ///< Sea­borgium
	Bh = 107, ///< Bohr­ium
	Hs = 108, ///< Has­sium
	Mt = 109, ///< Meit­nerium
	Ds = 110, ///< Darm­stadtium
	Rg = 111, ///< Roent­genium
	Cn = 112, ///< Coper­nicium
	Nh = 113, ///< Nihon­ium
	Fl = 114, ///< Flerov­ium
	Mc = 115, ///< Moscov­ium
	Lv = 116, ///< Liver­morium
	Ts = 117, ///< Tenness­ine
	Og = 118, ///< Oga­nesson

	Ce = 58, ///< Cerium
	Pr = 59, ///< Praseo­dymium
	Nd = 60, ///< Neo­dymium
	Pm = 61, ///< Prome­thium
	Sm = 62, ///< Sama­rium
	Eu = 63, ///< Europ­ium
	Gd = 64, ///< Gadolin­ium
	Tb = 65, ///< Ter­bium
	Dy = 66, ///< Dyspro­sium
	Ho = 67, ///< Hol­mium
	Er = 68, ///< Erbium
	Tm = 69, ///< Thulium
	Yb = 70, ///< Ytter­bium
	Lu = 71, ///< Lute­tium

	Th = 90,  ///< Thor­ium
	Pa = 91,  ///< Protac­tinium
	U = 92,   ///< Ura­nium
	Np = 93,  ///< Neptu­nium
	Pu = 94,  ///< Pluto­nium
	Am = 95,  ///< Ameri­cium
	Cm = 96,  ///< Curium
	Bk = 97,  ///< Berkel­ium
	Cf = 98,  ///< Califor­nium
	Es = 99,  ///< Einstei­nium
	Fm = 100, ///< Fer­mium
	Md = 101, ///< Mende­levium
	No = 102, ///< Nobel­ium
	Lr = 103, ///< Lawren­cium

	D = 119, ///< Deuterium
};

// --------------------------------------------------------------------

/// An enum used to select the desired radius for an atom.
/// All values are collected from the wikipedia pages on atom radii

enum class radius_type
{
	calculated, ///< Calculated radius from theoretical models
	empirical,  ///< Empirically measured covalent radii

	/// @deprecated It is a bit unclear where these values came from. So, better not use them
	covalent_empirical,

	single_bond, ///< Bond length for a single covalent bond calculated using statistically analysis
	double_bond, ///< Bond length for a double covalent bond calculated using statistically analysis
	triple_bond, ///< Bond length for a triple covalent bond calculated using statistically analysis

	van_der_waals, ///< Radius of an imaginary hard sphere representing the distance of closest approach for another atom

	type_count ///< Number of radii
};

/// @brief The number of radii per element which can be requested from atom_type_info
constexpr std::size_t kRadiusTypeCount = static_cast<std::size_t>(radius_type::type_count);

/// An enum used to select either the effective or the crystal radius of an ion.
/// See explanation on Wikipedia: https://en.wikipedia.org/wiki/Ionic_radius

enum class ionic_radius_type
{
	effective, ///< Based on distance between ions in a crystal structure as determined by X-ray crystallography
	crystal    ///< Calculated ion radius based on a function of ionic charge and spin
};

/// Requests for an unknown radius value return kNA
constexpr float kNA = std::numeric_limits<float>::quiet_NaN();

/// A struct holding the known information for all elements defined in atom_type

struct atom_type_info
{
	/// The type as an atom_type
	atom_type type;

	/// The official name for this element
	std::string name;

	/// The official symbol for this element
	std::string symbol;

	/// The weight of this element
	float weight;

	/// A flag indicating whether the element is a metal
	bool metal;

	/// Array containing all known radii for this element. A value of kNA is
	/// stored for unknown values
	std::array<float, kRadiusTypeCount> radii;
};

/// Array of atom_type_info struct for each of the defined elements in atom_type

extern CIFPP_EXPORT const atom_type_info kKnownAtoms[];

// --------------------------------------------------------------------
// AtomTypeTraits

/// A traits class to access information for known elements

class atom_type_traits
{
  public:
	/// Constructor taking an atom_type \a a
	atom_type_traits(atom_type a);

	/// Constructor based on the element as a string in \a symbol
	atom_type_traits(const std::string &symbol);

	[[nodiscard]] atom_type type() const { return m_info->type; }       ///< Returns the atom_type
	[[nodiscard]] std::string name() const { return m_info->name; }     ///< Returns the name of the element
	[[nodiscard]] std::string symbol() const { return m_info->symbol; } ///< Returns the symbol of the element
	[[nodiscard]] float weight() const { return m_info->weight; }       ///< Returns the average weight of the element

	[[nodiscard]] bool is_metal() const { return m_info->metal; } ///< Returns true if the element is a metal

	/// Return true if the symbol in \a symbol actually exists in the list of known elements in atom_type
	static bool is_element(const std::string &symbol);

	/// Return true if the symbol in \a symbol exists and is a metal
	static bool is_metal(const std::string &symbol);

	/// @brief Return the radius for the element, use \a type to select which radius to return
	/// @param type The selector for which radius to return
	/// @return The requested radius or kNA if not known (or applicable)
	[[nodiscard]] float radius(radius_type type = radius_type::single_bond) const
	{
		if (type >= radius_type::type_count)
			throw std::invalid_argument("invalid radius requested");
		return m_info->radii[static_cast<std::size_t>(type)] / 100.f;
	}

	/// \brief Return the radius for a charged version of this atom in a solid crystal
	///
	/// \param charge  The charge of the ion
	/// \return        The radius of the ion
	[[nodiscard]] float crystal_ionic_radius(int charge) const;

	/// \brief Return the radius for a charged version of this atom in a non-solid environment
	///
	/// \param charge  The charge of the ion
	/// \return        The radius of the ion
	[[nodiscard]] float effective_ionic_radius(int charge) const;

	/// \brief Return the radius for a charged version of this atom, returns the effective radius by default
	///
	/// \param charge  The charge of the ion
	/// \param type    The requested ion radius type
	/// \return        The radius of the ion
	[[nodiscard]] float ionic_radius(int charge, ionic_radius_type type = ionic_radius_type::effective) const
	{
		return type == ionic_radius_type::effective ? effective_ionic_radius(charge) : crystal_ionic_radius(charge);
	}

	/**
	 * @brief data type encapsulating the scattering factors
	 * in a simplified form (only a and b).
	 */
	struct SFData
	{
		/** @cond */
		double a[6], b[6];
		/** @endcond */
	};

	/// @brief to get the Cval and Siva scattering factor values, use this constant as charge:
	static constexpr int kWKSFVal = -99;

	/// @brief Return the Waasmaier & Kirfel scattering factor values for the element
	///
	/// The coefficients from Waasmaier & Kirfel (1995), Acta Cryst. A51, 416-431.
	///
	/// @param charge The charge for which the structure values should be returned, use kWSKFVal to return the *Cval* and *Siva* values
	/// @return The scattering factors as a SFData struct
	[[nodiscard]] const SFData &wksf(int charge = 0) const;

	/// @brief Return the electron scattering factor values for the element
	///
	/// @return The scattering factors as a SFData struct
	[[nodiscard]] const SFData &elsf() const;

	/// Clipper doesn't like atoms with charges that do not have a scattering factor. And
	/// rightly so, but we need to know in advance if this is the case
	[[nodiscard]] bool has_sf(int charge) const;

  private:
	const struct atom_type_info *m_info;
};

} // namespace cif
