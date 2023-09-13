/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2020-2022 NKI/AVL, Netherlands Cancer Institute
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

#include "cif++/atom_type.hpp"
#include "cif++/datablock.hpp"
#include "cif++/exports.hpp"
#include "cif++/point.hpp"
#include "cif++/utilities.hpp"

#include <map>
#include <set>
#include <tuple>
#include <vector>

/// \file compound.hpp
/// This file contains the definition for the class compound, encapsulating
/// the information found for compounds in the CCD.
///
/// The data is loaded by default from a file called `components.cif`. This file
/// is located using load_resource. (See documentation on cif::load_resource for more information)
///
/// But if the CCP4 environment is available at runtime, the compound information
/// may also be generated from the CCP4 monomer library.
///
/// Note that the information in CCP4 and CCD is not equal.
///

/// See also :doc:`/compound` for more information.

namespace cif
{

// --------------------------------------------------------------------

class compound;
struct compound_atom;
class compound_factory_impl;

/// \brief The bond type or bond order as defined in the CCD, possible values taken from the mmcif_pdbx file
enum class bond_type
{
	sing, ///< single bond
	doub, ///< double bond
	trip, ///< triple bond
	quad, ///< quadruple bond
	arom, ///< aromatic bond
	poly, ///< polymeric bond
	delo, ///< delocalized double bond
	pi,   ///< pi bond
};

/// @brief return the string representation of @a bondType
std::string bond_type_to_string(bond_type bondType);

/// @brief return the cif::bond_type for the string representation @a bondType
bond_type parse_bond_type_from_string(const std::string &bondType);

/// \brief The possible stereo config values for a compound_atom.
///
/// As the site https://psiberg.com/r-s-nomenclature/ states:
///
/// > RS nomenclature is currently the preferred system for assigning absolute
/// > configuration to chiral molecules. The letters R and S come from the Latin
/// > words ‘Rectus‘ and ‘Sinister‘ meaning ‘right’ and ‘left’. Molecules that
/// > rotate the plane of polarized light to right are referred to as ‘R isomers’
/// > and the molecules that rotate the plane of polarized light to left are
/// > referred to ‘S isomers’.
enum class stereo_config_type : uint8_t
{
	N = 'N', ///< Not polarizing
	R = 'R', ///< Rectus
	S = 'S'  ///< Sinister
};

/// @brief return the string representation of @a stereo_config
std::string to_string(stereo_config_type stereo_config);

/// @brief return the cif::stereo_config_type for the string representation @a stereo_config
stereo_config_type parse_stereo_config_from_string(const std::string &stereo_config);

/// --------------------------------------------------------------------
/// \brief struct containing information about an atom in a chemical compound.
/// This is a subset of the available information. Contact the author if you need more fields.

struct compound_atom
{
	std::string id;                                           ///< Identifier for each atom in the chemical component
	atom_type type_symbol;                                    ///< The element type for each atom in the chemical component.
	int charge = 0;                                           ///< The formal charge assigned to each atom in the chemical component.
	bool aromatic = false;                                    ///< Defines atoms in an aromatic moiety
	bool leaving_atom = false;                                ///< Flags atoms with "leaving" capability
	stereo_config_type stereo_config = stereo_config_type::N; ///< Defines the stereochemical configuration of the chiral center atom.
	float x,                                                  ///< The x component of the coordinates for each atom specified as orthogonal angstroms.
		y,                                                    ///< The y component of the coordinates for each atom specified as orthogonal angstroms.
		z;                                                    ///< The z component of the coordinates for each atom specified as orthogonal angstroms.

	/// Return the location of the atom as a point
	point get_location() const
	{
		return { x, y, z };
	}
};

/// --------------------------------------------------------------------
/// \brief struct containing information about the bonds

struct compound_bond
{
	std::string atom_id[2];    ///< The ID's of the two atoms that define the bond.
	bond_type type;            ///< The bond order of the chemical bond associated with the specified atoms.
	bool aromatic = false,     ///< Defines aromatic bonds.
		stereo_config = false; ///< Defines stereochemical bonds.
};

/// --------------------------------------------------------------------
/// \brief a class that contains information about a chemical compound.
/// This information is derived from the CDD by default.
///
/// To create compounds, you use the factory method. You can add your own
/// compound definitions by calling the addExtraComponents function and
/// pass it a valid CCD formatted file.

class compound
{
  public:
	// accessors

	std::string id() const { return m_id; }                   ///< Return the alphanumeric code for the chemical component.
	std::string name() const { return m_name; }               ///< Return the name of the chemical component.
	std::string type() const { return m_type; }               ///< Return the type of monomer.
	std::string formula() const { return m_formula; }         ///< Return the chemical formula of the chemical component.
	float formula_weight() const { return m_formula_weight; } ///< Return the formula mass of the chemical component in Daltons.
	int formal_charge() const { return m_formal_charge; }     ///< Return the formal charge on the chemical component.

	/// The group record is only available in CCP4 monomer library files.
	/// For CCD entries this value will always contain 'non-polymer'
	std::string group() const { return m_group; }

	const std::vector<compound_atom> &atoms() const { return m_atoms; } ///< Return the list of atoms for this compound
	const std::vector<compound_bond> &bonds() const { return m_bonds; } ///< Return the list of bonds for this compound

	compound_atom get_atom_by_atom_id(const std::string &atom_id) const; ///< Return the atom with id @a atom_id

	bool atoms_bonded(const std::string &atomId_1, const std::string &atomId_2) const; ///< Return true if @a atomId_1 is bonded to @a atomId_2
	float bond_length(const std::string &atomId_1, const std::string &atomId_2) const; ///< Return the bond length between @a atomId_1 and @a atomId_2

	bool is_water() const ///< Return if the compound is actually a water
	{
		return m_id == "HOH" or m_id == "H2O" or m_id == "WAT";
	}

  private:
	friend class compound_factory_impl;
	friend class CCD_compound_factory_impl;
	friend class CCP4_compound_factory_impl;

	compound(cif::datablock &db);
	compound(cif::datablock &db, const std::string &id, const std::string &name, const std::string &type, const std::string &group);

	std::string m_id;
	std::string m_name;
	std::string m_type;
	std::string m_group;
	std::string m_formula;
	float m_formula_weight = 0;
	int m_formal_charge = 0;
	std::vector<compound_atom> m_atoms;
	std::vector<compound_bond> m_bonds;
};

// --------------------------------------------------------------------
// Factory class for compound and Link objects

/// Use the compound_factory singleton instance to create compound objects

class compound_factory
{
  public:
	/// \brief Initialise a singleton instance.
	///
	/// If you have a multithreaded application and want to have different
	/// compounds in each thread (e.g. a web service processing user requests
	/// with different sets of compounds) you can set the \a useThreadLocalInstanceOnly
	/// flag to true.

	static void init(bool useThreadLocalInstanceOnly);

	/// Return the singleton instance. If initialized with local threads, this is the
	/// instance for the current thread.
	static compound_factory &instance();

	/// Delete and reset the singleton instance. If initialized with local threads, this is the
	/// instance for the current thread.
	static void clear();

	/// Set the default dictionary file to @a inDictFile
	void set_default_dictionary(const std::filesystem::path &inDictFile);

	/// Override any previously loaded dictionary with @a inDictFile
	void push_dictionary(const std::filesystem::path &inDictFile);

	/// Remove the last pushed dictionary
	void pop_dictionary();

	/// Return whether @a res_name is a valid and known peptide
	bool is_known_peptide(const std::string &res_name) const;

	/// Return whether @a res_name is a valid and known base
	bool is_known_base(const std::string &res_name) const;

	/// \brief Create the compound object for \a id
	///
	/// This will create the compound instance for \a id if it doesn't exist already.
	/// The result is owned by this factory and should not be deleted by the user.
	/// \param id	The compound ID, a three letter code usually
	/// \result		The compound, or nullptr if it could not be created (missing info)
	const compound *create(std::string id);

	~compound_factory();

	CIFPP_EXPORT static const std::map<std::string, char> kAAMap, ///< Globally accessible static list of the default amino acids
		kBaseMap;                                                 ///< Globally accessible static list of the default bases

  private:
	compound_factory();

	compound_factory(const compound_factory &) = delete;
	compound_factory &operator=(const compound_factory &) = delete;

	static std::unique_ptr<compound_factory> s_instance;
	static thread_local std::unique_ptr<compound_factory> tl_instance;
	static bool s_use_thread_local_instance;

	std::shared_ptr<compound_factory_impl> m_impl;
};

} // namespace cif
