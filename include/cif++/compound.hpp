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

/// \file This file contains the definition for the class compound, encapsulating
/// the information found for compounds in the CCD.

#include <map>
#include <set>
#include <tuple>
#include <vector>

#include <cif++.hpp>
#include <cif++/atom_type.hpp>

namespace cif
{

// --------------------------------------------------------------------

class compound;
struct compound_atom;
class compound_factory_impl;

/// \brief The bond type as defined in the CCD, possible values taken from the mmcif_pdbx file
enum class bond_type
{
	sing, // 'single bond'
	doub, // 'double bond'
	trip, // 'triple bond'
	quad, // 'quadruple bond'
	arom, // 'aromatic bond'
	poly, // 'polymeric bond'
	delo, // 'delocalized double bond'
	pi,   // 'pi bond'
};

std::string to_string(bond_type bondType);
bond_type from_string(const std::string &bondType);

/// --------------------------------------------------------------------
/// \brief struct containing information about an atom in a chemical compound.
/// This is a subset of the available information. Contact the author if you need more fields.

struct compound_atom
{
	std::string id;
	atom_type type_symbol;
	int charge = 0;
	bool aromatic = false;
	bool leaving_atom = false;
	bool stereo_config = false;
	float x, y, z;
};

/// --------------------------------------------------------------------
/// \brief struct containing information about the bonds

struct compound_bond
{
	std::string atom_id[2];
	bond_type type;
	bool aromatic = false, stereo_config = false;
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

	std::string id() const { return m_id; }
	std::string name() const { return m_name; }
	std::string type() const { return m_type; }
	std::string group() const { return m_group; }
	std::string formula() const { return m_formula; }
	float formula_weight() const { return m_formula_weight; }
	int formal_charge() const { return m_formal_charge; }

	const std::vector<compound_atom> &atoms() const { return m_atoms; }
	const std::vector<compound_bond> &bonds() const { return m_bonds; }

	compound_atom get_atom_by_atom_id(const std::string &atom_id) const;

	bool atoms_bonded(const std::string &atomId_1, const std::string &atomId_2) const;
	// float atomBondValue(const std::string &atomId_1, const std::string &atomId_2) const;
	// float bondAngle(const std::string &atomId_1, const std::string &atomId_2, const std::string &atomId_3) const;
	// float chiralVolume(const std::string &centreID) const;

	bool is_water() const
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
	static compound_factory &instance();
	static void clear();

	void set_default_dictionary(const std::filesystem::path &inDictFile);
	void push_dictionary(const std::filesystem::path &inDictFile);
	void pop_dictionary();

	bool is_known_peptide(const std::string &res_name) const;
	bool is_known_base(const std::string &res_name) const;

	/// \brief Create the compound object for \a id
	///
	/// This will create the compound instance for \a id if it doesn't exist already.
	/// The result is owned by this factory and should not be deleted by the user.
	/// \param id	The compound ID, a three letter code usually
	/// \result		The compound, or nullptr if it could not be created (missing info)
	const compound *create(std::string id);

	~compound_factory();

	static const std::map<std::string, char> kAAMap, kBaseMap;

  private:
	compound_factory();

	compound_factory(const compound_factory &) = delete;
	compound_factory &operator=(const compound_factory &) = delete;

	static std::unique_ptr<compound_factory> s_instance;
	static thread_local std::unique_ptr<compound_factory> tl_instance;
	static bool s_use_thread_local_instance;

	std::shared_ptr<compound_factory_impl> m_impl;
};

} // namespace pdbx
