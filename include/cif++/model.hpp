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

#include "cif++/atom_type.hpp"
#include "cif++/datablock.hpp"
#include "cif++/point.hpp"
#include "cif++/row.hpp"

#include <memory>

#if __cpp_lib_format
# include <format>
# include <utility>
#endif

/** @file model.hpp
 *
 * This file contains code to work with models of molecules.
 *
 * The classes available encapsulate the real world concepts of
 * atoms, residues, monomers, polymers and everything is then
 * bound together in a structure.
 *
 * This code is not finished yet, ideally it would be a high
 * level interface to manipulate macro molecular structures
 * and an attempt has been made to start work on this. But
 * there's still a lot that needs to be implemented.
 *
 * However, the code that is here is still useful in
 * manipulating the underlying mmCIF data model.
 *
 */

namespace cif::mm
{

class atom;
class residue;
class monomer;
class polymer;
class structure;

// --------------------------------------------------------------------

/**
 * @brief The class atom encapsulates the data in _atom_site and
 * _atom_site_anisotrop
 *
 * The class atom is a kind of flyweight class. It can be copied
 * with low overhead. All data is stored in the underlying mmCIF
 * categories but some very often used items are cached in the
 * impl.
 *
 * It is also possible to have symmetry copies of atoms. They
 * share the same data in the cif::category but their location
 * differs by using a symmetry operator.
 */

class atom
{
  private:
	/** @cond */
	struct atom_impl : public std::enable_shared_from_this<atom_impl>
	{
		atom_impl(const datablock &db, std::string_view id)
			: m_db(db)
			, m_cat(db["atom_site"])
			, m_id(id)
		{
			auto r = row();
			if (r)
				tie(m_location.m_x, m_location.m_y, m_location.m_z) = r.get("Cartn_x", "Cartn_y", "Cartn_z");
		}

		// constructor for a symmetry copy of an atom
		atom_impl(const atom_impl &impl, const point &loc, const std::string &sym_op)
			: atom_impl(impl)
		{
			m_location = loc;
			m_symop = sym_op;
		}

		atom_impl(const atom_impl &i) = default;

		[[nodiscard]] int compare(const atom_impl &b) const;

		// bool getAnisoU(float anisou[6]) const;

		[[nodiscard]] int get_charge() const;

		void moveTo(const point &p);

		// const compound *compound() const;

		[[nodiscard]] std::string get_property(std::string_view name) const;
		[[nodiscard]] int get_property_int(std::string_view name) const;
		[[nodiscard]] float get_property_float(std::string_view name) const;

		void set_property(const std::string_view name, const std::string &value);

		row_handle row()
		{
			return m_cat[{ { .name = "id", .value = m_id } }];
		}

		[[nodiscard]] const row_handle row() const
		{
			return m_cat[{ { .name = "id", .value = m_id } }];
		}

		row_handle row_aniso()
		{
			row_handle result{};
			auto cat = m_db.get("atom_site_anisotrop");
			if (cat)
				result = cat->operator[]({ { .name = "id", .value = m_id } });
			return result;
		}

		[[nodiscard]] const row_handle row_aniso() const
		{
			row_handle result{};
			auto cat = m_db.get("atom_site_anisotrop");
			if (cat)
				result = cat->operator[]({ { .name = "id", .value = m_id } });
			return result;
		}

		const datablock &m_db;
		const category &m_cat;
		std::string m_id;
		point m_location;
		std::string m_symop = "1_555";
	};
	/** @endcond */

  public:
	/**
	 * @brief Construct a new, empty atom object
	 */
	atom() = default;

	/**
	 * @brief Construct a new atom object using @a impl as impl
	 *
	 * @param impl The implementation objectt
	 */
	atom(std::shared_ptr<atom_impl> impl)
		: m_impl(std::move(impl))
	{
	}

	/**
	 * @brief Copy construct a new atom object
	 */
	atom(const atom &rhs) // NOLINT(modernize-use-equals-default)
		: m_impl(rhs.m_impl)
	{
	}

	/**
	 * @brief Construct a new atom object based on a cif::row
	 *
	 * @param db The datablock where the _atom_site category resides
	 * @param row The row containing the data for this atom
	 */
	atom(const datablock &db, const row_handle &row)
		: atom(std::make_shared<atom_impl>(db, row["id"].as<std::string>()))
	{
	}

	/**
	 * @brief A special constructor to create symmetry copies
	 *
	 * @param rhs The original atom to copy
	 * @param symmmetry_location The symmetry location
	 * @param symmetry_operation The symmetry operator used
	 */
	atom(const atom &rhs, const point &symmmetry_location, const std::string &symmetry_operation)
		: atom(std::make_shared<atom_impl>(*rhs.m_impl, symmmetry_location, symmetry_operation))
	{
	}

	/// \brief To quickly test if the atom has data
	explicit operator bool() const { return m_impl.operator bool(); }

	/// \brief Copy assignement operator
	atom &operator=(const atom &rhs) = default;

	/// \brief Return the item named @a name in the _atom_site category for this atom
	[[nodiscard]] std::string get_property(std::string_view name) const
	{
		if (not m_impl)
			throw std::logic_error("Error trying to fetch a property from an uninitialized atom");
		return m_impl->get_property(name);
	}

	/// \brief Return the item named @a name in the _atom_site category for this atom cast to an int
	[[nodiscard]] int get_property_int(std::string_view name) const
	{
		if (not m_impl)
			throw std::logic_error("Error trying to fetch a property from an uninitialized atom");
		return m_impl->get_property_int(name);
	}

	/// \brief Return the item named @a name in the _atom_site category for this atom cast to a float
	[[nodiscard]] float get_property_float(std::string_view name) const
	{
		if (not m_impl)
			throw std::logic_error("Error trying to fetch a property from an uninitialized atom");
		return m_impl->get_property_float(name);
	}

	/// \brief Set value for the item named @a name in the _atom_site category to @a value
	void set_property(const std::string_view name, const std::string &value)
	{
		if (not m_impl)
			throw std::logic_error("Error trying to modify an uninitialized atom");
		m_impl->set_property(name, value);
	}

	/// \brief Set value for the item named @a name in the _atom_site category to @a value
	template <typename T>
	void set_property(const std::string_view name, const T &value)
		requires(std::is_arithmetic_v<T>)
	{
		set_property(name, std::to_string(value));
	}

	/** Return the ID of the _atom_site record.
	 *
	 * @note Although I've never seen anything other than integers,
	 * the standard says this should be a string and so we use that.
	 */
	[[nodiscard]] const std::string &id() const { return impl().m_id; }

	/// \brief Return the type of the atom
	[[nodiscard]] cif::atom_type get_type() const { return atom_type_traits(get_property("type_symbol")).type(); }

	/// \brief Return the cached location of this atom
	[[nodiscard]] point get_location() const { return impl().m_location; }

	/// \brief Set the location of this atom, will set both the cached data as well as the data in the underlying _atom_site category
	void set_location(point p)
	{
		if (not m_impl)
			throw std::logic_error("Error trying to modify an uninitialized atom");
		m_impl->moveTo(p);
	}

	/// \brief Translate the position of this atom by \a t
	void translate(point t)
	{
		set_location(get_location() + t);
	}

	/// \brief Rotate the position of this atom by \a q
	void rotate(quaternion q)
	{
		auto loc = get_location();
		loc.rotate(q);
		set_location(loc);
	}

	/// \brief rotate the coordinates of this atom by \a q around point \a p
	void rotate(quaternion q, point p)
	{
		auto loc = get_location();
		loc.rotate(q, p);
		set_location(loc);
	}

	/// \brief Translate and rotate the position of this atom by \a t and \a q
	void translate_and_rotate(point t, quaternion q)
	{
		auto loc = get_location();
		loc += t;
		loc.rotate(q);
		set_location(loc);
	}

	/// \brief Translate, rotate and translate again the coordinates this atom by \a t1 , \a q and \a t2
	void translate_rotate_and_translate(point t1, quaternion q, point t2)
	{
		auto loc = get_location();
		loc += t1;
		loc.rotate(q);
		loc += t2;
		set_location(loc);
	}

	/// for direct access to underlying data, be careful!
	[[nodiscard]] const row_handle get_row() const { return impl().row(); }

	/// for direct access to underlying data, be careful!
	[[nodiscard]] const row_handle get_row_aniso() const { return impl().row_aniso(); }

	/// Return if the atom is actually a symmetry copy or the original one
	[[nodiscard]] bool is_symmetry_copy() const { return impl().m_symop != "1_555"; }

	/// Return the symmetry operator used
	[[nodiscard]] std::string symmetry() const { return impl().m_symop; }

	/// Return true if this atom is part of a water molecule
	[[nodiscard]] bool is_water() const
	{
		auto comp_id = get_label_comp_id();
		return comp_id == "HOH" or comp_id == "H2O" or comp_id == "WAT";
	}

	/// Return the charge
	[[nodiscard]] int get_charge() const { return impl().get_charge(); }

	/// Return the occupancy
	[[nodiscard]] float get_occupancy() const { return get_property_float("occupancy"); }

	// specifications

	[[nodiscard]] std::string get_label_asym_id() const { return get_property("label_asym_id"); }     ///< Return the label_asym_id property
	[[nodiscard]] int get_label_seq_id() const { return get_property_int("label_seq_id"); }           ///< Return the label_seq_id property
	[[nodiscard]] std::string get_label_atom_id() const { return get_property("label_atom_id"); }     ///< Return the label_atom_id property
	[[nodiscard]] std::string get_label_alt_id() const { return get_property("label_alt_id"); }       ///< Return the label_alt_id property
	[[nodiscard]] std::string get_label_comp_id() const { return get_property("label_comp_id"); }     ///< Return the label_comp_id property
	[[nodiscard]] std::string get_label_entity_id() const { return get_property("label_entity_id"); } ///< Return the label_entity_id property

	[[nodiscard]] std::string get_auth_asym_id() const { return get_property("auth_asym_id"); }      ///< Return the auth_asym_id property
	[[nodiscard]] std::string get_auth_seq_id() const { return get_property("auth_seq_id"); }        ///< Return the auth_seq_id property
	[[nodiscard]] std::string get_auth_atom_id() const { return get_property("auth_atom_id"); }      ///< Return the auth_atom_id property
	[[nodiscard]] std::string get_auth_alt_id() const { return get_property("pdbx_auth_alt_id"); }   ///< Return the auth_alt_id property
	[[nodiscard]] std::string get_auth_comp_id() const { return get_property("auth_comp_id"); }      ///< Return the auth_comp_id property
	[[nodiscard]] std::string get_pdb_ins_code() const { return get_property("pdbx_PDB_ins_code"); } ///< Return the pdb_ins_code property

	/// Return true if this atom is an alternate
	[[nodiscard]] bool is_alternate() const
	{
		if (auto alt_id = get_label_alt_id(); alt_id.empty() or alt_id == ".")
			return false;
		return true;
	}

	/// Convenience method to return a string that might be ID in PDB space
	[[nodiscard]] std::string pdb_id() const
	{
		return get_label_comp_id() + '_' + get_auth_asym_id() + '_' + get_auth_seq_id() + get_pdb_ins_code();
	}

	/// Compare two atoms
	bool operator==(const atom &rhs) const
	{
		if (m_impl == rhs.m_impl)
			return true;

		if (not(m_impl and rhs.m_impl))
			return false;

		return &m_impl->m_db == &rhs.m_impl->m_db and m_impl->m_id == rhs.m_impl->m_id;
	}

	/// Compare two atoms
	bool operator!=(const atom &rhs) const
	{
		return not operator==(rhs);
	}

	/// Is this atom a backbone atom
	[[nodiscard]] bool is_back_bone() const
	{
		auto atomID = get_label_atom_id();
		return atomID == "N" or atomID == "O" or atomID == "C" or atomID == "CA";
	}

	/// swap
	void swap(atom &b)
	{
		std::swap(m_impl, b.m_impl);
	}

	/// Compare this atom with @a b
	[[nodiscard]] int compare(const atom &b) const { return impl().compare(*b.m_impl); }

	/// Should this atom sort before @a rhs
	bool operator<(const atom &rhs) const
	{
		return compare(rhs) < 0;
	}

	/// Write the atom to std::ostream @a os
	friend std::ostream &operator<<(std::ostream &os, const atom &atom);

  private:
	friend class structure;

	[[nodiscard]] const atom_impl &impl() const
	{
		if (not m_impl)
			throw std::runtime_error("Uninitialized atom, not found?");
		return *m_impl;
	}

	std::shared_ptr<atom_impl> m_impl;
};

/** swap */
inline void swap(atom &a, atom &b)
{
	a.swap(b);
}

/** Calculate the distance between atoms @a and @a b in ångström */
inline float distance(const atom &a, const atom &b)
{
	return distance(a.get_location(), b.get_location());
}

/** Calculate the square of the distance between atoms @a and @a b in ångström
 *
 * @note Use this whenever possible instead of simply using distance since
 * this function does not have to calculate a square root which is expensive.
 */
inline float distance_squared(const atom &a, const atom &b)
{
	return distance_squared(a.get_location(), b.get_location());
}

// --------------------------------------------------------------------

/**
 * @brief The entity types that can be found in a mmCIF file
 *
 */
enum class EntityType
{
	Polymer,    ///< entity is a polymer
	NonPolymer, ///< entity is not a polymer
	Macrolide,  ///< entity is a macrolide
	Water,      ///< water in the solvent model
	Branched    ///< entity is branched
};

// --------------------------------------------------------------------

/**
 * @brief The class residue is a collection of atoms forming a molecule
 *
 * This class is used to store ligand e.g. Derived classes are monomer
 * and sugar.
 */

class residue
{
  public:
	friend class structure;

	/**
	 * @brief Construct a new residue object based on key items
	 */
	residue(structure &structure, std::string compoundID, std::string asymID, int seqID,
		std::string authAsymID, std::string authSeqID, std::string pdbInsCode)
		: m_structure(&structure)
		, m_compound_id(std::move(compoundID))
		, m_asym_id(std::move(asymID))
		, m_seq_id(seqID)
		, m_pdb_strand_id(std::move(authAsymID))
		, m_pdb_seq_num(std::move(authSeqID))
		, m_pdb_ins_code(std::move(pdbInsCode))
	{
	}

	/** Construct a new residue in structure with the atoms in @a atoms */
	residue(structure &structure, const std::vector<atom> &atoms);

	/** @cond */
	residue(const residue &rhs) = default;
	residue(residue &&rhs)
	{
		swap(*this, rhs);
	}

	residue &operator=(residue rhs)
	{
		swap(*this, rhs);
		return *this;
	}

	friend void swap(residue &a, residue &b) noexcept
	{
		if (&a != &b)
		{
			std::swap(a.m_structure, b.m_structure);
			std::swap(a.m_asym_id, b.m_asym_id);
			std::swap(a.m_seq_id, b.m_seq_id);
			std::swap(a.m_pdb_ins_code, b.m_pdb_ins_code);
			std::swap(a.m_atoms, b.m_atoms);
		}
	}

	virtual ~residue() = default;
	/** @endcond */

	/** Return the entity_id of this residue */
	[[nodiscard]] std::string get_entity_id() const;

	/** Return the entity type of this residue */
	[[nodiscard]] EntityType entity_type() const;

	[[nodiscard]] const std::string &get_asym_id() const { return m_asym_id; } ///< Return the asym_id
	[[nodiscard]] int get_seq_id() const { return m_seq_id; }                  ///< Return the seq_id

	[[nodiscard]] const std::string get_pdb_strand_id() const { return m_pdb_strand_id; } ///< Return the pdb_strand_id
	[[nodiscard]] const std::string get_pdb_seq_num() const { return m_pdb_seq_num; }     ///< Return the pdb_seq_num
	[[nodiscard]] std::string get_pdb_ins_code() const { return m_pdb_ins_code; }         ///< Return the pdb_ins_code

	[[nodiscard]] const std::string &get_compound_id() const { return m_compound_id; } ///< Return the compound_id
	void set_compound_id(const std::string &id) { m_compound_id = id; }                ///< Set the compound_id to @a id

	/** Return the structure this residue belongs to */
	[[nodiscard]] structure *get_structure() const { return m_structure; }

	/** Return a list of the atoms for this residue */
	std::vector<atom> &atoms()
	{
		return m_atoms;
	}

	/** Return a const list of the atoms for this residue */
	[[nodiscard]] const std::vector<atom> &atoms() const
	{
		return m_atoms;
	}

	/** Add atom @a atom to the atoms in this residue */
	void add_atom(atom &atom);

	/// \brief Unique atoms returns only the atoms without alternates and the first of each alternate atom id.
	[[nodiscard]] std::vector<atom> unique_atoms() const;

	/// \brief Return the atom with atom_id @a atomID
	[[nodiscard]] atom get_atom_by_atom_id(const std::string &atomID) const;

	/// \brief Return the atom with atom_id @a atomID and alternate_id @a altID
	[[nodiscard]] atom get_atom_by_atom_id(const std::string &atomID, const std::string &altID) const;

	/// \brief Return the list of atoms having ID \a atomID
	///
	/// This includes all alternate atoms with this ID
	/// whereas get_atom_by_atom_id only returns the first unique atom
	[[nodiscard]] std::vector<atom> get_atoms_by_id(const std::string &atomID) const;

	/// \brief Is this residue a single entity?
	[[nodiscard]] bool is_entity() const;

	/// \brief Is this residue a water molecule?
	[[nodiscard]] bool is_water() const { return m_compound_id == "HOH"; }

	/// \brief Return true if this residue has alternate atoms
	[[nodiscard]] bool has_alternate_atoms() const;

	/// \brief Return true if this residue has alternate atoms for the atom \a atomID
	[[nodiscard]] bool has_alternate_atoms_for(const std::string &atomID) const;

	/// \brief Return the list of unique alt ID's present in this residue
	[[nodiscard]] std::set<std::string> get_alternate_ids() const;

	/// \brief Return the list of unique atom ID's
	[[nodiscard]] std::set<std::string> get_atom_ids() const;

	/// \brief Return a tuple containing the center location and the radius for the atoms of this residue
	[[nodiscard]] std::tuple<point, float> center_and_radius() const;

	/// \brief Write the residue @a res to the std::ostream @a os
	friend std::ostream &operator<<(std::ostream &os, const residue &res);

	/// \brief Return true if this residue is equal to @a rhs
	bool operator==(const residue &rhs) const
	{
		return this == &rhs or (m_structure == rhs.m_structure and
								   m_seq_id == rhs.m_seq_id and
								   m_asym_id == rhs.m_asym_id and
								   m_compound_id == rhs.m_compound_id and
								   m_pdb_seq_num == rhs.m_pdb_seq_num);
	}

	/// @brief Create a new atom and add it to the list
	/// @return newly created atom
	virtual atom create_new_atom(atom_type inType, const std::string &inAtomID, point inLocation);

  protected:
	/** @cond */
	residue() = default;

	structure *m_structure = nullptr;
	std::string m_compound_id, m_asym_id;
	int m_seq_id = 0;
	std::string m_pdb_strand_id, m_pdb_seq_num, m_pdb_ins_code;
	std::vector<atom> m_atoms;
	/** @endcond */
};

// --------------------------------------------------------------------

/**
 * @brief a monomer models a single residue in a protein chain
 *
 */

class monomer : public residue
{
  public:
	/// \brief constructor with actual values
	monomer(const polymer &polymer, std::size_t index, int seqID, const std::string &authSeqID,
		const std::string &pdbInsCode, const std::string &compoundID);

	/// \brief Copy constructor
	monomer(const monomer &rhs) = default;

	/// \brief Move constructor
	monomer(monomer &&rhs)
	{
		swap(*this, rhs);
	}
	monomer &operator=(monomer rhs)
	{
		swap(*this, rhs);
		return *this;
	}

	friend void swap(monomer &a, monomer &b) noexcept
	{
		assert(a.m_polymer == b.m_polymer);
		std::swap(a.m_index, b.m_index);
		swap(static_cast<residue &>(a), static_cast<residue &>(b));
	}

	[[nodiscard]] bool is_first_in_chain() const; ///< Return if this residue is the first residue in the chain
	[[nodiscard]] bool is_last_in_chain() const;  ///< Return if this residue is the last residue in the chain

	[[nodiscard]] const monomer &prev() const; // Return previous monomer in polymer
	[[nodiscard]] const monomer &next() const; // Return next monomer in polymer

	// convenience
	[[nodiscard]] bool has_alpha() const; ///< Return if a alpha value can be calculated (depends on location in chain)
	[[nodiscard]] bool has_kappa() const; ///< Return if a kappa value can be calculated (depends on location in chain)

	// Assuming this is really an amino acid...

	[[nodiscard]] float phi() const;   ///< Return the phi value for this residue
	[[nodiscard]] float psi() const;   ///< Return the psi value for this residue
	[[nodiscard]] float alpha() const; ///< Return the alpha value for this residue
	[[nodiscard]] float kappa() const; ///< Return the kappa value for this residue
	[[nodiscard]] float tco() const;   ///< Return the tco value for this residue
	[[nodiscard]] float omega() const; ///< Return the omega value for this residue

	// torsion angles
	[[nodiscard]] std::size_t nr_of_chis() const; ///< Return how many torsion angles can be calculated
	[[nodiscard]] float chi(std::size_t i) const; ///< Return torsion angle @a i

	[[nodiscard]] bool is_cis() const; ///< Return true if this residue is in a cis conformation

	/// \brief Returns true if the four atoms C, CA, N and O are present
	[[nodiscard]] bool is_complete() const;

	/// \brief Returns true if any of the backbone atoms has an alternate
	[[nodiscard]] bool has_alternate_backbone_atoms() const;

	[[nodiscard]] atom CAlpha() const { return get_atom_by_atom_id("CA"); } ///< Return the CAlpha atom
	[[nodiscard]] atom C() const { return get_atom_by_atom_id("C"); }       ///< Return the C atom
	[[nodiscard]] atom N() const { return get_atom_by_atom_id("N"); }       ///< Return the N atom
	[[nodiscard]] atom O() const { return get_atom_by_atom_id("O"); }       ///< Return the O atom
	[[nodiscard]] atom H() const { return get_atom_by_atom_id("H"); }       ///< Return the H atom

	/// \brief Return true if this monomer is bonded to monomer @a rhs
	[[nodiscard]] bool is_bonded_to(const monomer &rhs) const
	{
		return this != &rhs and are_bonded(*this, rhs);
	}

	/**
	 * @brief Return true if the distance between the CA atoms of the
	 * two monomers @a a and @a b are within the expected range with
	 * an error margin of @a errorMargin.
	 *
	 * The expected distance is 3.0 ångström for a cis conformation
	 * and 3.8 ångström for trans.
	 */
	static bool are_bonded(const monomer &a, const monomer &b, float errorMargin = 0.5f);

	/// \brief Return true if the bond between @a a and @a b is cis
	static bool is_cis(const monomer &a, const monomer &b);

	/// \brief Return the omega angle between @a a and @a b
	static float omega(const monomer &a, const monomer &b);

	/// \brief Return the chiral volume, only for LEU and VAL
	[[nodiscard]] float chiral_volume() const;

	/// \brief Compare this monomer with \a rhs
	bool operator==(const monomer &rhs) const
	{
		return m_polymer == rhs.m_polymer and m_index == rhs.m_index;
	}

	atom create_new_atom(atom_type inType, const std::string &inAtomID, point inLocation) override;

  private:
	const polymer *m_polymer;
	std::size_t m_index{};
};

// --------------------------------------------------------------------

/**
 * @brief A polymer is simply a list of monomers
 *
 */
class polymer : public std::vector<monomer>
{
  public:
	/// \brief Constructor
	polymer(structure &s, std::string entityID, std::string asymID, std::string auth_asym_id);

	polymer(const polymer &) = delete;
	polymer &operator=(const polymer &) = delete;

	[[nodiscard]] structure *get_structure() const { return m_structure; } ///< Return the structure

	[[nodiscard]] std::string get_asym_id() const { return m_asym_id; }             ///< Return the asym_id
	[[nodiscard]] std::string get_pdb_strand_id() const { return m_pdb_strand_id; } ///< Return the PDB chain ID, actually
	[[nodiscard]] std::string get_entity_id() const { return m_entity_id; }         ///< Return the entity_id

  private:
	structure *m_structure;
	std::string m_entity_id;
	std::string m_asym_id;
	std::string m_pdb_strand_id;
};

// --------------------------------------------------------------------
// sugar and branch, to describe glycosylation sites

class branch;

/**
 * @brief A sugar is a residue that is part of a glycosylation site
 *
 */

class sugar : public residue
{
  public:
	/// \brief constructor
	sugar(branch &branch, const std::string &compoundID,
		const std::string &asymID, int authSeqID);

	/** @cond */
	sugar(const sugar &rhs) = default;

	sugar(sugar &&rhs)
	{
		swap(*this, rhs);
	}

	sugar &operator=(sugar rhs)
	{
		swap(*this, rhs);
		return *this;
	}

	friend void swap(sugar &a, sugar &b) noexcept
	{
		assert(a.m_branch == b.m_branch);
		std::swap(a.m_link, b.m_link);
		swap(static_cast<residue &>(a), static_cast<residue &>(b));
	}

	/** @endcond */

	/**
	 * @brief Return the sugar number in the glycosylation tree
	 *
	 * To store the sugar number, the auth_seq_id item has been overloaded
	 * in the specification. But since a sugar number should be, ehm, a number
	 * and auth_seq_id is specified to contain a string, we do a check here
	 * to see if it really is a number.
	 *
	 * @return The sugar number
	 */
	[[nodiscard]] int num() const
	{
		int result;
		auto r = std::from_chars(m_pdb_seq_num.data(), m_pdb_seq_num.data() + m_pdb_seq_num.length(), result);
		if (r.ec != std::errc{})
			throw std::runtime_error("The auth_seq_id should be a number for a sugar");
		return result;
	}

	/// \brief Return the name of this sugar
	[[nodiscard]] std::string name() const;

	/// \brief Return the atom the C1 is linked to
	[[nodiscard]] atom get_link() const { return m_link; }

	/// \brief Set the link atom C1 is linked to to @a link
	void set_link(atom link) { m_link = link; }

	/// \brief Return the sugar number of the sugar linked to C1
	[[nodiscard]] std::size_t get_link_nr() const
	{
		std::size_t result = 0;
		if (m_link)
			result = m_link.get_property_int("auth_seq_id");
		return result;
	}

	/// \brief Construct an atom based on the info in @a atom_info and add it to this sugar
	atom add_atom(row_initializer atom_info);

  private:
	branch *m_branch;
	atom m_link;
};

/**
 * @brief A branch is a list of sugars
 *
 * A list is how it is stored, but a branch is like a branch in a tree,
 * with potentially lots of sub branches. Each sugar is linked to a sugar
 * up in the branch with its (almost always) C1 atom.
 *
 */
class branch : public std::vector<sugar>
{
  public:
	/// \brief constructor
	branch(structure &structure, std::string asym_id, std::string entity_id);

	branch(const branch &) = delete;
	branch &operator=(const branch &) = delete;

	/** @cond */
	branch(branch &&) = default;
	branch &operator=(branch &&) = default;
	/** @endcond */

	/// \brief Update the link atoms in all sugars in this branch
	void link_atoms();

	/// \brief Return the name of the branch
	[[nodiscard]] std::string name() const;

	/// \brief Return the weight of the branch based on the formulae of the sugars
	[[nodiscard]] float weight() const;

	[[nodiscard]] std::string get_asym_id() const { return m_asym_id; }     ///< Return the asym_id
	[[nodiscard]] std::string get_entity_id() const { return m_entity_id; } ///< Return the entity_id

	[[nodiscard]] structure &get_structure() { return *m_structure; }       ///< Return the structure
	[[nodiscard]] structure &get_structure() const { return *m_structure; } ///< Return the structure

	/// \brief Return a reference to the sugar with number @a num
	[[nodiscard]] sugar &get_sugar_by_num(int nr);

	/// \brief Return a const reference to the sugar with number @a num
	[[nodiscard]] const sugar &get_sugar_by_num(int nr) const
	{
		return const_cast<branch *>(this)->get_sugar_by_num(nr);
	}

	/// \brief Construct a new sugar with compound ID @a compound_id in this branch
	/// and return a reference to the newly created sugar. Use this to create a first
	/// sugar in a branch.
	sugar &construct_sugar(const std::string &compound_id);

	/// \brief Construct a new sugar with compound ID @a compound_id in this branch
	/// and return a reference to the newly created sugar. The newly created sugar
	/// will be connected to an already created sugar in the branch using the
	/// information in @a atom_id, @a linked_sugar_nr and @a linked_atom_id
	sugar &construct_sugar(const std::string &compound_id, const std::string &atom_id,
		int linked_sugar_nr, const std::string &linked_atom_id);

  private:
	friend sugar;

	[[nodiscard]] std::string name(const sugar &s) const;

	structure *m_structure;
	std::string m_asym_id, m_entity_id;
};

/** @brief Enumeration for controlling atom selection based on occupancy. */
enum class occupancy_policy
{
	/** @brief Include all atoms regardless of their occupancy factor. */
	ALL = 0,

	/** @brief Select only alternate atoms with the maximum occupancy factor.
	 * If multiple atoms have the same maximum occupancy, choose the one with the minimum B-factor.
	 * If multiple atoms share both the maximum occupancy and the minimum B-factor, select the first encountered atom.
	 */
	MAX = 1,

	/** @brief Select only alternate atoms with the minimum occupancy factor.
	 * Similar to MAX, if multiple atoms have the same minimum occupancy, choose the one with the minimum B-factor.
	 * If multiple atoms share both the minimum occupancy and the minimum B-factor, select the first encountered atom.
	 */
	MIN = 2,

	/** @brief Exclude all atoms with an occupancy factor greater than zero. */
	UNOCCUPIED = 3
};

struct structure_open_options
{
	bool skip_hydrogen = false;                              ///< Do not include hydrogen atoms in the structure object
	bool skip_hetatom = false;                               ///< Do not include HET atoms in the structure object
	bool skip_water = false;                                 ///< Do not include water atoms in the structure object
	occupancy_policy occupancy_mode = occupancy_policy::ALL; ///< By default, the occupancy policy is set to occupancy_policy::ALL
	std::vector<std::string> asyms;                          ///< The asyms to load, if empty load all
	std::optional<float> min_b_factor;                       ///< Only load atoms with at least this b_factor
	std::optional<float> max_b_factor;                       ///< Only load atoms with at most this b_factor
};

// --------------------------------------------------------------------

/**
 * @brief A structure is the combination of polymers, ligand and sugar branches found
 * in the mmCIF file. This will always contain one model, the first model is taken
 * if not otherwise specified.
 *
 */
class structure
{
  public:
	/// \brief Read the structure from cif::file @a p
	structure(file &p, std::size_t modelNr = 1, structure_open_options options = {});

	/// \brief Load the structure from already parsed mmCIF data in @a db
	structure(datablock &db, std::size_t modelNr = 1, structure_open_options options = {});

	/** @cond */
	structure(structure &&s) = default;
	/** @endcond */

	// structures cannot be copied.

	structure(const structure &) = delete;
	structure &operator=(const structure &) = delete;
	~structure() = default;

	/// \brief Return the model number
	[[nodiscard]] std::size_t get_model_nr() const { return m_model_nr; }

	/// \brief Return a list of all the atoms in this structure
	[[nodiscard]] const std::vector<atom> &atoms() const { return m_atoms; }

	[[nodiscard]] EntityType get_entity_type_for_entity_id(const std::string entityID) const; ///< Return the entity type for the entity with id @a entity_id
	[[nodiscard]] EntityType get_entity_type_for_asym_id(const std::string asymID) const;     ///< Return the entity type for the asym with id @a asym_id

	[[nodiscard]] const std::list<polymer> &polymers() const { return m_polymers; } ///< Return the list of polymers
	[[nodiscard]] std::list<polymer> &polymers() { return m_polymers; }             ///< Return the list of polymers

	[[nodiscard]] polymer &get_polymer_by_asym_id(const std::string &asymID);            ///< Return the polymer having asym ID @a asymID
	[[nodiscard]] const polymer &get_polymer_by_asym_id(const std::string &asymID) const ///< Return the polymer having asym ID @a asymID
	{
		return const_cast<structure *>(this)->get_polymer_by_asym_id(asymID);
	}

	[[nodiscard]] const std::list<branch> &branches() const { return m_branches; } ///< Return the list of all branches
	[[nodiscard]] std::list<branch> &branches() { return m_branches; }             ///< Return the list of all branches

	[[nodiscard]] branch &get_branch_by_asym_id(const std::string &asymID);             ///< Return the branch having asym ID @a asymID
	[[nodiscard]] const branch &get_branch_by_asym_id(const std::string &asymID) const; ///< Return the branch having asym ID @a asymID

	[[nodiscard]] const std::vector<residue> &non_polymers() const { return m_non_polymers; } ///< Return the list of non-polymers, actually the list of ligands

	[[nodiscard]] bool has_atom_id(const std::string &id) const;    ///< Return true if an atom with ID @a id exists in this structure
	[[nodiscard]] atom get_atom_by_id(const std::string &id) const; ///< Return the atom with ID @a id

	/// \brief Return the atom identified by the label_ values specified
	[[nodiscard]] atom get_atom_by_label(const std::string &atomID, const std::string &asymID,
		const std::string &compID, int seqID, const std::string &altID = "");

	/// \brief Return the atom closest to point \a p
	[[nodiscard]] atom get_atom_by_position(point p) const;

	/// \brief Return the atom closest to point \a p with atom type \a type in a residue of type \a res_type
	[[nodiscard]] atom get_atom_by_position_and_type(point p, std::string_view type, std::string_view res_type) const;

	/// \brief Create a non-poly residue based on atoms already present in this structure.
	residue &create_residue(const std::vector<atom> &atoms);

	/// \brief Get a non-poly residue for an asym with id \a asymID
	[[nodiscard]] residue &get_residue(const std::string &asymID)
	{
		return get_residue(asymID, 0, "");
	}

	/// \brief Get a non-poly residue for an asym with id \a asymID
	[[nodiscard]] const residue &get_residue(const std::string &asymID) const
	{
		return get_residue(asymID, 0, "");
	}

	/// \brief Get a residue for an asym with id \a asymID seq id \a seqID and authSeqID \a authSeqID
	[[nodiscard]] residue &get_residue(const std::string &asymID, int seqID, const std::string &authSeqID);

	/// \brief Get a the single residue for an asym with id \a asymID seq id \a seqID and authSeqID \a authSeqID
	[[nodiscard]] const residue &get_residue(const std::string &asymID, int seqID, const std::string &authSeqID) const
	{
		return const_cast<structure *>(this)->get_residue(asymID, seqID, authSeqID);
	}

	/// \brief Get a residue for an asym with id \a asymID, compound id \a compID, seq id \a seqID and authSeqID \a authSeqID
	[[nodiscard]] residue &get_residue(const std::string &asymID, const std::string &compID, int seqID, const std::string &authSeqID);

	/// \brief Get a residue for an asym with id \a asymID, compound id \a compID, seq id \a seqID and authSeqID \a authSeqID
	[[nodiscard]] const residue &get_residue(const std::string &asymID, const std::string &compID, int seqID, const std::string &authSeqID) const
	{
		return const_cast<structure *>(this)->get_residue(asymID, compID, seqID, authSeqID);
	}

	/// \brief Get a the residue for atom \a atom
	[[nodiscard]] residue &get_residue(const atom &atom)
	{
		return get_residue(atom.get_label_asym_id(), atom.get_label_comp_id(), atom.get_label_seq_id(), atom.get_auth_seq_id());
	}

	/// \brief Get a the residue for atom \a atom
	[[nodiscard]] const residue &get_residue(const atom &atom) const
	{
		return get_residue(atom.get_label_asym_id(), atom.get_label_comp_id(), atom.get_label_seq_id(), atom.get_auth_seq_id());
	}

	// Actions. Originally a lot more actions were expected here

	/// \brief Remove atom @a a
	void remove_atom(atom &a)
	{
		remove_atom(a, true);
	}

	void swap_atoms(atom a1, atom a2); ///< swap the labels for these atoms
	void move_atom(atom a, point p);   ///< move atom to a new location

	/**
	 * @brief Change residue @a res to a new compound ID optionally
	 * remapping atoms.
	 *
	 * A new chem_comp entry as well as an entity is created if needed and
	 * if the list of @a remappedAtoms is not empty it is used to remap.
	 *
	 * The array in @a remappedAtoms contains tuples of strings, both
	 * strings contain an atom_id. The first is the one in the current
	 * residue and the second is the atom_id that should be used instead.
	 * If the second string is empty, the atom is removed from the residue.
	 *
	 * @param res
	 * @param newcompound
	 * @param remappedAtoms
	 */
	void change_residue(residue &res, const std::string &newcompound,
		const std::vector<std::tuple<std::string, std::string>> &remappedAtoms);

	/// \brief Remove a residue, can be monomer or nonpoly
	///
	/// \param asym_id     The asym ID
	/// \param seq_id      The sequence ID
	/// \param auth_seq_id The auth sequence ID
	void remove_residue(const std::string &asym_id, int seq_id, const std::string &auth_seq_id);

	/// \brief Create a new non-polymer entity, returns new ID
	/// \param mon_id	The mon_id for the new nonpoly, must be an existing and known compound from CCD
	/// \return			The ID of the created entity
	std::string create_non_poly_entity(const std::string &mon_id);

	/// \brief Create a new NonPolymer struct_asym with atoms constructed from \a atoms, returns asym_id.
	/// This method assumes you are copying data from one cif file to another.
	///
	/// \param entity_id	The entity ID of the new nonpoly
	/// \param atoms		The array of atom_site rows containing the data.
	/// \return				The newly create asym ID
	std::string create_non_poly(const std::string &entity_id, const std::vector<atom> &atoms);

	/// \brief Create a new NonPolymer struct_asym with atoms constructed from info in \a atom_info, returns asym_id.
	/// This method creates new atom records filled with info from the info.
	///
	/// \param entity_id	The entity ID of the new nonpoly
	/// \param atoms		The array of sets of item data containing the data for the atoms.
	/// \return				The newly create asym ID
	std::string create_non_poly(const std::string &entity_id, std::vector<row_initializer> atoms);

	/// \brief Create a new NonPolymer struct_asym for a compound of type \a compound_id, returns asym_id.
	/// This method creates new atom records filled with info from the CCD compound info.
	///
	/// \param compound_id	 The compound ID of the new nonpoly
	/// \param skip_hydrogen Do not create hydrogen atoms when true
	/// \return				 The newly create asym ID
	std::string create_non_poly(const std::string &compound_id, bool skip_hydrogen);

	/// \brief Create a new water with atom constructed from info in \a atom_info
	/// This method creates a new atom record filled with info from the info.
	///
	/// \param atom			The set of item data containing the data for the atoms.
	void create_water(row_initializer atom);

	/// \brief Create a link, a struct_conn record for two atoms.
	///
	/// \param a1			Atom 1
	/// \param a2			Atom 2
	/// \param link_type	The struct_conn_type ID for the link
	/// \param role			The pdbx_role field value
	/// \return 			The ID of the struct_conn record created

	std::string create_link(atom a1, atom a2, const std::string &link_type, const std::string &role);

	/// \brief Create a new and empty (sugar) branch
	branch &create_branch();

	// /// \brief Create a new (sugar) branch with one first NAG containing atoms constructed from \a atoms
	// branch &create_branch(std::vector<row_initializer> atoms);

	// /// \brief Extend an existing (sugar) branch identified by \a asymID with one sugar containing atoms constructed from \a atom_info
	// ///
	// /// \param asym_id      The asym id of the branch to extend
	// /// \param atom_info    Array containing the info for the atoms to construct for the new sugar
	// /// \param link_sugar   The sugar to link to, note: this is the sugar number (1 based)
	// /// \param link_atom    The atom id of the atom linked in the sugar
	// branch &extend_branch(const std::string &asym_id, std::vector<row_initializer> atom_info,
	// 	int link_sugar, const std::string &link_atom);

	/// \brief Remove \a branch
	void remove_branch(branch &branch);

	/// \brief Remove residue \a res
	///
	/// \param res         The residue to remove
	void remove_residue(residue &res);

	/// \brief Translate the coordinates of all atoms in the structure by \a t
	void translate(point t);

	/// \brief Rotate the coordinates of all atoms in the structure by \a q
	void rotate(quaternion t);

	/// \brief Translate and rotate the coordinates of all atoms in the structure by \a t and \a q
	void translate_and_rotate(point t, quaternion q);

	/// \brief Translate, rotate and translate again the coordinates of all atoms in the structure by \a t1 , \a q and \a t2
	void translate_rotate_and_translate(point t1, quaternion q, point t2);

	/// \brief Remove all categories that have no rows left
	void cleanup_empty_categories();

	/// \brief Direct access to underlying data
	[[nodiscard]] category &get_category(std::string_view name) const
	{
		return m_db[name];
	}

	/// \brief Direct access to underlying data
	[[nodiscard]] datablock &get_datablock() const
	{
		return m_db;
	}

	/// \brief Check if all atoms are part of either a polymer, a branch or one of the non-polymer residues
	void validate_atoms() const;

	/// \brief emplace a newly created atom using @a args
	template <typename... Args>
	atom &emplace_atom(Args &...args)
	{
		return emplace_atom(atom{ std::forward<Args>(args)... });
	}

	/// \brief emplace the moved atom @a atom
	atom &emplace_atom(atom &&atom);

	/// \brief Reorder atom_site atoms based on 'natural' ordering
	void reorder_atoms();

  private:
	friend polymer;
	friend residue;

	void load_atoms_for_model(structure_open_options options);

	std::string insert_compound(const std::string &compoundID, bool is_entity);

	std::string create_entity_for_branch(branch &branch);

	void load_data();

	void remove_atom(atom &a, bool removeFromResidue);
	void remove_sugar(sugar &sugar);

	datablock &m_db;
	std::size_t m_model_nr;
	std::vector<atom> m_atoms;
	std::vector<std::size_t> m_atom_index;
	std::list<polymer> m_polymers;
	std::list<branch> m_branches;
	std::vector<residue> m_non_polymers;
};

} // namespace cif::mm
