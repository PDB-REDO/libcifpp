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

#include <numeric>

#if __cpp_lib_format
#include <format>
#endif

namespace cif::mm
{

class atom;
class residue;
class monomer;
class polymer;
class structure;

// --------------------------------------------------------------------

class atom
{
  private:
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

		void prefetch();

		int compare(const atom_impl &b) const;

		// bool getAnisoU(float anisou[6]) const;

		int get_charge() const;

		void moveTo(const point &p);

		// const compound *compound() const;

		std::string get_property(std::string_view name) const;
		int get_property_int(std::string_view name) const;
		float get_property_float(std::string_view name) const;

		void set_property(const std::string_view name, const std::string &value);

		row_handle row()
		{
			return m_cat[{{"id", m_id}}];
		}

		const row_handle row() const
		{
			return m_cat[{{"id", m_id}}];
		}

		row_handle row_aniso()
		{
			auto cat = m_db.get("atom_site_anisotrop");
			return cat ? cat->operator[]({ {"id", m_id} }) : row_handle{};
		}

		const row_handle row_aniso() const
		{
			auto cat = m_db.get("atom_site_anisotrop");
			return cat ? cat->operator[]({ {"id", m_id} }) : row_handle{};
		}

		const datablock &m_db;
		const category &m_cat;
		std::string m_id;
		point m_location;
		std::string m_symop = "1_555";
	};

  public:
	atom() {}

	atom(std::shared_ptr<atom_impl> impl)
		: m_impl(impl)
	{
	}

	atom(const atom &rhs)
		: m_impl(rhs.m_impl)
	{
	}

	atom(const datablock &db, const row_handle &row)
		: atom(std::make_shared<atom_impl>(db, row["id"].as<std::string>()))
	{
	}

	// a special constructor to create symmetry copies
	atom(const atom &rhs, const point &symmmetry_location, const std::string &symmetry_operation)
		: atom(std::make_shared<atom_impl>(*rhs.m_impl, symmmetry_location, symmetry_operation))
	{
	}

	explicit operator bool() const { return (bool)m_impl; }

	// // return a copy of this atom, with data copied instead of referenced
	// atom clone() const
	// {
	// 	auto copy = std::make_shared<atom_impl>(*m_impl);
	// 	copy->mClone = true;
	// 	return atom(copy);
	// }

	atom &operator=(const atom &rhs) = default;

	// template <typename T>
	// T get_property(const std::string_view name) const;

	std::string get_property(std::string_view name) const
	{
		if (not m_impl)
			throw std::logic_error("Error trying to fetch a property from an uninitialized atom");
		return m_impl->get_property(name);
	}

	int get_property_int(std::string_view name) const
	{
		if (not m_impl)
			throw std::logic_error("Error trying to fetch a property from an uninitialized atom");
		return m_impl->get_property_int(name);
	}

	float get_property_float(std::string_view name) const
	{
		if (not m_impl)
			throw std::logic_error("Error trying to fetch a property from an uninitialized atom");
		return m_impl->get_property_float(name);
	}

	void set_property(const std::string_view name, const std::string &value)
	{
		if (not m_impl)
			throw std::logic_error("Error trying to modify an uninitialized atom");
		m_impl->set_property(name, value);
	}

	template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	void set_property(const std::string_view name, const T &value)
	{
		set_property(name, std::to_string(value));
	}

	const std::string &id() const { return impl().m_id; }

	cif::atom_type get_type() const { return atom_type_traits(get_property("type_symbol")).type(); }

	point get_location() const { return impl().m_location; }
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

	// for direct access to underlying data, be careful!
	const row_handle get_row() const { return impl().row(); }
	const row_handle get_row_aniso() const { return impl().row_aniso(); }

	bool is_symmetry_copy() const { return impl().m_symop != "1_555"; }
	std::string symmetry() const { return impl().m_symop; }

	// const compound &compound() const;

	bool is_water() const
	{
		auto comp_id = get_label_comp_id();
		return comp_id == "HOH" or comp_id == "H2O" or comp_id == "WAT";
	}

	int get_charge() const { return impl().get_charge(); }

	// float uIso() const;
	// bool getAnisoU(float anisou[6]) const { return impl().getAnisoU(anisou); }
	
	float get_occupancy() const { return get_property_float("occupancy"); }

	// specifications

	std::string get_label_asym_id() const { return get_property("label_asym_id"); }
	int get_label_seq_id() const { return get_property_int("label_seq_id"); }
	std::string get_label_atom_id() const { return get_property("label_atom_id"); }
	std::string get_label_alt_id() const { return get_property("label_alt_id"); }
	std::string get_label_comp_id() const { return get_property("label_comp_id"); }
	std::string get_label_entity_id() const { return get_property("label_entity_id"); }

	std::string get_auth_asym_id() const { return get_property("auth_asym_id"); }
	std::string get_auth_seq_id() const { return get_property("auth_seq_id"); }
	std::string get_auth_atom_id() const { return get_property("auth_atom_id"); }
	std::string get_auth_alt_id() const { return get_property("auth_alt_id"); }
	std::string get_auth_comp_id() const { return get_property("auth_comp_id"); }
	std::string get_pdb_ins_code() const { return get_property("pdbx_PDB_ins_code"); }

	bool is_alternate() const { return not get_label_alt_id().empty(); }

	// std::string labelID() const; // label_comp_id + '_' + label_asym_id + '_' + label_seq_id
	
	std::string pdb_id() const
	{
		return get_label_comp_id() + '_' + get_auth_asym_id() + '_' + get_auth_seq_id() + get_pdb_ins_code();
	}

	bool operator==(const atom &rhs) const
	{
		if (m_impl == rhs.m_impl)
			return true;

		if (not(m_impl and rhs.m_impl))
			return false;

		return &m_impl->m_db == &rhs.m_impl->m_db and m_impl->m_id == rhs.m_impl->m_id;
	}

	bool operator!=(const atom &rhs) const
	{
		return not operator==(rhs);
	}

	// // access data in compound for this atom

	// convenience routine
	bool is_back_bone() const
	{
		auto atomID = get_label_atom_id();
		return atomID == "N" or atomID == "O" or atomID == "C" or atomID == "CA";
	}

	void swap(atom &b)
	{
		std::swap(m_impl, b.m_impl);
	}

	int compare(const atom &b) const { return impl().compare(*b.m_impl); }

	bool operator<(const atom &rhs) const
	{
		return compare(rhs) < 0;
	}

	friend std::ostream &operator<<(std::ostream &os, const atom &atom);

  private:
	friend class structure;

	const atom_impl &impl() const
	{
		if (not m_impl)
			throw std::runtime_error("Uninitialized atom, not found?");
		return *m_impl;
	}

	std::shared_ptr<atom_impl> m_impl;
};

inline void swap(atom &a, atom &b)
{
	a.swap(b);
}

inline float distance(const atom &a, const atom &b)
{
	return distance(a.get_location(), b.get_location());
}

inline float distance_squared(const atom &a, const atom &b)
{
	return distance_squared(a.get_location(), b.get_location());
}

// --------------------------------------------------------------------

enum class EntityType
{
	polymer,
	NonPolymer,
	Macrolide,
	Water,
	Branched
};

// --------------------------------------------------------------------

class residue
{
  public:
	friend class structure;

	// constructor
	residue(structure &structure, const std::string &compoundID,
		const std::string &asymID, int seqID,
		const std::string &authAsymID, const std::string &authSeqID,
		const std::string &pdbInsCode)
		: m_structure(&structure)
		, m_compound_id(compoundID)
		, m_asym_id(asymID)
		, m_seq_id(seqID)
		, m_auth_asym_id(authAsymID)
		, m_auth_seq_id(authSeqID)
		, m_pdb_ins_code(pdbInsCode)
	{
	}

	residue(structure &structure, const std::vector<atom> &atoms);

	residue(const residue &rhs) = delete;
	residue &operator=(const residue &rhs) = delete;

	residue(residue &&rhs) = default;
	residue &operator=(residue &&rhs) = default;

	virtual ~residue() = default;

	std::string get_entity_id() const;

	EntityType entity_type() const;

	const std::string &get_asym_id() const { return m_asym_id; }
	int get_seq_id() const { return m_seq_id; }

	const std::string get_auth_asym_id() const { return m_auth_asym_id; }
	const std::string get_auth_seq_id() const { return m_auth_seq_id; }
	std::string get_pdb_ins_code() const { return m_pdb_ins_code; }

	const std::string &get_compound_id() const { return m_compound_id; }
	void set_compound_id(const std::string &id) { m_compound_id = id; }

	structure *get_structure() const { return m_structure; }

	// const compound &compound() const;

	std::vector<atom> &atoms()
	{
		return m_atoms;
	}

	const std::vector<atom> &atoms() const
	{
		return m_atoms;
	}

	void add_atom(atom &atom);

	/// \brief Unique atoms returns only the atoms without alternates and the first of each alternate atom id.
	std::vector<atom> unique_atoms() const;

	/// \brief The alt ID used for the unique atoms
	std::string unique_alt_id() const;

	atom get_atom_by_atom_id(const std::string &atomID) const;

	// Is this residue a single entity?
	bool is_entity() const;
	bool is_water() const { return m_compound_id == "HOH"; }
	// bool empty() const { return m_structure == nullptr; }

	bool has_alternate_atoms() const;

	/// \brief Return the list of unique alt ID's present in this residue
	std::set<std::string> get_alternate_ids() const;

	/// \brief Return the list of unique atom ID's
	std::set<std::string> get_atom_ids() const;

	/// \brief Return the list of atoms having ID \a atomID
	std::vector<atom> get_atoms_by_id(const std::string &atomID) const;

	// some routines for 3d work
	std::tuple<point, float> center_and_radius() const;

	friend std::ostream &operator<<(std::ostream &os, const residue &res);

	bool operator==(const residue &rhs) const
	{
		return this == &rhs or (m_structure == rhs.m_structure and
								   m_seq_id == rhs.m_seq_id and
								   m_asym_id == rhs.m_asym_id and
								   m_compound_id == rhs.m_compound_id and
								   m_auth_seq_id == rhs.m_auth_seq_id);
	}

  protected:
	residue() {}

	structure *m_structure = nullptr;
	std::string m_compound_id, m_asym_id;
	int m_seq_id = 0;
	std::string m_auth_asym_id, m_auth_seq_id, m_pdb_ins_code;
	std::vector<atom> m_atoms;
};

// --------------------------------------------------------------------
// a monomer models a single residue in a protein chain

class monomer : public residue
{
  public:
	//	monomer();
	monomer(const monomer &rhs) = delete;
	monomer &operator=(const monomer &rhs) = delete;

	monomer(monomer &&rhs);
	monomer &operator=(monomer &&rhs);

	monomer(const polymer &polymer, size_t index, int seqID, const std::string &authSeqID,
		const std::string &pdbInsCode, const std::string &compoundID);

	bool is_first_in_chain() const;
	bool is_last_in_chain() const;

	// convenience
	bool has_alpha() const;
	bool has_kappa() const;

	// Assuming this is really an amino acid...

	float phi() const;
	float psi() const;
	float alpha() const;
	float kappa() const;
	float tco() const;
	float omega() const;

	// torsion angles
	size_t nr_of_chis() const;
	float chi(size_t i) const;

	bool is_cis() const;

	/// \brief Returns true if the four atoms C, CA, N and O are present
	bool is_complete() const;

	/// \brief Returns true if any of the backbone atoms has an alternate
	bool has_alternate_backbone_atoms() const;

	atom CAlpha() const { return get_atom_by_atom_id("CA"); }
	atom C() const { return get_atom_by_atom_id("C"); }
	atom N() const { return get_atom_by_atom_id("N"); }
	atom O() const { return get_atom_by_atom_id("O"); }
	atom H() const { return get_atom_by_atom_id("H"); }

	bool is_bonded_to(const monomer &rhs) const
	{
		return this != &rhs and are_bonded(*this, rhs);
	}

	static bool are_bonded(const monomer &a, const monomer &b, float errorMargin = 0.5f);
	static bool is_cis(const monomer &a, const monomer &b);
	static float omega(const monomer &a, const monomer &b);

	// for LEU and VAL
	float chiral_volume() const;

	bool operator==(const monomer &rhs) const
	{
		return m_polymer == rhs.m_polymer and m_index == rhs.m_index;
	}

  private:
	const polymer *m_polymer;
	size_t m_index;
};

// --------------------------------------------------------------------

class polymer : public std::vector<monomer>
{
  public:
	polymer(structure &s, const std::string &entityID, const std::string &asymID, const std::string &auth_asym_id);

	polymer(const polymer &) = delete;
	polymer &operator=(const polymer &) = delete;

	// monomer &getBySeqID(int seqID);
	// const monomer &getBySeqID(int seqID) const;

	structure *get_structure() const { return m_structure; }

	std::string get_asym_id() const { return m_asym_id; }
	std::string get_auth_asym_id() const { return m_auth_asym_id; }	// The PDB chain ID, actually
	std::string get_entity_id() const { return m_entity_id; }

	// int Distance(const monomer &a, const monomer &b) const;

  private:
	structure *m_structure;
	std::string m_entity_id;
	std::string m_asym_id;
	std::string m_auth_asym_id;
};

// --------------------------------------------------------------------
// sugar and branch, to describe glycosylation sites

class branch;

class sugar : public residue
{
  public:
	sugar(branch &branch, const std::string &compoundID,
		const std::string &asymID, int authSeqID);

	sugar(sugar &&rhs);
	sugar &operator=(sugar &&rhs);

	int num() const {
		int result;
		auto r = std::from_chars(m_auth_seq_id.data(), m_auth_seq_id.data() + m_auth_seq_id.length(), result);
		if (r.ec != std::errc())
			throw std::runtime_error("The auth_seq_id should be a number for a sugar");
		return result;
	}
	std::string name() const;

	/// \brief Return the atom the C1 is linked to
	atom get_link() const { return m_link; }
	void set_link(atom link) { m_link = link; }

	size_t get_link_nr() const
	{
		size_t result = 0;
		if (m_link)
			result = m_link.get_property_int("auth_seq_id");
		return result;
	}

	cif::mm::atom add_atom(row_initializer atom_info);

  private:
	branch *m_branch;
	atom m_link;
};

class branch : public std::vector<sugar>
{
  public:
	branch(structure &structure, const std::string &asym_id, const std::string &entity_id);

	branch(const branch &) = delete;
	branch &operator=(const branch &) = delete;

	branch(branch &&) = default;
	branch &operator=(branch &&) = default;

	void link_atoms();

	std::string name() const;
	float weight() const;
	std::string get_asym_id() const { return m_asym_id; }
	std::string get_entity_id() const { return m_entity_id; }

	structure &get_structure() { return *m_structure; }
	structure &get_structure() const { return *m_structure; }

	sugar &get_sugar_by_num(int nr);

	const sugar &get_sugar_by_num(int nr) const
	{
		return const_cast<branch *>(this)->get_sugar_by_num(nr);
	}

	sugar &construct_sugar(const std::string &compound_id);
	sugar &construct_sugar(const std::string &compound_id, const std::string &atom_id,
		int linked_sugar_nr, const std::string &linked_atom_id);

  private:
	friend sugar;

	std::string name(const sugar &s) const;

	structure *m_structure;
	std::string m_asym_id, m_entity_id;
};

// --------------------------------------------------------------------

enum class StructureOpenOptions
{
	SkipHydrogen = 1 << 0
};

constexpr inline bool operator&(StructureOpenOptions a, StructureOpenOptions b)
{
	return static_cast<int>(a) bitand static_cast<int>(b);
}

// --------------------------------------------------------------------

class structure
{
  public:
	structure(file &p, size_t modelNr = 1, StructureOpenOptions options = {});

	structure(datablock &db, size_t modelNr = 1, StructureOpenOptions options = {});

	structure(structure &&s) = default;

	// Create a read-only clone of the current structure (for multithreaded calculations that move atoms)
	// NOTE: removed, simply create a new structure for each thread
	structure(const structure &) = delete;

	structure &operator=(const structure &) = delete;
	// Structure &operator=(Structure &&s) = default;

	~structure() = default;

	size_t get_model_nr() const { return m_model_nr; }

	const std::vector<atom> &atoms() const { return m_atoms; }
	// std::vector<atom> &atoms() { return m_atoms; }

	EntityType get_entity_type_for_entity_id(const std::string entityID) const;
	EntityType get_entity_type_for_asym_id(const std::string asymID) const;

	// std::vector<atom> waters() const;

	const std::list<polymer> &polymers() const { return m_polymers; }
	std::list<polymer> &polymers() { return m_polymers; }

	polymer &get_polymer_by_asym_id(const std::string &asymID);

	const polymer &get_polymer_by_asym_id(const std::string &asymID) const
	{
		return const_cast<structure *>(this)->get_polymer_by_asym_id(asymID);
	}

	const std::list<branch> &branches() const { return m_branches; }
	std::list<branch> &branches() { return m_branches; }

	branch &get_branch_by_asym_id(const std::string &asymID);
	const branch &get_branch_by_asym_id(const std::string &asymID) const;

	const std::vector<residue> &non_polymers() const { return m_non_polymers; }

	bool has_atom_id(const std::string &id) const;
	atom get_atom_by_id(const std::string &id) const;
	// atom getAtomByLocation(point pt, float maxDistance) const;

	atom get_atom_by_label(const std::string &atomID, const std::string &asymID,
		const std::string &compID, int seqID, const std::string &altID = "");

	// /// \brief Return the atom closest to point \a p
	atom get_atom_by_position(point p) const;

	/// \brief Return the atom closest to point \a p with atom type \a type in a residue of type \a res_type
	atom get_atom_by_position_and_type(point p, std::string_view type, std::string_view res_type) const;

	/// \brief Create a non-poly residue based on atoms already present in this structure.
	residue &create_residue(const std::vector<atom> &atoms);

	/// \brief Get a non-poly residue for an asym with id \a asymID
	residue &get_residue(const std::string &asymID)
	{
		return get_residue(asymID, 0, "");
	}

	/// \brief Get a non-poly residue for an asym with id \a asymID
	const residue &get_residue(const std::string &asymID) const
	{
		return get_residue(asymID, 0, "");
	}

	/// \brief Get a residue for an asym with id \a asymID seq id \a seqID and authSeqID \a authSeqID
	residue &get_residue(const std::string &asymID, int seqID, const std::string &authSeqID);

	/// \brief Get a the single residue for an asym with id \a asymID seq id \a seqID and authSeqID \a authSeqID
	const residue &get_residue(const std::string &asymID, int seqID, const std::string &authSeqID) const
	{
		return const_cast<structure *>(this)->get_residue(asymID, seqID, authSeqID);
	}

	/// \brief Get a residue for an asym with id \a asymID, compound id \a compID, seq id \a seqID and authSeqID \a authSeqID
	residue &get_residue(const std::string &asymID, const std::string &compID, int seqID, const std::string &authSeqID);

	/// \brief Get a residue for an asym with id \a asymID, compound id \a compID, seq id \a seqID and authSeqID \a authSeqID
	const residue &get_residue(const std::string &asymID, const std::string &compID, int seqID, const std::string &authSeqID) const
	{
		return const_cast<structure *>(this)->get_residue(asymID, compID, seqID, authSeqID);
	}

	/// \brief Get a the residue for atom \a atom
	residue &get_residue(const atom &atom)
	{
		return get_residue(atom.get_label_asym_id(), atom.get_label_comp_id(), atom.get_label_seq_id(), atom.get_auth_seq_id());
	}

	/// \brief Get a the residue for atom \a atom
	const residue &get_residue(const atom &atom) const
	{
		return get_residue(atom.get_label_asym_id(), atom.get_label_comp_id(), atom.get_label_seq_id(), atom.get_auth_seq_id());
	}

	// Actions
	void remove_atom(atom &a)
	{
		remove_atom(a, true);
	}

	void swap_atoms(atom a1, atom a2); // swap the labels for these atoms
	void move_atom(atom a, point p);   // move atom to a new location
	void change_residue(residue &res, const std::string &newcompound,
		const std::vector<std::tuple<std::string, std::string>> &remappedAtoms);

	/// \brief Remove a residue, can be monomer or nonpoly
	///
	/// \param asym_id     The asym ID
	/// \param seq_id      The sequence ID
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

	/// \brief Create a new water with atom constructed from info in \a atom_info
	/// This method creates a new atom record filled with info from the info.
	///
	/// \param atom			The set of item data containing the data for the atoms.
	void create_water(row_initializer atom);

	/// \brief Create a new and empty (sugar) branch
	branch &create_branch();

	/// \brief Create a new (sugar) branch with one first NAG containing atoms constructed from \a atoms
	branch &create_branch(std::vector<row_initializer> atoms);

	/// \brief Extend an existing (sugar) branch identified by \a asymID with one sugar containing atoms constructed from \a atom_info
	///
	/// \param asym_id      The asym id of the branch to extend
	/// \param atom_info    Array containing the info for the atoms to construct for the new sugar
	/// \param link_sugar   The sugar to link to, note: this is the sugar number (1 based)
	/// \param link_atom    The atom id of the atom linked in the sugar
	branch &extend_branch(const std::string &asym_id, std::vector<row_initializer> atom_info,
		int link_sugar, const std::string &link_atom);

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

	void cleanup_empty_categories();

	/// \brief Direct access to underlying data
	category &get_category(std::string_view name) const
	{
		return m_db[name];
	}

	datablock &get_datablock() const
	{
		return m_db;
	}

	void validate_atoms() const;

	// TODO: make this protected?

	void load_atoms_for_model(StructureOpenOptions options);

	template <typename... Args>
	atom &emplace_atom(Args&... args)
	{
		return emplace_atom(atom{ std::forward<Args>(args)... });
	}

	atom &emplace_atom(atom &&atom);

  private:
	friend polymer;
	friend residue;

	std::string insert_compound(const std::string &compoundID, bool is_entity);

	std::string create_entity_for_branch(branch &branch);

	void load_data();

	void remove_atom(atom &a, bool removeFromResidue);
	void remove_sugar(sugar &sugar);

	datablock &m_db;
	size_t m_model_nr;
	std::vector<atom> m_atoms;
	std::vector<size_t> m_atom_index;
	std::list<polymer> m_polymers;
	std::list<branch> m_branches;
	std::vector<residue> m_non_polymers;
};

// --------------------------------------------------------------------

/// \brief Reconstruct all missing categories for an assumed PDBx file.
/// Some people believe that simply dumping some atom records is enough.
/// \param db The cif::datablock that hopefully contains some valid data
void reconstruct_pdbx(datablock &db);

} // namespace cif::mm
