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

#include <numeric>

#if __cpp_lib_format
#include <format>
#endif

/*
    To modify a structure, you will have to use actions.

    The currently supported actions are:

//	- Move atom to new location
    - Remove atom
//	- Add new atom that was formerly missing
//	- Add alternate residue
    -

*/

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
		atom_impl(datablock &db, std::string_view id, row_handle row)
			: m_db(db)
			, m_id(id)
			, m_row(row)
		{
			prefetch();
		}

		// constructor for a symmetry copy of an atom
		atom_impl(const atom_impl &impl, const point &loc, const std::string &sym_op);

		atom_impl(const atom_impl &i) = default;

		void prefetch();

		int compare(const atom_impl &b) const;

		// bool getAnisoU(float anisou[6]) const;

		// int charge() const;

		void moveTo(const point &p)
		{
			if (m_symop != "1_555")
				throw std::runtime_error("Moving symmetry copy");

#if __cpp_lib_format
			m_row.assign("Cartn_x", std::format("{:.3f}", p.getX()), false, false);
			m_row.assign("Cartn_y", std::format("{:.3f}", p.getY()), false, false);
			m_row.assign("Cartn_z", std::format("{:.3f}", p.getZ()), false, false);
#else
			m_row.assign("Cartn_x", format("%.3f", p.m_x).str(), false, false);
			m_row.assign("Cartn_y", format("%.3f", p.m_y).str(), false, false);
			m_row.assign("Cartn_z", format("%.3f", p.m_z).str(), false, false);
#endif
			m_location = p;
		}

		// const compound *compound() const;

		std::string get_property(std::string_view name) const
		{
			return m_row[name].as<std::string>();
			// 	for (auto &&[tag, ref] : mCachedRefs)
			// 	{
			// 		if (tag == name)
			// 			return ref.as<std::string>();
			// 	}

			// 	mCachedRefs.emplace_back(name, const_cast<Row &>(mRow)[name]);
			// return std::get<1>(mCachedRefs.back()).as<std::string>();
		}

		int get_property_int(std::string_view name) const
		{
			int result = 0;
			if (not m_row[name].empty())
			{
				auto s = get_property(name);

				std::from_chars_result r = std::from_chars(s.data(), s.data() + s.length(), result);
				if (r.ec != std::errc() and VERBOSE > 0)
					std::cerr << "Error converting " << s << " to number for property " << name << std::endl;
			}
			return result;
		}

		void set_property(const std::string_view name, const std::string &value)
		{
			m_row.assign(name, value, true, true);
		}

		// const datablock &m_db;
		// std::string mID;
		// atom_type mType;

		// std::string mAtomID;
		// std::string mCompID;
		// std::string m_asym_id;
		// int m_seq_id;
		// std::string mAltID;
		// std::string m_auth_seq_id;

		// point mLocation;
		// row_handle mRow;

		// // mutable std::vector<std::tuple<std::string, detail::ItemReference>> mCachedRefs;

		// mutable const compound *mcompound = nullptr;

		// bool mSymmetryCopy = false;
		// bool mClone = false;

		const datablock &m_db;
		std::string m_id;
		row_handle m_row;
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

	atom(datablock &db, row_handle &row)
		: atom(std::make_shared<atom_impl>(db, row["id"].as<std::string>(), row))
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
	// AtomType type() const { return impl().mType; }

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

	// // for direct access to underlying data, be careful!
	// const row_handle getRow() const { return impl().mRow; }
	// const row_handle getRowAniso() const;

	// bool isSymmetryCopy() const { return impl().mSymmetryCopy; }
	// std::string symmetry() const { return impl().mSymmetryOperator; }

	// const compound &compound() const;
	// bool isWater() const { return impl().mCompID == "HOH" or impl().mCompID == "H2O" or impl().mCompID == "WAT"; }
	// int charge() const;

	// float uIso() const;
	// bool getAnisoU(float anisou[6]) const { return impl().getAnisoU(anisou); }
	// float occupancy() const;

	// specifications

	std::string get_label_asym_id() const { return get_property("label_asym_id"); }
	int get_label_seq_id() const { return get_property_int("label_seq_id"); }
	std::string get_label_atom_id() const { return get_property("label_atom_id"); }
	std::string get_label_alt_id() const { return get_property("label_alt_id"); }
	std::string get_label_comp_id() const { return get_property("label_comp_id"); }
	std::string get_label_entity_id() const { return get_property("label_entity_id"); }

	std::string get_auth_asym_id() const { return get_property("auth_asym_id"); }
	std::string get_auth_seq_id() const { return get_property("auth_seq_id"); }
	std::string get_pdb_ins_code() const { return get_property("pdbx_PDB_ins_code"); }

	// const std::string &labelAtomID() const { return impl().mAtomID; }
	// const std::string &get_label_comp_id() const { return impl().mCompID; }
	// const std::string &get_label_asym_id() const { return impl().m_asym_id; }
	// std::string labelEntityID() const;
	// int get_label_seq_id() const { return impl().m_seq_id; }
	// const std::string &labelAltID() const { return impl().mAltID; }

	bool is_alternate() const { return not get_label_alt_id().empty(); }

	// std::string authAtomID() const;
	// std::string authCompID() const;
	// std::string authAsymID() const;
	// const std::string &authSeqID() const { return impl().m_auth_seq_id; }
	// std::string pdbxAuthInsCode() const;
	// std::string pdbxAuthAltID() const;

	// std::string labelID() const; // label_comp_id + '_' + label_asym_id + '_' + label_seq_id
	// std::string pdbID() const;   // auth_comp_id + '_' + auth_asym_id + '_' + auth_seq_id + pdbx_PDB_ins_code

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

	// // convenience routine
	// bool isBackBone() const
	// {
	// 	auto atomID = labelAtomID();
	// 	return atomID == "N" or atomID == "O" or atomID == "C" or atomID == "CA";
	// }

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

	// /// \brief Synchronize data with underlying cif data
	// void sync()
	// {
	// 	if (m_impl)
	// 		m_impl->prefetch();
	// }

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

// template <>
// inline std::string atom::get_property<std::string>(const std::string_view name) const
// {
// 	return get_property(name);
// }

// template <>
// inline int atom::get_property<int>(const std::string_view name) const
// {
// 	auto v = impl().get_property(name);
// 	return v.empty() ? 0 : stoi(v);
// }

// template <>
// inline float atom::get_property<float>(const std::string_view name) const
// {
// 	return stof(impl().get_property(name));
// }

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
	residue(const structure &structure, const std::string &compoundID,
		const std::string &asymID, int seqID, const std::string &authSeqID)
		: m_structure(&structure)
		, m_compound_id(compoundID)
		, m_asym_id(asymID)
		, m_seq_id(seqID)
		, m_auth_seq_id(authSeqID)
	{
	}

	residue(const residue &rhs) = delete;
	residue &operator=(const residue &rhs) = delete;

	residue(residue &&rhs) = default;
	residue &operator=(residue &&rhs) = default;

	virtual ~residue() = default;

	const std::string &get_asym_id() const { return m_asym_id; }
	int get_seq_id() const { return m_seq_id; }

	const std::string get_auth_asym_id() const
	{
		return m_atoms.empty() ? m_asym_id : m_atoms.front().get_auth_asym_id();
	}
	const std::string get_auth_seq_id() const { return m_auth_seq_id; }

	const std::string &get_compound_id() const { return m_compound_id; }
	void set_compound_id(const std::string &id) { m_compound_id = id; }

	const structure *get_structure() const { return m_structure; }

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

	// /// \brief Unique atoms returns only the atoms without alternates and the first of each alternate atom id.
	// std::vector<atom> unique_atoms() const;

	// /// \brief The alt ID used for the unique atoms
	// std::string unique_alt_id() const;

	atom get_atom_by_atom_id(const std::string &atomID) const;

	// const std::string &asymID() const { return m_asym_id; }
	// int seqID() const { return m_seq_id; }
	std::string get_entity_id() const;

	EntityType entity_type() const;

	// std::string authAsymID() const;
	// std::string authSeqID() const;
	// std::string authInsCode() const;

	// // return a human readable PDB-like auth id (chain+seqnr+iCode)
	// std::string authID() const;

	// // similar for mmCIF space
	// std::string labelID() const;

	// Is this residue a single entity?
	bool is_entity() const;

	// bool isWater() const { return m_compound_id == "HOH"; }

	// bool empty() const { return m_structure == nullptr; }

	// bool hasAlternateAtoms() const;

	// /// \brief Return the list of unique alt ID's present in this residue
	// std::set<std::string> getAlternateIDs() const;

	// /// \brief Return the list of unique atom ID's
	// std::set<std::string> getAtomIDs() const;

	// /// \brief Return the list of atoms having ID \a atomID
	// std::vector<atom> getAtomsByID(const std::string &atomID) const;

	// // some routines for 3d work
	// std::tuple<point, float> centerAndRadius() const;

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

	const structure *m_structure = nullptr;
	std::string m_compound_id, m_asym_id;
	int m_seq_id = 0;
	std::string m_auth_seq_id;
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
		const std::string &compoundID);

	// bool is_first_in_chain() const;
	// bool is_last_in_chain() const;

	// // convenience
	// bool has_alpha() const;
	// bool has_kappa() const;

	// // Assuming this is really an amino acid...

	// float phi() const;
	// float psi() const;
	// float alpha() const;
	// float kappa() const;
	// float tco() const;
	// float omega() const;

	// // torsion angles
	// size_t nrOfChis() const;
	// float chi(size_t i) const;

	// bool isCis() const;

	// /// \brief Returns true if the four atoms C, CA, N and O are present
	// bool isComplete() const;

	// /// \brief Returns true if any of the backbone atoms has an alternate
	// bool hasAlternateBackboneAtoms() const;

	// atom CAlpha() const { return get_atom_by_atom_id("CA"); }
	// atom C() const { return get_atom_by_atom_id("C"); }
	// atom N() const { return get_atom_by_atom_id("N"); }
	// atom O() const { return get_atom_by_atom_id("O"); }
	// atom H() const { return get_atom_by_atom_id("H"); }

	// bool isBondedTo(const monomer &rhs) const
	// {
	// 	return this != &rhs and areBonded(*this, rhs);
	// }

	// static bool areBonded(const monomer &a, const monomer &b, float errorMargin = 0.5f);
	// static bool isCis(const monomer &a, const monomer &b);
	// static float omega(const monomer &a, const monomer &b);

	// // for LEU and VAL
	// float chiralVolume() const;

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
	polymer(const structure &s, const std::string &entityID, const std::string &asymID);

	polymer(const polymer &) = delete;
	polymer &operator=(const polymer &) = delete;

	// monomer &getBySeqID(int seqID);
	// const monomer &getBySeqID(int seqID) const;

	const structure *get_structure() const { return m_structure; }

	std::string get_asym_id() const { return m_asym_id; }
	std::string get_entity_id() const { return m_entity_id; }

	// std::string chainID() const;

	// int Distance(const monomer &a, const monomer &b) const;

  private:
	const structure *m_structure;
	std::string m_entity_id;
	std::string m_asym_id;
};

// --------------------------------------------------------------------
// sugar and branch, to describe glycosylation sites

class branch;

class sugar : public residue
{
  public:
	sugar(const branch &branch, const std::string &compoundID,
		const std::string &asymID, int authSeqID);

	sugar(sugar &&rhs);
	sugar &operator=(sugar &&rhs);

	int num() const { return std::stoi(m_auth_seq_id); }
	std::string name() const;

	/// \brief Return the atom the C1 is linked to
	atom get_link() const { return m_link; }
	void set_link(atom link) { m_link = link; }

	size_t get_link_nr() const
	{
		return m_link ? std::stoi(m_link.get_auth_seq_id()) : 0;
	}

  private:
	const branch *m_branch;
	atom m_link;
};

class branch : public std::vector<sugar>
{
  public:
	branch(structure &structure, const std::string &asymID);

	void link_atoms();

	std::string name() const;
	float weight() const;
	std::string get_asym_id() const { return m_asym_id; }

	structure &get_structure() { return *m_structure; }
	const structure &get_structure() const { return *m_structure; }

	sugar &getSugarByNum(int nr);
	const sugar &getSugarByNum(int nr) const;

  private:
	friend sugar;

	std::string name(const sugar &s) const;

	structure *m_structure;
	std::string m_asym_id;
};

// // --------------------------------------------------------------------
// // file is a reference to the data stored in e.g. the cif file.
// // This object is not copyable.

// class File : public file
// {
//   public:
// 	File() {}

// 	// File(const std::filesystem::path &path)
// 	// {
// 	// 	load(path);
// 	// }

// 	// File(const char *data, size_t length)
// 	// {
// 	// 	load(data, length);
// 	// }

// 	File(const File &) = delete;
// 	File &operator=(const File &) = delete;

// 	// void load(const std::filesystem::path &p) override;
// 	// void save(const std::filesystem::path &p) override;

// 	// using file::load;
// 	// using file::save;

// 	datablock &data() { return front(); }
// };

// --------------------------------------------------------------------

enum class StructureOpenOptions
{
	SkipHydrogen = 1 << 0
};

inline bool operator&(StructureOpenOptions a, StructureOpenOptions b)
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

	const std::vector<atom> &atoms() const { return m_atoms; }
	// std::vector<atom> &atoms() { return m_atoms; }

	EntityType get_entity_type_for_entity_id(const std::string entityID) const;
	EntityType get_entity_type_for_asym_id(const std::string asymID) const;

	// std::vector<atom> waters() const;

	const std::list<polymer> &polymers() const { return m_polymers; }
	std::list<polymer> &polymers() { return m_polymers; }

	// polymer &getPolymerByAsymID(const std::string &asymID);

	// const polymer &getPolymerByAsymID(const std::string &asymID) const
	// {
	// 	return const_cast<structure *>(this)->getPolymerByAsymID(asymID);
	// }

	const std::list<branch> &branches() const { return m_branches; }
	std::list<branch> &branches() { return m_branches; }

	branch &get_branch_by_asym_id(const std::string &asymID);
	const branch &get_branch_by_asym_id(const std::string &asymID) const;

	const std::vector<residue> &non_polymers() const { return m_non_polymers; }

	atom get_atom_by_id(const std::string &id) const;
	// atom getAtomByLocation(point pt, float maxDistance) const;

	// atom getAtomByLabel(const std::string &atomID, const std::string &asymID,
	// 	const std::string &compID, int seqID, const std::string &altID = "");

	// /// \brief Return the atom closest to point \a p
	// atom getAtomByPosition(point p) const;

	// /// \brief Return the atom closest to point \a p with atom type \a type in a residue of type \a res_type
	// atom getAtomByPositionAndType(point p, std::string_view type, std::string_view res_type) const;

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
	void remove_residue(const std::string &asym_id, int seq_id, const std::string &auth_seq_id)
	{
		remove_residue(get_residue(asym_id, seq_id, auth_seq_id));
	}

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
	std::string create_non_poly(const std::string &entity_id, std::vector<std::vector<item>> &atom_info);

	/// \brief Create a new (sugar) branch with one first NAG containing atoms constructed from \a nag_atom_info
	branch &create_branch(std::vector<std::vector<item>> &nag_atom_info);

	/// \brief Extend an existing (sugar) branch identified by \a asymID with one sugar containing atoms constructed from \a atom_info
	///
	/// \param asym_id      The asym id of the branch to extend
	/// \param atom_info    Array containing the info for the atoms to construct for the new sugar
	/// \param link_sugar   The sugar to link to, note: this is the sugar number (1 based)
	/// \param link_atom    The atom id of the atom linked in the sugar
	branch &extend_branch(const std::string &asym_id, std::vector<std::vector<item>> &atom_info,
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

	// /// \brief Direct access to underlying data
	// category &category(std::string_view name) const
	// {
	// 	return m_db[name];
	// }

	datablock &get_datablock() const
	{
		return m_db;
	}

	void validate_atoms() const;

  private:
	friend polymer;
	friend residue;

	std::string insert_compound(const std::string &compoundID, bool is_entity);

	std::string create_entity_for_branch(branch &branch);

	void loadData();

	void load_atoms_for_model(StructureOpenOptions options);

	template <typename... Args>
	atom &emplace_atom(Args... args)
	{
		return emplace_atom(atom{ std::forward<Args>(args)... });
	}

	atom &emplace_atom(atom &&atom);

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

} // namespace cif::mm