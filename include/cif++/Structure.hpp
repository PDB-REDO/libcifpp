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

#include "cif++/AtomType.hpp"
#include "cif++/Cif++.hpp"
#include "cif++/Compound.hpp"
#include "cif++/Point.hpp"

/*
    To modify a structure, you will have to use actions.

    The currently supported actions are:

//	- Move atom to new location
    - Remove atom
//	- Add new atom that was formerly missing
//	- Add alternate Residue
    -

*/

namespace mmcif
{

class Atom;
class Residue;
class Monomer;
class Polymer;
class Structure;
class File;

// --------------------------------------------------------------------

class Atom
{
  private:
	struct AtomImpl : public std::enable_shared_from_this<AtomImpl>
	{
		AtomImpl(cif::Datablock &db, const std::string &id, cif::Row row);

		// constructor for a symmetry copy of an atom
		AtomImpl(const AtomImpl &impl, const Point &loc, const std::string &sym_op);

		AtomImpl(const AtomImpl &i) = default;

		void prefetch();

		int compare(const AtomImpl &b) const;

		bool getAnisoU(float anisou[6]) const;

		int charge() const;

		void moveTo(const Point &p);

		const Compound &comp() const;

		const std::string get_property(const std::string_view name) const;
		void set_property(const std::string_view name, const std::string &value);

		const cif::Datablock &mDb;
		std::string mID;
		AtomType mType;

		std::string mAtomID;
		std::string mCompID;
		std::string mAsymID;
		int mSeqID;
		std::string mAltID;
		std::string mAuthSeqID;

		Point mLocation;
		int mRefcount;
		cif::Row mRow;

		mutable std::vector<std::tuple<std::string, cif::detail::ItemReference>> mCachedRefs;

		mutable const Compound *mCompound = nullptr;

		bool mSymmetryCopy = false;
		bool mClone = false;

		std::string mSymmetryOperator = "1_555";
	};

  public:
	Atom() {}

	Atom(std::shared_ptr<AtomImpl> impl)
		: mImpl(impl)
	{
	}

	Atom(const Atom &rhs)
		: mImpl(rhs.mImpl)
	{
	}

	Atom(cif::Datablock &db, cif::Row &row);

	// a special constructor to create symmetry copies
	Atom(const Atom &rhs, const Point &symmmetry_location, const std::string &symmetry_operation);

	explicit operator bool() const { return (bool)mImpl; }

	// return a copy of this atom, with data copied instead of referenced
	Atom clone() const
	{
		auto copy = std::make_shared<AtomImpl>(*mImpl);
		copy->mClone = true;
		return Atom(copy);
	}

	Atom &operator=(const Atom &rhs) = default;

	template <typename T>
	T get_property(const std::string_view name) const;

	void set_property(const std::string_view name, const std::string &value)
	{
		if (not mImpl)
			throw std::logic_error("Error trying to modify an uninitialized atom");
		mImpl->set_property(name, value);
	}

	template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	void property(const std::string_view name, const T &value)
	{
		set_property(name, std::to_string(value));
	}

	const std::string &id() const { return impl().mID; }
	AtomType type() const { return impl().mType; }

	Point location() const { return impl().mLocation; }
	void location(Point p)
	{
		if (not mImpl)
			throw std::logic_error("Error trying to modify an uninitialized atom");
		mImpl->moveTo(p);
	}

	/// \brief Translate the position of this atom by \a t
	void translate(Point t);

	/// \brief Rotate the position of this atom by \a q
	void rotate(Quaternion q);

	/// \brief Translate and rotate the position of this atom by \a t and \a q
	void translateAndRotate(Point t, Quaternion q);

	/// \brief Translate, rotate and translate again the coordinates this atom by \a t1 , \a q and \a t2
	void translateRotateAndTranslate(Point t1, Quaternion q, Point t2);

	// for direct access to underlying data, be careful!
	const cif::Row getRow() const { return impl().mRow; }
	const cif::Row getRowAniso() const;

	bool isSymmetryCopy() const { return impl().mSymmetryCopy; }
	std::string symmetry() const { return impl().mSymmetryOperator; }

	const Compound &comp() const { return impl().comp(); }
	bool isWater() const { return impl().mCompID == "HOH" or impl().mCompID == "H2O" or impl().mCompID == "WAT"; }
	int charge() const;

	float uIso() const;
	bool getAnisoU(float anisou[6]) const { return impl().getAnisoU(anisou); }
	float occupancy() const;

	// specifications
	const std::string &labelAtomID() const { return impl().mAtomID; }
	const std::string &labelCompID() const { return impl().mCompID; }
	const std::string &labelAsymID() const { return impl().mAsymID; }
	std::string labelEntityID() const;
	int labelSeqID() const { return impl().mSeqID; }
	const std::string &labelAltID() const { return impl().mAltID; }
	bool isAlternate() const { return not impl().mAltID.empty(); }

	std::string authAtomID() const;
	std::string authCompID() const;
	std::string authAsymID() const;
	const std::string &authSeqID() const { return impl().mAuthSeqID; }
	std::string pdbxAuthInsCode() const;
	std::string pdbxAuthAltID() const;

	std::string labelID() const; // label_comp_id + '_' + label_asym_id + '_' + label_seq_id
	std::string pdbID() const;   // auth_comp_id + '_' + auth_asym_id + '_' + auth_seq_id + pdbx_PDB_ins_code

	bool operator==(const Atom &rhs) const;
	bool operator!=(const Atom &rhs) const
	{
		return not operator==(rhs);
	}

	// access data in compound for this atom

	// convenience routine
	bool isBackBone() const
	{
		auto atomID = labelAtomID();
		return atomID == "N" or atomID == "O" or atomID == "C" or atomID == "CA";
	}

	void swap(Atom &b)
	{
		std::swap(mImpl, b.mImpl);
	}

	int compare(const Atom &b) const { return impl().compare(*b.mImpl); }

	bool operator<(const Atom &rhs) const
	{
		return compare(rhs) < 0;
	}

	friend std::ostream &operator<<(std::ostream &os, const Atom &atom);

	/// \brief Synchronize data with underlying cif data
	void sync()
	{
		if (mImpl)
			mImpl->prefetch();
	}

  private:
	friend class Structure;

	const AtomImpl &impl() const
	{
		if (not mImpl)
			throw std::runtime_error("Uninitialized atom, not found?");
		return *mImpl;
	}

	std::shared_ptr<AtomImpl> mImpl;
};

template <>
inline std::string Atom::get_property<std::string>(const std::string_view name) const
{
	return impl().get_property(name);
}

template <>
inline int Atom::get_property<int>(const std::string_view name) const
{
	auto v = impl().get_property(name);
	return v.empty() ? 0 : stoi(v);
}

template <>
inline float Atom::get_property<float>(const std::string_view name) const
{
	return stof(impl().get_property(name));
}

inline void swap(mmcif::Atom &a, mmcif::Atom &b)
{
	a.swap(b);
}

inline double Distance(const Atom &a, const Atom &b)
{
	return Distance(a.location(), b.location());
}

inline double DistanceSquared(const Atom &a, const Atom &b)
{
	return DistanceSquared(a.location(), b.location());
}

typedef std::vector<Atom> AtomView;

// --------------------------------------------------------------------

enum class EntityType
{
	Polymer, NonPolymer, Macrolide, Water, Branched
};

// --------------------------------------------------------------------

class Residue
{
  public:
	// constructor
	Residue(const Structure &structure, const std::string &compoundID,
		const std::string &asymID, int seqID, const std::string &authSeqID)
		: mStructure(&structure)
		, mCompoundID(compoundID)
		, mAsymID(asymID)
		, mSeqID(seqID)
		, mAuthSeqID(authSeqID)
	{
	}

	Residue(const Residue &rhs) = delete;
	Residue &operator=(const Residue &rhs) = delete;

	Residue(Residue &&rhs);
	Residue &operator=(Residue &&rhs);

	virtual ~Residue();

	const Compound &compound() const;

	AtomView &atoms();
	const AtomView &atoms() const;

	void addAtom(Atom &atom);

	/// \brief Unique atoms returns only the atoms without alternates and the first of each alternate atom id.
	AtomView unique_atoms() const;

	/// \brief The alt ID used for the unique atoms
	std::string unique_alt_id() const;

	Atom atomByID(const std::string &atomID) const;

	const std::string &compoundID() const { return mCompoundID; }
	void setCompoundID(const std::string &id) { mCompoundID = id; }

	const std::string &asymID() const { return mAsymID; }
	int seqID() const { return mSeqID; }
	std::string entityID() const;

	EntityType entityType() const;

	std::string authAsymID() const;
	std::string authSeqID() const;
	std::string authInsCode() const;

	// return a human readable PDB-like auth id (chain+seqnr+iCode)
	std::string authID() const;

	// similar for mmCIF space
	std::string labelID() const;

	// Is this residue a single entity?
	bool isEntity() const;

	bool isWater() const { return mCompoundID == "HOH"; }

	const Structure &structure() const { return *mStructure; }

	bool empty() const { return mStructure == nullptr; }

	bool hasAlternateAtoms() const;

	/// \brief Return the list of unique alt ID's present in this residue
	std::set<std::string> getAlternateIDs() const;

	/// \brief Return the list of unique atom ID's
	std::set<std::string> getAtomIDs() const;

	/// \brief Return the list of atoms having ID \a atomID
	AtomView getAtomsByID(const std::string &atomID) const;

	// some routines for 3d work
	std::tuple<Point, float> centerAndRadius() const;

	friend std::ostream &operator<<(std::ostream &os, const Residue &res);

	friend Structure;

	bool operator==(const Residue &rhs) const
	{
		return mStructure == rhs.mStructure and mCompoundID == rhs.mCompoundID and mAsymID == rhs.mAsymID and mSeqID == rhs.mSeqID and mAuthSeqID == rhs.mAuthSeqID;
	}

  protected:
	Residue() {}

	friend class Polymer;

	const Structure *mStructure = nullptr;
	std::string mCompoundID, mAsymID;
	int mSeqID = 0;
	std::string mAuthSeqID;
	AtomView mAtoms;
};

// --------------------------------------------------------------------
// a monomer models a single Residue in a protein chain

class Monomer : public Residue
{
  public:
	//	Monomer();
	Monomer(const Monomer &rhs) = delete;
	Monomer &operator=(const Monomer &rhs) = delete;

	Monomer(Monomer &&rhs);
	Monomer &operator=(Monomer &&rhs);

	Monomer(const Polymer &polymer, size_t index, int seqID, const std::string &authSeqID,
		const std::string &compoundID);

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
	size_t nrOfChis() const;
	float chi(size_t i) const;

	bool isCis() const;

	/// \brief Returns true if the four atoms C, CA, N and O are present
	bool isComplete() const;

	/// \brief Returns true if any of the backbone atoms has an alternate
	bool hasAlternateBackboneAtoms() const;

	Atom CAlpha() const { return atomByID("CA"); }
	Atom C() const { return atomByID("C"); }
	Atom N() const { return atomByID("N"); }
	Atom O() const { return atomByID("O"); }
	Atom H() const { return atomByID("H"); }

	bool isBondedTo(const Monomer &rhs) const
	{
		return this != &rhs and areBonded(*this, rhs);
	}

	static bool areBonded(const Monomer &a, const Monomer &b, float errorMargin = 0.5f);
	static bool isCis(const Monomer &a, const Monomer &b);
	static float omega(const Monomer &a, const Monomer &b);

	// for LEU and VAL
	float chiralVolume() const;

	bool operator==(const Monomer &rhs) const
	{
		return mPolymer == rhs.mPolymer and mIndex == rhs.mIndex;
	}

  private:
	const Polymer *mPolymer;
	size_t mIndex;
};

// --------------------------------------------------------------------

class Polymer : public std::vector<Monomer>
{
  public:
	Polymer(const Structure &s, const std::string &entityID, const std::string &asymID);

	Polymer(const Polymer &) = delete;
	Polymer &operator=(const Polymer &) = delete;

	//	Polymer(Polymer&& rhs) = delete;
	//	Polymer& operator=(Polymer&& rhs) = de;

	Monomer &getBySeqID(int seqID);
	const Monomer &getBySeqID(int seqID) const;

	Structure *structure() const { return mStructure; }

	std::string asymID() const { return mAsymID; }
	std::string entityID() const { return mEntityID; }

	std::string chainID() const;

	int Distance(const Monomer &a, const Monomer &b) const;

  private:
	Structure *mStructure;
	std::string mEntityID;
	std::string mAsymID;
	cif::RowSet mPolySeq;
};

// --------------------------------------------------------------------
// Sugar and Branch, to describe glycosylation sites

class Branch;

class Sugar : public Residue
{
  public:
	Sugar(const Branch &branch, const std::string &compoundID,
		const std::string &asymID, int authSeqID);

	int num() const { return std::stoi(mAuthSeqID); }
	std::string name() const;

	/// \brief Return the atom the C1 is linked to
	Atom getLink() const { return mLink; }
	void setLink(Atom link) { mLink = link; }

  private:

	const Branch &mBranch;
	Atom mLink;
};

class Branch : public std::vector<Sugar>
{
  public:
	Branch(Structure &structure, const std::string &asymID);

	void linkAtoms();

	std::string name() const;
	float weight() const;
	std::string asymID() const { return mAsymID; }

	Structure &structure() { return *mStructure; }
	const Structure &structure() const { return *mStructure; }

	Sugar &getSugarByNum(int nr);
	const Sugar &getSugarByNum(int nr) const;

  private:
	friend Sugar;

	std::string name(const Sugar &s) const;

	Structure *mStructure;
	std::string mAsymID;
};

// --------------------------------------------------------------------
// file is a reference to the data stored in e.g. the cif file.
// This object is not copyable.

class File : public cif::File
{
  public:
	File() {}

	File(const std::filesystem::path &path)
	{
		load(path);
	}

	File(const char *data, size_t length)
	{
		load(data, length);
	}

	File(const File &) = delete;
	File &operator=(const File &) = delete;

	void load(const std::filesystem::path &p) override;
	void save(const std::filesystem::path &p) override;

	void load(std::istream &is) override;

	using cif::File::load;
	using cif::File::save;

	cif::Datablock &data() { return front(); }
};

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

class Structure
{
  public:
	Structure(cif::File &p, size_t modelNr = 1, StructureOpenOptions options = {})
		: Structure(p.front(), modelNr, options)
	{
	}

	Structure(cif::Datablock &db, size_t modelNr = 1, StructureOpenOptions options = {});

	// Create a read-only clone of the current structure (for multithreaded calculations that move atoms)
	Structure(const Structure &);

	Structure &operator=(const Structure &) = delete;
	~Structure();

	const AtomView &atoms() const { return mAtoms; }
	// AtomView &atoms() { return mAtoms; }

	EntityType getEntityTypeForEntityID(const std::string entityID) const;
	EntityType getEntityTypeForAsymID(const std::string asymID) const;

	AtomView waters() const;

	const std::list<Polymer> &polymers() const { return mPolymers; }
	std::list<Polymer> &polymers() { return mPolymers; }

	Polymer &getPolymerByAsymID(const std::string &asymID);

	const Polymer &getPolymerByAsymID(const std::string &asymID) const
	{
		return const_cast<Structure *>(this)->getPolymerByAsymID(asymID);
	}

	const std::list<Branch> &branches() const { return mBranches; }
	std::list<Branch> &branches() { return mBranches; }

	Branch &getBranchByAsymID(const std::string &asymID);
	const Branch &getBranchByAsymID(const std::string &asymID) const;

	const std::vector<Residue> &nonPolymers() const { return mNonPolymers; }
	
	Atom getAtomByID(const std::string &id) const;
	// Atom getAtomByLocation(Point pt, float maxDistance) const;

	Atom getAtomByLabel(const std::string &atomID, const std::string &asymID,
		const std::string &compID, int seqID, const std::string &altID = "");

	/// \brief Return the atom closest to point \a p
	Atom getAtomByPosition(Point p) const;

	/// \brief Return the atom closest to point \a p with atom type \a type in a residue of type \a res_type
	Atom getAtomByPositionAndType(Point p, std::string_view type, std::string_view res_type) const;

	/// \brief Get a non-poly residue for an asym with id \a asymID
	Residue &getResidue(const std::string &asymID)
	{
		return getResidue(asymID, 0, "");
	}

	/// \brief Get a non-poly residue for an asym with id \a asymID
	const Residue &getResidue(const std::string &asymID) const
	{
		return getResidue(asymID, 0, "");
	}

	/// \brief Get a residue for an asym with id \a asymID seq id \a seqID and authSeqID \a authSeqID
	Residue &getResidue(const std::string &asymID, int seqID, const std::string &authSeqID);

	/// \brief Get a the single residue for an asym with id \a asymID seq id \a seqID and authSeqID \a authSeqID
	const Residue &getResidue(const std::string &asymID, int seqID, const std::string &authSeqID) const
	{
		return const_cast<Structure *>(this)->getResidue(asymID, seqID, authSeqID);
	}

	/// \brief Get a residue for an asym with id \a asymID, compound id \a compID, seq id \a seqID and authSeqID \a authSeqID
	Residue &getResidue(const std::string &asymID, const std::string &compID, int seqID, const std::string &authSeqID);

	/// \brief Get a residue for an asym with id \a asymID, compound id \a compID, seq id \a seqID and authSeqID \a authSeqID
	const Residue &getResidue(const std::string &asymID, const std::string &compID, int seqID, const std::string &authSeqID) const
	{
		return const_cast<Structure *>(this)->getResidue(asymID, compID, seqID, authSeqID);
	}

	/// \brief Get a the residue for atom \a atom
	Residue &getResidue(const mmcif::Atom &atom)
	{
		return getResidue(atom.labelAsymID(), atom.labelCompID(), atom.labelSeqID(), atom.authSeqID());
	}

	/// \brief Get a the residue for atom \a atom
	const Residue &getResidue(const mmcif::Atom &atom) const
	{
		return getResidue(atom.labelAsymID(), atom.labelCompID(), atom.labelSeqID(), atom.authSeqID());
	}

	// Actions
	void removeAtom(Atom &a);
	void swapAtoms(Atom a1, Atom a2); // swap the labels for these atoms
	void moveAtom(Atom a, Point p);  // move atom to a new location
	void changeResidue(Residue &res, const std::string &newCompound,
		const std::vector<std::tuple<std::string, std::string>> &remappedAtoms);

	/// \brief Remove a residue, can be monomer or nonpoly
	///
	/// \param asym_id     The asym ID
	/// \param seq_id      The sequence ID
	void removeResidue(const std::string &asym_id, int seq_id, const std::string &auth_seq_id)
	{
		removeResidue(getResidue(asym_id, seq_id, auth_seq_id));
	}

	/// \brief Remove residue \a res, can be monomer or nonpoly
	void removeResidue(Residue &res);
	
	/// \brief Create a new non-polymer entity, returns new ID
	/// \param mon_id	The mon_id for the new nonpoly, must be an existing and known compound from CCD
	/// \return			The ID of the created entity
	std::string createNonPolyEntity(const std::string &mon_id);

	/// \brief Create a new NonPolymer struct_asym with atoms constructed from \a atoms, returns asym_id.
	/// This method assumes you are copying data from one cif file to another.
	///
	/// \param entity_id	The entity ID of the new nonpoly
	/// \param atoms		The array of atom_site rows containing the data.
	/// \return				The newly create asym ID
	std::string createNonpoly(const std::string &entity_id, const std::vector<mmcif::Atom> &atoms);

	/// \brief Create a new NonPolymer struct_asym with atoms constructed from info in \a atom_info, returns asym_id.
	/// This method creates new atom records filled with info from the info.
	///
	/// \param entity_id	The entity ID of the new nonpoly
	/// \param atoms		The array of sets of cif::item data containing the data for the atoms.
	/// \return				The newly create asym ID
	std::string createNonpoly(const std::string &entity_id, std::vector<std::vector<cif::Item>> &atom_info);

	/// \brief Create a new (sugar) branch with one first NAG containing atoms constructed from \a nag_atom_info
	Branch &createBranch(std::vector<std::vector<cif::Item>> &nag_atom_info);

	/// \brief Extend an existing (sugar) branch identified by \a asymID with one sugar containing atoms constructed from \a atom_info
	///
	/// \param asym_id      The asym id of the branch to extend
	/// \param atom_info    Array containing the info for the atoms to construct for the new sugar
	/// \param link_sugar   The sugar to link to, note: this is the sugar number (1 based)
	/// \param link_atom    The atom id of the atom linked in the sugar
	Branch &extendBranch(const std::string &asym_id, std::vector<std::vector<cif::Item>> &atom_info,
		int link_sugar, const std::string &link_atom);

	/// \brief Remove \a branch
	void removeBranch(Branch &branch);

	/// \brief Translate the coordinates of all atoms in the structure by \a t
	void translate(Point t);

	/// \brief Rotate the coordinates of all atoms in the structure by \a q
	void rotate(Quaternion t);

	/// \brief Translate and rotate the coordinates of all atoms in the structure by \a t and \a q
	void translateAndRotate(Point t, Quaternion q);

	/// \brief Translate, rotate and translate again the coordinates of all atoms in the structure by \a t1 , \a q and \a t2
	void translateRotateAndTranslate(Point t1, Quaternion q, Point t2);

	const std::vector<Residue> &getNonPolymers() const { return mNonPolymers; }

	void cleanupEmptyCategories();

	/// \brief Direct access to underlying data
	cif::Category &category(std::string_view name) const
	{
		return mDb[name];
	}

	cif::Datablock &datablock() const
	{
		return mDb;
	}

  private:
	friend Polymer;
	friend Residue;

	std::string insertCompound(const std::string &compoundID, bool isEntity);

	std::string createEntityForBranch(Branch &branch);

	void loadData();

	void loadAtomsForModel(StructureOpenOptions options);

	template<typename... Args>
	Atom& emplace_atom(Args ...args)
	{
		return emplace_atom(Atom{std::forward<Args>(args)...});
	}

	Atom &emplace_atom(Atom &&atom);

	cif::Datablock &mDb;
	size_t mModelNr;
	AtomView mAtoms;
	std::vector<size_t> mAtomIndex;
	std::list<Polymer> mPolymers;
	std::list<Branch> mBranches;
	std::vector<Residue> mNonPolymers;
};

} // namespace mmcif
