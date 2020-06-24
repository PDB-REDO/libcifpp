// Lib for working with structures as contained in mmCIF and PDB files

#pragma once

#include <numeric>

#include <boost/math/quaternion.hpp>

#include <boost/any.hpp>

#include <clipper/core/coords.h>

#include "cif++/AtomType.h"
#include "cif++/Point.h"
#include "cif++/Compound.h"
#include "cif++/Cif++.h"

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

struct SecondaryStructure;

// --------------------------------------------------------------------
// We do not want to introduce a dependency on cif++ here, we might want
// to change the backend storage in the future.
// So, in order to access the data we use properties based on boost::any
// Eventually this should be moved to std::variant, but that's only when
// c++17 is acceptable.

struct Property
{
	Property() {}
	Property(const std::string& name, const boost::any& value)
		: name(name), value(value) {}
	
	std::string name;
	boost::any value;
};

typedef std::vector<Property>	PropertyList;

// --------------------------------------------------------------------

class Atom
{
  public:
//	Atom(const structure& s, const std::string& id);
	Atom(struct AtomImpl* impl);
	Atom(const Atom& rhs);

	~Atom();
	
	// return a copy of this atom, with data copied instead of referenced
	Atom clone() const;
	
	Atom& operator=(const Atom& rhs);

	const std::string& id() const;
	AtomType type() const;

	Point location() const;
	void location(Point p);
	
	Atom symmetryCopy(const Point& d, const clipper::RTop_orth& rt);
	bool isSymmetryCopy() const;
	std::string symmetry() const;
	const clipper::RTop_orth& symop() const;
	
	const Compound& comp() const;
	bool isWater() const;
	int charge() const;

	float uIso() const;
	bool getAnisoU(float anisou[6]) const;
	float occupancy() const;
	
	template<typename T>
	T property(const std::string& name) const;
	
	template<typename T>
	void property(const std::string& name, const T& value);
	
	// specifications
	std::string labelAtomId() const;
	std::string labelCompId() const;
	std::string labelAsymId() const;
	int labelSeqId() const;
	std::string labelAltId() const;
	
	std::string authAtomId() const;
	std::string authCompId() const;
	std::string authAsymId() const;
	std::string authSeqId() const;
	std::string pdbxAuthInsCode() const;
	std::string pdbxAuthAltId() const;
	
	std::string labelID() const;// label_comp_id + '_' + label_asym_id + '_' + label_seq_id
	std::string pdbID() const;	// auth_comp_id + '_' + auth_asym_id + '_' + auth_seq_id + pdbx_PDB_ins_code
	
	bool operator==(const Atom& rhs) const;

	// get clipper format Atom
	clipper::Atom toClipper() const;

	// Radius calculation based on integrating the density until perc of electrons is found
	void calculateRadius(float resHigh, float resLow, float perc);
	float radius() const;
	
	// access data in compound for this atom

	// the energy-type field
	std::string energyType() const;
	
	// convenience routine
	bool isBackBone() const
	{
		return labelAtomId() == "N" or labelAtomId() == "O" or
			labelAtomId() == "C" or labelAtomId() == "CA";
	}
	
	void swap(Atom& b)
	{
		std::swap(mImpl, b.mImpl);
	}

	int compare(const Atom& b) const;

  private:
	friend class Structure;
	void setID(int id);

 	struct AtomImpl*			mImpl;
};

inline void swap(mmcif::Atom& a, mmcif::Atom& b)
{
	a.swap(b);
}

inline double Distance(const Atom& a, const Atom& b)
{
	return Distance(a.location(), b.location());
}

inline double DistanceSquared(const Atom& a, const Atom& b)
{
	return DistanceSquared(a.location(), b.location());
}

typedef std::vector<Atom> AtomView;

// --------------------------------------------------------------------

class Residue
{
  public:
	// constructors should be private, but that's not possible for now (needed in emplace)
	// constructor for waters
	Residue(const Structure& structure, const std::string& compoundID,
		const std::string& asymID, const std::string& authSeqID);

	Residue(const Structure& structure, const std::string& compoundID,
		const std::string& asymID, int seqID = 0);

	Residue(const Residue& rhs) = delete;
	Residue& operator=(const Residue& rhs) = delete;

	Residue(Residue&& rhs);
	Residue& operator=(Residue&& rhs);

	virtual ~Residue();
	
	const Compound&		compound() const;
	const AtomView&		atoms() const;
	
	Atom				atomByID(const std::string& atomID) const;

	const std::string&	compoundID() const	{ return mCompoundID; }
	const std::string&	asymID() const		{ return mAsymID; }
	int					seqID() const		{ return mSeqID; }
	
	std::string			authSeqID() const;
	std::string			authInsCode() const;
	
	// return a human readable PDB-like auth id (chain+seqnr+iCode)
	std::string			authID() const;

	// similar for mmCIF space
	std::string			labelID() const;
	
	// Is this residue a single entity?
	bool isEntity() const;
	
	bool isWater() const					{ return mCompoundID == "HOH"; }
	
	const Structure& structure() const		{ return *mStructure; }

	bool empty() const						{ return mStructure == nullptr; }
	
	// some routines for 3d work
	std::tuple<Point,float> centerAndRadius() const;

  protected:

	Residue() {}

	friend class Polymer;

	const Structure* mStructure = nullptr;
	std::string	mCompoundID, mAsymID;
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
	Monomer(const Monomer& rhs) = delete;
	Monomer& operator=(const Monomer& rhs) = delete;
	
	Monomer(Monomer&& rhs);
	Monomer& operator=(Monomer&& rhs);

//	Monomer(const Polymer& polymer, uint32_t index);
	Monomer(const Polymer& polymer, uint32_t index, int seqID,
		const std::string& compoundID);

	// Assuming this is really an amino acid...
	
	float phi() const;
	float psi() const;
	float alpha() const;
	float kappa() const;
	float tco() const;

    // torsion angles
    size_t nrOfChis() const;
    float chi(size_t i) const;
	
	bool isCis() const;

	/// \brief Returns true if the four atoms C, CA, N and O are present
	bool isComplete() const;

	Atom CAlpha() const		{ return atomByID("CA"); }
	Atom C() const			{ return atomByID("C"); }
	Atom N() const			{ return atomByID("N"); }
	Atom O() const			{ return atomByID("O"); }
	Atom H() const			{ return atomByID("H"); }

	bool isBondedTo(const Monomer& rhs) const
	{
		return this != &rhs and areBonded(*this, rhs);
	}

	static bool areBonded(const Monomer& a, const Monomer& b, float errorMargin = 0.5f);
	static bool isCis(const Monomer& a, const Monomer& b);
	static float omega(const Monomer& a, const Monomer& b);

    // for LEU and VAL
    float chiralVolume() const;
	
  private:
	const Polymer*	mPolymer;
	uint32_t		mIndex;
};

// --------------------------------------------------------------------

class Polymer : public std::vector<Monomer>
{
  public:
	Polymer(const Structure& s, const std::string& entityID, const std::string& asymID);

	Polymer(const Polymer&) = delete;
	Polymer& operator=(const Polymer&) = delete;
	
//	Polymer(Polymer&& rhs) = delete;
//	Polymer& operator=(Polymer&& rhs) = de;

	Monomer& getBySeqID(int seqID);
	const Monomer& getBySeqID(int seqID) const;

	Structure* structure() const	{ return mStructure; }
	
	std::string asymID() const		{ return mAsymID; }
	std::string entityID() const	{ return mEntityID; }
	
	std::string chainID() const;
	
	int Distance(const Monomer& a, const Monomer& b) const;

  private:

	Structure*					mStructure;
	std::string					mEntityID;
	std::string					mAsymID;
	cif::RowSet					mPolySeq;
};

// --------------------------------------------------------------------
// file is a reference to the data stored in e.g. the cif file.
// This object is not copyable.

class File : public std::enable_shared_from_this<File>
{
  public:
	File();
	File(const std::string& path);
	~File();

	File(const File&) = delete;
	File& operator=(const File&) = delete;

	void load(const std::string& path);
	void save(const std::string& path);
	
	Structure* model(size_t nr = 1);

	struct FileImpl& impl() const						{ return *mImpl; }

	cif::Datablock& data();
	cif::File& file();

  private:

	struct FileImpl*	mImpl;
};

// --------------------------------------------------------------------

class Structure
{
  public:
	Structure(File& p, uint32_t modelNr = 1);
	Structure& operator=(const Structure&) = delete;
	~Structure();

	// Create a read-only clone of the current structure (for multithreaded calculations that move atoms)
	Structure(const Structure&);

	File& getFile() const;

	const AtomView& atoms() const							{ return mAtoms; }
	AtomView waters() const;
	
	const std::list<Polymer>& polymers() const				{ return mPolymers; }
	const std::vector<Residue>& nonPolymers() const			{ return mNonPolymers; }

	Atom getAtomById(std::string id) const;
	// Atom getAtomByLocation(Point pt, float maxDistance) const;
	
	Atom getAtomByLabel(const std::string& atomId, const std::string& asymId,
		const std::string& compId, int seqId, const std::string& altId = "");
	
	// Atom getAtomByAuth(const std::string& atomId, const std::string& asymId,
	// 	const std::string& compId, int seqId, const std::string& altId = "",
	// 	const std::string& pdbxAuthInsCode = "");
	
	// map between auth and label locations
	
	std::tuple<std::string,int,std::string> MapAuthToLabel(const std::string& asymId,
		const std::string& seqId, const std::string& compId, const std::string& insCode = "");
	
	std::tuple<std::string,std::string,std::string,std::string> MapLabelToAuth(
		const std::string& asymId, int seqId, const std::string& compId);

	// returns chain, seqnr, icode
	std::tuple<char,int,char> MapLabelToAuth(
		const std::string& asymId, int seqId) const;

	// returns chain,seqnr,comp,iCode
	std::tuple<std::string,int,std::string,std::string> MapLabelToPDB(
		const std::string& asymId, int seqId, const std::string& compId,
		const std::string& authSeqID) const;

	std::tuple<std::string,int,std::string> MapPDBToLabel(
		const std::string& asymId, int seqId, const std::string& compId, const std::string& iCode) const;
	
	// Actions
	void removeAtom(Atom& a);
	void swapAtoms(Atom& a1, Atom& a2);	// swap the labels for these atoms
	void moveAtom(Atom& a, Point p);	// move atom to a new location
	void changeResidue(const Residue& res, const std::string& newCompound,
		const std::vector<std::tuple<std::string,std::string>>& remappedAtoms);

	/// To sort the atoms in order of model > asym-id > res-id > atom-id
	/// Will asssign new atom_id's to all atoms. Be carefull
	void sortAtoms();
	
	// iterator for all residues
	
	class residue_iterator : public std::iterator<std::forward_iterator_tag, const Residue>
	{
	  public:
		typedef std::iterator<std::forward_iterator_tag, const Residue>	baseType;
		typedef typename baseType::pointer								pointer;
		typedef typename baseType::reference							reference;
		
		typedef std::list<Polymer>::const_iterator						poly_iterator;
		
		residue_iterator(const Structure* s, poly_iterator polyIter, size_t polyResIndex, size_t nonPolyIndex);
		
		reference operator*();
		pointer operator->();
		
		residue_iterator& operator++();
		residue_iterator operator++(int);
		
		bool operator==(const residue_iterator& rhs) const;
		bool operator!=(const residue_iterator& rhs) const;
		
	  private:
		const Structure&	mStructure;
		poly_iterator		mPolyIter;
		size_t				mPolyResIndex;
		size_t				mNonPolyIndex;
	};
	
	class residue_view
	{
	  public:
		residue_view(const Structure* s) : mStructure(s) {}
		residue_view(const residue_view& rhs) : mStructure(rhs.mStructure) {}
		residue_view& operator=(residue_view& rhs)
		{
			mStructure = rhs.mStructure;
			return *this;
		}
		
		residue_iterator begin() const		{ return residue_iterator(mStructure, mStructure->mPolymers.begin(), 0, 0); }
		residue_iterator end() const		{ return residue_iterator(mStructure, mStructure->mPolymers.end(), 0, mStructure->mNonPolymers.size()); }
		size_t size() const
		{
			size_t ps = std::accumulate(mStructure->mPolymers.begin(), mStructure->mPolymers.end(), 0UL, [](size_t s, auto& p) { return s + p.size(); });
			return ps + mStructure->mNonPolymers.size();
		}

	  private:
		const Structure* mStructure;
	};
	
	residue_view residues() const			{ return residue_view(this); }
	
  private:
	friend Polymer;
	friend Residue;
	friend residue_view;
	friend residue_iterator;

	cif::Category& category(const char* name) const;
	cif::Datablock& datablock() const;

	void insertCompound(const std::string& compoundID, bool isEntity);
	
	void loadData();
	void updateAtomIndex();
	
	File&					mFile;
	uint32_t					mModelNr;
	AtomView				mAtoms;
	std::vector<size_t>		mAtomIndex;
	std::list<Polymer>		mPolymers;
	std::vector<Residue>	mNonPolymers;
};

}
