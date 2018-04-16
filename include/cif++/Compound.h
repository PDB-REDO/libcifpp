// Lib for working with structures as contained in mmCIF and PDB files

#pragma once

#include <set>
#include <tuple>
#include <vector>
#include <map>

#include "cif++/AtomType.h"

namespace libcif
{

// --------------------------------------------------------------------
// The chemical composition of the structure in an mmCIF file is 
// defined in the class composition. A compositon consists of
// entities. Each Entity can be either a polymer, a non-polymer
// a macrolide or a water molecule.
// Entities themselves are made up of compounds. And compounds
// contain CompoundAtom records for each atom.

class Composition;
class Entity;
class Compound;
struct CompoundAtom;

// --------------------------------------------------------------------
// struct containing information about an atom in a chemical compound
// This information comes from the CCP4 monomer library. 

struct CompoundAtom
{
	std::string	id;
	AtomType	typeSymbol;
	std::string	typeEnergy;
	float		partialCharge;
};

// --------------------------------------------------------------------
// struct containing information about the bonds
// This information comes from the CCP4 monomer library. 

enum CompoundBondType { singleBond, doubleBond, tripleBond, delocalizedBond };

struct CompoundBond
{
	std::string			atomID[2];
	CompoundBondType	type;
	float				distance;
	float				esd;
};

// --------------------------------------------------------------------
// struct containing information about the bond-angles
// This information comes from the CCP4 monomer library. 

struct CompoundAngle
{
	std::string			atomID[3];
	float				angle;
	float				esd;
};

// --------------------------------------------------------------------
// struct containing information about a chiral centre
// This information comes from the CCP4 monomer library. 

enum ChiralVolumeSign { negativ, positiv, both };

struct ChiralCentre
{
	std::string			id;
	std::string			atomIDCentre;
	std::string			atomID[3];
	ChiralVolumeSign	volumeSign;
};

// --------------------------------------------------------------------
// a class that contains information about a chemical compound.
// This information is derived from the ccp4 monomer library by default.
// To create compounds, you'd best use the factory method.

class Compound
{
  public:

	Compound(const std::string& id, const std::string& name,
		const std::string& group, std::vector<CompoundAtom>&& atoms,
		std::vector<CompoundBond>&& bonds, std::vector<CompoundAngle>&& angles,
		std::vector<ChiralCentre>&& chiralCentres)
		: mId(id), mName(name), mGroup(group)
		, mAtoms(std::move(atoms)), mBonds(std::move(bonds))
		, mAngles(std::move(angles))
		, mChiralCentres(std::move(chiralCentres))
	{
	}

	// factory method, create a Compound based on the three letter code
	// (for amino acids) or the one-letter code (for bases) or the
	// code as it is known in the CCP4 monomer library.

	static const Compound* create(const std::string& id);

	// this second factory method can create a Compound even if it is not
	// recorded in the library. It will take the values from the CCP4 lib
	// unless the value passed to this function is not empty.
	static const Compound* create(const std::string& id, const std::string& name,
		const std::string& type, const std::string& formula);
	
	// add an additional path to the monomer library.
	static void addMonomerLibraryPath(const std::string& dir);

	// accessors
	std::string id() const						{ return mId; }
	std::string	name() const					{ return mName; }
	std::string	type() const;
	std::string group() const					{ return mGroup; }
	std::vector<CompoundAtom> atoms() const		{ return mAtoms; }
	std::vector<CompoundBond> bonds() const		{ return mBonds; }
	std::vector<CompoundAngle> angles() const	{ return mAngles; }
	
	CompoundAtom getAtomById(const std::string& atomId) const;
	
	bool atomsBonded(const std::string& atomId_1, const std::string& atomId_2) const;
	float atomBondValue(const std::string& atomId_1, const std::string& atomId_2) const;

	std::string formula() const;
	float formulaWeight() const;
	int charge() const;
	bool isWater() const;
	bool isSugar() const;

	std::vector<ChiralCentre> chiralCentres() const			{ return mChiralCentres; }
	
	std::vector<std::string>	isomers() const;
	bool isIsomerOf(const Compound& c) const;
	std::vector<std::tuple<std::string,std::string>> mapToIsomer(const Compound& c) const;

  private:

	~Compound();

//	Entity&						mEntity;
	std::string					mId;
	std::string					mName;
	std::string					mGroup;
	std::vector<CompoundAtom>	mAtoms;
	std::vector<CompoundBond>	mBonds;
	std::vector<CompoundAngle>	mAngles;
	std::vector<ChiralCentre>	mChiralCentres;
};

// --------------------------------------------------------------------
// an Entity. This is a base class for PolymerEntity and NonPolyEntity
// The latter can be either a regular non-polymer (residue), a macrolide or
// water.

class Entity
{
  public:
	Entity(const std::string& id, const std::string& type, const std::string& description);
	virtual ~Entity();

	std::string id() const;
	std::string	type() const;
	std::string description() const;

	virtual float formulaWeight() const = 0;
	
  private:
	std::string		mId;
	std::string		mType;
	std::string		mDescription;
};

// --------------------------------------------------------------------
// A polymer Entity

class PolymerEntity : public Entity
{
  public:
	PolymerEntity(const std::string& id, const std::string& description);
	~PolymerEntity();
	
	std::string		seqOneLetterCode(bool cannonical) const;
	std::string		pdbxStrandId() const;
	virtual float	formulaWeight() const;
	
	class monomer
	{
	  public:
		friend class PolymerEntity;
		
		size_t				num() const;					// sequence number
		bool				hetero() const;					// whether this position contains alternate Compounds
		const Compound&		comp(size_t altNr) const;		// the chemical Compound of this monomer
		
	  private:
		monomer*	mNext;
		monomer*	mAlt;
		size_t		mNum;
		Compound*	mComp;
	};
	
	class iterator : public std::iterator<std::forward_iterator_tag, const monomer>
	{
	  public:
		typedef std::iterator<std::forward_iterator_tag, const monomer>	baseType;
		typedef baseType::reference									reference;
		typedef baseType::pointer										pointer;
		
		iterator(monomer* monomer = nullptr)
			: mCursor(monomer) {}

		iterator(const iterator& rhs)
			: mCursor(rhs.mCursor)
		{
		}

		iterator& operator=(const iterator& rhs)
		{
			mCursor = rhs.mCursor;
			return *this;
		}
		
		reference	operator*()			{ return *mCursor; }
		pointer		operator->()		{ return mCursor; }
		
		iterator&	operator++()		{ mCursor = mCursor->mNext; return *this; }
		iterator	operator++(int)
		{
			iterator tmp(*this);
			operator++();
			return tmp;
		}
		
		bool		operator==(const iterator& rhs) const		{ return mCursor == rhs.mCursor; }
		bool		operator!=(const iterator& rhs) const		{ return mCursor != rhs.mCursor; }

	  private:
		monomer*	mCursor;
	};
	
	iterator begin() const		{ return iterator(mSeq); }
	iterator end() const		{ return iterator(); }
	
	const monomer& operator[](size_t index) const;

  private:
	Entity&		mEntity;
	monomer*	mSeq;
};

// --------------------------------------------------------------------
// nonPoly Entity 

class NonPolyEntity : public Entity
{
  public:
	NonPolyEntity(const std::string& id, const std::string& type, const std::string& description);
	~NonPolyEntity();
	
	Compound&		comp() const;
	virtual float	formulaWeight() const;

  private:
 	Compound*	mCompound;
};

}
