// Lib for working with structures as contained in mmCIF and PDB files

#pragma once

#include <boost/filesystem/operations.hpp>
#include <boost/math/quaternion.hpp>

#include <boost/any.hpp>

#include <clipper/core/coords.h>

#include "cif++/AtomType.h"
#include "cif++/Point.h"
#include "cif++/Compound.h"

/*
	To modify a structure, you will have to use actions.
	
	The currently supported actions are:
	
//	- Move atom to new location
	- Remove atom
//	- Add new atom that was formerly missing
//	- Add alternate Residue
	- 
	
	Other important design principles:
	
	- all objects here are references to the actual data. Not models of
	  the data itself. That means that if you copy an atom, you copy the
	  reference to an atom in the structure. You're not creating a new
	  atom. This may sound obvious, but it is not if you are used to
	  copy semantics in the C++ world.
	
	
*/

// forward declaration
namespace cif
{
	class Datablock;
	class File;
};

namespace libcif
{

class Atom;
class Residue;
class Monomer;
class Polymer;
class Structure;
class File;

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
	Atom(const File& f, const std::string& id);
	Atom(const Atom& rhs);

	~Atom();
	
	Atom& operator=(const Atom& rhs);

	std::string id() const;
	AtomType type() const;

	Point location() const;

	const Compound& comp() const;
	const Entity& ent() const;
	bool isWater() const;
	int charge() const;

	float uIso() const;
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
	int authSeqId() const;
	std::string pdbxAuthInsCode() const;
	std::string authAltId() const;
	
	bool operator==(const Atom& rhs) const;

	const File& getFile() const;
	
	// get clipper format Atom
	clipper::Atom toClipper() const;

  private:
 	struct AtomImpl*			mImpl;
};

typedef std::vector<Atom> AtomView;

// --------------------------------------------------------------------

class Residue : public std::enable_shared_from_this<Residue>
{
  public:
	Residue(const Compound& cmp) : mCompound(cmp) {}

	const Compound&		comp() const		{ return mCompound; }
	virtual AtomView	atoms();

  private:
	const Compound&		mCompound;
};

//// --------------------------------------------------------------------
//// a monomer models a single Residue in a protein chain 
//
//class monomer : public Residue
//{
//  public:
//	monomer(polymer& polymer, size_t seqId, const std::string& compId,
//		const std::string& altId);
//
//	int num() const								{ return mNum; }
////	polymer& getPolymer();
//
////	std::vector<monomer_ptr> alternates();
//
//  private:
//	polymer_ptr	mPolymer;
//	int			mNum;
//};
//
//// --------------------------------------------------------------------
//
//class polymer : public std::enable_shared_from_this<polymer>
//{
//  public:
//	polymer(const polymerEntity& pe, const std::string& asymId);
//	
//	struct iterator : public std::iterator<std::random_access_iterator_tag, monomer>
//	{
//		typedef std::iterator<std::bidirectional_iterator_tag, monomer>	base_type;
//		typedef base_type::reference									reference;
//		typedef base_type::pointer										pointer;
//		
//		iterator(polymer& list, uint32 index);
//		iterator(iterator&& rhs);
//		iterator(const iterator& rhs);
//		iterator& operator=(const iterator& rhs);
//		iterator& operator=(iterator&& rhs);
//		
//		reference	operator*();
//		pointer		operator->();
//		
//		iterator&	operator++();
//		iterator	operator++(int);
//		
//		iterator&	operator--();
//		iterator	operator--(int);
//
//		bool		operator==(const iterator& rhs) const;
//		bool		operator!=(const iterator& rhs) const;
//	};
//
//	iterator begin();
//	iterator end();
//
//  private:
//	polymer_entity				mEntity;
//	std::string					mAsymId;
//	std::vector<Residue_ptr>	mMonomers;
//};

// --------------------------------------------------------------------
// file is a reference to the data stored in e.g. the cif file.
// This object is not copyable.

class File : public std::enable_shared_from_this<File>
{
  public:
	File();
	File(boost::filesystem::path p);
	~File();

	File(const File&) = delete;
	File& operator=(const File&) = delete;

	void load(boost::filesystem::path p);
	void save(boost::filesystem::path p);
	
	Structure* model(size_t nr = 1);

	struct FileImpl& impl() const						{ return *mImpl; }

	std::vector<const Entity*> entities();

	cif::Datablock& data();
	cif::File& file();

  private:

	struct FileImpl*	mImpl;
};

// --------------------------------------------------------------------

class Structure
{
  public:
	Structure(File& p, uint32 modelNr = 1);
	Structure(const Structure&);
	Structure& operator=(const Structure&);
	~Structure();

	File& getFile() const;

	AtomView atoms() const;
	AtomView waters() const;

	Atom getAtomById(std::string id) const;
	Atom getAtomByLocation(Point pt, float maxDistance) const;
	
	Atom getAtomForLabel(const std::string& atomId, const std::string& asymId,
		const std::string& compId, int seqId, const std::string& altId = "");
	
	Atom getAtomForAuth(const std::string& atomId, const std::string& asymId,
		const std::string& compId, int seqId, const std::string& altId = "",
		const std::string& pdbxAuthInsCode = "");
	
	// map between auth and label locations
	
	std::tuple<std::string,int,std::string> MapAuthToLabel(const std::string& asymId,
		const std::string& seqId, const std::string& compId, const std::string& insCode = "");
	
	std::tuple<std::string,std::string,std::string,std::string> MapLabelToAuth(
		const std::string& asymId, int seqId, const std::string& compId);

	// returns chain, seqnr
	std::tuple<std::string,std::string> MapLabelToAuth(
		const std::string& asymId, int seqId);

	// returns chain,seqnr,comp,iCode
	std::tuple<std::string,int,std::string,std::string> MapLabelToPDB(
		const std::string& asymId, int seqId, const std::string& compId);

	std::tuple<std::string,int,std::string,std::string> MapPDBToLabel(
		const std::string& asymId, int seqId, const std::string& compId, const std::string& iCode);
	
	// Actions
	void removeAtom(Atom& a);

  private:
	struct StructureImpl*	mImpl;
};

}
