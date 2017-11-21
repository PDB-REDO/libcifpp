// Lib for working with structures as contained in mmCIF and PDB files

#pragma once

#include <boost/filesystem/operations.hpp>
#include <boost/math/quaternion.hpp>

#include <boost/any.hpp>

#include "libcif/atom_type.h"
#include "libcif/point.h"
#include "libcif/compound.h"

/*
	To modify a structure, you will have to use actions.
	
	The currently supported actions are:
	
//	- Move atom to new location
	- Remove atom
//	- Add new atom that was formerly missing
//	- Add alternate residue
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
	class datablock;
};


namespace libcif
{

class atom;
class residue;
class monomer;
class polymer;
class structure;
class file;

// --------------------------------------------------------------------
// We do not want to introduce a dependency on cif++ here, we might want
// to change the backend storage in the future.
// So, in order to access the data we use properties based on boost::any
// Eventually this should be moved to std::variant, but that's only when
// c++17 is acceptable.

struct property
{
	property() {}
	property(const std::string& name, const boost::any& value)
		: name(name), value(value) {}
	
	std::string name;
	boost::any value;
};

typedef std::vector<property>	property_list;

// --------------------------------------------------------------------

class atom
{
  public:
//	atom(const structure& s, const std::string& id);
	atom(struct atom_impl* impl);
	atom(const file& f, const std::string& id);
	atom(const atom& rhs);

	~atom();
	
	atom& operator=(const atom& rhs);

	std::string id() const;
	atom_type type() const;

	point location() const;

	const compound& comp() const;
	const entity& ent() const;
	bool is_water() const;
	int charge() const;
	
	boost::any property(const std::string& name) const;
	void property(const std::string& name, const boost::any& value);
	
	// specifications
	std::string label_atom_id() const;
	std::string label_comp_id() const;
	std::string label_asym_id() const;
	int label_seq_id() const;
	std::string label_alt_id() const;
	
	std::string auth_atom_id() const;
	std::string auth_comp_id() const;
	std::string auth_asym_id() const;
	int auth_seq_id() const;
	std::string pdbx_auth_ins_code() const;
	std::string auth_alt_id() const;
	
	bool operator==(const atom& rhs) const;

	const file& get_file() const;

  private:
 	struct atom_impl*			m_impl;
};

typedef std::vector<atom> atom_view;

// --------------------------------------------------------------------

class residue : public std::enable_shared_from_this<residue>
{
  public:
	residue(const compound& cmp) : m_compound(cmp) {}

	const compound&		comp() const		{ return m_compound; }
	virtual atom_view	atoms();

  private:
	const compound&		m_compound;
};

//// --------------------------------------------------------------------
//// a monomer models a single residue in a protein chain 
//
//class monomer : public residue
//{
//  public:
//	monomer(polymer& polymer, size_t seq_id, const std::string& comp_id,
//		const std::string& alt_id);
//
//	int num() const								{ return m_num; }
////	polymer& get_polymer();
//
////	std::vector<monomer_ptr> alternates();
//
//  private:
//	polymer_ptr	m_polymer;
//	int			m_num;
//};
//
//// --------------------------------------------------------------------
//
//class polymer : public std::enable_shared_from_this<polymer>
//{
//  public:
//	polymer(const polymer_entity& pe, const std::string& asym_id);
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
//	polymer_entity				m_entity;
//	std::string					m_asym_id;
//	std::vector<residue_ptr>	m_monomers;
//};

// --------------------------------------------------------------------
// file is a reference to the data stored in e.g. the cif file.
// This object is not copyable.

class file : public std::enable_shared_from_this<file>
{
  public:
	file();
	file(boost::filesystem::path p);
	~file();

	file(const file&) = delete;
	file& operator=(const file&) = delete;

	void load(boost::filesystem::path p);
	void save(boost::filesystem::path p);
	
	structure* model(size_t nr = 1);

	struct file_impl& impl() const						{ return *m_impl; }

	std::vector<const entity*> entities();

	cif::datablock& data();

  private:

	struct file_impl*	m_impl;
};

// --------------------------------------------------------------------

class structure
{
  public:
	structure(file& p, uint32 model_nr = 1);
	structure(const structure&);
	structure& operator=(const structure&);
	~structure();

	file& get_file() const;

	atom_view atoms() const;
	atom_view waters() const;

	atom get_atom_by_id(std::string id) const;
	atom get_atom_by_location(point pt, float max_distance) const;
	
	atom get_atom_for_label(const std::string& atom_id, const std::string& asym_id,
		const std::string& comp_id, int seq_id, const std::string& alt_id = "");
	
	atom get_atom_for_auth(const std::string& atom_id, const std::string& asym_id,
		const std::string& comp_id, int seq_id, const std::string& alt_id = "",
		const std::string& pdbx_auth_ins_code = "");
	
	// map between auth and label locations
	
	std::tuple<std::string,int,std::string> MapAuthToLabel(const std::string& asym_id,
		const std::string& seq_id, const std::string& comp_id, const std::string& ins_code = "");
	
	std::tuple<std::string,std::string,std::string,std::string> MapLabelToAuth(
		const std::string& asym_id, int seq_id, const std::string& comp_id);

	// returns chain, seqnr
	std::tuple<std::string,std::string> MapLabelToAuth(
		const std::string& asym_id, int seq_id);

	// returns chain,seqnr,comp,iCode
	std::tuple<std::string,int,std::string,std::string> MapLabelToPDB(
		const std::string& asym_id, int seq_id, const std::string& comp_id);

	std::tuple<std::string,int,std::string,std::string> MapPDBToLabel(
		const std::string& asym_id, int seq_id, const std::string& comp_id, const std::string& iCode);
	
	// Actions
	void remove_atom(atom& a);

  private:
	friend class action;

	struct structure_impl*	m_impl;
};

}
