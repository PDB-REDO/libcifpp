// Lib for working with structures as contained in mmCIF and PDB files

#pragma once

#include <set>
#include <tuple>
#include <vector>
#include <map>

#include "libcif/atom_type.h"

namespace libcif
{

// --------------------------------------------------------------------
// The chemical composition of the structure in an mmCIF file is 
// defined in the class composition. A compositon consists of
// entities. Each entity can be either a polymer, a non-polymer
// a macrolide or a water molecule.
// Entities themselves are made up of compounds. And compounds
// contain comp_atom records for each atom.

class composition;
class entity;
class compound;
struct comp_atom;

// --------------------------------------------------------------------
// struct containing information about an atom in a chemical compound
// This information comes from the CCP4 monomer library. 

struct comp_atom
{
	std::string	id;
	atom_type	type_symbol;
	std::string	type_energy;
	float		partial_charge;
};

// --------------------------------------------------------------------
// a class that contains information about a chemical compound.
// This information is derived from the ccp4 monomer library by default.
// To create compounds, you'd best use the factory method.

class compound
{
  public:

	compound(const std::string& id, const std::string& name,
		const std::string& group, std::vector<comp_atom>&& atoms,
		std::map<std::tuple<std::string,std::string>,float>&& bonds)
		: m_id(id), m_name(name), m_group(group)
		, m_atoms(std::move(atoms)), m_bonds(std::move(bonds))
	{
	}

	~compound();

	// factory method, create a compound based on the three letter code
	// (for amino acids) or the one-letter code (for bases) or the
	// code as it is known in the CCP4 monomer library.

	static const compound* create(const std::string& id);

	// this second factory method can create a compound even if it is not
	// recorded in the library. It will take the values from the CCP4 lib
	// unless the value passed to this function is not empty.
	static const compound* create(const std::string& id, const std::string& name,
		const std::string& type, const std::string& formula);
	
	// add an additional path to the monomer library.
	static void add_monomer_library_path(const std::string& dir);

	// accessors
	std::string id() const					{ return m_id; }
	std::string	name() const				{ return m_name; }
	std::string	type() const;
//	std::string group() const				{ return m_group; }
	std::vector<comp_atom> atoms() const	{ return m_atoms; }
	
	comp_atom get_atom_by_id(const std::string& atom_id) const;
	
	bool atoms_bonded(const std::string& atom_id_1, const std::string& atom_id_2) const;
	float atom_bond_value(const std::string& atom_id_1, const std::string& atom_id_2) const;

	std::string formula() const;
	float formula_weight() const;
	int charge() const;
	bool is_water() const;

  private:
//	entity&					m_entity;
	std::string				m_id;
	std::string				m_name;
	std::string				m_group;
	std::vector<comp_atom>	m_atoms;
	std::map<std::tuple<std::string,std::string>,float>	m_bonds;
};

// --------------------------------------------------------------------
// an entity. This is a base class for polymer_entity and non_poly_entity
// The latter can be either a regular non-polymer (residue), a macrolide or
// water.

class entity
{
  public:
	entity(const std::string& id, const std::string& type, const std::string& description);
	virtual ~entity();

	std::string id() const;
	std::string	type() const;
	std::string description() const;

	virtual float formula_weight() const = 0;
	
  private:
	std::string		m_id;
	std::string		m_type;
	std::string		m_description;
};

// --------------------------------------------------------------------
// A polymer entity

class polymer_entity : public entity
{
  public:
	polymer_entity(const std::string& id, const std::string& description);
	~polymer_entity();
	
	std::string		seq_one_letter_code(bool cannonical) const;
	std::string		pdbx_strand_id() const;
	virtual float	formula_weight() const;
	
	class monomer
	{
	  public:
		friend class polymer_entity;
		
		size_t				num() const;					// sequence number
		bool				hetero() const;					// whether this position contains alternate compounds
		const compound&		comp(size_t alt_nr) const;		// the chemical compound of this monomer
		
	  private:
		monomer*	m_next;
		monomer*	m_alt;
		size_t		m_num;
		compound*	m_comp;
	};
	
	class iterator : public std::iterator<std::forward_iterator_tag, const monomer>
	{
	  public:
		typedef std::iterator<std::forward_iterator_tag, const monomer>	base_type;
		typedef base_type::reference									reference;
		typedef base_type::pointer										pointer;
		
		iterator(monomer* monomer = nullptr)
			: m_cursor(monomer) {}

		iterator(const iterator& rhs)
			: m_cursor(rhs.m_cursor)
		{
		}

		iterator& operator=(const iterator& rhs)
		{
			m_cursor = rhs.m_cursor;
			return *this;
		}
		
		reference	operator*()			{ return *m_cursor; }
		pointer		operator->()		{ return m_cursor; }
		
		iterator&	operator++()		{ m_cursor = m_cursor->m_next; return *this; }
		iterator	operator++(int)
		{
			iterator tmp(*this);
			operator++();
			return tmp;
		}
		
		bool		operator==(const iterator& rhs) const		{ return m_cursor == rhs.m_cursor; }
		bool		operator!=(const iterator& rhs) const		{ return m_cursor != rhs.m_cursor; }

	  private:
		monomer*	m_cursor;
	};
	
	iterator begin() const		{ return iterator(m_seq); }
	iterator end() const		{ return iterator(); }
	
	const monomer& operator[](size_t index) const;

  private:
	entity&		m_entity;
	monomer*	m_seq;
};

// --------------------------------------------------------------------
// non_poly entity 

class non_poly_entity : public entity
{
  public:
	non_poly_entity(const std::string& id, const std::string& type, const std::string& description);
	~non_poly_entity();
	
	compound&		comp() const;
	virtual float	formula_weight() const;

  private:
 	compound*	m_compound;
};

}
