// Lib for working with structures as contained in mmCIF and PDB files

#include "libcif/config.h"

#include <map>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include "libcif/compound.h"
#include "libcif/cif++.h"

using namespace std;
namespace ba = boost::algorithm;
namespace fs = boost::filesystem;

namespace libcif
{

class compound_factory
{
  public:
	
	static compound_factory& instance();
	const compound* create(string id);

  private:
	compound_factory();
	~compound_factory();
	
	static compound_factory* sInstance;

	fs::path m_clibd_mon;
	vector<compound*> m_compounds;
};


// --------------------------------------------------------------------
// compound

string compound::formula() const
{
	string result;
	
	map<string,uint32> atoms;
	float charge_sum = 0;

	for (auto r: m_atoms)
	{
		atoms[atom_type_traits(r.type_symbol).symbol()] += 1;
		charge_sum += r.partial_charge;
	}
	
	auto c = atoms.find("C");
	if (c != atoms.end())
	{
		result = "C";
		
		if (c->second > 1)
			result += to_string(c->second);
		
		atoms.erase(c);
		
		auto h = atoms.find("H");
		if (h != atoms.end())
		{
			result += " H";
			if (h->second > 1)
				result += to_string(h->second);
			
			atoms.erase(h);
		}
	}
	
	for (auto a: atoms)
	{
		if (not result.empty())
			result += ' ';
		
		result += a.first;
		if (a.second > 1)
			result += to_string(a.second);	
	}

	int charge = lrint(charge_sum);
	if (charge != 0)
		result += ' ' + to_string(charge);

	return result;
}

int compound::charge() const
{
	float result = 0;

	for (auto r: m_atoms)
		result += r.partial_charge;

	return lrint(result);
}

string compound::type() const
{
	string result;
	
	// known groups are (counted from ccp4 monomer dictionary)

	//	D-pyranose
	//	DNA
	//	L-PEPTIDE LINKING
	//	L-SACCHARIDE
	//	L-peptide
	//	L-pyranose
	//	M-peptide
	//	NON-POLYMER
	//	P-peptide
	//	RNA
	//	furanose
	//	non-polymer
	//	non_polymer
	//	peptide
	//	pyranose
	//	saccharide
	
	if (cif::iequals(m_id, "gly"))
		result = "peptide linking";
	else if (cif::iequals(m_group, "l-peptide") or cif::iequals(m_group, "L-peptide linking") or cif::iequals(m_group, "peptide"))
		result = "L-peptide linking";
	else if (cif::iequals(m_group, "DNA"))
		result = "DNA linking";
	else if (cif::iequals(m_group, "RNA"))
		result = "RNA linking";
	
	return result;
}

bool compound::is_water() const
{
	return m_id == "HOH" or m_id == "H2O";
}

comp_atom compound::get_atom_by_id(const string& atom_id) const
{
	comp_atom result;
	for (auto& a: m_atoms)
	{
		if (a.id == atom_id)
		{
			result = a;
			break;
		}
	}

	if (result.id != atom_id)	
		throw out_of_range("No atom " + atom_id + " in compound " + m_id);
	
	return result;
}

const compound* compound::create(const string& id)
{
	return compound_factory::instance().create(id);
}

// --------------------------------------------------------------------
// a factory class to generate compounds

compound_factory* compound_factory::sInstance = nullptr;

compound_factory::compound_factory()
{
	const char* clibd_mon = getenv("CLIBD_MON");
	if (clibd_mon == nullptr)
		throw runtime_error("Cannot locate peptide list, please souce the CCP4 environment");
	m_clibd_mon = clibd_mon;
}

compound_factory::~compound_factory()
{
}

compound_factory& compound_factory::instance()
{
	if (sInstance == nullptr)
		sInstance = new compound_factory();
	return *sInstance;
}

// id is the three letter code
const compound* compound_factory::create(std::string id)
{
	ba::to_upper(id);

	compound* result = nullptr;
	
	for (auto cmp: m_compounds)
	{
		if (cmp->id() == id)
		{
			result = cmp;
			break;
		}
	}
	
	if (result == nullptr)
	{
		fs::path resFile = m_clibd_mon / ba::to_lower_copy(id.substr(0, 1)) / (id + ".cif");
		fs::ifstream file(resFile);
		if (file.is_open())
		{
			cif::file cf;
			
			try
			{
				cf.load(file);
			}
			catch (const exception& ex)
			{
				cerr << "Error while loading " << resFile << endl;
				throw ex;
			}
			
			auto& list = cf["comp_list"];
			auto row = list["chem_comp"][cif::key("id") == id];
			
			string name, group;
			uint32 number_atoms_all, number_atoms_nh;
			cif::tie(name, group, number_atoms_all, number_atoms_nh) =
				row.get("name", "group", "number_atoms_all", "number_atoms_nh");
	
			ba::trim(name);
			ba::trim(group);
			
			auto& comp_atoms = cf["comp_" + id]["chem_comp_atom"];
			
			vector<comp_atom> atoms;
			for (auto row: comp_atoms)
			{
				string id, symbol, energy;
				float charge;
				
				cif::tie(id, symbol, energy, charge) = row.get("atom_id", "type_symbol", "type_energy", "partial_charge");
				
				atoms.push_back({
					id, atom_type_traits(symbol).type(), energy, charge
				});
			}

			auto& comp_bonds = cf["comp_" + id]["chem_comp_bond"];
			
			map<tuple<string,string>,float> bonds;
			for (auto row: comp_bonds)
			{
				string atom_id_1, atom_id_2, type;
				
				cif::tie(atom_id_1, atom_id_2, type) = row.get("atom_id_1", "atom_id_2", "type");
				
				float value = 0;
				if (type == "single")		value = 1;
				else if (type == "double")	value = 2;
				else if (type == "triple")	value = 3;
				else if (type == "deloc" or type == "aromat")
											value = 1.5;
				else
				{
					cerr << "Unimplemented chem_comp_bond.type " << type << " in file " << resFile << endl;
					value = 1.0;
				}
				
				bonds[make_tuple(atom_id_1, atom_id_2)] = value;
			}
			
			result = new compound(id, name, group, move(atoms), move(bonds));
			m_compounds.push_back(result);
		}
	}
	
	return result;
}

bool compound::atoms_bonded(const string& atom_id_1, const string& atom_id_2) const
{
	return m_bonds.count(make_tuple(atom_id_1, atom_id_2)) or m_bonds.count(make_tuple(atom_id_2, atom_id_1));
}

float compound::atom_bond_value(const string& atom_id_1, const string& atom_id_2) const
{
	auto i = m_bonds.find(make_tuple(atom_id_1, atom_id_2));
	if (i == m_bonds.end())
		i = m_bonds.find(make_tuple(atom_id_2, atom_id_1));
	
	return i == m_bonds.end() ? 0 : i->second;
}

}
