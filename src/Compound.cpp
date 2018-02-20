// Lib for working with structures as contained in mmCIF and PDB files

#include "cif++/Config.h"

#include <map>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include "cif++/Cif++.h"
#include "cif++/Compound.h"

using namespace std;
namespace ba = boost::algorithm;
namespace fs = boost::filesystem;

namespace libcif
{

class CompoundFactory
{
  public:
	
	static CompoundFactory& instance();
	const Compound* create(string id);

  private:
	CompoundFactory();
	~CompoundFactory();
	
	static CompoundFactory* sInstance;

	fs::path mClibdMon;
	vector<Compound*> mCompounds;
};


// --------------------------------------------------------------------
// Compound

string Compound::formula() const
{
	string result;
	
	map<string,uint32> atoms;
	float chargeSum = 0;

	for (auto r: mAtoms)
	{
		atoms[AtomTypeTraits(r.typeSymbol).symbol()] += 1;
		chargeSum += r.partialCharge;
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

	int charge = lrint(chargeSum);
	if (charge != 0)
		result += ' ' + to_string(charge);

	return result;
}

int Compound::charge() const
{
	float result = 0;

	for (auto r: mAtoms)
		result += r.partialCharge;

	return lrint(result);
}

string Compound::type() const
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
	
	if (cif::iequals(mId, "gly"))
		result = "peptide linking";
	else if (cif::iequals(mGroup, "l-peptide") or cif::iequals(mGroup, "L-peptide linking") or cif::iequals(mGroup, "peptide"))
		result = "L-peptide linking";
	else if (cif::iequals(mGroup, "DNA"))
		result = "DNA linking";
	else if (cif::iequals(mGroup, "RNA"))
		result = "RNA linking";
	
	return result;
}

bool Compound::isWater() const
{
	return mId == "HOH" or mId == "H2O";
}

CompoundAtom Compound::getAtomById(const string& atomId) const
{
	CompoundAtom result;
	for (auto& a: mAtoms)
	{
		if (a.id == atomId)
		{
			result = a;
			break;
		}
	}

	if (result.id != atomId)	
		throw out_of_range("No atom " + atomId + " in Compound " + mId);
	
	return result;
}

const Compound* Compound::create(const string& id)
{
	return CompoundFactory::instance().create(id);
}

// --------------------------------------------------------------------
// a factory class to generate compounds

CompoundFactory* CompoundFactory::sInstance = nullptr;

CompoundFactory::CompoundFactory()
{
	const char* clibdMon = getenv("CLIBD_MON");
	if (clibdMon == nullptr)
		throw runtime_error("Cannot locate peptide list, please souce the CCP4 environment");
	mClibdMon = clibdMon;
}

CompoundFactory::~CompoundFactory()
{
}

CompoundFactory& CompoundFactory::instance()
{
	if (sInstance == nullptr)
		sInstance = new CompoundFactory();
	return *sInstance;
}

// id is the three letter code
const Compound* CompoundFactory::create(std::string id)
{
	ba::to_upper(id);

	Compound* result = nullptr;
	
	for (auto cmp: mCompounds)
	{
		if (cmp->id() == id)
		{
			result = cmp;
			break;
		}
	}
	
	if (result == nullptr)
	{
		fs::path resFile = mClibdMon / ba::to_lower_copy(id.substr(0, 1)) / (id + ".cif");
		fs::ifstream file(resFile);
		if (file.is_open())
		{
			cif::File cf;
			
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
			auto row = list["chem_comp"][cif::Key("id") == id];
			
			string name, group;
			uint32 numberAtomsAll, numberAtomsNh;
			cif::tie(name, group, numberAtomsAll, numberAtomsNh) =
				row.get("name", "group", "number_atoms_all", "number_atoms_nh");
	
			ba::trim(name);
			ba::trim(group);
			
			auto& compoundAtoms = cf["comp_" + id]["chem_comp_atom"];
			
			vector<CompoundAtom> atoms;
			for (auto row: compoundAtoms)
			{
				string id, symbol, energy;
				float charge;
				
				cif::tie(id, symbol, energy, charge) = row.get("atom_id", "type_symbol", "type_energy", "partial_charge");
				
				atoms.push_back({
					id, AtomTypeTraits(symbol).type(), energy, charge
				});
			}

			auto& compBonds = cf["comp_" + id]["chem_comp_bond"];
			
			vector<CompoundBond> bonds;
			for (auto row: compBonds)
			{
				CompoundBond b;
				string type, aromatic;
				
				cif::tie(b.atomID[0], b.atomID[1], type, b.distance, aromatic) =
					row.get("atom_id_1", "atom_id_2", "type", "distance", "aromatic");
				
				using cif::iequals;
				
				if (iequals(type, "single"))		b.type = singleBond;
				else if (iequals(type, "double"))	b.type = doubleBond;
				else if (iequals(type, "triple"))	b.type = tripleBond;
				else if (iequals(type, "deloc") or iequals(type, "aromat") or iequals(type, "aromatic"))
													b.type = delocalizedBond;
				else
				{
					if (VERBOSE)
						cerr << "Unimplemented chem_comp_bond.type " << type << " in file " << resFile << endl;
					b.type = singleBond;
				}
				
				bonds.push_back(b);
			}

			auto& compChir = cf["comp_" + id]["chem_comp_chir"];
			
			vector<ChiralCentre> chiralCentres;
			for (auto row: compChir)
			{
				ChiralCentre cc;
				string volumeSign;
				
				cif::tie(cc.id, cc.atomIDCentre, cc.atomID[0],
					cc.atomID[1], cc.atomID[2], volumeSign) = 
					row.get("id", "atom_id_centre", "atom_id_1",
						"atom_id_2", "atom_id_3", "volume_sign");
				
				if (volumeSign == "negativ")
					cc.volumeSign = negativ;
				else if (volumeSign == "positiv")
					cc.volumeSign = positiv;
				else if (volumeSign == "both")
					cc.volumeSign = both;
				else
				{
					if (VERBOSE)
						cerr << "Unimplemented chem_comp_chir.volume_sign " << volumeSign << " in file " << resFile << endl;
					continue;
				}
				
				chiralCentres.push_back(cc);
			}
			
			result = new Compound(id, name, group, move(atoms), move(bonds),
				move(chiralCentres));
			mCompounds.push_back(result);
		}
	}
	
	return result;
}

bool Compound::atomsBonded(const string& atomId_1, const string& atomId_2) const
{
	auto i = find_if(mBonds.begin(), mBonds.end(),
		[&](const CompoundBond& b)
		{
			return (b.atomID[0] == atomId_1 and b.atomID[1] == atomId_2)
				or (b.atomID[0] == atomId_2 and b.atomID[1] == atomId_1);
		});
	
	return i != mBonds.end();
}

float Compound::atomBondValue(const string& atomId_1, const string& atomId_2) const
{
	auto i = find_if(mBonds.begin(), mBonds.end(),
		[&](const CompoundBond& b)
		{
			return (b.atomID[0] == atomId_1 and b.atomID[1] == atomId_2)
				or (b.atomID[0] == atomId_2 and b.atomID[1] == atomId_1);
		});
	
	return i != mBonds.end() ? i->distance : 0;
}

}
