// Lib for working with structures as contained in mmCIF and PDB files

#include "cif++/Config.h"

#include <map>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/thread.hpp>

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

	const Compound* get(string id);
	const Compound* create(string id);

  private:
	CompoundFactory();
	~CompoundFactory();
	
	static CompoundFactory* sInstance;

	fs::path mClibdMon;
	vector<Compound*> mCompounds;
	boost::shared_mutex mMutex;
};

// --------------------------------------------------------------------
// Compound helper classes

struct CompoundAtomLess
{
	bool operator()(const CompoundAtom& a, const CompoundAtom& b) const
	{
		int d = a.id.compare(b.id);
		if (d == 0)
			d = a.typeSymbol - b.typeSymbol;
		return d < 0;
	}
};

struct CompoundBondLess
{
	bool operator()(const CompoundBond& a, const CompoundBond& b) const
	{
		int d = a.atomID[0].compare(b.atomID[0]);
		if (d == 0)
			d = a.atomID[1].compare(b.atomID[1]);
		if (d == 0)
			d = a.type - b.type;
		return d < 0;
	}
};

// --------------------------------------------------------------------
// Brute force comparison of two structures, when they are isomers the
// mapping between the atoms of both is returned

struct Node
{
	string									id;
	AtomType								symbol;
	vector<tuple<size_t,CompoundBondType>>	links;
};

// Check to see if the nodes a[iA] and b[iB] are the start of a similar sub structure
bool SubStructuresAreIsomeric(
	const vector<Node>& a, const vector<Node>& b, size_t iA, size_t iB,
	vector<bool> visitedA, vector<bool> visitedB, vector<tuple<string,string>>& outMapping)
{
	bool result = false;

	auto& na = a[iA];
	auto& nb = b[iB];
	size_t N = na.links.size();
	
	if (na.symbol == nb.symbol and nb.links.size() == N)
	{
		result = true;
		
		visitedA[iA] = true;
		visitedB[iB] = true;
		
		// we now have two sets of links to compare. 
		// To keep code clean, first create the list of possible permutations
		
		vector<size_t> ilb(N);
		iota(ilb.begin(), ilb.end(), 0);
		
		for (;;)
		{
			result = true;
			vector<tuple<string,string>> m = outMapping;
			
			for (size_t i = 0; result and i < N; ++i)
			{
				size_t lA, lB;
				CompoundBondType typeA, typeB;

				tie(lA, typeA) = na.links[i];
				tie(lB, typeB) = nb.links[ilb[i]];
				
				if (typeA != typeB or visitedA[lA] != visitedB[lB])
					result = false;
				else if (not visitedA[lA])
					result = SubStructuresAreIsomeric(a, b, lA, lB, visitedA, visitedB, m);
			}
			
			if (result)
			{
				swap(m, outMapping);
				break;
			}

			if (not next_permutation(ilb.begin(), ilb.end()))
				break;
		}
		
		if (result)
			outMapping.emplace_back(na.id, nb.id);
	}
	
	return result;
}

bool StructuresAreIsomeric(vector<CompoundAtom> atomsA, const vector<CompoundBond>& bondsA,
	vector<CompoundAtom> atomsB, const vector<CompoundBond>& bondsB,
	vector<tuple<string,string>>& outMapping)
{
	assert(atomsA.size() == atomsB.size());
	assert(bondsA.size() == bondsB.size());

	vector<Node> a, b;
	map<string,size_t> ma, mb;
	
	for (auto& atomA: atomsA)
	{
		ma[atomA.id] = a.size();
		a.push_back({atomA.id, atomA.typeSymbol});
	}
	
	for (auto& bondA: bondsA)
	{
		size_t atom1 = ma.at(bondA.atomID[0]);
		size_t atom2 = ma.at(bondA.atomID[1]);
		
		a[atom1].links.emplace_back(atom2, bondA.type);
		a[atom2].links.emplace_back(atom1, bondA.type);
	}

	for (auto& atomB: atomsB)
	{
		mb[atomB.id] = b.size();
		b.push_back({atomB.id, atomB.typeSymbol});
	}

	for (auto& bondB: bondsB)
	{
		size_t atom1 = mb.at(bondB.atomID[0]);
		size_t atom2 = mb.at(bondB.atomID[1]);
		
		b[atom1].links.emplace_back(atom2, bondB.type);
		b[atom2].links.emplace_back(atom1, bondB.type);
	}

	size_t N = atomsA.size();
	
	bool result = false;

	// try each atom in B to see if it can be traced to be similar to A starting at zero
	for (size_t ib = 0; ib < N; ++ib)
	{
		vector<bool> va(N, false), vb(N, false);

		if (SubStructuresAreIsomeric(a, b, 0, ib, va, vb, outMapping))
		{
			result = true;
			break;
		}
	}
	
	return result;
}

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
	auto result = CompoundFactory::instance().get(id);
	if (result == nullptr)
		result = CompoundFactory::instance().create(id);
	return result;
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

const Compound* CompoundFactory::get(std::string id)
{
	boost::shared_lock<boost::shared_mutex> lock(mMutex);

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
	
	return result;
}

const Compound* CompoundFactory::create(std::string id)
{
	boost::upgrade_lock<boost::shared_mutex> lock(mMutex);

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

		if (not fs::exists(resFile) and	(id == "COM" or id == "CON" or "PRN")) 		// seriously...
			mClibdMon / ba::to_lower_copy(id.substr(0, 1)) / (id + '_' + id + ".cif");

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
			sort(atoms.begin(), atoms.end(), CompoundAtomLess());

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
				
				if (b.atomID[0] > b.atomID[1])
					swap(b.atomID[0], b.atomID[1]);
				
				bonds.push_back(b);
			}
			sort(bonds.begin(), bonds.end(), CompoundBondLess());

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

			boost::upgrade_to_unique_lock<boost::shared_mutex> uniqueLock(lock);
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

bool Compound::isIsomerOf(const Compound& c) const
{
	bool result = false;
	
	for (;;)
	{
		// easy tests first
		if (mId == c.mId)
		{
			result = true;
			break;
		}
		
		if (mAtoms.size() != c.mAtoms.size())
			break;
		
		if (mBonds.size() != c.mBonds.size())
			break;
		
		if (mChiralCentres.size() != c.mChiralCentres.size())
			break;

		// same number of atoms of each type?
		map<AtomType,int> aTypeCount, bTypeCount;
		
		bool sameAtomNames = true;
		for (size_t i = 0; i < mAtoms.size(); ++i)
		{
			auto& a = mAtoms[i];
			auto& b = c.mAtoms[i];
			
			aTypeCount[a.typeSymbol] += 1;
			bTypeCount[b.typeSymbol] += 1;
			
			if (a.id != b.id or a.typeSymbol != b.typeSymbol)
				sameAtomNames = false;
		}
		
		if (not sameAtomNames and aTypeCount != bTypeCount)
			break;
		
		bool sameBonds = sameAtomNames;
		for (size_t i = 0; sameBonds and i < mBonds.size(); ++i)
		{
			sameBonds =
				mBonds[i].atomID[0] == c.mBonds[i].atomID[0] and
				mBonds[i].atomID[1] == c.mBonds[i].atomID[1] and
				mBonds[i].type == c.mBonds[i].type;
		}
		
		if (sameBonds)
		{
			result = true;
			break;
		}
		
		// implement rest of tests
		
		vector<tuple<string,string>> mapping;
		result = StructuresAreIsomeric(mAtoms, mBonds, c.mAtoms, c.mBonds, mapping);
		
		if (VERBOSE and result)
		{
			for (auto& m: mapping)
				cerr << "  " << get<0>(m) << " => " << get<1>(m) << endl;
		}
		
		break;	
	}
	
	return result;
}

}
