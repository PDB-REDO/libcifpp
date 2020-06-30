// Lib for working with structures as contained in mmCIF and PDB files

#include "cif++/Config.h"

#include <map>
#include <numeric>
#include <shared_mutex>

#include <boost/algorithm/string.hpp>

#if __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#elif __has_include(<boost/filesystem.hpp>)
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#endif

#include "cif++/Cif++.h"
#include "cif++/Point.h"
#include "cif++/Compound.h"
#include "cif++/CifUtils.h"

using namespace std;
namespace ba = boost::algorithm;

namespace mmcif
{

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
// Isomers in the ccp4 dictionary

struct IsomerSet
{
	vector<string>	compounds;
};

struct IsomerSets
{
	vector<IsomerSet>	isomers;
};

class IsomerDB
{
  public:
	
	static IsomerDB& instance();

	size_t count(const string& compound) const;
	const vector<string>& operator[](const string& compound) const;
	
  private:
	
	IsomerDB();
	
	IsomerSets			mData;
};

IsomerDB::IsomerDB()
{
	auto isomers = cif::rsrc_loader::load("isomers.txt");
	if (not isomers)
		throw runtime_error("Missing isomers.txt resource");
	
	struct membuf : public streambuf
	{
		membuf(char* data, size_t length)		{ this->setg(data, data, data + length); }
	} buffer(const_cast<char*>(isomers.data()), isomers.size());
	
	istream is(&buffer);
// #else
// 	cerr << "resource support was not compiled in, falling back to a local file" << endl;
// 	fs::path isomersFile = "isomers.txt";
// 	if (not fs::exists(isomersFile))
// 		isomersFile = fs::path(cif::get_executable_path()).parent_path() / "isomers.txt";

// 	fs::ifstream is(isomersFile);
// 	if (not is.is_open())
// 		throw runtime_error("Could not open the file isomers.txt");
// #endif

	string line;

	while (getline(is, line))
	{
		IsomerSet compounds;
		ba::split(compounds.compounds, line, ba::is_any_of(":"));

		if (not compounds.compounds.empty())
			mData.isomers.emplace_back(std::move(compounds));
	}
}

IsomerDB& IsomerDB::instance()
{
	static IsomerDB sInstance;
	return sInstance;
}

size_t IsomerDB::count(const string& compound) const
{
	size_t n = 0;
	for (auto& d: mData.isomers)
	{
		if (find(d.compounds.begin(), d.compounds.end(), compound) != d.compounds.end())
			++n;
	}
	return n;
}

const vector<string>& IsomerDB::operator[](const string& compound) const
{
	for (auto& d: mData.isomers)
	{
		if (find(d.compounds.begin(), d.compounds.end(), compound) != d.compounds.end())
			return d.compounds;
	}
	
	throw runtime_error("No isomer set found containing " + compound);
}

// --------------------------------------------------------------------
// Brute force comparison of two structures, when they are isomers the
// mapping between the atoms of both is returned
// This is not an optimal solution, but it works good enough for now

struct Node
{
	string							id;
	AtomType						symbol;
	vector<tuple<size_t,BondType>>	links;
	size_t							hydrogens = 0;
};

// Check to see if the nodes a[iA] and b[iB] are the start of a similar sub structure
bool SubStructuresAreIsomeric(
	const vector<Node>& a, const vector<Node>& b, size_t iA, size_t iB,
	vector<bool> visitedA, vector<bool> visitedB, vector<tuple<string,string>>& outMapping)
{
	auto& na = a[iA];
	auto& nb = b[iB];
	size_t N = na.links.size();
	
	assert(na.symbol == nb.symbol);
	assert(nb.links.size() == N);

	// we're optimistic today
	bool result = true;
	
	visitedA[iA] = true;
	visitedB[iB] = true;
	
	// we now have two sets of links to compare. 
	// Compare each permutation of the second set of indices with the first
	
	vector<size_t> ilb(N);
	iota(ilb.begin(), ilb.end(), 0);
	
	for (;;)
	{
		result = true;
		vector<tuple<string,string>> m;
		
		for (size_t i = 0; result and i < N; ++i)
		{
			size_t lA, lB;
			BondType typeA, typeB;

			tie(lA, typeA) = na.links[i];			assert(lA < a.size());
			tie(lB, typeB) = nb.links[ilb[i]];		assert(lB < b.size());
			
			if (typeA != typeB or visitedA[lA] != visitedB[lB])
			{
				result = false;
				break;
			}

			auto& la = a[lA];
			auto& lb = b[lB];
			
			if (la.symbol != lb.symbol or la.hydrogens != lb.hydrogens or la.links.size() != lb.links.size())
			{
				result = false;
				break;
			}

			if (not visitedA[lA])
				result = SubStructuresAreIsomeric(a, b, lA, lB, visitedA, visitedB, m);
		}
		
		if (result)
		{
			outMapping.insert(outMapping.end(), m.begin(), m.end());
			break;
		}

		if (not next_permutation(ilb.begin(), ilb.end()))
			break;
	}
	
	if (result and na.id != nb.id)
		outMapping.emplace_back(na.id, nb.id);
	
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
		
		if (a[atom2].symbol == H)
			a[atom1].hydrogens += 1;
		else
			a[atom1].links.emplace_back(atom2, bondA.type);
			
		if (a[atom1].symbol == H)
			a[atom2].hydrogens += 1;
		else
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
		
		if (b[atom2].symbol == H)
			b[atom1].hydrogens += 1;
		else
			b[atom1].links.emplace_back(atom2, bondB.type);
		
		if (b[atom1].symbol == H)
			b[atom2].hydrogens += 1;
		else
			b[atom2].links.emplace_back(atom1, bondB.type);
	}

	size_t N = atomsA.size();
	size_t ia = 0;
	
	bool result = false;

	// try each atom in B to see if it can be traced to be similar to A starting at zero
	for (size_t ib = 0; ib < N; ++ib)
	{
		if (b[ib].symbol != a[ia].symbol or a[ia].hydrogens != b[ib].hydrogens or a[ia].links.size() != b[ib].links.size())
			continue;
		
		vector<bool> va(N, false), vb(N, false);

		if (SubStructuresAreIsomeric(a, b, ia, ib, va, vb, outMapping))
		{
			result = true;
			break;
		}
	}
	
	return result;
}

// --------------------------------------------------------------------
// Compound

Compound::Compound(const std::string& file, const std::string& id,
	const std::string& name, const std::string& group)
	: mID(id), mName(name), mGroup(group)
{
	try
	{
		mCF.load(file);
		
		// locate the datablock
		auto& db = mCF["comp_" + id];
		
		auto& compoundAtoms = db["chem_comp_atom"];
		
		for (auto row: compoundAtoms)
		{
			string id, symbol, energy;
			float charge;
			
			cif::tie(id, symbol, energy, charge) = row.get("atom_id", "type_symbol", "type_energy", "partial_charge");
			
			mAtoms.push_back({
				id, AtomTypeTraits(symbol).type(), energy, charge
			});
		}
		sort(mAtoms.begin(), mAtoms.end(), CompoundAtomLess());

		auto& compBonds = db["chem_comp_bond"];
		
		for (auto row: compBonds)
		{
			CompoundBond b;
			string type, aromatic;
			
			cif::tie(b.atomID[0], b.atomID[1], type, b.distance, b.esd) =
				row.get("atom_id_1", "atom_id_2", "type", "value_dist", "value_dist_esd");
			
			using cif::iequals;
			
			if (iequals(type, "single") or iequals(type, "sing"))		b.type = singleBond;
			else if (iequals(type, "double") or iequals(type, "doub"))	b.type = doubleBond;
			else if (iequals(type, "triple") or iequals(type, "trip"))	b.type = tripleBond;
			else if (iequals(type, "deloc") or iequals(type, "aromat") or iequals(type, "aromatic"))
				b.type = delocalizedBond;
			else
			{
				if (cif::VERBOSE)
					cerr << "Unimplemented chem_comp_bond.type " << type << " in " << id << endl;
				b.type = singleBond;
			}
			
			if (b.atomID[0] > b.atomID[1])
				swap(b.atomID[0], b.atomID[1]);
			
			mBonds.push_back(b);
		}
		sort(mBonds.begin(), mBonds.end(), CompoundBondLess());

		for (auto row: db["chem_comp_angle"])
		{
			CompoundAngle a;
			
			cif::tie(a.atomID[0], a.atomID[1], a.atomID[2], a.angle, a.esd) =
				row.get("atom_id_1", "atom_id_2", "atom_id_3", "value_angle", "value_angle_esd");
			
			mAngles.push_back(a);
		}

		for (auto row: db["chem_comp_tor"])
		{
			CompoundTorsion a;
			
			cif::tie(a.atomID[0], a.atomID[1], a.atomID[2], a.atomID[3], a.angle, a.esd, a.period) =
				row.get("atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4", "value_angle", "value_angle_esd", "period");
			
			mTorsions.push_back(a);
		}

		for (auto row: db["chem_comp_chir"])
		{
			CompoundChiralCentre cc;
			string volumeSign;
			
			cif::tie(cc.id, cc.atomIDCentre, cc.atomID[0],
				cc.atomID[1], cc.atomID[2], volumeSign) = 
				row.get("id", "atom_id_centre", "atom_id_1",
					"atom_id_2", "atom_id_3", "volume_sign");
			
			if (volumeSign == "negativ" or volumeSign == "negative")
				cc.volumeSign = negativ;
			else if (volumeSign == "positiv" or volumeSign == "positive")
				cc.volumeSign = positiv;
			else if (volumeSign == "both")
				cc.volumeSign = both;
			else
			{
				if (cif::VERBOSE)
					cerr << "Unimplemented chem_comp_chir.volume_sign " << volumeSign << " in " << id << endl;
				continue;
			}
			
			mChiralCentres.push_back(cc);
		}

		auto& compPlanes = db["chem_comp_plane_atom"];
		
		for (auto row: compPlanes)
		{
			string atom_id, plane_id;
			float esd;	
			
			cif::tie(atom_id, plane_id, esd) = row.get("atom_id", "plane_id", "dist_esd");

			auto i = find_if(mPlanes.begin(), mPlanes.end(), [&](auto& p) { return p.id == plane_id;});
			if (i == mPlanes.end())
				mPlanes.emplace_back(CompoundPlane{ plane_id, { atom_id }, esd });
			else
				i->atomID.push_back(atom_id);
		}
	}
	catch (const exception& ex)
	{
		cerr << "Error loading ccp4 file for " << id << " from file " << file << endl;
		throw;
	}
}

string Compound::formula() const
{
	string result;
	
	map<string,uint32_t> atoms;
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

float Compound::formulaWeight() const
{
	float result = 0;

	for (auto r: mAtoms)
		result += AtomTypeTraits(r.typeSymbol).weight();

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
	
	if (cif::iequals(mID, "gly"))
		result = "peptide linking";
	else if (cif::iequals(mGroup, "l-peptide") or cif::iequals(mGroup, "L-peptide linking") or cif::iequals(mGroup, "peptide"))
		result = "L-peptide linking";
	else if (cif::iequals(mGroup, "DNA"))
		result = "DNA linking";
	else if (cif::iequals(mGroup, "RNA"))
		result = "RNA linking";
//	else
//		result = mGroup;
	
	return result;
}

bool Compound::isWater() const
{
	return mID == "HOH" or mID == "H2O";
}

bool Compound::isSugar() const
{
	return cif::iequals(mGroup, "furanose") or cif::iequals(mGroup, "pyranose");
}

CompoundAtom Compound::getAtomByID(const string& atomID) const
{
	CompoundAtom result = {};
	for (auto& a: mAtoms)
	{
		if (a.id == atomID)
		{
			result = a;
			break;
		}
	}

	if (result.id != atomID)	
		throw out_of_range("No atom " + atomID + " in Compound " + mID);
	
	return result;
}

const Compound* Compound::create(const string& id)
{
	auto result = CompoundFactory::instance().get(id);
	if (result == nullptr)
		result = CompoundFactory::instance().create(id);
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
		if (mID == c.mID)
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
		
		if (cif::VERBOSE and result)
		{
			for (auto& m: mapping)
				cerr << "  " << get<0>(m) << " => " << get<1>(m) << endl;
		}
		
		break;	
	}
	
	return result;
}

vector<tuple<string,string>> Compound::mapToIsomer(const Compound& c) const
{
	vector<tuple<string,string>> result;
	
	bool check = StructuresAreIsomeric(mAtoms, mBonds, c.mAtoms, c.mBonds, result);
	if (not check)
		throw runtime_error("Compounds " + id() + " and " + c.id() + " are not isomers in call to mapToIsomer");
	
	return result;
}

vector<string> Compound::isomers() const
{
	vector<string> result;
	
	auto& db = IsomerDB::instance();
	if (db.count(mID))
	{
		result = db[mID];
		
		auto i = find(result.begin(), result.end(), mID);
		assert(i != result.end());
		
		result.erase(i);
	}
	
	return result;
}

float Compound::bondAngle(const string& atomId_1, const string& atomId_2, const string& atomId_3) const
{
	float result = nanf("1");
	
	for (auto& a: mAngles)
	{
		if (not (a.atomID[1] == atomId_2 and
				 ((a.atomID[0] == atomId_1 and a.atomID[2] == atomId_3) or
				  (a.atomID[2] == atomId_1 and a.atomID[0] == atomId_3))))
			continue;
		
		result = a.angle;
		break;
	}
	
	return result;
}

//static float calcC(float a, float b, float alpha)
//{
//	float f = b * sin(alpha * kPI / 180);
//	float d = sqrt(b * b - f * f);
//	float e = a - d;
//	float c = sqrt(f * f + e * e);
//	
//	return c;
//}

float Compound::chiralVolume(const string& centreID) const
{
	float result = 0;
	
	for (auto& cv: mChiralCentres)
	{
		if (cv.id != centreID)
			continue;
		
		// calculate the expected chiral volume

		// the edges
		
		float a = atomBondValue(cv.atomIDCentre, cv.atomID[0]);
		float b = atomBondValue(cv.atomIDCentre, cv.atomID[1]);
		float c = atomBondValue(cv.atomIDCentre, cv.atomID[2]);
		
		// the angles for the top of the tetrahedron
		
		float alpha = bondAngle(cv.atomID[0], cv.atomIDCentre, cv.atomID[1]);
		float beta  = bondAngle(cv.atomID[1], cv.atomIDCentre, cv.atomID[2]);
		float gamma = bondAngle(cv.atomID[2], cv.atomIDCentre, cv.atomID[0]);
		
		float cosa = cos(alpha * kPI / 180);
		float cosb = cos(beta * kPI / 180);
		float cosc = cos(gamma * kPI / 180);
		
		result = (a * b * c * sqrt(1 + 2 * cosa * cosb * cosc - (cosa * cosa) - (cosb * cosb) - (cosc * cosc))) / 6;
		
		if (cv.volumeSign == negativ)
			result = -result;
		
		break;
	}
	
	return result;
}

// --------------------------------------------------------------------

Link::Link(cif::Datablock& db)
{
	mID = db.getName();
	
	auto& linkBonds = db["chem_link_bond"];
	
	for (auto row: linkBonds)
	{
		LinkBond b;
		string type, aromatic;
		
		cif::tie(b.atom[0].compID, b.atom[0].atomID,
				b.atom[1].compID, b.atom[1].atomID, type, b.distance, b.esd) =
			row.get("atom_1_comp_id", "atom_id_1", "atom_2_comp_id", "atom_id_2",
				"type", "value_dist", "value_dist_esd");
		
		using cif::iequals;
		
		if (iequals(type, "single") or iequals(type, "sing"))		b.type = singleBond;
		else if (iequals(type, "double") or iequals(type, "doub"))	b.type = doubleBond;
		else if (iequals(type, "triple") or iequals(type, "trip"))	b.type = tripleBond;
		else if (iequals(type, "deloc") or iequals(type, "aromat") or iequals(type, "aromatic"))
											b.type = delocalizedBond;
		else
		{
			if (cif::VERBOSE)
				cerr << "Unimplemented chem_link_bond.type " << type << " in " << mID << endl;
			b.type = singleBond;
		}
		
//		if (b.atom[0] > b.atom[1])
//			swap(b.atom[0], b.atom[1]);
		
		mBonds.push_back(b);
	}
//	sort(mBonds.begin(), mBonds.end(), LinkBondLess());

	auto& linkAngles = db["chem_link_angle"];
	for (auto row: linkAngles)
	{
		LinkAngle a;
		
		cif::tie(a.atom[0].compID, a.atom[0].atomID, a.atom[1].compID, a.atom[1].atomID,
				a.atom[2].compID, a.atom[2].atomID, a.angle, a.esd) =
			row.get("atom_1_comp_id", "atom_id_1", "atom_2_comp_id", "atom_id_2",
				"atom_3_comp_id", "atom_id_3", "value_angle", "value_angle_esd");
		
		mAngles.push_back(a);
	}

	for (auto row: db["chem_link_tor"])
	{
		LinkTorsion a;
		
		cif::tie(a.atom[0].compID, a.atom[0].atomID, a.atom[1].compID, a.atom[1].atomID,
				a.atom[2].compID, a.atom[2].atomID, a.atom[3].compID, a.atom[3].atomID,
				a.angle, a.esd, a.period) =
			row.get("atom_1_comp_id", "atom_id_1", "atom_2_comp_id", "atom_id_2",
				"atom_3_comp_id", "atom_id_3", "atom_4_comp_id", "atom_id_4",
				"value_angle", "value_angle_esd", "period");
		
		mTorsions.push_back(a);
	}

	auto& linkChir = db["chem_link_chir"];
	for (auto row: linkChir)
	{
		LinkChiralCentre cc;
		string volumeSign;
		
		cif::tie(cc.id, cc.atomCentre.compID, cc.atomCentre.atomID,
				cc.atom[0].compID, cc.atom[0].atomID,
				cc.atom[1].compID, cc.atom[1].atomID,
				cc.atom[2].compID, cc.atom[2].atomID,
				volumeSign) = 
			row.get("id", "atom_centre_comp_id", "atom_id_centre",
				"atom_1_comp_id", "atom_id_1", "atom_2_comp_id", "atom_id_2",
				"atom_3_comp_id", "atom_id_3", "volume_sign");
		
		if (volumeSign == "negativ" or volumeSign == "negative")
			cc.volumeSign = negativ;
		else if (volumeSign == "positiv" or volumeSign == "positive")
			cc.volumeSign = positiv;
		else if (volumeSign == "both")
			cc.volumeSign = both;
		else
		{
			if (cif::VERBOSE)
				cerr << "Unimplemented chem_link_chir.volume_sign " << volumeSign << " in " << mID << endl;
			continue;
		}
		
		mChiralCentres.push_back(cc);
	}

	auto& linkPlanes = db["chem_link_plane"];
	for (auto row: linkPlanes)
	{
		int compID;
		string atomID, planeID;
		float esd;	
		
		cif::tie(planeID, compID, atomID, esd) = row.get("plane_id", "atom_comp_id", "atom_id", "dist_esd");

		auto i = find_if(mPlanes.begin(), mPlanes.end(), [&](auto& p) { return p.id == planeID;});
		if (i == mPlanes.end())
		{
			vector<LinkAtom> atoms{LinkAtom{ compID, atomID }};
			mPlanes.emplace_back(LinkPlane{ planeID, move(atoms), esd });
		}
		else
			i->atoms.push_back({ compID, atomID });
	}
}

const Link& Link::create(const std::string& id)
{
	auto result = CompoundFactory::instance().getLink(id);
	if (result == nullptr)
		result = CompoundFactory::instance().createLink(id);
	
	if (result == nullptr)
		throw runtime_error("Link with id " + id + " not found");
	
	return *result;
}

Link::~Link()
{
}

float Link::atomBondValue(const LinkAtom& atom1, const LinkAtom& atom2) const
{
	auto i = find_if(mBonds.begin(), mBonds.end(),
		[&](auto& b)
		{
			return (b.atom[0] == atom1 and b.atom[1] == atom2)
				or (b.atom[0] == atom2 and b.atom[1] == atom1);
		});
	
	return i != mBonds.end() ? i->distance : 0;
}

float Link::bondAngle(const LinkAtom& atom1, const LinkAtom& atom2, const LinkAtom& atom3) const
{
	float result = nanf("1");
	
	for (auto& a: mAngles)
	{
		if (not (a.atom[1] == atom2 and
				 ((a.atom[0] == atom1 and a.atom[2] == atom3) or
				  (a.atom[2] == atom1 and a.atom[0] == atom3))))
			continue;
		
		result = a.angle;
		break;
	}
	
	return result;
}

float Link::chiralVolume(const string& centreID) const
{
	float result = 0;
	
	for (auto& cv: mChiralCentres)
	{
		if (cv.id != centreID)
			continue;
		
		// calculate the expected chiral volume

		// the edges
		
		float a = atomBondValue(cv.atomCentre, cv.atom[0]);
		float b = atomBondValue(cv.atomCentre, cv.atom[1]);
		float c = atomBondValue(cv.atomCentre, cv.atom[2]);
		
		// the angles for the top of the tetrahedron
		
		float alpha = bondAngle(cv.atom[0], cv.atomCentre, cv.atom[1]);
		float beta  = bondAngle(cv.atom[1], cv.atomCentre, cv.atom[2]);
		float gamma = bondAngle(cv.atom[2], cv.atomCentre, cv.atom[0]);
		
		float cosa = cos(alpha * kPI / 180);
		float cosb = cos(beta * kPI / 180);
		float cosc = cos(gamma * kPI / 180);
		
		result = (a * b * c * sqrt(1 + 2 * cosa * cosb * cosc - (cosa * cosa) - (cosb * cosb) - (cosc * cosc))) / 6;
		
		if (cv.volumeSign == negativ)
			result = -result;
		
		break;
	}
	
	return result;
}

// --------------------------------------------------------------------
// a factory class to generate compounds

const map<string,char> kAAMap{
	{ "ALA", 'A' },
	{ "ARG", 'R' },
	{ "ASN", 'N' },
	{ "ASP", 'D' },
	{ "CYS", 'C' },
	{ "GLN", 'Q' },
	{ "GLU", 'E' },
	{ "GLY", 'G' },
	{ "HIS", 'H' },
	{ "ILE", 'I' },
	{ "LEU", 'L' },
	{ "LYS", 'K' },
	{ "MET", 'M' },
	{ "PHE", 'F' },
	{ "PRO", 'P' },
	{ "SER", 'S' },
	{ "THR", 'T' },
	{ "TRP", 'W' },
	{ "TYR", 'Y' },
	{ "VAL", 'V' },
	{ "GLX", 'Z' },
	{ "ASX", 'B' }
};

const map<string,char> kBaseMap{
	{ "A", 'A' },
	{ "C", 'C' },
	{ "G", 'G' },
	{ "T", 'T' },
	{ "U", 'U' },
	{ "DA", 'A' },
	{ "DC", 'C' },
	{ "DG", 'G' },
	{ "DT", 'T' }
};

// --------------------------------------------------------------------

class CompoundFactoryImpl
{
  public:
	CompoundFactoryImpl(const std::string& file, CompoundFactoryImpl* next);

	~CompoundFactoryImpl()
	{
		delete mNext;
	}

	Compound* get(string id);
	Compound* create(string id);

	const Link* getLink(std::string id);
	const Link* createLink(std::string id);
	
	CompoundFactoryImpl* pop()
	{
		auto result = mNext;
		mNext = nullptr;
		delete this;
		return result;
	}

	string unalias(const string& resName) const
	{
		string result = resName;
		
		auto& e = const_cast<cif::File&>(mFile)["comp_synonym_list"];
		
		for (auto& synonym: e["chem_comp_synonyms"])
		{
			if (ba::iequals(synonym["comp_alternative_id"].as<string>(), resName) == false)
				continue;
			
			result = synonym["comp_id"].as<string>();
			ba::trim(result);
			break;
		}

		if (result.empty() and mNext)
			result = mNext->unalias(resName);
		
		return result;
	}

	bool isKnownPeptide(const string& resName)
	{
		return mKnownPeptides.count(resName) or
			(mNext != nullptr and mNext->isKnownPeptide(resName));
	}

	bool isKnownBase(const string& resName)
	{
		return mKnownBases.count(resName) or
			(mNext != nullptr and mNext->isKnownBase(resName));
	}

  private:
	std::shared_timed_mutex		mMutex;

	std::string					mPath;
	vector<Compound*>			mCompounds;
	vector<const Link*>			mLinks;
	/*unordered_*/set<string>	mKnownPeptides;
	set<string>					mKnownBases;
	set<string>					mMissing;
	cif::File					mFile;
	cif::Category&				mChemComp;
	CompoundFactoryImpl*		mNext;
};

// --------------------------------------------------------------------

CompoundFactoryImpl::CompoundFactoryImpl(const std::string& file, CompoundFactoryImpl* next)
	: mPath(file), mFile(file), mChemComp(mFile.firstDatablock()["chem_comp"]), mNext(next)
{
	const std::regex peptideRx("(?:[lmp]-)?peptide", std::regex::icase);

	for (auto& chemComp: mChemComp)
	{
		string group, threeLetterCode;
		
		cif::tie(group, threeLetterCode) = chemComp.get("group", "three_letter_code");

		if (std::regex_match(group, peptideRx))
			mKnownPeptides.insert(threeLetterCode);
		else if (ba::iequals(group, "DNA") or ba::iequals(group, "RNA"))
			mKnownBases.insert(threeLetterCode);
	}
}

//struct CompoundFactoryImpl
//{
//	CompoundFactoryImpl();
//
//	Compound* get(string id);
//	Compound* create(string id);
//	
//	void pushDictionary(const boost::filesystem::path& dictFile);
//
//	// cif::File cf acts as owner of datablock
//	cif::Datablock& getDatablock(const string& name, cif::File& cf);
//	
//	deque<cif::File> mDictionaries;
//	fs::path mClibdMon;
//	vector<Compound*> mCompounds;
//	boost::shared_mutex mMutex;
//};
//
//CompoundFactoryImpl::CompoundFactoryImpl()
//{
//	const char* clibdMon = getenv("CLIBD_MON");
//	if (clibdMon == nullptr)
//		throw runtime_error("Cannot locate peptide list, please source the CCP4 environment");
//	mClibdMon = clibdMon;
//}
//
//void CompoundFactoryImpl::pushDictionary(const boost::filesystem::path& dictFile)
//{
//	boost::upgrade_lock<boost::shared_mutex> lock(mMutex);
//
//	mDictionaries.emplace_back(dictFile);
//}
//
//cif::Datablock& CompoundFactoryImpl::getDatablock(const string& name, cif::File& cf)
//{
//	cif::Datablock* result = nullptr;
//	
//	// first check all loaded dictionaries
//	
//	for (auto d: mDictionaries)
//	{
//		result = d.get(name);
//		
//		if (result != nullptr)
//			break;
//	}	
//	
//	if (not result)
//	{
//		fs::path resFile = mClibdMon / ba::to_lower_copy(name.substr(0, 1)) / (name + ".cif");
//	
//		if (not fs::exists(resFile) and	(name == "COM" or name == "CON" or "PRN")) 		// seriously...
//			mClibdMon / ba::to_lower_copy(name.substr(0, 1)) / (name + '_' + name + ".cif");
//	
//		fs::ifstream file(resFile);
//		if (file.is_open())
//		{
//			try
//			{
//				cf.load(file);
//			}
//			catch (...)
//			{
//				throw_with_nested(runtime_error("Error while loading " + resFile));
//			}
//			
//			result = cf.get(name);
//		}
//	}
//
//	if (result == nullptr)
//		throw runtime_error("Failed to load datablock " + name);	
//
//	return *result;
//}

Compound* CompoundFactoryImpl::get(string id)
{
	std::shared_lock lock(mMutex);

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
	
	if (result == nullptr and mNext != nullptr)
		result = mNext->get(id);
	
	return result;
}

Compound* CompoundFactoryImpl::create(std::string id)
{
	ba::to_upper(id);

	Compound* result = get(id);
	if (result == nullptr and mMissing.count(id) == 0)
	{
		std::unique_lock lock(mMutex);

		auto rs = mChemComp.find(cif::Key("three_letter_code") == id);
	
		if (not rs.empty())
		{
			auto row = rs.front();

			string name, group;
			uint32_t numberAtomsAll, numberAtomsNh;
			cif::tie(name, group, numberAtomsAll, numberAtomsNh) =
				row.get("name", "group", "number_atoms_all", "number_atoms_nh");
			
			ba::trim(name);
			ba::trim(group);

			if (mFile.get("comp_" + id) == nullptr)
			{
				auto clibd_mon = fs::path(getenv("CLIBD_MON"));
				
				fs::path resFile = clibd_mon / ba::to_lower_copy(id.substr(0, 1)) / (id + ".cif");
			
				if (not fs::exists(resFile) and	(id == "COM" or id == "CON" or "PRN")) 		// seriously...
					resFile = clibd_mon / ba::to_lower_copy(id.substr(0, 1)) / (id + '_' + id + ".cif");

				if (not fs::exists(resFile))
					mMissing.insert(id);
				else
				{
					mCompounds.push_back(new Compound(resFile, id, name, group));
					result = mCompounds.back();
				}
			}				
			else
			{
				mCompounds.push_back(new Compound(mPath, id, name, group));
				result = mCompounds.back();
			}
		}
		
		if (result == nullptr and mNext != nullptr)
			result = mNext->create(id);
	}
	
	return result;
}

const Link* CompoundFactoryImpl::getLink(string id)
{
	std::shared_lock lock(mMutex);

	ba::to_upper(id);

	const Link* result = nullptr;
	
	for (auto link: mLinks)
	{
		if (link->id() == id)
		{
			result = link;
			break;
		}
	}
	
	if (result == nullptr and mNext != nullptr)
		result = mNext->getLink(id);
	
	return result;
}

const Link* CompoundFactoryImpl::createLink(std::string id)
{
	ba::to_upper(id);

	const Link* result = getLink(id);

	if (result == nullptr)
	{
		std::unique_lock lock(mMutex);

		auto db = mFile.get("link_" + id);
		if (db != nullptr)
		{
			result = new Link(*db);
			mLinks.push_back(result);
		}
		
		if (result == nullptr and mNext != nullptr)
			result = mNext->createLink(id);
	}
	
	return result;
}

// --------------------------------------------------------------------

CompoundFactory* CompoundFactory::sInstance;

CompoundFactory::CompoundFactory()
	: mImpl(nullptr)
{
	const char* clibdMon = getenv("CLIBD_MON");
	if (clibdMon == nullptr)
		throw runtime_error("Cannot locate peptide list, please source the CCP4 environment");
	
	fs::path db = fs::path(clibdMon) / "list" / "mon_lib_list.cif";

	pushDictionary(db.string());
}

CompoundFactory::~CompoundFactory()
{
	delete mImpl;
}

CompoundFactory& CompoundFactory::instance()
{
	if (sInstance == nullptr)
		sInstance = new CompoundFactory();
	return *sInstance;
}

void CompoundFactory::pushDictionary(const string& inDictFile)
{
	if (not fs::exists(inDictFile))
		throw runtime_error("file not found: " + inDictFile);

//	ifstream file(inDictFile);
//	if (not file.is_open())
//		throw runtime_error("Could not open peptide list " + inDictFile);

	try
	{
		mImpl = new CompoundFactoryImpl(inDictFile, mImpl);
	}
	catch (const exception& ex)
	{
		cerr << "Error loading dictionary " << inDictFile << endl;
		throw;
	}
}

void CompoundFactory::popDictionary()
{
	if (mImpl != nullptr)
		mImpl = mImpl->pop();
}

// id is the three letter code
const Compound* CompoundFactory::get(std::string id)
{
	return mImpl->get(id);
}

const Compound* CompoundFactory::create(std::string id)
{
	return mImpl->create(id);
}

const Link* CompoundFactory::getLink(std::string id)
{
	return mImpl->getLink(id);
}

const Link* CompoundFactory::createLink(std::string id)
{
	return mImpl->createLink(id);
}

bool CompoundFactory::isKnownPeptide(const string& resName) const
{
	return mImpl->isKnownPeptide(resName);
}

bool CompoundFactory::isKnownBase(const string& resName) const
{
	return mImpl->isKnownBase(resName);
}

string CompoundFactory::unalias(const string& resName) const
{
	return mImpl->unalias(resName);
}

}
