// copyright

#include "cif++/Config.h"

#include "cif++/Cif++.h"
#include "cif++/BondMap.h"
#include "cif++/Compound.h"
#include "cif++/CifUtils.h"

using namespace std;

namespace mmcif
{

// --------------------------------------------------------------------

BondMap::BondMap(const Structure& p)
{
	auto atoms = p.atoms();
	dim = atoms.size();

//	bond = vector<bool>(dim * (dim - 1), false);

	for (auto& atom: atoms)
	{
		size_t ix = index.size();
		index[atom.id()] = ix;
	};
	
	auto bindAtoms = [this](const string& a, const string& b)
	{
		uint32 ixa = index[a];
		uint32 ixb = index[b];
		
		bond.insert(key(ixa, ixb));
	};

	auto linkAtoms = [this,&bindAtoms](const string& a, const string& b)
	{
		bindAtoms(a, b);

		link[a].insert(b);
		link[b].insert(a);
	};

	cif::Datablock& db = p.getFile().data();

	// collect all compounds first
	set<string> compounds;
	for (auto c: db["chem_comp"])
		compounds.insert(c["id"].as<string>());
	
	// make sure we also have all residues in the polyseq
	for (auto m: db["entity_poly_seq"])
	{
		string c = m["mon_id"].as<string>();
		if (compounds.count(c))
			continue;
		
		if (VERBOSE)
			cerr << "Warning: mon_id " << c << " is missing in the chem_comp category" << endl;
		compounds.insert(c);
	}

	// first link all residues in a polyseq
	
	string lastAsymID;
	int lastSeqID = 0;
	for (auto r: db["pdbx_poly_seq_scheme"])
	{
		string asymId;
		int seqId;

		cif::tie(asymId, seqId) = r.get("asym_id", "seq_id");

		if (asymId != lastAsymID)		// first in a new sequece
		{
			lastAsymID = asymId;
			lastSeqID = seqId;
			continue;
		}
		
		auto c = db["atom_site"].find(cif::Key("label_asym_id") == asymId and cif::Key("label_seq_id") == lastSeqID and cif::Key("label_atom_id") == "C");
		if (c.size() != 1 and VERBOSE > 1)
			cerr << "Unexpected number (" << c.size() << ") of atoms with atom ID C in asym_id " << asymId << " with seq id " << lastSeqID << endl;
		
		auto n = db["atom_site"].find(cif::Key("label_asym_id") == asymId and cif::Key("label_seq_id") == seqId and cif::Key("label_atom_id") == "N");
		if (n.size() != 1 and VERBOSE > 1)
			cerr << "Unexpected number (" << n.size() << ") of atoms with atom ID N in asym_id " << asymId << " with seq id " << seqId << endl;
		
		if (not (c.empty() or n.empty()))
			bindAtoms(c.front()["id"].as<string>(), n.front()["id"].as<string>());
		
		lastSeqID = seqId;
	}

	for (auto l: db["struct_conn"])
	{
		string asym1, asym2, atomId1, atomId2;
		int seqId1, seqId2;
		cif::tie(asym1, asym2, atomId1, atomId2, seqId1, seqId2) =
			l.get("ptnr1_label_asym_id", "ptnr2_label_asym_id",
				  "ptnr1_label_atom_id", "ptnr2_label_atom_id",
				  "ptnr1_label_seq_id", "ptnr2_label_seq_id");
		
		auto a = 
			l["ptnr1_label_seq_id"].empty() ?
				db["atom_site"].find(cif::Key("label_asym_id") == asym1 and cif::Key("label_atom_id") == atomId1) :
				db["atom_site"].find(cif::Key("label_asym_id") == asym1 and cif::Key("label_seq_id") == seqId1 and cif::Key("label_atom_id") == atomId1);
		
		if (a.size() != 1 and VERBOSE > 1)
			cerr << "Unexpected number (" << a.size() << ") of atoms for link with asym_id " << asym1 << " seq_id " << seqId1 << " atom_id " << atomId1 << endl;
		
		auto b =
			l["ptnr2_label_seq_id"].empty() ?
				db["atom_site"].find(cif::Key("label_asym_id") == asym2 and cif::Key("label_atom_id") == atomId2) :
				db["atom_site"].find(cif::Key("label_asym_id") == asym2 and cif::Key("label_seq_id") == seqId2 and cif::Key("label_atom_id") == atomId2);

		if (b.size() != 1 and VERBOSE > 1)
			cerr << "Unexpected number (" << b.size() << ") of atoms for link with asym_id " << asym2 << " seq_id " << seqId2 << " atom_id " << atomId2 << endl;
		
		if (not (a.empty() or b.empty()))
			linkAtoms(a.front()["id"].as<string>(), b.front()["id"].as<string>());
	}

	// then link all atoms in the compounds
	
	cif::Progress progress(compounds.size(), "Creating bond map");

	for (auto c: compounds)
	{
		auto* compound = mmcif::Compound::create(c);
		if (not compound)
		{
			if (VERBOSE)
				cerr << "Missing compound information for " << c << endl;
			continue;
		}
		
		if (compound->isWater())
		{
			if (VERBOSE)
				cerr << "skipping water in bond map calculation" << endl;
			continue;
		}
		
		// loop over poly_seq_scheme
		for (auto r: db["pdbx_poly_seq_scheme"].find(cif::Key("mon_id") == c))
		{
			string asymId;
			int seqId;
			cif::tie(asymId, seqId) = r.get("asym_id", "seq_id");
			
			vector<Atom> rAtoms;
			copy_if(atoms.begin(), atoms.end(), back_inserter(rAtoms),
				[&](auto& a) { return a.labelAsymId() == asymId and a.labelSeqId() == seqId; });
			
			for (uint32 i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (uint32 j = i + 1; j < rAtoms.size(); ++j)
				{
					if (compound->atomsBonded(rAtoms[i].labelAtomId(), rAtoms[j].labelAtomId()))
						bindAtoms(rAtoms[i].id(), rAtoms[j].id());
				}
			}
		}

		// loop over pdbx_nonpoly_scheme
		for (auto r: db["pdbx_nonpoly_scheme"].find(cif::Key("mon_id") == c))
		{
			string asymId;
			cif::tie(asymId) = r.get("asym_id");
			
			vector<Atom> rAtoms;
			copy_if(atoms.begin(), atoms.end(), back_inserter(rAtoms),
				[&](auto& a) { return a.labelAsymId() == asymId; });
//			for (auto a: db["atom_site"].find(cif::Key("label_asym_id") == asymId))
//				rAtoms.push_back(p.getAtomById(a["id"].as<string>()));
			
			for (uint32 i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (uint32 j = i + 1; j < rAtoms.size(); ++j)
				{
					if (compound->atomsBonded(rAtoms[i].labelAtomId(), rAtoms[j].labelAtomId()))
					{
						uint32 ixa = index[rAtoms[i].id()];
						uint32 ixb = index[rAtoms[j].id()];
						
						bond.insert(key(ixa, ixb));
					}
				}
			}
		}
	}
	
	
	// start by creating an index for single bonds

//cout << "Maken van b1_2 voor " << bond.size() << " bindingen" << endl;
	
	multimap<uint32,uint32> b1_2;
	for (auto& bk: bond)
	{
		uint32 a, b;
		tie(a, b) = dekey(bk);
		
		b1_2.insert({ a, b });
		b1_2.insert({ b, a });
	}
	
//cout << "Afmeting b1_2: " << b1_2.size() << endl;

	multimap<uint32,uint32> b1_3;
	for (uint32 i = 0; i < dim; ++i)
	{
		auto a = b1_2.equal_range(i);
		
		vector<uint32> s;
		for (auto j = a.first; j != a.second; ++j)
			s.push_back(j->second);

		for (size_t si1 = 0; si1 + 1 < s.size(); ++si1)
		{
			for (size_t si2 = si1 + 1; si2 < s.size(); ++si2)
			{
				uint32 x = s[si1];
				uint32 y = s[si2];
				
				if (isBonded(x, y))
					continue;
				
				b1_3.insert({ x, y });
				b1_3.insert({ y, x });
			}
		}
	}

//cout << "Afmeting b1_3: " << b1_3.size() << endl;

	for (uint32 i = 0; i < dim; ++i)
	{
		auto a1 = b1_2.equal_range(i);
		auto a2 = b1_3.equal_range(i);
		
		for (auto ai1 = a1.first; ai1 != a1.second; ++ai1)
		{
			for (auto ai2 = a2.first; ai2 != a2.second; ++ai2)
			{
				uint32 b1 = ai1->second;
				uint32 b2 = ai2->second;
				
				if (isBonded(b1, b2))
					continue;
				
				bond_1_4.insert(key(b1, b2));
			}
		}
	}
//cout << "Afmeting b1_4: " << bond_1_4.size() << endl;
}

vector<string> BondMap::linked(const Atom& a) const
{
	auto i = link.find(a.id());
	
	vector<string> result;
	
	if (i != link.end())
		result = vector<string>(i->second.begin(), i->second.end());

	return result;
}


}
