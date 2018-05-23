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
	: dim(0)
{
	auto atoms = p.atoms();
	dim = atoms.size();
	
	bond = vector<bool>(dim * (dim - 1), false);

	for (auto& atom: atoms)
	{
		size_t ix = index.size();
		index[atom.id()] = ix;
	};
	
	auto bindAtoms = [this](const string& a, const string& b)
	{
		size_t ixa = index[a];
		size_t ixb = index[b];
		
		if (ixb < ixa)
			swap(ixa, ixb);
		
		size_t ix = ixb + ixa * dim - ixa * (ixa + 1) / 2;
		
		assert(ix < bond.size());
		bond[ix] = true;
		
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
			bindAtoms(a.front()["id"].as<string>(), b.front()["id"].as<string>());
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
			for (auto a: db["atom_site"].find(cif::Key("label_asym_id") == asymId and cif::Key("label_seq_id") == seqId))
				rAtoms.push_back(p.getAtomById(a["id"].as<string>()));
			
			for (size_t i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (size_t j = i + 1; j < rAtoms.size(); ++j)
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
			for (auto a: db["atom_site"].find(cif::Key("label_asym_id") == asymId))
				rAtoms.push_back(p.getAtomById(a["id"].as<string>()));
			
			for (size_t i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (size_t j = i + 1; j < rAtoms.size(); ++j)
				{
					if (compound->atomsBonded(rAtoms[i].labelAtomId(), rAtoms[j].labelAtomId()))
					{
						size_t ixa = index[rAtoms[i].id()];
						size_t ixb = index[rAtoms[j].id()];
						
						if (ixb < ixa)
							swap(ixa, ixb);
						
						size_t ix = ixb + ixa * dim - ixa * (ixa + 1) / 2;
						
						assert(ix < bond.size());
						bond[ix] = true;
					}
				}
			}
		}
		
		progress.consumed(1);
	}
}

//bool BondMap::operator()(const Atom& a, const Atom& b) const
//{
//	size_t ixa = index.at(a.id());
//	size_t ixb = index.at(b.id());
//	
//	if (ixb < ixa)
//		swap(ixa, ixb);
//	
//	size_t ix = ixb + ixa * dim - ixa * (ixa + 1) / 2;
//	
//	assert(ix < bond.size());
//	return bond[ix];
//}

bool BondMap::isBonded(size_t ixa, size_t ixb) const
{
	if (ixa == ixb)
		return false;
	
	if (ixa > ixb)
		swap(ixa, ixb);
	
	size_t ix = ixb + ixa * dim - ixa * (ixa + 1) / 2;
	
	assert(ix < bond.size());
	return bond[ix];
}

bool BondMap::is1_4(const Atom& a, const Atom& b) const
{
	size_t ixa = index.at(a.id());
	size_t ixb = index.at(b.id());
	
	if (ixb < ixa)
		swap(ixa, ixb);
	
	bool result = false;
	
	for (size_t ia = 0; result == false and ia + 1 < dim; ++ia)
	{
		if (ia == ixa or ia == ixb)
			continue;
		
		for (size_t ib = ia + 1; result == false and ib < dim; ++ib)
		{
			if (ib == ixa or ib == ixb)
				continue;
			
			result = isBonded(ixa, ia) and isBonded(ia, ib) and isBonded(ib, ixb);
		}
	}
	
	return result;
}

}
