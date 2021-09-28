/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <fstream>
#include <algorithm>
#include <mutex>

#include "cif++/Cif++.hpp"
#include "cif++/Compound.hpp"
#include "cif++/CifUtils.hpp"
#include "cif++/BondMap.hpp"

namespace mmcif
{

namespace
{

union IDType
{
	IDType() : id_n(0){}
	IDType(const IDType& rhs) : id_n(rhs.id_n) {}
	IDType(const std::string& s)
		: IDType()
	{
		assert(s.length() <= 4);
		if (s.length() > 4)
			throw BondMapException("Atom ID '" + s + "' is too long");
		std::copy(s.begin(), s.end(), id_s);
	}

	IDType& operator=(const IDType& rhs)
	{
		id_n = rhs.id_n;
		return *this;
	}

	IDType& operator=(const std::string& s)
	{
		id_n = 0;
		assert(s.length() <= 4);
		if (s.length() > 4)
			throw BondMapException("Atom ID '" + s + "' is too long");
		std::copy(s.begin(), s.end(), id_s);
		return *this;
	}

	bool operator<(const IDType& rhs) const
	{
		return id_n < rhs.id_n;
	}

	bool operator<=(const IDType& rhs) const
	{
		return id_n <= rhs.id_n;
	}

	bool operator==(const IDType& rhs) const
	{
		return id_n == rhs.id_n;
	}

	bool operator!=(const IDType& rhs) const
	{
		return id_n != rhs.id_n;
	}

	char		id_s[4];
	uint32_t	id_n;
};

static_assert(sizeof(IDType) == 4, "atom_id_type should be 4 bytes");	
}

// // --------------------------------------------------------------------

// void createBondInfoFile(const fs::path& components, const fs::path& infofile)
// {
// 	std::ofstream outfile(infofile.string() + ".tmp", std::ios::binary);
// 	if (not outfile.is_open())
// 		throw BondMapException("Could not create bond info file " + infofile.string() + ".tmp");

// 	cif::File infile(components);

// 	std::set<atom_id_type> atomIDs;
// 	std::vector<atom_id_type> compoundIDs;

// 	for (auto& db: infile)
// 	{
// 		auto chem_comp_bond = db.get("chem_comp_bond");
// 		if (not chem_comp_bond)
// 		{
// 			if (cif::VERBOSE > 1)
// 				std::cerr << "Missing chem_comp_bond category in data block " << db.getName() << std::endl;
// 			continue;
// 		}

// 		for (const auto& [atom_id_1, atom_id_2]: chem_comp_bond->rows<std::string,std::string>({"atom_id_1", "atom_id_2"}))
// 		{
// 			atomIDs.insert(atom_id_1);
// 			atomIDs.insert(atom_id_2);
// 		}

// 		compoundIDs.push_back({ db.getName() });
// 	}

// 	if (cif::VERBOSE)
// 		std::cout << "Number of unique atom names is " << atomIDs.size() << std::endl
// 				  << "Number of unique residue names is " << compoundIDs.size() << std::endl;

// 	CompoundBondInfoFileHeader header = {};
// 	header.indexEntries = compoundIDs.size();
// 	header.atomEntries = atomIDs.size();

// 	outfile << header;
	
// 	for (auto atomID: atomIDs)
// 		outfile << atomID;

// 	auto dataOffset = outfile.tellp();

// 	std::vector<CompoundBondInfo> entries;
// 	entries.reserve(compoundIDs.size());

// 	std::map<atom_id_type, uint16_t> atomIDMap;
// 	for (auto& atomID: atomIDs)
// 		atomIDMap[atomID] = atomIDMap.size();

// 	for (auto& db: infile)
// 	{
// 		auto chem_comp_bond = db.get("chem_comp_bond");
// 		if (not chem_comp_bond)
// 			continue;

// 		std::set<uint16_t> bondedAtoms;

// 		for (const auto& [atom_id_1, atom_id_2]: chem_comp_bond->rows<std::string,std::string>({"atom_id_1", "atom_id_2"}))
// 		{
// 			bondedAtoms.insert(atomIDMap[atom_id_1]);
// 			bondedAtoms.insert(atomIDMap[atom_id_2]);
// 		}

// 		std::map<uint16_t, int32_t> bondedAtomMap;
// 		for (auto id: bondedAtoms)
// 			bondedAtomMap[id] = static_cast<int32_t>(bondedAtomMap.size());
		
// 		CompoundBondInfo info = {
// 			db.getName(),
// 			static_cast<uint32_t>(bondedAtomMap.size()),
// 			outfile.tellp() - dataOffset
// 		};

// 		entries.push_back(info);

// 		// An now first write the array of atom ID's in this compound
// 		for (uint16_t id: bondedAtoms)
// 			write(outfile, id);

// 		// And then the symmetric matrix with bonds
// 		size_t N = bondedAtoms.size();
// 		size_t M = (N * (N - 1)) / 2;

// 		size_t K = M / 8;
// 		if (M % 8)
// 			K += 1;
		
// 		std::vector<uint8_t> m(K);

// 		for (const auto& [atom_id_1, atom_id_2]: chem_comp_bond->rows<std::string,std::string>({"atom_id_1", "atom_id_2"}))
// 		{
// 			auto a = bondedAtomMap[atomIDMap[atom_id_1]];
// 			auto b = bondedAtomMap[atomIDMap[atom_id_2]];

// 			assert(a != b);
// 			assert((int)b < (int)N);

// 			if (a > b)
// 				std::swap(a, b);
			
// 			size_t ix = ((b - 1) * b) / 2 + a;
// 			assert(ix < M);

// 			auto Bix = ix / 8;
// 			auto bix = ix % 8;

// 			m[Bix] |= 1 << bix;
// 		}

// 		outfile.write(reinterpret_cast<char*>(m.data()), m.size());
// 	}

// 	header.dataSize = outfile.tellp() - dataOffset;

// 	std::sort(entries.begin(), entries.end(), [](CompoundBondInfo& a, CompoundBondInfo& b)
// 	{
// 		return a.id < b.id;
// 	});

// 	for (auto& info: entries)
// 		outfile << info;

// 	outfile.seekp(0);
// 	outfile << header;

// 	// compress
// 	outfile.close();

// 	std::ifstream in(infofile.string() + ".tmp", std::ios::binary);
// 	std::ofstream out(infofile, std::ios::binary);

// 	{
// 		io::filtering_stream<io::output> os;
// 		os.push(io::gzip_compressor());
// 		os.push(out);
// 		io::copy(in, os);
// 	}

// 	in.close();
// 	out.close();

// 	fs::remove(infofile.string() + ".tmp");
// }

// --------------------------------------------------------------------

struct CompoundBondInfo
{
	IDType	mID;
	std::set<std::tuple<uint32_t,uint32_t>> mBonded;

	bool bonded(uint32_t a1, uint32_t a2) const
	{
		return mBonded.count({ a1, a2 }) > 0;
	}

};

// --------------------------------------------------------------------

class CompoundBondMap
{
  public:

	static CompoundBondMap &instance()
	{
		static std::unique_ptr<CompoundBondMap> s_instance(new CompoundBondMap);
		return *s_instance;
	}

	bool bonded(const std::string& compoundID, const std::string& atomID1, const std::string& atomID2);

  private:

	CompoundBondMap() {}

	uint32_t getAtomID(const std::string& atomID)
	{
		IDType id(atomID);

		uint32_t result;

		auto i = mAtomIDIndex.find(id);
		if (i == mAtomIDIndex.end())
		{
			result = uint32_t(mAtomIDIndex.size());
			mAtomIDIndex[id] = result;
		}
		else
			result = i->second;
		
		return result;
	}

	std::map<IDType,uint32_t> mAtomIDIndex;
	std::vector<CompoundBondInfo> mCompounds;
	std::mutex mMutex;
};

bool CompoundBondMap::bonded(const std::string &compoundID, const std::string& atomID1, const std::string& atomID2)
{
	std::lock_guard lock(mMutex);

	using namespace std::literals;

	IDType id(compoundID);
	uint32_t a1 = getAtomID(atomID1);
	uint32_t a2 = getAtomID(atomID2);
	if (a1 > a2)
		std::swap(a1, a2);
	
	for (auto &bi: mCompounds)
	{
		if (bi.mID != id)
			continue;
		
		return bi.bonded(a1, a2);
	}

	bool result = false;

	// not found in our cache, calculate
	CompoundBondInfo bondInfo{ id };

	auto compound = mmcif::CompoundFactory::instance().create(compoundID);
	if (not compound)
		std::cerr << "Missing compound bond info for " << compoundID << std::endl;
	else
	{
		for (auto &atom: compound->bonds())
		{
			uint32_t ca1 = getAtomID(atom.atomID[0]);
			uint32_t ca2 = getAtomID(atom.atomID[1]);
			if (ca1 > ca2)
				std::swap(ca1, ca2);
			
			bondInfo.mBonded.insert({ca1, ca2});
			result = result or (a1 == ca1 and a2 == ca2);
		}
	}

	mCompounds.push_back(bondInfo);

	return result;
}

// --------------------------------------------------------------------

BondMap::BondMap(const Structure& p)
{
	auto& compoundBondInfo = CompoundBondMap::instance();

	auto atoms = p.atoms();
	dim = uint32_t(atoms.size());

//	bond = std::vector<bool>(dim * (dim - 1), false);

	for (auto& atom: atoms)
		index[atom.id()] = uint32_t(index.size());
	
	auto bindAtoms = [this](const std::string& a, const std::string& b)
	{
		uint32_t ixa = index[a];
		uint32_t ixb = index[b];
		
		bond.insert(key(ixa, ixb));
	};

	auto linkAtoms = [this,&bindAtoms](const std::string& a, const std::string& b)
	{
		bindAtoms(a, b);

		link[a].insert(b);
		link[b].insert(a);
	};

	cif::Datablock& db = p.getFile().data();

	// collect all compounds first
	std::set<std::string> compounds;
	for (auto c: db["chem_comp"])
		compounds.insert(c["id"].as<std::string>());
	
	// make sure we also have all residues in the polyseq
	for (auto m: db["entity_poly_seq"])
	{
		std::string c = m["mon_id"].as<std::string>();
		if (compounds.count(c))
			continue;
		
		if (cif::VERBOSE > 1)
			std::cerr << "Warning: mon_id " << c << " is missing in the chem_comp category" << std::endl;
		compounds.insert(c);
	}

	cif::Progress progress(compounds.size(), "Creating bond map");

	// some helper indices to speed things up a bit
	std::map<std::tuple<std::string,int,std::string>,std::string> atomMapByAsymSeqAndAtom;
	for (auto& a: p.atoms())
	{
		auto key = make_tuple(a.labelAsymID(), a.labelSeqID(), a.labelAtomID());
		atomMapByAsymSeqAndAtom[key] = a.id();
	}
	
	// first link all residues in a polyseq
	
	std::string lastAsymID;
	int lastSeqID = 0;
	for (auto r: db["pdbx_poly_seq_scheme"])
	{
		std::string asymID;
		int seqID;

		cif::tie(asymID, seqID) = r.get("asym_id", "seq_id");

		if (asymID != lastAsymID)		// first in a new sequece
		{
			lastAsymID = asymID;
			lastSeqID = seqID;
			continue;
		}
		
		auto c = atomMapByAsymSeqAndAtom[make_tuple(asymID, lastSeqID, "C")];
		auto n = atomMapByAsymSeqAndAtom[make_tuple(asymID, seqID, "N")];

		if (not (c.empty() or n.empty()))
			bindAtoms(c, n);
		
		lastSeqID = seqID;
	}

	for (auto l: db["struct_conn"])
	{
		std::string asym1, asym2, atomId1, atomId2;
		int seqId1 = 0, seqId2 = 0;
		cif::tie(asym1, asym2, atomId1, atomId2, seqId1, seqId2) =
			l.get("ptnr1_label_asym_id", "ptnr2_label_asym_id",
				  "ptnr1_label_atom_id", "ptnr2_label_atom_id",
				  "ptnr1_label_seq_id", "ptnr2_label_seq_id");

		std::string a = atomMapByAsymSeqAndAtom[make_tuple(asym1, seqId1, atomId1)];
		std::string b = atomMapByAsymSeqAndAtom[make_tuple(asym2, seqId2, atomId2)];
			
		if (not (a.empty() or b.empty()))
			linkAtoms(a, b);
	}

	// then link all atoms in the compounds
	
	for (auto c: compounds)
	{
		if (c == "HOH" or c == "H2O" or c == "WAT")
		{
			if (cif::VERBOSE)
				std::cerr << "skipping water in bond map calculation" << std::endl;
			continue;
		}

		auto bonded = [c, &compoundBondInfo](const Atom& a, const Atom& b)
		{
			auto label_a = a.labelAtomID();
			auto label_b = b.labelAtomID();

			return compoundBondInfo.bonded(c, label_a, label_b);
		};

		// loop over poly_seq_scheme
		for (auto r: db["pdbx_poly_seq_scheme"].find(cif::Key("mon_id") == c))
		{
			std::string asymID;
			int seqID;
			cif::tie(asymID, seqID) = r.get("asym_id", "seq_id");
			
			std::vector<Atom> rAtoms;
			copy_if(atoms.begin(), atoms.end(), back_inserter(rAtoms),
				[&](auto& a) { return a.labelAsymID() == asymID and a.labelSeqID() == seqID; });
			
			for (uint32_t i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (uint32_t j = i + 1; j < rAtoms.size(); ++j)
				{
					if (bonded(rAtoms[i], rAtoms[j]))
						bindAtoms(rAtoms[i].id(), rAtoms[j].id());
				}
			}
		}

		// loop over pdbx_nonpoly_scheme
		for (auto r: db["pdbx_nonpoly_scheme"].find(cif::Key("mon_id") == c))
		{
			std::string asymID;
			cif::tie(asymID) = r.get("asym_id");
			
			std::vector<Atom> rAtoms;
			copy_if(atoms.begin(), atoms.end(), back_inserter(rAtoms),
				[&](auto& a) { return a.labelAsymID() == asymID; });
			
			for (uint32_t i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (uint32_t j = i + 1; j < rAtoms.size(); ++j)
				{
					if (bonded(rAtoms[i], rAtoms[j]))
					{
						uint32_t ixa = index[rAtoms[i].id()];
						uint32_t ixb = index[rAtoms[j].id()];
						
						bond.insert(key(ixa, ixb));
					}
				}
			}
		}

		// loop over pdbx_branch_scheme
		for (auto r: db["pdbx_branch_scheme"].find(cif::Key("mon_id") == c))
		{
			std::string asymID;
			cif::tie(asymID) = r.get("asym_id");
			
			std::vector<Atom> rAtoms;
			copy_if(atoms.begin(), atoms.end(), back_inserter(rAtoms),
				[&](auto& a) { return a.labelAsymID() == asymID; });
			
			for (uint32_t i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (uint32_t j = i + 1; j < rAtoms.size(); ++j)
				{
					if (bonded(rAtoms[i], rAtoms[j]))
					{
						uint32_t ixa = index[rAtoms[i].id()];
						uint32_t ixb = index[rAtoms[j].id()];
						
						bond.insert(key(ixa, ixb));
					}
				}
			}
		}
	}
	
	// start by creating an index for single bonds
	
	std::multimap<uint32_t,uint32_t> b1_2;
	for (auto& bk: bond)
	{
		uint32_t a, b;
		std::tie(a, b) = dekey(bk);
		
		b1_2.insert({ a, b });
		b1_2.insert({ b, a });
	}
	
	std::multimap<uint32_t,uint32_t> b1_3;
	for (uint32_t i = 0; i < dim; ++i)
	{
		auto a = b1_2.equal_range(i);
		
		std::vector<uint32_t> s;
		for (auto j = a.first; j != a.second; ++j)
			s.push_back(j->second);

		for (size_t si1 = 0; si1 + 1 < s.size(); ++si1)
		{
			for (size_t si2 = si1 + 1; si2 < s.size(); ++si2)
			{
				uint32_t x = s[si1];
				uint32_t y = s[si2];
				
				if (isBonded(x, y))
					continue;
				
				b1_3.insert({ x, y });
				b1_3.insert({ y, x });
			}
		}
	}

	for (uint32_t i = 0; i < dim; ++i)
	{
		auto a1 = b1_2.equal_range(i);
		auto a2 = b1_3.equal_range(i);
		
		for (auto ai1 = a1.first; ai1 != a1.second; ++ai1)
		{
			for (auto ai2 = a2.first; ai2 != a2.second; ++ai2)
			{
				uint32_t b1 = ai1->second;
				uint32_t b2 = ai2->second;
				
				if (isBonded(b1, b2))
					continue;
				
				bond_1_4.insert(key(b1, b2));
			}
		}
	}
}

std::vector<std::string> BondMap::linked(const Atom& a) const
{
	auto i = link.find(a.id());
	
	std::vector<std::string> result;
	
	if (i != link.end())
		result = std::vector<std::string>(i->second.begin(), i->second.end());

	return result;
}

std::vector<std::string> BondMap::atomIDsForCompound(const std::string& compoundID)
{
	std::vector<std::string> result;

	auto* compound = mmcif::CompoundFactory::instance().create(compoundID);

	if (compound == nullptr)
		throw BondMapException("Missing bond information for compound " + compoundID);

	for (auto& compAtom: compound->atoms())
		result.push_back(compAtom.id);

	return result;
}

}
