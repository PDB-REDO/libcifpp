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

#if __has_include("Config.hpp")
#include "Config.hpp"
#endif

#include <map>
#include <mutex>
#include <numeric>
#include <shared_mutex>

#include <boost/algorithm/string.hpp>

#include <filesystem>
#include <fstream>

#include "cif++/Cif++.hpp"
#include "cif++/CifUtils.hpp"
#include "cif++/Compound.hpp"
#include "cif++/Point.hpp"
#include "cif++/CifParser.hpp"

namespace ba = boost::algorithm;
namespace fs = std::filesystem;

namespace mmcif
{

// --------------------------------------------------------------------

std::string to_string(BondType bondType)
{
	switch (bondType)
	{
		case BondType::sing:	return "sing";
		case BondType::doub:	return "doub";
		case BondType::trip:	return "trip";
		case BondType::quad:	return "quad";
		case BondType::arom:	return "arom";
		case BondType::poly:	return "poly";
		case BondType::delo:	return "delo";
		case BondType::pi:		return "pi";
	}
}

BondType from_string(const std::string& bondType)
{
	if (cif::iequals(bondType, "sing"))	return BondType::sing;
	if (cif::iequals(bondType, "doub"))	return BondType::doub;
	if (cif::iequals(bondType, "trip"))	return BondType::trip;
	if (cif::iequals(bondType, "quad"))	return BondType::quad;
	if (cif::iequals(bondType, "arom"))	return BondType::arom;
	if (cif::iequals(bondType, "poly"))	return BondType::poly;
	if (cif::iequals(bondType, "delo"))	return BondType::delo;
	if (cif::iequals(bondType, "pi"))	return BondType::pi;
	throw std::invalid_argument("Invalid bondType: " + bondType);
}

// --------------------------------------------------------------------
// Compound helper classes

struct CompoundAtomLess
{
	bool operator()(const CompoundAtom &a, const CompoundAtom &b) const
	{
		int d = a.id.compare(b.id);
		if (d == 0)
			d = a.typeSymbol - b.typeSymbol;
		return d < 0;
	}
};

struct CompoundBondLess
{
	bool operator()(const CompoundBond &a, const CompoundBond &b) const
	{
		int d = a.atomID[0].compare(b.atomID[0]);
		if (d == 0)
			d = a.atomID[1].compare(b.atomID[1]);
		if (d == 0)
			d = static_cast<int>(a.type) - static_cast<int>(b.type);
		return d < 0;
	}
};

// --------------------------------------------------------------------
// Compound

Compound::Compound(cif::Datablock &db)
{
	auto &chemComp = db["chem_comp"];

	if (chemComp.size() != 1)
		throw std::runtime_error("Invalid compound file, chem_comp should contain a single row");

	cif::tie(mID, mName, mType, mFormula, mFormulaWeight, mFormalCharge) =
		chemComp.front().get("id", "name", "type", "formula", "formula_weight", "pdbx_formal_charge");

	auto &chemCompAtom = db["chem_comp_atom"];
	for (auto row: chemCompAtom)
	{
		CompoundAtom atom;
		std::string typeSymbol;
		cif::tie(atom.id, typeSymbol, atom.charge, atom.aromatic, atom.leavingAtom, atom.stereoConfig, atom.x, atom.y, atom.z) =
			row.get("atom_id", "type_symbol", "charge", "pdbx_aromatic_flag", "pdbx_leaving_atom_flag", "pdbx_stereo_config",
			"model_Cartn_x", "model_Cartn_y", "model_Cartn_z");
		atom.typeSymbol = AtomTypeTraits(typeSymbol).type();
		mAtoms.push_back(std::move(atom));
	}

	auto &chemCompBond = db["chem_comp_bond"];
	for (auto row: chemCompBond)
	{
		CompoundBond bond;
		std::string valueOrder;
		cif::tie(bond.atomID[0], bond.atomID[1], valueOrder, bond.aromatic, bond.stereoConfig)
			= row.get("atom_id_1", "atom_id_2", "value_order", "pdbx_aromatic_flag", "pdbx_stereo_config");
		bond.type = from_string(valueOrder);
		mBonds.push_back(std::move(bond));
	}
}

Compound::Compound(cif::Datablock &db, const std::string &id, const std::string &name, const std::string &type)
	: mID(id)
	, mName(name)
	, mType(type)
{
	auto &chemCompAtom = db["chem_comp_atom"];
	for (auto row: chemCompAtom)
	{
		CompoundAtom atom;
		std::string typeSymbol;
		cif::tie(atom.id, typeSymbol, atom.charge, atom.x, atom.y, atom.z) =
			row.get("atom_id", "type_symbol", "charge", "x", "y", "z");
		atom.typeSymbol = AtomTypeTraits(typeSymbol).type();

		mFormalCharge += atom.charge;
		mFormulaWeight += AtomTypeTraits(atom.typeSymbol).weight();

		mAtoms.push_back(std::move(atom));
	}

	auto &chemCompBond = db["chem_comp_bond"];
	for (auto row: chemCompBond)
	{
		CompoundBond bond;
		std::string type;
		cif::tie(bond.atomID[0], bond.atomID[1], type, bond.aromatic)
			= row.get("atom_id_1", "atom_id_2", "type", "aromatic");

		using cif::iequals;

			 if (iequals(type, "single"))	bond.type = BondType::sing;
		else if (iequals(type, "double"))	bond.type = BondType::doub;
		else if (iequals(type, "triple"))	bond.type = BondType::trip;
		else if (iequals(type, "deloc") or iequals(type, "aromat") or iequals(type, "aromatic"))
			bond.type = BondType::delo;
		else
		{
			if (cif::VERBOSE)
				std::cerr << "Unimplemented chem_comp_bond.type " << type << " in " << id << std::endl;
			bond.type = BondType::sing;
		}
		mBonds.push_back(std::move(bond));
	}
}

CompoundAtom Compound::getAtomByID(const std::string &atomID) const
{
	CompoundAtom result = {};
	for (auto &a : mAtoms)
	{
		if (a.id == atomID)
		{
			result = a;
			break;
		}
	}

	if (result.id != atomID)
		throw std::out_of_range("No atom " + atomID + " in Compound " + mID);

	return result;
}

const Compound *Compound::create(const std::string &id)
{
	auto result = CompoundFactory::instance().get(id);
	if (result == nullptr)
		result = CompoundFactory::instance().create(id);
	return result;
}

bool Compound::atomsBonded(const std::string &atomId_1, const std::string &atomId_2) const
{
	auto i = find_if(mBonds.begin(), mBonds.end(),
		[&](const CompoundBond &b) {
			return (b.atomID[0] == atomId_1 and b.atomID[1] == atomId_2) or (b.atomID[0] == atomId_2 and b.atomID[1] == atomId_1);
		});

	return i != mBonds.end();
}

// float Compound::atomBondValue(const std::string &atomId_1, const std::string &atomId_2) const
// {
// 	auto i = find_if(mBonds.begin(), mBonds.end(),
// 		[&](const CompoundBond &b) {
// 			return (b.atomID[0] == atomId_1 and b.atomID[1] == atomId_2) or (b.atomID[0] == atomId_2 and b.atomID[1] == atomId_1);
// 		});

// 	return i != mBonds.end() ? i->distance : 0;
// }


// float Compound::bondAngle(const std::string &atomId_1, const std::string &atomId_2, const std::string &atomId_3) const
// {
// 	float result = nanf("1");

// 	for (auto &a : mAngles)
// 	{
// 		if (not(a.atomID[1] == atomId_2 and
// 				((a.atomID[0] == atomId_1 and a.atomID[2] == atomId_3) or
// 					(a.atomID[2] == atomId_1 and a.atomID[0] == atomId_3))))
// 			continue;

// 		result = a.angle;
// 		break;
// 	}

// 	return result;
// }

//static float calcC(float a, float b, float alpha)
//{
//	float f = b * sin(alpha * kPI / 180);
//	float d = sqrt(b * b - f * f);
//	float e = a - d;
//	float c = sqrt(f * f + e * e);
//
//	return c;
//}

// float Compound::chiralVolume(const std::string &centreID) const
// {
// 	float result = 0;

// 	for (auto &cv : mChiralCentres)
// 	{
// 		if (cv.id != centreID)
// 			continue;

// 		// calculate the expected chiral volume

// 		// the edges

// 		float a = atomBondValue(cv.atomIDCentre, cv.atomID[0]);
// 		float b = atomBondValue(cv.atomIDCentre, cv.atomID[1]);
// 		float c = atomBondValue(cv.atomIDCentre, cv.atomID[2]);

// 		// the angles for the top of the tetrahedron

// 		float alpha = bondAngle(cv.atomID[0], cv.atomIDCentre, cv.atomID[1]);
// 		float beta = bondAngle(cv.atomID[1], cv.atomIDCentre, cv.atomID[2]);
// 		float gamma = bondAngle(cv.atomID[2], cv.atomIDCentre, cv.atomID[0]);

// 		float cosa = cos(alpha * kPI / 180);
// 		float cosb = cos(beta * kPI / 180);
// 		float cosc = cos(gamma * kPI / 180);

// 		result = (a * b * c * sqrt(1 + 2 * cosa * cosb * cosc - (cosa * cosa) - (cosb * cosb) - (cosc * cosc))) / 6;

// 		if (cv.volumeSign == negativ)
// 			result = -result;

// 		break;
// 	}

// 	return result;
// }

// --------------------------------------------------------------------
// a factory class to generate compounds

const std::map<std::string, char> kAAMap{
	{"ALA", 'A'},
	{"ARG", 'R'},
	{"ASN", 'N'},
	{"ASP", 'D'},
	{"CYS", 'C'},
	{"GLN", 'Q'},
	{"GLU", 'E'},
	{"GLY", 'G'},
	{"HIS", 'H'},
	{"ILE", 'I'},
	{"LEU", 'L'},
	{"LYS", 'K'},
	{"MET", 'M'},
	{"PHE", 'F'},
	{"PRO", 'P'},
	{"SER", 'S'},
	{"THR", 'T'},
	{"TRP", 'W'},
	{"TYR", 'Y'},
	{"VAL", 'V'},
	{"GLX", 'Z'},
	{"ASX", 'B'}};

const std::map<std::string, char> kBaseMap{
	{"A", 'A'},
	{"C", 'C'},
	{"G", 'G'},
	{"T", 'T'},
	{"U", 'U'},
	{"DA", 'A'},
	{"DC", 'C'},
	{"DG", 'G'},
	{"DT", 'T'}};

// --------------------------------------------------------------------

class CompoundFactoryImpl
{
  public:
	CompoundFactoryImpl();

	CompoundFactoryImpl(const std::string &file, CompoundFactoryImpl *next);

	virtual ~CompoundFactoryImpl()
	{
		delete mNext;
	}

	Compound *get(std::string id);
	virtual Compound *create(std::string id);

	CompoundFactoryImpl *pop()
	{
		auto result = mNext;
		mNext = nullptr;
		delete this;
		return result;
	}

	bool isKnownPeptide(const std::string &resName)
	{
		return mKnownPeptides.count(resName) or
			   (mNext != nullptr and mNext->isKnownPeptide(resName));
	}

	bool isKnownBase(const std::string &resName)
	{
		return mKnownBases.count(resName) or
			   (mNext != nullptr and mNext->isKnownBase(resName));
	}

  protected:
	std::shared_timed_mutex mMutex;

	std::vector<Compound *> mCompounds;
	std::set<std::string> mKnownPeptides;
	std::set<std::string> mKnownBases;
	std::set<std::string> mMissing;
	CompoundFactoryImpl *mNext = nullptr;
};

// --------------------------------------------------------------------

CompoundFactoryImpl::CompoundFactoryImpl()
{
	for (const auto &[key, value] : kAAMap)
		mKnownPeptides.insert(key);

	for (const auto &[key, value] : kBaseMap)
		mKnownBases.insert(key);
}

CompoundFactoryImpl::CompoundFactoryImpl(const std::string &file, CompoundFactoryImpl *next)
	: mNext(next)
{
	cif::File cifFile(file);

	auto compList = cifFile.get("comp_list");
	if (compList)	// So this is a CCP4 restraints file, special handling
	{
		auto &chemComp = (*compList)["chem_comp"];

		for (const auto &[id, name, group]: chemComp.rows<std::string,std::string,std::string>({"id", "name", "group"}))
		{
			std::string type;
			
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
			
			if (cif::iequals(id, "gly"))
				type = "peptide linking";
			else if (cif::iequals(group, "l-peptide") or cif::iequals(group, "L-peptide linking") or cif::iequals(group, "peptide"))
				type = "L-peptide linking";
			else if (cif::iequals(group, "DNA"))
				type = "DNA linking";
			else if (cif::iequals(group, "RNA"))
				type = "RNA linking";
			else
				type = "non-polymer";
			
			auto &db = cifFile["comp_" + id];

			mCompounds.push_back(new Compound(db, id, name, type));
		}
	}
	else
	{
		// A CCD components file, validate it first
		cifFile.loadDictionary("mmcif_pdbx_v50");

		if (not cifFile.isValid())
			throw std::runtime_error("Invalid compound file");

		for (auto &db: cifFile)
			mCompounds.push_back(new Compound(db));
	}
}

Compound *CompoundFactoryImpl::get(std::string id)
{
	std::shared_lock lock(mMutex);

	ba::to_upper(id);

	Compound *result = nullptr;

	for (auto cmp : mCompounds)
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

Compound *CompoundFactoryImpl::create(std::string id)
{
	ba::to_upper(id);

	Compound *result = get(id);
	// if (result == nullptr and mMissing.count(id) == 0 and not mFile.empty())
	// {
	// 	std::unique_lock lock(mMutex);

	// 	auto &cat = mFile.firstDatablock()["chem_comp"];

	// 	auto rs = cat.find(cif::Key("three_letter_code") == id);

	// 	if (not rs.empty())
	// 	{
	// 		auto row = rs.front();

	// 		std::string name, group;
	// 		uint32_t numberAtomsAll, numberAtomsNh;
	// 		cif::tie(name, group, numberAtomsAll, numberAtomsNh) =
	// 			row.get("name", "group", "number_atoms_all", "number_atoms_nh");

	// 		ba::trim(name);
	// 		ba::trim(group);

	// 		if (mFile.get("comp_" + id) == nullptr)
	// 		{
	// 			auto clibd_mon = fs::path(getenv("CLIBD_MON"));

	// 			fs::path resFile = clibd_mon / ba::to_lower_copy(id.substr(0, 1)) / (id + ".cif");

	// 			if (not fs::exists(resFile) and (id == "COM" or id == "CON" or "PRN")) // seriously...
	// 				resFile = clibd_mon / ba::to_lower_copy(id.substr(0, 1)) / (id + '_' + id + ".cif");

	// 			if (not fs::exists(resFile))
	// 				mMissing.insert(id);
	// 			else
	// 			{
	// 				mCompounds.push_back(new Compound(resFile.string(), id, name, group));
	// 				result = mCompounds.back();
	// 			}
	// 		}
	// 		else
	// 		{
	// 			mCompounds.push_back(new Compound(mPath, id, name, group));
	// 			result = mCompounds.back();
	// 		}
	// 	}

	// 	if (result == nullptr and mNext != nullptr)
	// 		result = mNext->create(id);
	// }

	return result;
}

// --------------------------------------------------------------------
// Version for the default compounds, based on the cached components.cif file from CCD

class CCDCompoundFactoryImpl : public CompoundFactoryImpl
{
  public:
	CCDCompoundFactoryImpl() {}

	Compound *create(std::string id) override;

	cif::DatablockIndex mIndex;
};

Compound *CCDCompoundFactoryImpl::create(std::string id)
{
	ba::to_upper(id);

	Compound *result = get(id);

	if (result)
		return result;

	auto ccd = cif::loadResource("components.cif");
	if (not ccd)
		throw std::runtime_error("Could not locate the CCD components.cif file, please make sure the software is installed properly and/or use the update-dictionary-script to fetch the data.");

	cif::File file;

	if (mIndex.empty())
	{
		if (cif::VERBOSE > 1)
		{
			std::cout << "Creating component index " << "...";
			std::cout.flush();
		}

		cif::Parser parser(*ccd, file, false);
		mIndex = parser.indexDatablocks();

		if (cif::VERBOSE > 1)
			std::cout << " done" << std::endl;
		
		// reload the resource, perhaps this should be improved...
		ccd = cif::loadResource("components.cif");
	}

	if (cif::VERBOSE > 1)
	{
		std::cout << "Loading component " << id << "...";
		std::cout.flush();
	}

	cif::Parser parser(*ccd, file, false);
	parser.parseSingleDatablock(id, mIndex);

	if (cif::VERBOSE > 1)
		std::cout << " done" << std::endl;

	if (not file.empty())
	{
		auto &db = file.firstDatablock();
		if (db.getName() == id)
		{
			result = new Compound(db);

			std::shared_lock lock(mMutex);
			mCompounds.push_back(result);
		}
	}

	if (result == nullptr and cif::VERBOSE > 1)
		std::cerr << "Could not locate compound " << id << " in the CCD components file" << std::endl;

	return result;
}

// --------------------------------------------------------------------

CompoundFactory *CompoundFactory::sInstance;
thread_local std::unique_ptr<CompoundFactory> CompoundFactory::tlInstance;
bool CompoundFactory::sUseThreadLocalInstance;

void CompoundFactory::init(bool useThreadLocalInstanceOnly)
{
	sUseThreadLocalInstance = useThreadLocalInstanceOnly;
}

CompoundFactory::CompoundFactory()
	: mImpl(new CCDCompoundFactoryImpl)
{
}

CompoundFactory::~CompoundFactory()
{
	delete mImpl;
}

CompoundFactory &CompoundFactory::instance()
{
	if (sUseThreadLocalInstance)
	{
		if (not tlInstance)
			tlInstance.reset(new CompoundFactory());
		return *tlInstance;
	}
	else
	{
		if (sInstance == nullptr)
			sInstance = new CompoundFactory();
		return *sInstance;
	}
}

void CompoundFactory::clear()
{
	if (sUseThreadLocalInstance)
		tlInstance.reset(nullptr);
	else
	{
		delete sInstance;
		sInstance = nullptr;
	}
}

void CompoundFactory::pushDictionary(const std::string &inDictFile)
{
	if (not fs::exists(inDictFile))
		throw std::runtime_error("file not found: " + inDictFile);

	//	ifstream file(inDictFile);
	//	if (not file.is_open())
	//		throw std::runtime_error("Could not open peptide list " + inDictFile);

	try
	{
		mImpl = new CompoundFactoryImpl(inDictFile, mImpl);
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Error loading dictionary " << inDictFile << std::endl;
		throw;
	}
}

void CompoundFactory::popDictionary()
{
	if (mImpl != nullptr)
		mImpl = mImpl->pop();
}

// id is the three letter code
const Compound *CompoundFactory::get(std::string id)
{
	return mImpl->get(id);
}

const Compound *CompoundFactory::create(std::string id)
{
	return mImpl->create(id);
}

bool CompoundFactory::isKnownPeptide(const std::string &resName) const
{
	return mImpl->isKnownPeptide(resName);
}

bool CompoundFactory::isKnownBase(const std::string &resName) const
{
	return mImpl->isKnownBase(resName);
}

} // namespace mmcif
