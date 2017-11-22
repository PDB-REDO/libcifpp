#include "cif++/Config.h"

#include <set>
#include <map>
#include <unordered_set>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>

#include "cif++/Cif++.h"
#include "cif++/PeptideDB.h"

using namespace std;
namespace fs = boost::filesystem;
namespace ba = boost::algorithm;

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

struct PeptideDBImpl
{
	PeptideDBImpl(istream& data, PeptideDBImpl* next);

	~PeptideDBImpl()
	{
		delete mNext;
	}

	/*unordered_*/set<string>	mKnownPeptides;
	set<string>					mKnownBases;
	cif::File					mFile;
	cif::Category&				mChemComp;
	PeptideDBImpl*				mNext;

	string nameFor(const string& resName) const
	{
		string result;
		
		for (auto& chemComp: mChemComp)
		{
			if (ba::iequals(chemComp["three_letter_code"].as<string>(), resName) == false)
				continue;
			
			result = chemComp["name"].as<string>();
			ba::trim(result);
			break;
		}
		
		if (result.empty() and mNext)
			result = mNext->nameFor(resName);
		
		return result;
	}

	string formulaFor(string resName) const;

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
};

PeptideDBImpl::PeptideDBImpl(istream& data, PeptideDBImpl* next)
	: mFile(data), mChemComp(mFile.firstDatablock()["chem_comp"]), mNext(next)
{
	for (auto& chemComp: mChemComp)
	{
		string group, threeLetterCode;
		
		cif::tie(group, threeLetterCode) = chemComp.get("group", "three_letter_code");
		
		if (group == "peptide" or group == "M-peptide" or group == "P-peptide")
			mKnownPeptides.insert(threeLetterCode);
		else if (group == "DNA" or group == "RNA")
			mKnownBases.insert(threeLetterCode);
	}
}

string PeptideDBImpl::formulaFor(string res) const
{
	string result;

	ba::to_upper(res);
	
	for (auto& db: mFile)
	{
		if (db.getName() != "comp_" + res)
			continue;
		
		auto& cat = db["chem_comp_atom"];
				
		map<string,uint32> atoms;
		for (auto r: cat)
			atoms[r["type_symbol"].as<string>()] += 1;
		
		for (auto a: atoms)
		{
			if (not result.empty())
				result += ' ';
			
			result += a.first;
			if (a.second > 1)
				result += to_string(a.second);	
		}		
	}
	
	if (result.empty())
	{
		if (mNext != nullptr)
			result = mNext->formulaFor(res);
		else
		{
			const char* clibdMon = getenv("CLIBD_MON");
			if (clibdMon == nullptr)
				throw runtime_error("Cannot locate peptide list, please souce the CCP4 environment");
			
			fs::path resFile = fs::path(clibdMon) / ba::to_lower_copy(res.substr(0, 1)) / (res + ".cif");
			if (fs::exists(resFile))
			{
				fs::ifstream file(resFile);
				if (file.is_open())
				{
					try
					{
						cif::File cf(file);
						
						auto& cat = cf["comp_" + res]["chem_comp_atom"];
						
						map<string,uint32> atoms;
						for (auto r: cat)
							atoms[r["type_symbol"].as<string>()] += 1;
						
						for (auto a: atoms)
						{
							if (not result.empty())
								result += ' ';
							
							result += a.first;
							if (a.second > 1)
								result += to_string(a.second);	
						}
					}
					catch (exception& ex)
					{
						if (VERBOSE)
							cerr << ex.what();
						result.clear();
					}
				}
			}
		}
	}

	return result;
}

// --------------------------------------------------------------------

PeptideDB* PeptideDB::sInstance;

PeptideDB& PeptideDB::Instance()
{
	if (sInstance == nullptr)
		sInstance = new PeptideDB();
	return *sInstance;
}

PeptideDB::PeptideDB()
{
	const char* clibdMon = getenv("CLIBD_MON");
	if (clibdMon == nullptr)
		throw runtime_error("Cannot locate peptide list, please souce the CCP4 environment");
	
	fs::path db = fs::path(clibdMon) / "list" / "mon_lib_list.cif";

	pushDictionary(db);
	
	sInstance = this;
}

void PeptideDB::pushDictionary(boost::filesystem::path dict)
{
	if (not fs::exists(dict))
		throw runtime_error("file not found: " + dict.string());

	fs::ifstream file(dict);
	if (not file.is_open())
		throw runtime_error("Could not open peptide list " + dict.string());

	mImpl = new PeptideDBImpl(file, mImpl);
}

void PeptideDB::popDictionary()
{
	if (mImpl != nullptr)
	{
		auto i = mImpl;
		mImpl = i->mNext;
		i->mNext = nullptr;
		delete i;
	}
}

PeptideDB::~PeptideDB()
{
	delete mImpl;
}

bool PeptideDB::isKnownPeptide(const string& resName) const
{
	return mImpl->mKnownPeptides.count(resName) > 0;
}

bool PeptideDB::isKnownBase(const string& resName) const
{
	return mImpl->mKnownBases.count(resName) > 0;
}

string PeptideDB::nameForResidue(const string& resName) const
{
	return mImpl->nameFor(resName);
}

string PeptideDB::formulaForResidue(const string& resName) const
{
	return mImpl->formulaFor(resName);
}

string PeptideDB::unalias(const string& resName) const
{
	return mImpl->unalias(resName);
}
