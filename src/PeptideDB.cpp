#include "libpr.h"

#include <set>
#include <map>
#include <unordered_set>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>

#include "cif++.h"

#include "peptidedb.h"

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
		delete m_next;
	}

	/*unordered_*/set<string>	m_known_peptides;
	set<string>					m_known_bases;
	cif::file					m_file;
	cif::category&				m_chem_comp;
	PeptideDBImpl*				m_next;

	string name_for(const string& res_name) const
	{
		string result;
		
		for (auto& chem_comp: m_chem_comp)
		{
			if (ba::iequals(chem_comp["three_letter_code"].as<string>(), res_name) == false)
				continue;
			
			result = chem_comp["name"].as<string>();
			ba::trim(result);
			break;
		}
		
		if (result.empty() and m_next)
			result = m_next->name_for(res_name);
		
		return result;
	}

	string formula_for(string res_name) const;

	string unalias(const string& res_name) const
	{
		string result = res_name;
		
		auto& e = const_cast<cif::file&>(m_file)["comp_synonym_list"];
		
		for (auto& synonym: e["chem_comp_synonyms"])
		{
			if (ba::iequals(synonym["comp_alternative_id"].as<string>(), res_name) == false)
				continue;
			
			result = synonym["comp_id"].as<string>();
			ba::trim(result);
			break;
		}

		if (result.empty() and m_next)
			result = m_next->unalias(res_name);
		
		return result;
	}
};

PeptideDBImpl::PeptideDBImpl(istream& data, PeptideDBImpl* next)
	: m_file(data), m_chem_comp(m_file.first_datablock()["chem_comp"]), m_next(next)
{
	for (auto& chem_comp: m_chem_comp)
	{
		string group, three_letter_code;
		
		cif::tie(group, three_letter_code) = chem_comp.get("group", "three_letter_code");
		
		if (group == "peptide" or group == "M-peptide" or group == "P-peptide")
			m_known_peptides.insert(three_letter_code);
		else if (group == "DNA" or group == "RNA")
			m_known_bases.insert(three_letter_code);
	}
}

string PeptideDBImpl::formula_for(string res) const
{
	string result;

	ba::to_upper(res);
	
	for (auto& db: m_file)
	{
		if (db.name() != "comp_" + res)
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
		if (m_next != nullptr)
			result = m_next->formula_for(res);
		else
		{
			const char* clibd_mon = getenv("CLIBD_MON");
			if (clibd_mon == nullptr)
				throw runtime_error("Cannot locate peptide list, please souce the CCP4 environment");
			
			fs::path resFile = fs::path(clibd_mon) / ba::to_lower_copy(res.substr(0, 1)) / (res + ".cif");
			if (fs::exists(resFile))
			{
				fs::ifstream file(resFile);
				if (file.is_open())
				{
					try
					{
						cif::file cf(file);
						
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
	const char* clibd_mon = getenv("CLIBD_MON");
	if (clibd_mon == nullptr)
		throw runtime_error("Cannot locate peptide list, please souce the CCP4 environment");
	
	fs::path db = fs::path(clibd_mon) / "list" / "mon_lib_list.cif";

	PushDictionary(db);
	
	sInstance = this;
}

void PeptideDB::PushDictionary(boost::filesystem::path dict)
{
	if (not fs::exists(dict))
		throw runtime_error("file not found: " + dict.string());

	fs::ifstream file(dict);
	if (not file.is_open())
		throw runtime_error("Could not open peptide list " + dict.string());

	mImpl = new PeptideDBImpl(file, mImpl);
}

void PeptideDB::PopDictionary()
{
	if (mImpl != nullptr)
	{
		auto i = mImpl;
		mImpl = i->m_next;
		i->m_next = nullptr;
		delete i;
	}
}

PeptideDB::~PeptideDB()
{
	delete mImpl;
}

bool PeptideDB::IsKnownPeptide(const string& res_name) const
{
	return mImpl->m_known_peptides.count(res_name) > 0;
}

bool PeptideDB::IsKnownBase(const string& res_name) const
{
	return mImpl->m_known_bases.count(res_name) > 0;
}

string PeptideDB::GetNameForResidue(const string& res_name) const
{
	return mImpl->name_for(res_name);
}

string PeptideDB::GetFormulaForResidue(const string& res_name) const
{
	return mImpl->formula_for(res_name);
}

string PeptideDB::Unalias(const string& res_name) const
{
	return mImpl->unalias(res_name);
}
