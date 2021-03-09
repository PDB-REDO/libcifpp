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

#include <map>
#include <set>
#include <regex>

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/concepts.hpp>    // output_filter
#include <boost/iostreams/operations.hpp>  // put

#include "cif++/Cif2PDB.hpp"
#include "cif++/AtomType.hpp"
#include "cif++/Compound.hpp"

namespace ba = boost::algorithm;
namespace io = boost::iostreams;

using cif::Datablock;
using cif::Category;
using cif::Row;

// --------------------------------------------------------------------
// FillOutLineFilter is used to make sure all lines in PDB files
// are at filled out with spaces to be as much as 80 characters wide.

class FillOutLineFilter : public io::output_filter
{
 public:
	FillOutLineFilter()
		: mLineCount(0), mColumnCount(0) {}
	
    template<typename Sink>
    bool put(Sink& dest, int c)
    {
    	bool result = true;
    	
    	if (c == '\n')
    	{
    		for (int i = mColumnCount; result and i < 80; ++i)
    			result = io::put(dest, ' ');
    	}
    	
    	if (result)
    		result = io::put(dest, c);
    	
    	if (result)
    	{
    		if (c == '\n')
    		{
    			mColumnCount = 0;
    			++mLineCount;
    		}
    		else
    			++mColumnCount;
    	}
    	
    	return result;
    }

    template<typename Sink>
    void close(Sink&)
    {
    	mLineCount = 0;
    	mColumnCount = 0;
    }
    
    int GetLineCount() const 		{ return mLineCount; }
	
  private:
	int	mLineCount;
	int mColumnCount;
};

// --------------------------------------------------------------------
// conversion routines between cif and pdb format

std::string cif2pdbDate(const std::string& d)
{
	const std::regex rx(R"((\d{4})-(\d{2})(?:-(\d{2}))?)");
	const char* kMonths[12] = {
		"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
	};

	std::smatch m;
	std::string result;
	
	if (std::regex_match(d, m, rx))
	{
		int year = std::stoi(m[1].str());
		int month = std::stoi(m[2].str());
		
		if (m[3].matched)
			result = (boost::format("%02.2d-%3.3s-%02.2d") % stoi(m[3].str()) % kMonths[month - 1] % (year % 100)).str();
		else
			result = (boost::format("%3.3s-%02.2d") % kMonths[month - 1] % (year % 100)).str();
	}
	
	return result;
}

std::string cif2pdbAuth(std::string name)
{
	const std::regex rx(R"(([^,]+), (\S+))");
	
	std::smatch m;
	if (std::regex_match(name, m, rx))
		name = m[2].str() + m[1].str();
		
	return name;
}

std::string cif2pdbSymmetry(std::string s)
{
	auto i = s.rfind('_');
	if (i != std::string::npos)
		s.erase(i, 1);
	return s;
}

std::string cif2pdbAtomName(std::string name, std::string resName, Datablock& db)
{
	if (name.length() < 4)
	{
		for (auto r: db["atom_site"].find(cif::Key("label_atom_id") == name and cif::Key("label_comp_id") == resName))
		{
			std::string element = r["type_symbol"].as<std::string>();
			
			if (element.length() == 1 or not cif::iequals(name, element))
				name.insert(name.begin(), ' ');
			
			break;
		}
	}	
	
	return name;
}

enum SoftwareType { eRefinement, eDataScaling, eDataExtraction, eDataReduction, ePhasing };

std::string cifSoftware(Datablock& db, SoftwareType sw)
{
	std::string result = "NULL";
	
	try
	{
		switch (sw)
		{
			case eRefinement:		result = db["computing"][cif::Key("entry_id") == db.getName()]["structure_refinement"].as<std::string>();	break;
			case eDataScaling:		result = db["computing"][cif::Key("entry_id") == db.getName()]["pdbx_data_reduction_ds"].as<std::string>(); break;
			case eDataReduction:	result = db["computing"][cif::Key("entry_id") == db.getName()]["pdbx_data_reduction_ii"].as<std::string>(); break;
			default: break;
		}
		
		if (result.empty() or result == "NULL")
		{
			auto& software = db["software"];
		
			Row r;
			
			switch (sw)
			{
				case eRefinement:		r = software[cif::Key("classification") == "refinement"];		break;
				case eDataScaling:		r = software[cif::Key("classification") == "data scaling"];		break;
				case eDataExtraction:	r = software[cif::Key("classification") == "data extraction"];	break;
				case eDataReduction:	r = software[cif::Key("classification") == "data reduction"];	break;
				case ePhasing:			r = software[cif::Key("classification") == "phasing"];			break;
			}
			
			result = r["name"].as<std::string>() + " " + r["version"].as<std::string>();
		}

		ba::trim(result);
		ba::to_upper(result);
		
		if (result.empty())
			result = "NULL";
	}
	catch (...) {}
	
	return result;
}

// Map asym ID's back to PDB Chain ID's
std::vector<std::string> MapAsymIDs2ChainIDs(const std::vector<std::string>& asymIDs, Datablock& db)
{
	std::set<std::string> result;
	
	for (auto asym: asymIDs)
	{
		for (auto r: db["pdbx_poly_seq_scheme"].find(cif::Key("asym_id") == asym))
		{
			result.insert(r["pdb_strand_id"].as<std::string>());
			break;
		}
		
		for (auto r: db["pdbx_nonpoly_scheme"].find(cif::Key("asym_id") == asym))
		{
			result.insert(r["pdb_strand_id"].as<std::string>());
			break;
		}
	}
	
	return { result.begin(), result.end() };
}

// support for wrapping text using a 'continuation marker'
int WriteContinuedLine(std::ostream& pdbFile, std::string header, int& count, int cLen, std::string text, int lStart = 0)
{
	if (lStart == 0)
	{
		if (cLen == 0)
			lStart = header.length() + 1;
		else
			lStart = header.length() + cLen;
	}

	int maxLength = 80 - lStart - 1;

	std::vector<std::string> lines = cif::wordWrap(text, maxLength);

	for (auto& line: lines)
	{
		ba::to_upper(line);

		pdbFile << header;
		
		if (++count <= 1 or cLen == 0)
		{
			pdbFile << std::string(lStart - header.length(), ' ');
			if (count == 1)
				lStart = header.length() + cLen + 1;
		}
		else
			pdbFile << std::fixed << std::setw(cLen) << std::right << count << ' ';
		
		pdbFile << line << std::endl;
	}
	
	return lines.size();

}

int WriteOneContinuedLine(std::ostream& pdbFile, std::string header, int cLen, std::string line, int lStart = 0)
{
	int count = 0;
	return WriteContinuedLine(pdbFile, header, count, cLen, line, lStart);
}

int WriteCitation(std::ostream& pdbFile, Datablock& db, Row r, int reference)
{
	int result = 0;
	
	std::string s1;
	
	if (reference > 0)
	{
		pdbFile << "REMARK   1 REFERENCE " << std::to_string(reference) << std::endl;
		result = 1;
		s1 = "REMARK   1  ";
	}
	else
		s1 = "JRNL        ";
	
	std::string id, title, pubname, volume, astm, country, issn, csd, publ, pmid, doi, pageFirst, pageLast, year;
	
	cif::tie(id, title, pubname, volume, astm, country, issn, csd, publ, pmid, doi, pageFirst, pageLast, year) =
		r.get("id", "title", "journal_abbrev", "journal_volume", "journal_id_ASTM", "country", "journal_id_ISSN", 
			  "journal_id_CSD", "book_publisher", "pdbx_database_id_PubMed", "pdbx_database_id_DOI",
			  "page_first", "page_last", "year");
	
	std::vector<std::string> authors;
	for (auto r1: db["citation_author"].find(cif::Key("citation_id") == id))
		authors.push_back(cif2pdbAuth(r1["name"].as<std::string>()));

	if (not authors.empty())
		result += WriteOneContinuedLine(pdbFile, s1 + "AUTH", 2, ba::join(authors, ","), 19);
	
	result += WriteOneContinuedLine(pdbFile, s1 + "TITL", 2, title, 19);
	
	if (not pubname.empty())
	{
		ba::to_upper(pubname);
	
		const std::string kRefHeader = s1 + "REF %2.2d %-28.28s  %2.2s%4.4d %5.5d %4.4d";
		pdbFile << (boost::format(kRefHeader)
					% ""	// continuation
					% pubname
					% (volume.empty() ? "" : "V.")
					% volume
					% pageFirst
					% year).str()
				<< std::endl;
		++result;
	}
	
	if (not issn.empty())
	{
		const std::string kRefHeader = s1 + "REFN                   ISSN %-25.25s";
		pdbFile << (boost::format(kRefHeader) % issn).str() << std::endl;
		++result;
	}

//		if (not issn.empty() or astm.empty())
//		{
////    0         1         2         3         4         5         6         7         8
////    HEADER    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxDDDDDDDDD   IIII
//const char kRefHeader[] =
//         "REMARK   1  REFN    %4.4s %-6.6s  %2.2s %-25.25s";
//			
//			pdbFile << (boost::format(kRefHeader)
//						% (astm.empty() ? "" : "ASTN")
//						% astm
//						% country
//						% issn).str()
//					<< std::endl;
//		}

	if (not pmid.empty())
	{
		const std::string kPMID = s1 + "PMID   %-60.60s ";
		pdbFile << (boost::format(kPMID) % pmid).str() << std::endl;
		++result;
	}

	if (not doi.empty())
	{
		const std::string kDOI = s1 + "DOI    %-60.60s ";
		pdbFile << (boost::format(kDOI) % doi).str() << std::endl;
		++result;
	}
	
	return result;
}

void WriteHeaderLines(std::ostream& pdbFile, Datablock& db)
{
	//    0         1         2         3         4         5         6         7         8
	//    HEADER    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxDDDDDDDDD   IIII
	const char kHeader[] =
		 "HEADER    %1$-40.40s"                           "%2$-9.9s""   %3$-4.4s";

	// HEADER
	
	std::string keywords;
	auto& cat1 = db["struct_keywords"];
	
	for (auto r: cat1)
	{
		keywords = r["pdbx_keywords"].as<std::string>();
		break;
	}
	
	std::string date;
	for (auto r: db["pdbx_database_status"])
	{
		date = r["recvd_initial_deposition_date"].as<std::string>();
		if (date.empty())
			continue;
		date = cif2pdbDate(date);
		break;
	}
	
	if (date.empty())
	{
		for (auto r: db["database_PDB_rev"])
		{
			date = r["date_original"].as<std::string>();
			if (date.empty())
				continue;
			date = cif2pdbDate(date);
			break;
		}
	}
	
	pdbFile << (boost::format(kHeader) % keywords % date % db.getName()).str() << std::endl;
	
	// TODO: implement
	// OBSLTE (skip for now)
	
	// TITLE
	for (auto r: db["struct"])
	{
		std::string title = r["title"].as<std::string>();
		ba::trim(title);
		WriteOneContinuedLine(pdbFile, "TITLE   ", 2, title);
		break;
	}
	
	// COMPND
	using namespace std::placeholders;
	
	int molID = 0;
	std::vector<std::string> cmpnd;
	
	for (auto r: db["entity"])
	{
		if (r["type"] != "polymer")
			continue;
		
		std::string entityID = r["id"].as<std::string>();
		
		++molID;
		cmpnd.push_back("MOL_ID: " + std::to_string(molID));
		
		std::string molecule = r["pdbx_description"].as<std::string>();
		cmpnd.push_back("MOLECULE: " + molecule);

		auto poly = db["entity_poly"].find(cif::Key("entity_id") == entityID);
		if (not poly.empty())
		{
			std::string chains = poly.front()["pdbx_strand_id"].as<std::string>();
			ba::replace_all(chains, ",", ", ");
			cmpnd.push_back("CHAIN: " + chains);
		}

		std::string fragment = r["pdbx_fragment"].as<std::string>();
		if (not fragment.empty())
			cmpnd.push_back("FRAGMENT: " + fragment);
		
		for (auto sr: db["entity_name_com"].find(cif::Key("entity_id") == entityID))
		{
			std::string syn = sr["name"].as<std::string>();
			if (not syn.empty())
				cmpnd.push_back("SYNONYM: " + syn);
		}
		
		std::string mutation = r["pdbx_mutation"].as<std::string>();
		if (not mutation.empty())
			cmpnd.push_back("MUTATION: " + mutation);
		
		std::string ec = r["pdbx_ec"].as<std::string>();
		if (not ec.empty())
			cmpnd.push_back("EC: " + ec);
		
		if (r["src_method"] == "man" or r["src_method"] == "syn")
			cmpnd.push_back("ENGINEERED: YES");

		std::string details = r["details"].as<std::string>();
		if (not details.empty())
			cmpnd.push_back("OTHER_DETAILS: " + details);
	}

	WriteOneContinuedLine(pdbFile, "COMPND ", 3, ba::join(cmpnd, ";\n"));

	// SOURCE
	
	molID = 0;
	std::vector<std::string> source;
	
	for (auto r: db["entity"])
	{
		if (r["type"] != "polymer")
			continue;
		
		std::string entityID = r["id"].as<std::string>();
		
		++molID;
		source.push_back("MOL_ID: " + std::to_string(molID));
		
		if (r["src_method"] == "syn")
			source.push_back("SYNTHETIC: YES");
		
		auto& gen = db["entity_src_gen"];
		const std::pair<const char*,const char*> kGenSourceMapping[] = {
			{ "gene_src_common_name",			"ORGANISM_COMMON" },
			{ "pdbx_gene_src_gene",				"GENE" },
			{ "gene_src_strain",				"STRAIN" },
			{ "pdbx_gene_src_cell_line",		"CELL_LINE" },
			{ "pdbx_gene_src_organelle",		"ORGANELLE" },
			{ "pdbx_gene_src_cellular_location","CELLULAR_LOCATION" },
			{ "pdbx_gene_src_scientific_name",	"ORGANISM_SCIENTIFIC" },
			{ "pdbx_gene_src_ncbi_taxonomy_id",	"ORGANISM_TAXID" },
			{ "pdbx_host_org_scientific_name",	"EXPRESSION_SYSTEM" },
			{ "pdbx_host_org_ncbi_taxonomy_id",	"EXPRESSION_SYSTEM_TAXID" },
			{ "pdbx_host_org_strain",			"EXPRESSION_SYSTEM_STRAIN" },
			{ "pdbx_host_org_variant",			"EXPRESSION_SYSTEM_VARIANT" },
			{ "pdbx_host_org_cellular_location","EXPRESSION_SYSTEM_CELLULAR_LOCATION" },
			{ "pdbx_host_org_vector_type",		"EXPRESSION_SYSTEM_VECTOR_TYPE" },
			{ "pdbx_host_org_vector",			"EXPRESSION_SYSTEM_VECTOR" },
			{ "pdbx_host_org_gene",				"EXPRESSION_SYSTEM_GENE" },
			{ "plasmid_name",					"EXPRESSION_SYSTEM_PLASMID" },
			{ "details",						"OTHER_DETAILS" }
		};

		for (auto gr: gen.find(cif::Key("entity_id") == entityID))
		{
			for (auto m: kGenSourceMapping)
			{
				std::string cname, sname;
				tie(cname, sname) = m;
				
				std::string s = gr[cname].as<std::string>();
				if (not s.empty())
					source.push_back(sname + ": " + s);
			}
		}
		
		auto& nat = db["entity_src_nat"];
		const std::pair<const char*, const char*> kNatSourceMapping[] = {
			{ "common_name",				"ORGANISM_COMMON" },
			{ "strain",						"STRAIN" },
			{ "pdbx_organism_scientific",	"ORGANISM_SCIENTIFIC" },
			{ "pdbx_ncbi_taxonomy_id",		"ORGANISM_TAXID" },
			{ "pdbx_cellular_location",		"CELLULAR_LOCATION" },
			{ "pdbx_plasmid_name",			"PLASMID" },
			{ "pdbx_organ",					"ORGAN" },
			{ "details",					"OTHER_DETAILS" }
		};
		
		for (auto nr: nat.find(cif::Key("entity_id") == entityID))
		{
			for (auto m: kNatSourceMapping)
			{
				std::string cname, sname;
				tie(cname, sname) = m;
				
				std::string s = nr[cname].as<std::string>();
				if (not s.empty())
					source.push_back(sname + ": " + s);
			}
		}
	}

	WriteOneContinuedLine(pdbFile, "SOURCE ", 3, ba::join(source, ";\n"));
	
	// KEYWDS
	
	keywords.clear();
	for (auto r: cat1)
	{
		if (not r["text"].empty())
			keywords += r["text"].as<std::string>();
		else
			keywords += r["pdbx_keywords"].as<std::string>();
	}
	
	if (not keywords.empty())
		WriteOneContinuedLine(pdbFile, "KEYWDS  ", 2, keywords);
	
	// EXPDTA
	
	auto& dbexpt = db["exptl"];
	if (not dbexpt.empty())
	{
		std::vector<std::string> method;
		for (auto r: dbexpt)
			method.push_back(r["method"].as<std::string>());
		if (not method.empty())
			WriteOneContinuedLine(pdbFile, "EXPDTA  ", 2, ba::join(method, "; "));
	}
	
	// NUMMDL
	// TODO...

	// MDLTYP
	// TODO...

	// AUTHOR
	std::vector<std::string> authors;
	for (auto r: db["audit_author"])
		authors.push_back(cif2pdbAuth(r["name"].as<std::string>()));
	if (not authors.empty())
		WriteOneContinuedLine(pdbFile, "AUTHOR  ", 2, ba::join(authors, ","));
}

void WriteTitle(std::ostream& pdbFile, Datablock& db)
{
	WriteHeaderLines(pdbFile, db);
	
	// REVDAT
	boost::format kRevDat("REVDAT %3.3d%2.2s %9.9s %4.4s    %1.1d      ");
	auto& cat2 = db["database_PDB_rev"];
	std::vector<Row> rev(cat2.begin(), cat2.end());
	sort(rev.begin(), rev.end(), [](Row a, Row b) -> bool { return a["num"].as<int>() > b["num"].as<int>(); });
	for (auto r: rev)
	{
		int revNum, modType;
		std::string date, replaces;
		
		cif::tie(revNum, modType, date, replaces) = r.get("num", "mod_type", "date", "replaces");
		
		date = cif2pdbDate(date);
		
		std::vector<std::string> types;
		
		for (auto r1: db["database_PDB_rev_record"].find(cif::Key("rev_num") == revNum))
			types.push_back(r1["type"].as<std::string>());
		
		int continuation = 0;
		do
		{
			std::string cs = ++continuation > 1 ? std::to_string(continuation) : std::string();
			
			pdbFile << (kRevDat % revNum % cs % date % db.getName() % modType).str();
			for (size_t i = 0; i < 4; ++i)
				pdbFile << (boost::format(" %-6.6s") % (i < types.size() ? types[i] : std::string())).str();
			pdbFile << std::endl;
			
			if (types.size() > 4)
				types.erase(types.begin(), types.begin() + 4);
			else
				types.clear();
		}
		while (types.empty() == false);
	}
	
	// SPRSDE
	// TODO...
	
	// JRNL
	for (auto r: db["citation"])
	{
		WriteCitation(pdbFile, db, r, 0);
		break;
	}
}

void WriteRemark1(std::ostream& pdbFile, Datablock& db)
{
	int reference = 0;
	
	for (auto r: db["citation"])
	{
		if (reference > 0)
		{
			if (reference == 1)
				pdbFile << "REMARK   1" << std::endl;

			WriteCitation(pdbFile, db, r, reference);
		}

		++reference;
	}
}

void WriteRemark2(std::ostream& pdbFile, Datablock& db)
{
	auto& refine = db["refine"];
	if (refine.empty())
	{
		pdbFile << "REMARK   2" << std::endl
				<< "REMARK   2 RESOLUTION. NOT APPLICABLE." << std::endl;
	}
	else
	{
		try
		{
			float resHigh = refine.front()["ls_d_res_high"].as<float>();
			
			boost::format kREMARK2("REMARK   2 RESOLUTION. %7.2f ANGSTROMS.");

			pdbFile << "REMARK   2" << std::endl
					<< (kREMARK2 % resHigh) << std::endl;
		}
		catch (...) { /* skip it */ }
	}
}

// --------------------------------------------------------------------
// Code to help format RERMARK 3 data

class FBase
{
  public:
	virtual ~FBase() {}

	virtual void out(std::ostream& os) = 0;

  protected:
	FBase(Row r, const char* f)
		: mRow(r), mField(f) {}
	FBase(Category& cat, cif::Condition&& cond, const char* f)
		: mField(f)
	{
		auto r = cat.find(std::move(cond));
		if (not r.empty())
			mRow = r.front();
	}
	
	Row mRow;
	const char* mField;
};

class Fi : public FBase
{
  public:
	Fi(Row r, const char* f) : FBase(r, f) {}
	Fi(Category& cat, cif::Condition&& cond, const char* f) : FBase(cat, std::move(cond), f) {}
	
	virtual void out(std::ostream& os)
	{
		std::string s = mRow[mField].as<std::string>();
		if (s.empty())
		{
			os << "NULL";
			if (os.width() > 4)
				os << std::string(os.width() - 4, ' ');
		}
		else
			os << stol(s);
	}
};

class Ff : public FBase
{
  public:
	Ff(Row r, const char* f) : FBase(r, f) {}
	Ff(Category& cat, cif::Condition&& cond, const char* f) : FBase(cat, std::move(cond), f) {}
	
	virtual void out(std::ostream& os)
	{
		std::string s = mRow[mField].as<std::string>();
		if (s.empty())
		{
			os << "NULL";
			if (os.width() > 4)
				os << std::string(os.width() - 4, ' ');
		}
		else
        {
            try
            {
                os << std::stod(s);
            }
            catch (const std::exception& ex)
            {
                std::cerr << "Failed to write '" << s << "' as a double, this indicates an error in the code for writing PDB files" << std::endl;
                os << s;
            }
        }
	}
};

class Fs : public FBase
{
  public:
	Fs(Row r, const char* f, int remarkNr = 3) : FBase(r, f), mNr(remarkNr) {}
	Fs(Category& cat, cif::Condition&& cond, const char* f, int remarkNr = 3) : FBase(cat, std::move(cond), f), mNr(remarkNr) {}
	
	virtual void out(std::ostream& os)
	{
		std::string s = mRow[mField].as<std::string>();
		size_t width = os.width();
		
		if (s.empty())
		{
			os << "NULL";
			if (os.width() > 4)
				os << std::string(width - 4, ' ');
		}
		else if (width == 0 or s.length() <= width)
			os << s;
		else
		{
			os << std::endl;
			
			std::stringstream ss;
			ss << "REMARK " << std::setw(3) << std::right << mNr << ' ';
			WriteOneContinuedLine(os, ss.str(), 0, s);
		}
	}
	
	int mNr = 3;
};

std::ostream& operator<<(std::ostream& os, FBase&& fld)
{
	fld.out(os);
	return os;
}

template<int N>
struct RM
{
	RM(const char* desc, int width = 0, int precision = 6) : mDesc(desc), mWidth(width), mPrecision(precision) {}
	const char* mDesc;
	int mWidth, mPrecision;
};

typedef RM<3> RM3;

template<int N>
std::ostream& operator<<(std::ostream& os, RM<N>&& rm)
{
	os << "REMARK " << std::setw(3) << std::right << N << " " << rm.mDesc << (rm.mWidth > 0 ? std::left : std::right) << std::fixed << std::setw(abs(rm.mWidth)) << std::setprecision(rm.mPrecision);
	return os;
}

struct SEP
{
	SEP(const char* txt, int width, int precision = 6) : mText(txt), mWidth(width), mPrecision(precision) {}
	const char* mText;
	int mWidth, mPrecision;
};

std::ostream& operator<<(std::ostream& os, SEP&& sep)
{
	os << sep.mText << (sep.mWidth > 0 ? std::left : std::right) << std::fixed << std::setw(abs(sep.mWidth)) << std::setprecision(sep.mPrecision);
	return os;
}

// --------------------------------------------------------------------

void WriteRemark3BusterTNT(std::ostream& pdbFile, Datablock& db)
{
	auto refine = db["refine"].front();
	auto ls_shell = db["refine_ls_shell"].front();
	auto hist = db["refine_hist"].front();
	auto reflns = db["reflns"].front();
	auto analyze = db["refine_analyze"].front();
	auto& ls_restr = db["refine_ls_restr"];
//	auto ls_restr_ncs = db["refine_ls_restr_ncs"].front();
//	auto pdbx_xplor_file = db["pdbx_xplor_file"].front();
//	auto pdbx_refine = db["pdbx_refine"].front();
	
	pdbFile	<< RM3("") << std::endl
			<< RM3(" DATA USED IN REFINEMENT.") << std::endl
			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_high") << std::endl
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_low") << std::endl
			<< RM3("  DATA CUTOFF            (SIGMA(F)) : ", 6, 3)	<< Ff(refine, "pdbx_ls_sigma_F") << std::endl
			<< RM3("  COMPLETENESS FOR RANGE        (%) : ", 6, 1)	<< Ff(refine, "ls_percent_reflns_obs") << std::endl
			<< RM3("  NUMBER OF REFLECTIONS             : ", 12, 6)	<< Fi(refine, "ls_number_reflns_obs") << std::endl
	
			<< RM3("") << std::endl
			<< RM3(" FIT TO DATA USED IN REFINEMENT.") << std::endl
			<< RM3("  CROSS-VALIDATION METHOD          : ")	<< Fs(refine, "pdbx_ls_cross_valid_method") << std::endl
			<< RM3("  FREE R VALUE TEST SET SELECTION  : ")	<< Fs(refine, "pdbx_R_Free_selection_details") << std::endl
			<< RM3("  R VALUE     (WORKING + TEST SET) : ", 7, 3)	<< Ff(refine, "ls_R_factor_obs") << std::endl
			<< RM3("  R VALUE            (WORKING SET) : ", 7, 3)	<< Ff(refine, "ls_R_factor_R_work") << std::endl
			<< RM3("  FREE R VALUE                     : ", 7, 3)	<< Ff(refine, "ls_R_factor_R_free") << std::endl
			<< RM3("  FREE R VALUE TEST SET SIZE   (%) : ", 7, 3)	<< Ff(refine, "ls_percent_reflns_R_free") << std::endl
			<< RM3("  FREE R VALUE TEST SET COUNT      : ", 12, 6)	<< Fi(refine, "ls_number_reflns_R_free") << std::endl
			<< RM3("  ESTIMATED ERROR OF FREE R VALUE  : ", 7, 3)	<< Ff(refine, "ls_R_factor_R_free_error") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" FIT IN THE HIGHEST RESOLUTION BIN.") << std::endl
			<< RM3("  TOTAL NUMBER OF BINS USED               : ", 12, 6)	<< Fi(ls_shell, "pdbx_total_number_of_bins_used") << std::endl
			<< RM3("  BIN RESOLUTION RANGE HIGH   (ANGSTROMS) : ", 5, 2)	<< Ff(ls_shell, "d_res_high") << std::endl
			<< RM3("  BIN RESOLUTION RANGE LOW    (ANGSTROMS) : ", 5, 2)	<< Ff(ls_shell, "d_res_low") << std::endl
			<< RM3("  BIN COMPLETENESS     (WORKING+TEST) (%) : ", 6, 2)	<< Ff(ls_shell, "percent_reflns_obs") << std::endl
			<< RM3("  REFLECTIONS IN BIN (WORKING + TEST SET) : ", 12, 6)	<< Fi(ls_shell, "number_reflns_all") << std::endl
			<< RM3("  BIN R VALUE        (WORKING + TEST SET) : ", 8, 4)	<< Ff(ls_shell, "R_factor_all") << std::endl
			<< RM3("  REFLECTIONS IN BIN        (WORKING SET) : ", 12, 6)	<< Fi(ls_shell, "number_reflns_R_work") << std::endl
			<< RM3("  BIN R VALUE               (WORKING SET) : ", 8, 4)	<< Ff(ls_shell, "R_factor_R_work") << std::endl
			<< RM3("  BIN FREE R VALUE                        : ", 8, 4)	<< Ff(ls_shell, "R_factor_R_free") << std::endl
			<< RM3("  BIN FREE R VALUE TEST SET SIZE      (%) : ", 6, 2)	<< Ff(ls_shell, "percent_reflns_R_free") << std::endl
			<< RM3("  BIN FREE R VALUE TEST SET COUNT         : ", 12, 7)	<< Fi(ls_shell, "number_reflns_R_free") << std::endl
			<< RM3("  ESTIMATED ERROR OF BIN FREE R VALUE     : ", 7, 3)	<< Ff(ls_shell, "R_factor_R_free_error") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.") << std::endl
			<< RM3("  PROTEIN ATOMS            : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_protein") << std::endl
			<< RM3("  NUCLEIC ACID ATOMS       : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_nucleic_acid") << std::endl
			<< RM3("  HETEROGEN ATOMS          : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_ligand") << std::endl
			<< RM3("  SOLVENT ATOMS            : ", 12, 6)	<< Fi(hist, "number_atoms_solvent") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" B VALUES.") << std::endl
//			<< RM3("  B VALUE TYPE                      : ")		<< Fs(refine, "pdbx_TLS_residual_ADP_flag") << std::endl
			<< RM3("  FROM WILSON PLOT           (A**2) : ", 7, 2)	<< Ff(reflns, "B_iso_Wilson_estimate") << std::endl
			<< RM3("  MEAN B VALUE      (OVERALL, A**2) : ", 7, 2)	<< Ff(refine, "B_iso_mean") << std::endl
			
			<< RM3("  OVERALL ANISOTROPIC B VALUE.") << std::endl
			<< RM3("   B11 (A**2) : ", -8, 5)	<< Ff(refine, "aniso_B[1][1]") << std::endl
			<< RM3("   B22 (A**2) : ", -8, 5)	<< Ff(refine, "aniso_B[2][2]") << std::endl
			<< RM3("   B33 (A**2) : ", -8, 5)	<< Ff(refine, "aniso_B[3][3]") << std::endl
			<< RM3("   B12 (A**2) : ", -8, 5)	<< Ff(refine, "aniso_B[1][2]") << std::endl
			<< RM3("   B13 (A**2) : ", -8, 5)	<< Ff(refine, "aniso_B[1][3]") << std::endl
			<< RM3("   B23 (A**2) : ", -8, 5)	<< Ff(refine, "aniso_B[2][3]") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" ESTIMATED COORDINATE ERROR.") << std::endl
			<< RM3("  ESD FROM LUZZATI PLOT                    (A) : ", 7, 3)	<< Ff(analyze, "Luzzati_coordinate_error_obs") << std::endl
			<< RM3("  DPI (BLOW EQ-10) BASED ON R VALUE        (A) : ", 5, 3)	<< Ff(refine, "pdbx_overall_SU_R_Blow_DPI") << std::endl
			<< RM3("  DPI (BLOW EQ-9) BASED ON FREE R VALUE    (A) : ", 5, 3)	<< Ff(refine, "pdbx_overall_SU_R_free_Blow_DPI") << std::endl
			<< RM3("  DPI (CRUICKSHANK) BASED ON R VALUE       (A) : ", 5, 3)	<< Ff(refine, "overall_SU_R_Cruickshank_DPI") << std::endl
			<< RM3("  DPI (CRUICKSHANK) BASED ON FREE R VALUE  (A) : ", 5, 3)	<< Ff(refine, "pdbx_overall_SU_R_free_Cruickshank_DPI") << std::endl
			
			<< RM3("") << std::endl
			<< RM3("  REFERENCES: BLOW, D. (2002) ACTA CRYST D58, 792-797") << std::endl
			<< RM3("              CRUICKSHANK, D.W.J. (1999) ACTA CRYST D55, 583-601") << std::endl

			<< RM3("") << std::endl
			<< RM3("  CORRELATION COEFFICIENTS.") << std::endl
			<< RM3("  CORRELATION COEFFICIENT FO-FC      : ", 5, 3)	<< Ff(refine, "correlation_coeff_Fo_to_Fc") << std::endl
			<< RM3("  CORRELATION COEFFICIENT FO-FC FREE : ", 5, 3)	<< Ff(refine, "correlation_coeff_Fo_to_Fc_free") << std::endl
			
			<< RM3("") << std::endl
			<< RM3("  NUMBER OF GEOMETRIC FUNCTION TERMS DEFINED : 15") << std::endl
			<< RM3("  TERM                          COUNT    WEIGHT   FUNCTION.") << std::endl
			<< RM3("   BOND LENGTHS              : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_bond_d", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_bond_d", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_bond_d", "pdbx_restraint_function") << std::endl
			<< RM3("   BOND ANGLES               : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_angle_deg", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_angle_deg", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_angle_deg", "pdbx_restraint_function") << std::endl
			<< RM3("   TORSION ANGLES            : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_dihedral_angle_d", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_dihedral_angle_d", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_dihedral_angle_d", "pdbx_restraint_function") << std::endl
			<< RM3("   TRIGONAL CARBON PLANES    : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_trig_c_planes", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_trig_c_planes", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_trig_c_planes", "pdbx_restraint_function") << std::endl
			<< RM3("   GENERAL PLANES            : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_gen_planes", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_gen_planes", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_gen_planes", "pdbx_restraint_function") << std::endl
			<< RM3("   ISOTROPIC THERMAL FACTORS : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_it", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_it", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_it", "pdbx_restraint_function") << std::endl
			<< RM3("   BAD NON-BONDED CONTACTS   : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_nbd", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_nbd", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_nbd", "pdbx_restraint_function") << std::endl
			<< RM3("   IMPROPER TORSIONS         : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_improper_torsion", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_improper_torsion", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_improper_torsion", "pdbx_restraint_function") << std::endl
			<< RM3("   PSEUDOROTATION ANGLES     : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_pseud_angle", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_pseud_angle", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_pseud_angle", "pdbx_restraint_function") << std::endl
			<< RM3("   CHIRAL IMPROPER TORSION   : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_chiral_improper_torsion", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_chiral_improper_torsion", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_chiral_improper_torsion", "pdbx_restraint_function") << std::endl
			<< RM3("   SUM OF OCCUPANCIES        : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_sum_occupancies", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_sum_occupancies", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_sum_occupancies", "pdbx_restraint_function") << std::endl
			<< RM3("   UTILITY DISTANCES         : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_utility_distance", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_utility_distance", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_utility_distance", "pdbx_restraint_function") << std::endl
			<< RM3("   UTILITY ANGLES            : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_utility_angle", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_utility_angle", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_utility_angle", "pdbx_restraint_function") << std::endl
			<< RM3("   UTILITY TORSION           : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_utility_torsion", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_utility_torsion", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_utility_torsion", "pdbx_restraint_function") << std::endl
			<< RM3("   IDEAL-DIST CONTACT TERM   : ", 7, 0)	<< Ff(ls_restr, cif::Key("type") == "t_ideal_dist_contact", "number")
										<< SEP("; ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_ideal_dist_contact", "weight")
										<< SEP("; ", 12)	<< Fs(ls_restr, cif::Key("type") == "t_ideal_dist_contact", "pdbx_restraint_function") << std::endl

			
			<< RM3("") << std::endl
			<< RM3(" RMS DEVIATIONS FROM IDEAL VALUES.") << std::endl
			<< RM3("  BOND LENGTHS                       (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "t_bond_d", "dev_ideal") << std::endl
			<< RM3("  BOND ANGLES                  (DEGREES) : ", 7, 2)	<< Ff(ls_restr, cif::Key("type") == "t_angle_deg", "dev_ideal") << std::endl
			<< RM3("  PEPTIDE OMEGA TORSION ANGLES (DEGREES) : ", 7, 2)	<< Ff(ls_restr, cif::Key("type") == "t_omega_torsion", "dev_ideal") << std::endl
			<< RM3("  OTHER TORSION ANGLES         (DEGREES) : ", 7, 2)	<< Ff(ls_restr, cif::Key("type") == "t_other_torsion", "dev_ideal") << std::endl;

	auto& tls = db["pdbx_refine_tls"];	

	pdbFile << RM3("") << std::endl
			<< RM3(" TLS DETAILS") << std::endl
			<< RM3("  NUMBER OF TLS GROUPS  : ") << (tls.size() ? std::to_string(tls.size()) : "NULL") << std::endl;
	
	for (auto t: tls)
	{
		std::string id = t["id"].as<std::string>();
		auto g = db["pdbx_refine_tls_group"][cif::Key("refine_tls_id") == id];
		
		pdbFile << RM3("") << std::endl
				<< RM3("  TLS GROUP : ") << id << std::endl
				<< RM3("   SELECTION: ") << Fs(g, "selection_details") << std::endl;
		
		pdbFile << RM3("   ORIGIN FOR THE GROUP (A):", -9, 4) << Ff(t, "origin_x")
				<< SEP("", -9, 4) << Ff(t, "origin_y")
				<< SEP("", -9, 4) << Ff(t, "origin_z") << std::endl
				<< RM3("   T TENSOR") << std::endl
				<< RM3("     T11:", -9, 4) << Ff(t, "T[1][1]") << SEP(" T22:", -9, 4) << Ff(t, "T[2][2]") << std::endl
				<< RM3("     T33:", -9, 4) << Ff(t, "T[3][3]") << SEP(" T12:", -9, 4) << Ff(t, "T[1][2]") << std::endl
				<< RM3("     T13:", -9, 4) << Ff(t, "T[1][3]") << SEP(" T23:", -9, 4) << Ff(t, "T[2][3]") << std::endl
				<< RM3("   L TENSOR") << std::endl
				<< RM3("     L11:", -9, 4) << Ff(t, "L[1][1]") << SEP(" L22:", -9, 4) << Ff(t, "L[2][2]") << std::endl
				<< RM3("     L33:", -9, 4) << Ff(t, "L[3][3]") << SEP(" L12:", -9, 4) << Ff(t, "L[1][2]") << std::endl
				<< RM3("     L13:", -9, 4) << Ff(t, "L[1][3]") << SEP(" L23:", -9, 4) << Ff(t, "L[2][3]") << std::endl
				<< RM3("   S TENSOR") << std::endl
				<< RM3("     S11:", -9, 4) << Ff(t, "S[1][1]") << SEP(" S12:", -9, 4) << Ff(t, "S[1][2]") << SEP(" S13:", -9, 4) << Ff(t, "S[1][3]") << std::endl
				<< RM3("     S21:", -9, 4) << Ff(t, "S[2][1]") << SEP(" S22:", -9, 4) << Ff(t, "S[2][2]") << SEP(" S23:", -9, 4) << Ff(t, "S[2][3]") << std::endl
				<< RM3("     S31:", -9, 4) << Ff(t, "S[3][1]") << SEP(" S32:", -9, 4) << Ff(t, "S[3][2]") << SEP(" S33:", -9, 4) << Ff(t, "S[3][3]") << std::endl;
	}
			
	pdbFile	<< RM3("") << std::endl;
}

// --------------------------------------------------------------------

void WriteRemark3CNS(std::ostream& pdbFile, Datablock& db)
{
	auto refine = db["refine"].front();
	auto ls_shell = db["refine_ls_shell"].front();
	auto hist = db["refine_hist"].front();
	auto reflns = db["reflns"].front();
	auto analyze = db["refine_analyze"].front();
	auto& ls_restr = db["refine_ls_restr"];
	auto ls_restr_ncs = db["refine_ls_restr_ncs"].front();
//	auto pdbx_xplor_file = db["pdbx_xplor_file"].front();
//	auto pdbx_refine = db["pdbx_refine"].front();
	
	pdbFile	<< RM3("") << std::endl
			<< RM3("REFINEMENT TARGET : ") << Fs(refine, "pdbx_stereochemistry_target_values") << std::endl
			<< RM3("") << std::endl
			<< RM3(" DATA USED IN REFINEMENT.") << std::endl
			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_high") << std::endl
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_low") << std::endl
			<< RM3("  DATA CUTOFF            (SIGMA(F)) : ", 6, 3)	<< Ff(refine, "pdbx_ls_sigma_F") << std::endl
			<< RM3("  DATA CUTOFF HIGH         (ABS(F)) : ", 6, 3)	<< Ff(refine, "pdbx_data_cutoff_high_absF") << std::endl
			<< RM3("  DATA CUTOFF LOW          (ABS(F)) : ", 7, 4)	<< Ff(refine, "pdbx_data_cutoff_low_absF") << std::endl
			<< RM3("  COMPLETENESS (WORKING+TEST)   (%) : ", 4, 1)	<< Ff(refine, "ls_percent_reflns_obs") << std::endl
			<< RM3("  NUMBER OF REFLECTIONS             : ", 12, 6)	<< Fi(refine, "ls_number_reflns_obs") << std::endl
	
			<< RM3("") << std::endl
			<< RM3(" FIT TO DATA USED IN REFINEMENT.") << std::endl
			<< RM3("  CROSS-VALIDATION METHOD          : ")	<< Fs(refine, "pdbx_ls_cross_valid_method") << std::endl
			<< RM3("  FREE R VALUE TEST SET SELECTION  : ")	<< Fs(refine, "pdbx_R_Free_selection_details") << std::endl
//			<< RM3("  R VALUE     (WORKING + TEST SET) : ", 7, 5)	<< Ff(refine, "ls_R_factor_obs") << std::endl
			<< RM3("  R VALUE            (WORKING SET) : ", 7, 3)	<< Ff(refine, "ls_R_factor_R_work") << std::endl
			<< RM3("  FREE R VALUE                     : ", 7, 3)	<< Ff(refine, "ls_R_factor_R_free") << std::endl
			<< RM3("  FREE R VALUE TEST SET SIZE   (%) : ", 7, 3)	<< Ff(refine, "ls_percent_reflns_R_free") << std::endl
			<< RM3("  FREE R VALUE TEST SET COUNT      : ", 12, 6)	<< Fi(refine, "ls_number_reflns_R_free") << std::endl
			<< RM3("  ESTIMATED ERROR OF FREE R VALUE  : ", 7, 3)	<< Ff(refine, "ls_R_factor_R_free_error") << std::endl
			
//			<< RM3("") << std::endl
//			<< RM3(" FIT/AGREEMENT OF MODEL WITH ALL DATA.") << std::endl
//			<< RM3("  R VALUE   (WORKING + TEST SET, NO CUTOFF) : ", 7, 3)	<< Ff(pdbx_refine, "R_factor_all_no_cutoff") << std::endl
//			<< RM3("  R VALUE          (WORKING SET, NO CUTOFF) : ", 7, 3)	<< Ff(pdbx_refine, "R_factor_obs_no_cutoff") << std::endl
//			<< RM3("  FREE R VALUE                  (NO CUTOFF) : ", 7, 3)	<< Ff(pdbx_refine, "free_R_factor_no_cutoff") << std::endl
//			<< RM3("  FREE R VALUE TEST SET SIZE (%, NO CUTOFF) : ", 7, 3)	<< Ff(pdbx_refine, "free_R_val_test_set_size_perc_no_cutoff") << std::endl
//			<< RM3("  FREE R VALUE TEST SET COUNT   (NO CUTOFF) : ", 12, 6)	<< Fi(pdbx_refine, "free_R_val_test_set_ct_no_cutoff") << std::endl
//			<< RM3("  TOTAL NUMBER OF REFLECTIONS   (NO CUTOFF) : ", 12, 6)	<< Fi(refine, "ls_number_reflns_all") << std::endl

			<< RM3("") << std::endl
			<< RM3(" FIT IN THE HIGHEST RESOLUTION BIN.") << std::endl
			<< RM3("  TOTAL NUMBER OF BINS USED           : ", 12, 6)	<< Fi(ls_shell, "pdbx_total_number_of_bins_used") << std::endl
			<< RM3("  BIN RESOLUTION RANGE HIGH       (A) : ", 5, 2)	<< Ff(ls_shell, "d_res_high") << std::endl
			<< RM3("  BIN RESOLUTION RANGE LOW        (A) : ", 5, 2)	<< Ff(ls_shell, "d_res_low") << std::endl
			<< RM3("  BIN COMPLETENESS (WORKING+TEST) (%) : ", 6, 2)	<< Ff(ls_shell, "percent_reflns_obs") << std::endl
			<< RM3("  REFLECTIONS IN BIN    (WORKING SET) : ", 12, 6)	<< Fi(ls_shell, "number_reflns_R_work") << std::endl
			<< RM3("  BIN R VALUE           (WORKING SET) : ", 8, 4)	<< Ff(ls_shell, "R_factor_R_work") << std::endl
			<< RM3("  BIN FREE R VALUE                    : ", 8, 4)	<< Ff(ls_shell, "R_factor_R_free") << std::endl
			<< RM3("  BIN FREE R VALUE TEST SET SIZE  (%) : ", 6, 2)	<< Ff(ls_shell, "percent_reflns_R_free") << std::endl
			<< RM3("  BIN FREE R VALUE TEST SET COUNT     : ", 12, 7)	<< Fi(ls_shell, "number_reflns_R_free") << std::endl
			<< RM3("  ESTIMATED ERROR OF BIN FREE R VALUE : ", 7, 3)	<< Ff(ls_shell, "R_factor_R_free_error") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.") << std::endl
			<< RM3("  PROTEIN ATOMS            : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_protein") << std::endl
			<< RM3("  NUCLEIC ACID ATOMS       : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_nucleic_acid") << std::endl
			<< RM3("  HETEROGEN ATOMS          : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_ligand") << std::endl
			<< RM3("  SOLVENT ATOMS            : ", 12, 6)	<< Fi(hist, "number_atoms_solvent") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" B VALUES.") << std::endl
			<< RM3("  B VALUE TYPE                      : ")		<< Fs(refine, "pdbx_TLS_residual_ADP_flag") << std::endl
			<< RM3("  FROM WILSON PLOT           (A**2) : ", 7, 2)	<< Ff(reflns, "B_iso_Wilson_estimate") << std::endl
			<< RM3("  MEAN B VALUE      (OVERALL, A**2) : ", 7, 2)	<< Ff(refine, "B_iso_mean") << std::endl
			
			<< RM3("  OVERALL ANISOTROPIC B VALUE.") << std::endl
			<< RM3("   B11 (A**2) : ", -8, 5)	<< Ff(refine, "aniso_B[1][1]") << std::endl
			<< RM3("   B22 (A**2) : ", -8, 5)	<< Ff(refine, "aniso_B[2][2]") << std::endl
			<< RM3("   B33 (A**2) : ", -8, 5)	<< Ff(refine, "aniso_B[3][3]") << std::endl
			<< RM3("   B12 (A**2) : ", -8, 5)	<< Ff(refine, "aniso_B[1][2]") << std::endl
			<< RM3("   B13 (A**2) : ", -8, 5)	<< Ff(refine, "aniso_B[1][3]") << std::endl
			<< RM3("   B23 (A**2) : ", -8, 5)	<< Ff(refine, "aniso_B[2][3]") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" ESTIMATED COORDINATE ERROR.") << std::endl
			<< RM3("  ESD FROM LUZZATI PLOT        (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_coordinate_error_obs") << std::endl
			<< RM3("  ESD FROM SIGMAA              (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_sigma_a_obs") << std::endl
			<< RM3("  LOW RESOLUTION CUTOFF        (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_d_res_low_obs") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" CROSS-VALIDATED ESTIMATED COORDINATE ERROR.") << std::endl
			<< RM3("  ESD FROM C-V LUZZATI PLOT    (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_coordinate_error_free") << std::endl
			<< RM3("  ESD FROM C-V SIGMAA          (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_sigma_a_free") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" RMS DEVIATIONS FROM IDEAL VALUES.") << std::endl
			<< RM3("  BOND LENGTHS                 (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "c_bond_d", "dev_ideal") << std::endl
			<< RM3("  BOND ANGLES            (DEGREES) : ", 7, 2)	<< Ff(ls_restr, cif::Key("type") == "c_angle_deg", "dev_ideal") << std::endl
			<< RM3("  DIHEDRAL ANGLES        (DEGREES) : ", 7, 2)	<< Ff(ls_restr, cif::Key("type") == "c_dihedral_angle_d", "dev_ideal") << std::endl
			<< RM3("  IMPROPER ANGLES        (DEGREES) : ", 7, 2)	<< Ff(ls_restr, cif::Key("type") == "c_improper_angle_d", "dev_ideal") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" ISOTROPIC THERMAL MODEL : ") << Fs(refine, "pdbx_isotropic_thermal_model") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" ISOTROPIC THERMAL FACTOR RESTRAINTS.    RMS    SIGMA") << std::endl
			<< RM3("  MAIN-CHAIN BOND              (A**2) : ", 7, 3) << Ff(ls_restr, cif::Key("type") == "c_mcbond_it", "dev_ideal") << SEP("; ", 7, 3) 
																	 << Ff(ls_restr, cif::Key("type") == "c_mcbond_it", "dev_ideal_target") << std::endl
			<< RM3("  MAIN-CHAIN ANGLE             (A**2) : ", 7, 3) << Ff(ls_restr, cif::Key("type") == "c_mcangle_it", "dev_ideal") << SEP("; ", 7, 3) 
																	 << Ff(ls_restr, cif::Key("type") == "c_mcangle_it", "dev_ideal_target") << std::endl
			<< RM3("  SIDE-CHAIN BOND              (A**2) : ", 7, 3) << Ff(ls_restr, cif::Key("type") == "c_scbond_it", "dev_ideal") << SEP("; ", 7, 3) 
																	 << Ff(ls_restr, cif::Key("type") == "c_scbond_it", "dev_ideal_target") << std::endl
			<< RM3("  SIDE-CHAIN ANGLE             (A**2) : ", 7, 3) << Ff(ls_restr, cif::Key("type") == "c_scangle_it", "dev_ideal") << SEP("; ", 7, 3) 
																	 << Ff(ls_restr, cif::Key("type") == "c_scangle_it", "dev_ideal_target") << std::endl

			<< RM3("") << std::endl
			<< RM3(" BULK SOLVENT MODELING.") << std::endl
			<< RM3("  METHOD USED        : ")				<< Fs(refine, "solvent_model_details") << std::endl
			<< RM3("  KSOL               : ", 5, 2)         << Ff(refine, "solvent_model_param_ksol") << std::endl
			<< RM3("  BSOL               : ", 5, 2)         << Ff(refine, "solvent_model_param_bsol") << std::endl

			<< RM3("") << std::endl
			<< RM3(" NCS MODEL : ")	<< Fs(ls_restr_ncs, "ncs_model_details") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" NCS RESTRAINTS.                         RMS   SIGMA/WEIGHT") << std::endl
			
			// TODO: using only group 1 here, should this be fixed???
			<< RM3("  GROUP  1  POSITIONAL            (A) : ", 4, 2) << Ff(ls_restr_ncs, "rms_dev_position") << SEP("; ", 6, 2) 
																		<< Ff(ls_restr_ncs, "weight_position") << SEP("; ", 6, 2) << std::endl
			<< RM3("  GROUP  1  B-FACTOR           (A**2) : ", 4, 2) << Ff(ls_restr_ncs, "rms_dev_B_iso") << SEP("; ", 6, 2) 
																		<< Ff(ls_restr_ncs, "weight_B_iso") << SEP("; ", 6, 2) << std::endl
			
			// TODO: using only files from serial_no 1 here
//			<< RM3("") << std::endl
//			<< RM3(" PARAMETER FILE   1  : ") << Fs(pdbx_xplor_file, "param_file") << std::endl
//			<< RM3(" TOPOLOGY FILE   1   : ") << Fs(pdbx_xplor_file, "topol_file") << std::endl
			
			<< RM3("") << std::endl;
}



// --------------------------------------------------------------------

void WriteRemark3Refmac(std::ostream& pdbFile, Datablock& db)
{
	auto refine = db["refine"].front();
	auto ls_shell = db["refine_ls_shell"].front();
	auto hist = db["refine_hist"].front();
	auto reflns = db["reflns"].front();
//	auto analyze = db["refine_analyze"].front();
	auto& ls_restr = db["refine_ls_restr"];
//	auto pdbx_xplor_file = db["pdbx_xplor_file"].front();
	
	auto c = [](const char* t) -> cif::Condition { return cif::Key("type") == t; };
	
	pdbFile	<< RM3("") << std::endl
			<< RM3("REFINEMENT TARGET : ") << Fs(refine, "pdbx_stereochemistry_target_values") << std::endl
			<< RM3("") << std::endl
			<< RM3(" DATA USED IN REFINEMENT.") << std::endl
			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_high") << std::endl
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_low") << std::endl
			<< RM3("  DATA CUTOFF            (SIGMA(F)) : ", 6, 3)	<< Ff(refine, "pdbx_ls_sigma_F") << std::endl
			<< RM3("  COMPLETENESS FOR RANGE        (%) : ", 5, 2)	<< Ff(refine, "ls_percent_reflns_obs") << std::endl
			<< RM3("  NUMBER OF REFLECTIONS             : ", 12, 6)	<< Fi(refine, "ls_number_reflns_obs") << std::endl
	
			<< RM3("") << std::endl
			<< RM3(" FIT TO DATA USED IN REFINEMENT.") << std::endl
			<< RM3("  CROSS-VALIDATION METHOD          : ")	<< Fs(refine, "pdbx_ls_cross_valid_method") << std::endl
			<< RM3("  FREE R VALUE TEST SET SELECTION  : ")	<< Fs(refine, "pdbx_R_Free_selection_details") << std::endl
			<< RM3("  R VALUE     (WORKING + TEST SET) : ", 7, 5)	<< Ff(refine, "ls_R_factor_obs") << std::endl
			<< RM3("  R VALUE            (WORKING SET) : ", 7, 5)	<< Ff(refine, "ls_R_factor_R_work") << std::endl
			<< RM3("  FREE R VALUE                     : ", 7, 5)	<< Ff(refine, "ls_R_factor_R_free") << std::endl
			<< RM3("  FREE R VALUE TEST SET SIZE   (%) : ", 7, 1)	<< Ff(refine, "ls_percent_reflns_R_free") << std::endl
			<< RM3("  FREE R VALUE TEST SET COUNT      : ", 12, 6)	<< Fi(refine, "ls_number_reflns_R_free") << std::endl
			<< RM3("  ESTIMATED ERROR OF FREE R VALUE  : ", 7, 3)	<< Ff(refine, "ls_R_factor_R_free_error") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" FIT IN THE HIGHEST RESOLUTION BIN.") << std::endl
			<< RM3("  TOTAL NUMBER OF BINS USED           : ")			<< Fi(ls_shell, "pdbx_total_number_of_bins_used") << std::endl
			<< RM3("  BIN RESOLUTION RANGE HIGH       (A) : ", 5, 3)	<< Ff(ls_shell, "d_res_high") << std::endl
			<< RM3("  BIN RESOLUTION RANGE LOW        (A) : ", 5, 3)	<< Ff(ls_shell, "d_res_low") << std::endl
			<< RM3("  REFLECTION IN BIN     (WORKING SET) : ")			<< Fi(ls_shell, "number_reflns_R_work") << std::endl
			<< RM3("  BIN COMPLETENESS (WORKING+TEST) (%) : ", 5, 2)	<< Ff(ls_shell, "percent_reflns_obs") << std::endl
			<< RM3("  BIN R VALUE           (WORKING SET) : ", 7, 3)	<< Ff(ls_shell, "R_factor_R_work") << std::endl
			<< RM3("  BIN FREE R VALUE SET COUNT          : ")			<< Fi(ls_shell, "number_reflns_R_free") << std::endl
			<< RM3("  BIN FREE R VALUE                    : ", 7, 3)	<< Ff(ls_shell, "R_factor_R_free") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.") << std::endl
			<< RM3("  PROTEIN ATOMS            : ")	<< Fi(hist, "pdbx_number_atoms_protein") << std::endl
			<< RM3("  NUCLEIC ACID ATOMS       : ")	<< Fi(hist, "pdbx_number_atoms_nucleic_acid") << std::endl
			<< RM3("  HETEROGEN ATOMS          : ")	<< Fi(hist, "pdbx_number_atoms_ligand") << std::endl
			<< RM3("  SOLVENT ATOMS            : ")	<< Fi(hist, "number_atoms_solvent") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" B VALUES.") << std::endl
			<< RM3("  B VALUE TYPE                      : ")	<< Fs(refine, "pdbx_TLS_residual_ADP_flag") << std::endl
			<< RM3("  FROM WILSON PLOT           (A**2) : ", 8, 3)	<< Ff(reflns, "B_iso_Wilson_estimate") << std::endl
			<< RM3("  MEAN B VALUE      (OVERALL, A**2) : ", 8, 3)	<< Ff(refine, "B_iso_mean") << std::endl
			
			<< RM3("  OVERALL ANISOTROPIC B VALUE.") << std::endl
			<< RM3("   B11 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][1]") << std::endl
			<< RM3("   B22 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[2][2]") << std::endl
			<< RM3("   B33 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[3][3]") << std::endl
			<< RM3("   B12 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][2]") << std::endl
			<< RM3("   B13 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][3]") << std::endl
			<< RM3("   B23 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[2][3]") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" ESTIMATED OVERALL COORDINATE ERROR.") << std::endl
			<< RM3("  ESU BASED ON R VALUE                            (A): ", 6, 3) << Ff(refine, "pdbx_overall_ESU_R") << std::endl
			<< RM3("  ESU BASED ON FREE R VALUE                       (A): ", 6, 3) << Ff(refine, "pdbx_overall_ESU_R_Free") << std::endl
			<< RM3("  ESU BASED ON MAXIMUM LIKELIHOOD                 (A): ", 6, 3) << Ff(refine, "overall_SU_ML") << std::endl
			<< RM3("  ESU FOR B VALUES BASED ON MAXIMUM LIKELIHOOD (A**2): ", 6, 3) << Ff(refine, "overall_SU_B") << std::endl

			<< RM3("") << std::endl
			<< RM3(" CORRELATION COEFFICIENTS.") << std::endl
			<< RM3("  CORRELATION COEFFICIENT FO-FC      : ", 6, 3) << Ff(refine, "correlation_coeff_Fo_to_Fc") << std::endl
			<< RM3("  CORRELATION COEFFICIENT FO-FC FREE : ", 6, 3) << Ff(refine, "correlation_coeff_Fo_to_Fc_free") << std::endl

			<< RM3("") << std::endl
			<< RM3(" RMS DEVIATIONS FROM IDEAL VALUES        COUNT    RMS    WEIGHT") << std::endl
			<< RM3("  BOND LENGTHS REFINED ATOMS        (A): ", -5) << Fi(ls_restr, c("r_bond_refined_d"             ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_bond_refined_d"             ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_bond_refined_d"             ), "dev_ideal_target") << std::endl
			<< RM3("  BOND LENGTHS OTHERS               (A): ", -5) << Fi(ls_restr, c("r_bond_other_d"               ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_bond_other_d"               ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_bond_other_d"               ), "dev_ideal_target") << std::endl
			<< RM3("  BOND ANGLES REFINED ATOMS   (DEGREES): ", -5) << Fi(ls_restr, c("r_angle_refined_deg"          ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_angle_refined_deg"          ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_angle_refined_deg"          ), "dev_ideal_target") << std::endl
			<< RM3("  BOND ANGLES OTHERS          (DEGREES): ", -5) << Fi(ls_restr, c("r_angle_other_deg"            ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_angle_other_deg"            ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_angle_other_deg"            ), "dev_ideal_target") << std::endl
			<< RM3("  TORSION ANGLES, PERIOD 1    (DEGREES): ", -5) << Fi(ls_restr, c("r_dihedral_angle_1_deg"       ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_dihedral_angle_1_deg"       ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_dihedral_angle_1_deg"       ), "dev_ideal_target") << std::endl
			<< RM3("  TORSION ANGLES, PERIOD 2    (DEGREES): ", -5) << Fi(ls_restr, c("r_dihedral_angle_2_deg"       ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_dihedral_angle_2_deg"       ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_dihedral_angle_2_deg"       ), "dev_ideal_target") << std::endl
			<< RM3("  TORSION ANGLES, PERIOD 3    (DEGREES): ", -5) << Fi(ls_restr, c("r_dihedral_angle_3_deg"       ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_dihedral_angle_3_deg"       ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_dihedral_angle_3_deg"       ), "dev_ideal_target") << std::endl
			<< RM3("  TORSION ANGLES, PERIOD 4    (DEGREES): ", -5) << Fi(ls_restr, c("r_dihedral_angle_4_deg"       ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_dihedral_angle_4_deg"       ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_dihedral_angle_4_deg"       ), "dev_ideal_target") << std::endl
			<< RM3("  CHIRAL-CENTER RESTRAINTS       (A**3): ", -5) << Fi(ls_restr, c("r_chiral_restr"               ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_chiral_restr"               ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_chiral_restr"               ), "dev_ideal_target") << std::endl
			<< RM3("  GENERAL PLANES REFINED ATOMS      (A): ", -5) << Fi(ls_restr, c("r_gen_planes_refined"         ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_gen_planes_refined"         ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_gen_planes_refined"         ), "dev_ideal_target") << std::endl
			<< RM3("  GENERAL PLANES OTHERS             (A): ", -5) << Fi(ls_restr, c("r_gen_planes_other"           ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_gen_planes_other"           ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_gen_planes_other"           ), "dev_ideal_target") << std::endl
			<< RM3("  NON-BONDED CONTACTS REFINED ATOMS (A): ", -5) << Fi(ls_restr, c("r_nbd_refined"                ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_nbd_refined"                ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_nbd_refined"                ), "dev_ideal_target") << std::endl
			<< RM3("  NON-BONDED CONTACTS OTHERS        (A): ", -5) << Fi(ls_restr, c("r_nbd_other"                  ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_nbd_other"                  ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_nbd_other"                  ), "dev_ideal_target") << std::endl
			<< RM3("  NON-BONDED TORSION REFINED ATOMS  (A): ", -5) << Fi(ls_restr, c("r_nbtor_refined"              ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_nbtor_refined"              ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_nbtor_refined"              ), "dev_ideal_target") << std::endl
			<< RM3("  NON-BONDED TORSION OTHERS         (A): ", -5) << Fi(ls_restr, c("r_nbtor_other"                ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_nbtor_other"                ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_nbtor_other"                ), "dev_ideal_target") << std::endl
			<< RM3("  H-BOND (X...Y) REFINED ATOMS      (A): ", -5) << Fi(ls_restr, c("r_xyhbond_nbd_refined"        ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_xyhbond_nbd_refined"        ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_xyhbond_nbd_refined"        ), "dev_ideal_target") << std::endl
			<< RM3("  H-BOND (X...Y) OTHERS             (A): ", -5) << Fi(ls_restr, c("r_xyhbond_nbd_other"          ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_xyhbond_nbd_other"          ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_xyhbond_nbd_other"          ), "dev_ideal_target") << std::endl
			<< RM3("  POTENTIAL METAL-ION REFINED ATOMS (A): ", -5) << Fi(ls_restr, c("r_metal_ion_refined"          ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_metal_ion_refined"          ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_metal_ion_refined"          ), "dev_ideal_target") << std::endl
			<< RM3("  POTENTIAL METAL-ION OTHERS        (A): ", -5) << Fi(ls_restr, c("r_metal_ion_other"            ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_metal_ion_other"            ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_metal_ion_other"            ), "dev_ideal_target") << std::endl
			<< RM3("  SYMMETRY VDW REFINED ATOMS        (A): ", -5) << Fi(ls_restr, c("r_symmetry_vdw_refined"       ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_symmetry_vdw_refined"       ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_symmetry_vdw_refined"       ), "dev_ideal_target") << std::endl
			<< RM3("  SYMMETRY VDW OTHERS               (A): ", -5) << Fi(ls_restr, c("r_symmetry_vdw_other"         ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_symmetry_vdw_other"         ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_symmetry_vdw_other"         ), "dev_ideal_target") << std::endl
			<< RM3("  SYMMETRY H-BOND REFINED ATOMS     (A): ", -5) << Fi(ls_restr, c("r_symmetry_hbond_refined"     ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_symmetry_hbond_refined"     ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_symmetry_hbond_refined"     ), "dev_ideal_target") << std::endl
			<< RM3("  SYMMETRY H-BOND OTHERS            (A): ", -5) << Fi(ls_restr, c("r_symmetry_hbond_other"       ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_symmetry_hbond_other"       ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_symmetry_hbond_other"       ), "dev_ideal_target") << std::endl
			<< RM3("  SYMMETRY METAL-ION REFINED ATOMS  (A): ", -5) << Fi(ls_restr, c("r_symmetry_metal_ion_refined" ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_symmetry_metal_ion_refined" ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_symmetry_metal_ion_refined" ), "dev_ideal_target") << std::endl
			<< RM3("  SYMMETRY METAL-ION OTHERS         (A): ", -5) << Fi(ls_restr, c("r_symmetry_metal_ion_other"   ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_symmetry_metal_ion_other"   ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_symmetry_metal_ion_other"   ), "dev_ideal_target") << std::endl

			<< RM3("") << std::endl
			<< RM3(" ISOTROPIC THERMAL FACTOR RESTRAINTS.     COUNT   RMS    WEIGHT") << std::endl
			<< RM3("  MAIN-CHAIN BOND REFINED ATOMS  (A**2): ", -5) << Fi(ls_restr, c("r_mcbond_it"   ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_mcbond_it"   ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_mcbond_it"   ), "dev_ideal_target") << std::endl
			<< RM3("  MAIN-CHAIN BOND OTHER ATOMS    (A**2): ", -5) << Fi(ls_restr, c("r_mcbond_other"   ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_mcbond_other"   ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_mcbond_other"   ), "dev_ideal_target") << std::endl
			<< RM3("  MAIN-CHAIN ANGLE REFINED ATOMS (A**2): ", -5) << Fi(ls_restr, c("r_mcangle_it"   ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_mcangle_it"   ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_mcangle_it"   ), "dev_ideal_target") << std::endl
			<< RM3("  MAIN-CHAIN ANGLE OTHER ATOMS   (A**2): ", -5) << Fi(ls_restr, c("r_mcangle_other"   ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_mcangle_other"   ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_mcangle_other"   ), "dev_ideal_target") << std::endl
			<< RM3("  SIDE-CHAIN BOND REFINED ATOMS  (A**2): ", -5) << Fi(ls_restr, c("r_scbond_it"   ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_scbond_it"   ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_scbond_it"   ), "dev_ideal_target") << std::endl
			<< RM3("  SIDE-CHAIN BOND OTHER ATOMS    (A**2): ", -5) << Fi(ls_restr, c("r_scbond_other"   ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_scbond_other"   ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_scbond_other"   ), "dev_ideal_target") << std::endl
			<< RM3("  SIDE-CHAIN ANGLE REFINED ATOMS (A**2): ", -5) << Fi(ls_restr, c("r_scangle_it"   ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_scangle_it"   ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_scangle_it"   ), "dev_ideal_target") << std::endl
			<< RM3("  SIDE-CHAIN ANGLE OTHER ATOMS   (A**2): ", -5) << Fi(ls_restr, c("r_scangle_other"   ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_scangle_other"   ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_scangle_other"   ), "dev_ideal_target") << std::endl
			<< RM3("  LONG RANGE B REFINED ATOMS     (A**2): ", -5) << Fi(ls_restr, c("r_long_range_B_refined"   ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_long_range_B_refined"   ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_long_range_B_refined"   ), "dev_ideal_target") << std::endl
			<< RM3("  LONG RANGE B OTHER ATOMS       (A**2): ", -5) << Fi(ls_restr, c("r_long_range_B_other"   ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_long_range_B_other"   ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_long_range_B_other"   ), "dev_ideal_target") << std::endl

			<< RM3("") << std::endl
			<< RM3(" ANISOTROPIC THERMAL FACTOR RESTRAINTS.   COUNT   RMS    WEIGHT") << std::endl
			<< RM3("  RIGID-BOND RESTRAINTS          (A**2): ", -5) << Fi(ls_restr, c("r_rigid_bond_restr"   ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_rigid_bond_restr"   ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_rigid_bond_restr"   ), "dev_ideal_target") << std::endl
			<< RM3("  SPHERICITY; FREE ATOMS         (A**2): ", -5) << Fi(ls_restr, c("r_sphericity_free"   ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_sphericity_free"   ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_sphericity_free"   ), "dev_ideal_target") << std::endl
			<< RM3("  SPHERICITY; BONDED ATOMS       (A**2): ", -5) << Fi(ls_restr, c("r_sphericity_bonded"   ), "number") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_sphericity_bonded"   ), "dev_ideal") << SEP(" ;", -6, 3)
																	<< Ff(ls_restr, c("r_sphericity_bonded"   ), "dev_ideal_target") << std::endl

			<< RM3("") << std::endl
			<< RM3(" NCS RESTRAINTS STATISTICS") << std::endl;
	
		auto& ncs_dom = db["struct_ncs_dom"];
		if (ncs_dom.empty())
			pdbFile << RM3("  NUMBER OF DIFFERENT NCS GROUPS : NULL") << std::endl;
		else
		{
			std::set<std::string> ncs_groups;
			for (auto i: ncs_dom)
				ncs_groups.insert(i["pdbx_ens_id"].as<std::string>());

			pdbFile << RM3("  NUMBER OF DIFFERENT NCS GROUPS : ") << ncs_groups.size() << std::endl;
			
			for (auto ens_id: ncs_groups)
			{
				auto lim = db["struct_ncs_dom_lim"].find(cif::Key("pdbx_ens_id") == ens_id);
				
				std::set<std::string> chains;
				std::set<int> component_ids;
				
				for (auto l: lim)
				{
					chains.insert(l["beg_auth_asym_id"].as<std::string>());
					component_ids.insert(l["pdbx_component_id"].as<int>());
				}
				
				pdbFile << RM3("") << std::endl
						<< RM3(" NCS GROUP NUMBER               : ") << ens_id << std::endl
						<< RM3("    CHAIN NAMES                    : ") << ba::join(chains, " ") << std::endl
						<< RM3("    NUMBER OF COMPONENTS NCS GROUP : ") << component_ids.size() << std::endl
						<< RM3("      COMPONENT C  SSSEQI  TO  C   SSSEQI   CODE") << std::endl;
				
				for (auto l: lim)
				{
					pdbFile << RM3("         ", -2)		<< Fi(l, "pdbx_component_id")
							<< SEP(" ", -5)				<< Fs(l, "beg_auth_asym_id")
							<< SEP("  ", -5)			<< Fi(l, "beg_auth_seq_id")
							<< SEP("   ", -5)			<< Fs(l, "end_auth_asym_id")
							<< SEP("   ", -5)			<< Fi(l, "end_auth_seq_id")
							<< SEP("  ", -5)			<< Fs(l, "pdbx_refine_code")
							<< std::endl;
				}
				
				pdbFile << RM3("                  GROUP CHAIN        COUNT   RMS     WEIGHT") << std::endl;
				for (auto l: db["refine_ls_restr_ncs"].find(cif::Key("pdbx_ens_id") == ens_id))
				{
					std::string type = l["pdbx_type"].as<std::string>();
					ba::to_upper(type);
					
					std::string unit;
					if (ba::ends_with(type, "POSITIONAL"))
						unit = "    (A): ";
					else if (ba::ends_with(type, "THERMAL"))
						unit = " (A**2): ";
					else
						unit = "       : ";
					
					pdbFile << RM3("  ", 18)			<< type
							<< SEP("", -2)				<< Fi(l, "pdbx_ens_id")
							<< SEP("    ", 1)			<< Fs(l, "pdbx_auth_asym_id")
							<< SEP(unit.c_str(), -6)	<< Fi(l, "pdbx_number")
							<< SEP(" ;", -6, 3)			<< Ff(l, "rms_dev_position")
							<< SEP(" ;", -6, 3)			<< Ff(l, "weight_position")
							<< std::endl;
				}
			}
		}
		
		// TODO: add twin information

//	{ R"(TWIN DETAILS)", "", {} },
//	{ R"(NUMBER OF TWIN DOMAINS)", "", {} },

	auto& tls = db["pdbx_refine_tls"];	

	pdbFile << RM3("") << std::endl
			<< RM3(" TLS DETAILS") << std::endl
			<< RM3("  NUMBER OF TLS GROUPS  : ") << (tls.size() ? std::to_string(tls.size()) : "NULL") << std::endl;
	
	for (auto t: tls)
	{
		std::string id = t["id"].as<std::string>();
		auto g = db["pdbx_refine_tls_group"].find(cif::Key("refine_tls_id") == id);
		
		pdbFile << RM3("") << std::endl
				<< RM3("  TLS GROUP : ") << id << std::endl
				<< RM3("   NUMBER OF COMPONENTS GROUP : ") << g.size() << std::endl
				<< RM3("   COMPONENTS        C SSSEQI   TO  C SSSEQI") << std::endl;
		
		for (auto gi: g)
		{
			pdbFile << RM3("   RESIDUE RANGE :   ") << Fs(gi, "beg_auth_asym_id")
					<< SEP("", -6)					<< Fi(gi, "beg_auth_seq_id")
					<< SEP("", -9)					<< Fs(gi, "end_auth_asym_id")
					<< SEP("", -6)					<< Fi(gi, "end_auth_seq_id")
					<< std::endl;
		}
		
		pdbFile << RM3("   ORIGIN FOR THE GROUP (A):", -9, 4) << Ff(t, "origin_x")
				<< SEP("", -9, 4) << Ff(t, "origin_y")
				<< SEP("", -9, 4) << Ff(t, "origin_z") << std::endl
				<< RM3("   T TENSOR") << std::endl
				<< RM3("     T11:", -9, 4) << Ff(t, "T[1][1]") << SEP(" T22:", -9, 4) << Ff(t, "T[2][2]") << std::endl
				<< RM3("     T33:", -9, 4) << Ff(t, "T[3][3]") << SEP(" T12:", -9, 4) << Ff(t, "T[1][2]") << std::endl
				<< RM3("     T13:", -9, 4) << Ff(t, "T[1][3]") << SEP(" T23:", -9, 4) << Ff(t, "T[2][3]") << std::endl
				<< RM3("   L TENSOR") << std::endl
				<< RM3("     L11:", -9, 4) << Ff(t, "L[1][1]") << SEP(" L22:", -9, 4) << Ff(t, "L[2][2]") << std::endl
				<< RM3("     L33:", -9, 4) << Ff(t, "L[3][3]") << SEP(" L12:", -9, 4) << Ff(t, "L[1][2]") << std::endl
				<< RM3("     L13:", -9, 4) << Ff(t, "L[1][3]") << SEP(" L23:", -9, 4) << Ff(t, "L[2][3]") << std::endl
				<< RM3("   S TENSOR") << std::endl
				<< RM3("     S11:", -9, 4) << Ff(t, "S[1][1]") << SEP(" S12:", -9, 4) << Ff(t, "S[1][2]") << SEP(" S13:", -9, 4) << Ff(t, "S[1][3]") << std::endl
				<< RM3("     S21:", -9, 4) << Ff(t, "S[2][1]") << SEP(" S22:", -9, 4) << Ff(t, "S[2][2]") << SEP(" S23:", -9, 4) << Ff(t, "S[2][3]") << std::endl
				<< RM3("     S31:", -9, 4) << Ff(t, "S[3][1]") << SEP(" S32:", -9, 4) << Ff(t, "S[3][2]") << SEP(" S33:", -9, 4) << Ff(t, "S[3][3]") << std::endl;
	}

	pdbFile << RM3("") << std::endl
			<< RM3(" BULK SOLVENT MODELLING.") << std::endl
			<< RM3("  METHOD USED : ") << Fs(refine, "solvent_model_details") << std::endl
			<< RM3("  PARAMETERS FOR MASK CALCULATION") << std::endl
			<< RM3("  VDW PROBE RADIUS   : ", 5, 2) << Ff(refine, "pdbx_solvent_vdw_probe_radii") << std::endl
			<< RM3("  ION PROBE RADIUS   : ", 5, 2) << Ff(refine, "pdbx_solvent_ion_probe_radii") << std::endl
			<< RM3("  SHRINKAGE RADIUS   : ", 5, 2) << Ff(refine, "pdbx_solvent_shrinkage_radii") << std::endl
			
			<< RM3("") << std::endl;
}

void WriteRemark3Shelxl(std::ostream& pdbFile, Datablock& db)
{
	auto refine = db["refine"].front();
//	auto ls_shell = db["refine_ls_shell"].front();
	auto refine_hist = db["refine_hist"].front();
//	auto reflns = db["reflns"].front();
	auto refine_analyze = db["refine_analyze"].front();
	auto& ls_restr = db["refine_ls_restr"];
//	auto pdbx_xplor_file = db["pdbx_xplor_file"].front();
	auto pdbx_refine = db["pdbx_refine"].front();
	
	auto c = [](const char* t) -> cif::Condition { return cif::Key("type") == t; };
	
	pdbFile	<< RM3("") << std::endl
			<< RM3(" DATA USED IN REFINEMENT.") << std::endl
			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_high") << std::endl
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_low") << std::endl
			<< RM3("  DATA CUTOFF            (SIGMA(F)) : ", 6, 3)	<< Ff(refine, "pdbx_ls_sigma_F") << std::endl
			<< RM3("  COMPLETENESS FOR RANGE        (%) : ", 5, 2)	<< Ff(refine, "ls_percent_reflns_obs") << std::endl
			<< RM3("  CROSS-VALIDATION METHOD           : ")	<< Fs(refine, "pdbx_ls_cross_valid_method") << std::endl
			<< RM3("  FREE R VALUE TEST SET SELECTION   : ")	<< Fs(refine, "pdbx_R_Free_selection_details") << std::endl
	
			<< RM3("") << std::endl
			<< RM3(" FIT TO DATA USED IN REFINEMENT (NO CUTOFF).") << std::endl
			<< RM3("  R VALUE   (WORKING + TEST SET, NO CUTOFF) : ", 7, 3)	<< Ff(pdbx_refine, "R_factor_all_no_cutoff") << std::endl
			<< RM3("  R VALUE          (WORKING SET, NO CUTOFF) : ", 7, 3)	<< Ff(pdbx_refine, "R_factor_obs_no_cutoff") << std::endl
			<< RM3("  FREE R VALUE                  (NO CUTOFF) : ", 7, 3)	<< Ff(pdbx_refine, "free_R_factor_no_cutoff") << std::endl
			<< RM3("  FREE R VALUE TEST SET SIZE (%, NO CUTOFF) : ", 7, 3)	<< Ff(pdbx_refine, "free_R_val_test_set_size_perc_no_cutoff") << std::endl
			<< RM3("  FREE R VALUE TEST SET COUNT   (NO CUTOFF) : ", 12, 6)	<< Fi(pdbx_refine, "free_R_val_test_set_ct_no_cutoff") << std::endl
			<< RM3("  TOTAL NUMBER OF REFLECTIONS   (NO CUTOFF) : ", 12, 6)	<< Fi(refine, "ls_number_reflns_all") << std::endl

			<< RM3("") << std::endl
			<< RM3(" FIT/AGREEMENT OF MODEL FOR DATA WITH F>4SIG(F).") << std::endl
			<< RM3("  R VALUE   (WORKING + TEST SET, F>4SIG(F)) : ", 7, 3)	<< Ff(pdbx_refine, "R_factor_all_4sig_cutoff") << std::endl
			<< RM3("  R VALUE          (WORKING SET, F>4SIG(F)) : ", 7, 3)	<< Ff(pdbx_refine, "R_factor_obs_4sig_cutoff") << std::endl
			<< RM3("  FREE R VALUE                  (F>4SIG(F)) : ", 7, 3)	<< Ff(pdbx_refine, "free_R_factor_4sig_cutoff") << std::endl
			<< RM3("  FREE R VALUE TEST SET SIZE (%, F>4SIG(F)) : ", 7, 3)	<< Ff(pdbx_refine, "free_R_val_test_set_size_perc_4sig_cutoff") << std::endl
			<< RM3("  FREE R VALUE TEST SET COUNT   (F>4SIG(F)) : ")		<< Fi(pdbx_refine, "free_R_val_test_set_ct_4sig_cutoff") << std::endl
			<< RM3("  TOTAL NUMBER OF REFLECTIONS   (F>4SIG(F)) : ")		<< Fi(pdbx_refine, "number_reflns_obs_4sig_cutoff") << std::endl

			<< RM3("") << std::endl
			<< RM3(" NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.") << std::endl
			<< RM3("  PROTEIN ATOMS      : ")								<< Fi(refine_hist, "pdbx_number_atoms_protein") << std::endl
			<< RM3("  NUCLEIC ACID ATOMS : ")								<< Fi(refine_hist, "pdbx_number_atoms_nucleic_acid") << std::endl
			<< RM3("  HETEROGEN ATOMS    : ")								<< Fi(refine_hist, "pdbx_number_atoms_ligand") << std::endl
			<< RM3("  SOLVENT ATOMS      : ")								<< Fi(refine_hist, "number_atoms_solvent") << std::endl

			<< RM3("") << std::endl
			<< RM3(" MODEL REFINEMENT.")  << std::endl
			<< RM3("  OCCUPANCY SUM OF NON-HYDROGEN ATOMS      : ", 7, 3)	<< Ff(refine_analyze, "occupancy_sum_non_hydrogen") << std::endl
			<< RM3("  OCCUPANCY SUM OF HYDROGEN ATOMS          : ", 7, 3)	<< Ff(refine_analyze, "occupancy_sum_hydrogen") << std::endl
			<< RM3("  NUMBER OF DISCRETELY DISORDERED RESIDUES : ")			<< Fi(refine_analyze, "number_disordered_residues") << std::endl
			<< RM3("  NUMBER OF LEAST-SQUARES PARAMETERS       : ")			<< Fi(refine, "ls_number_parameters") << std::endl
			<< RM3("  NUMBER OF RESTRAINTS                     : ")			<< Fi(refine, "ls_number_restraints") << std::endl

			<< RM3("") << std::endl
			<< RM3(" RMS DEVIATIONS FROM RESTRAINT TARGET VALUES.")  << std::endl
			<< RM3("  BOND LENGTHS                         (A) : ", 7, 3)	<< Ff(ls_restr, c("s_bond_d"), "dev_ideal") << std::endl
			<< RM3("  ANGLE DISTANCES                      (A) : ", 7, 3)	<< Ff(ls_restr, c("s_angle_d"), "dev_ideal") << std::endl
			<< RM3("  SIMILAR DISTANCES (NO TARGET VALUES) (A) : ", 7, 3)	<< Ff(ls_restr, c("s_similar_dist"), "dev_ideal") << std::endl
			<< RM3("  DISTANCES FROM RESTRAINT PLANES      (A) : ", 7, 3)	<< Ff(ls_restr, c("s_from_restr_planes"), "dev_ideal") << std::endl
			<< RM3("  ZERO CHIRAL VOLUMES               (A**3) : ", 7, 3)	<< Ff(ls_restr, c("s_zero_chiral_vol"), "dev_ideal") << std::endl
			<< RM3("  NON-ZERO CHIRAL VOLUMES           (A**3) : ", 7, 3)	<< Ff(ls_restr, c("s_non_zero_chiral_vol"), "dev_ideal") << std::endl
			<< RM3("  ANTI-BUMPING DISTANCE RESTRAINTS     (A) : ", 7, 3)	<< Ff(ls_restr, c("s_anti_bump_dis_restr"), "dev_ideal") << std::endl
			<< RM3("  RIGID-BOND ADP COMPONENTS         (A**2) : ", 7, 3)	<< Ff(ls_restr, c("s_rigid_bond_adp_cmpnt"), "dev_ideal") << std::endl
			<< RM3("  SIMILAR ADP COMPONENTS            (A**2) : ", 7, 3)	<< Ff(ls_restr, c("s_similar_adp_cmpnt"), "dev_ideal") << std::endl
			<< RM3("  APPROXIMATELY ISOTROPIC ADPS      (A**2) : ", 7, 3)	<< Ff(ls_restr, c("s_approx_iso_adps"), "dev_ideal") << std::endl

			<< RM3("") << std::endl
			<< RM3(" BULK SOLVENT MODELING.")  << std::endl
			<< RM3("  METHOD USED: ")										<< Fs(refine, "solvent_model_details") << std::endl

			<< RM3("") << std::endl
			<< RM3(" STEREOCHEMISTRY TARGET VALUES : ")						<< Fs(refine, "pdbx_stereochemistry_target_values") << std::endl
			<< RM3("  SPECIAL CASE: ")										<< Fs(refine, "pdbx_stereochem_target_val_spec_case") << std::endl
			
			<< RM3("") << std::endl;
}

void WriteRemark3Phenix(std::ostream& pdbFile, Datablock& db)
{
	auto refine = db["refine"].front();
//	auto ls_shell = db["refine_ls_shell"].front();
//	auto hist = db["refine_hist"].front();
	auto reflns = db["reflns"].front();
//	auto analyze = db["refine_analyze"].front();
	auto& ls_restr = db["refine_ls_restr"];
//	auto pdbx_xplor_file = db["pdbx_xplor_file"].front();
	auto pdbx_reflns_twin = db["pdbx_reflns_twin"].front();
	
	auto c = [](const char* t) -> cif::Condition { return cif::Key("type") == t; };

	pdbFile	<< RM3("") << std::endl
			<< RM3("   REFINEMENT TARGET : ") << Fs(refine, "pdbx_stereochemistry_target_values") << std::endl
			<< RM3("") << std::endl
			<< RM3(" DATA USED IN REFINEMENT.") << std::endl
			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_high") << std::endl
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_low") << std::endl
			<< RM3("  MIN(FOBS/SIGMA_FOBS)              : ", 6, 3)	<< Ff(refine, "pdbx_ls_sigma_F") << std::endl
			<< RM3("  COMPLETENESS FOR RANGE        (%) : ", 5, 2)	<< Ff(refine, "ls_percent_reflns_obs") << std::endl
			<< RM3("  NUMBER OF REFLECTIONS             : ", 12, 6)	<< Fi(refine, "ls_number_reflns_obs") << std::endl
			<< RM3("") << std::endl
			<< RM3(" FIT TO DATA USED IN REFINEMENT.") << std::endl
			<< RM3("  R VALUE     (WORKING + TEST SET) : ", 7, 5)	<< Ff(refine, "ls_R_factor_obs") << std::endl
			<< RM3("  R VALUE            (WORKING SET) : ", 7, 5)	<< Ff(refine, "ls_R_factor_R_work") << std::endl
			<< RM3("  FREE R VALUE                     : ", 7, 5)	<< Ff(refine, "ls_R_factor_R_free") << std::endl
			<< RM3("  FREE R VALUE TEST SET SIZE   (%) : ", 7, 3)	<< Ff(refine, "ls_percent_reflns_R_free") << std::endl
			<< RM3("  FREE R VALUE TEST SET COUNT      : ", 12, 6)	<< Fi(refine, "ls_number_reflns_R_free") << std::endl

			<< RM3("") << std::endl
			<< RM3(" FIT TO DATA USED IN REFINEMENT (IN BINS).") << std::endl
			<< RM3("  BIN  RESOLUTION RANGE  COMPL.    NWORK NFREE   RWORK  RFREE") << std::endl;
		
	int bin = 1;
	std::vector<Row> bins;
	for (auto r: db["refine_ls_shell"])
		bins.push_back(r);
//	reverse(bins.begin(), bins.end());
	try
	{
		sort(bins.begin(), bins.end(), [](Row a, Row b) -> bool { return a["d_res_high"].as<float>() > b["d_res_high"].as<float>(); });
	}
	catch (...) {}

	for (auto r: bins)
	{
		boost::format fmt("%3.3d %7.4f - %7.4f    %4.2f %8.8d %5.5d  %6.4f %6.4f");
		
		float d_res_low, d_res_high, percent_reflns_obs, R_factor_R_work, R_factor_R_free;
		int number_reflns_R_work, number_reflns_R_free;
		
		cif::tie(d_res_low, d_res_high, percent_reflns_obs, number_reflns_R_work,
				number_reflns_R_free, R_factor_R_work, R_factor_R_free) =
			r.get("d_res_low", "d_res_high", "percent_reflns_obs", "number_reflns_R_work",
				"number_reflns_R_free", "R_factor_R_work", "R_factor_R_free");
		
		percent_reflns_obs /= 100;

		pdbFile << RM3("  ") << fmt % bin++ % d_res_low % d_res_high % percent_reflns_obs % number_reflns_R_work %
			number_reflns_R_free % R_factor_R_work % R_factor_R_free << std::endl;
	}
	
	pdbFile << RM3("") << std::endl
			<< RM3(" BULK SOLVENT MODELLING.") << std::endl
			<< RM3("  METHOD USED        : ")				<< Fs(refine, "solvent_model_details") << std::endl
			<< RM3("  SOLVENT RADIUS     : ", 5, 2)			<< Ff(refine, "pdbx_solvent_vdw_probe_radii") << std::endl
			<< RM3("  SHRINKAGE RADIUS   : ", 5, 2)			<< Ff(refine, "pdbx_solvent_shrinkage_radii") << std::endl
			<< RM3("  K_SOL              : ", 5, 2)         << Ff(refine, "solvent_model_param_ksol") << std::endl
			<< RM3("  B_SOL              : ", 5, 2)         << Ff(refine, "solvent_model_param_bsol") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" ERROR ESTIMATES.") << std::endl
			<< RM3("  COORDINATE ERROR (MAXIMUM-LIKELIHOOD BASED)     : ", 6, 3)	<< Ff(refine, "overall_SU_ML") << std::endl
			<< RM3("  PHASE ERROR (DEGREES, MAXIMUM-LIKELIHOOD BASED) : ", 6, 3)	<< Ff(refine, "pdbx_overall_phase_error") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" B VALUES.") << std::endl
			<< RM3("  B VALUE TYPE                      : ")    					<< Fs(refine, "pdbx_TLS_residual_ADP_flag") << std::endl
			<< RM3("  FROM WILSON PLOT           (A**2) : ", 7, 4)					<< Ff(reflns, "B_iso_Wilson_estimate") << std::endl
			<< RM3("  MEAN B VALUE      (OVERALL, A**2) : ", 7, 4)					<< Ff(refine, "B_iso_mean") << std::endl
			<< RM3("  OVERALL ANISOTROPIC B VALUE.") << std::endl
			<< RM3("   B11 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][1]") << std::endl
			<< RM3("   B22 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[2][2]") << std::endl
			<< RM3("   B33 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[3][3]") << std::endl
			<< RM3("   B12 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][2]") << std::endl
			<< RM3("   B13 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][3]") << std::endl
			<< RM3("   B23 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[2][3]") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" TWINNING INFORMATION.") << std::endl
			<< RM3("  FRACTION: ") << Fs(pdbx_reflns_twin, "fraction") << std::endl
			<< RM3("  OPERATOR: ") << Fs(pdbx_reflns_twin, "operator") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" DEVIATIONS FROM IDEAL VALUES.") << std::endl
			<< RM3("                RMSD          COUNT") << std::endl
			<< RM3("  BOND      : ", -6, 3)	<< Ff(ls_restr, c("f_bond_d"             ), "dev_ideal") << SEP("        ", -7)
											<< Fi(ls_restr, c("f_bond_d"             ), "number")
											<< std::endl
			<< RM3("  ANGLE     : ", -6, 3)	<< Ff(ls_restr, c("f_angle_d"            ), "dev_ideal") << SEP("        ", -7)
											<< Fi(ls_restr, c("f_angle_d"            ), "number")
											<< std::endl
			<< RM3("  CHIRALITY : ", -6, 3)	<< Ff(ls_restr, c("f_chiral_restr"       ), "dev_ideal") << SEP("        ", -7)
											<< Fi(ls_restr, c("f_chiral_restr"       ), "number")
											<< std::endl
			<< RM3("  PLANARITY : ", -6, 3)	<< Ff(ls_restr, c("f_plane_restr"        ), "dev_ideal") << SEP("        ", -7)
											<< Fi(ls_restr, c("f_plane_restr"        ), "number")
											<< std::endl
			<< RM3("  DIHEDRAL  : ", -6, 3)	<< Ff(ls_restr, c("f_dihedral_angle_d"   ), "dev_ideal") << SEP("        ", -7)
											<< Fi(ls_restr, c("f_dihedral_angle_d"   ), "number")
											<< std::endl;

	auto& tls = db["pdbx_refine_tls"];

	pdbFile << RM3("") << std::endl
			<< RM3(" TLS DETAILS") << std::endl
			<< RM3("  NUMBER OF TLS GROUPS  : ") << (tls.size() ? std::to_string(tls.size()) : "NULL") << std::endl;
	
	for (auto t: tls)
	{
		std::string id = t["id"].as<std::string>();

		auto pdbx_refine_tls_group = db["pdbx_refine_tls_group"][cif::Key("refine_tls_id") == id];
		
		pdbFile << RM3("  TLS GROUP : ") << id << std::endl
				<< RM3("   SELECTION: ") << Fs(pdbx_refine_tls_group, "selection_details") << std::endl
				<< RM3("   ORIGIN FOR THE GROUP (A):", -9, 4) << Ff(t, "origin_x")
				<< SEP("", -9, 4) << Ff(t, "origin_y")
				<< SEP("", -9, 4) << Ff(t, "origin_z") << std::endl
				<< RM3("   T TENSOR") << std::endl
				<< RM3("     T11:", -9, 4) << Ff(t, "T[1][1]") << SEP(" T22:", -9, 4) << Ff(t, "T[2][2]") << std::endl
				<< RM3("     T33:", -9, 4) << Ff(t, "T[3][3]") << SEP(" T12:", -9, 4) << Ff(t, "T[1][2]") << std::endl
				<< RM3("     T13:", -9, 4) << Ff(t, "T[1][3]") << SEP(" T23:", -9, 4) << Ff(t, "T[2][3]") << std::endl
				<< RM3("   L TENSOR") << std::endl
				<< RM3("     L11:", -9, 4) << Ff(t, "L[1][1]") << SEP(" L22:", -9, 4) << Ff(t, "L[2][2]") << std::endl
				<< RM3("     L33:", -9, 4) << Ff(t, "L[3][3]") << SEP(" L12:", -9, 4) << Ff(t, "L[1][2]") << std::endl
				<< RM3("     L13:", -9, 4) << Ff(t, "L[1][3]") << SEP(" L23:", -9, 4) << Ff(t, "L[2][3]") << std::endl
				<< RM3("   S TENSOR") << std::endl
				<< RM3("     S11:", -9, 4) << Ff(t, "S[1][1]") << SEP(" S12:", -9, 4) << Ff(t, "S[1][2]") << SEP(" S13:", -9, 4) << Ff(t, "S[1][3]") << std::endl
				<< RM3("     S21:", -9, 4) << Ff(t, "S[2][1]") << SEP(" S22:", -9, 4) << Ff(t, "S[2][2]") << SEP(" S23:", -9, 4) << Ff(t, "S[2][3]") << std::endl
				<< RM3("     S31:", -9, 4) << Ff(t, "S[3][1]") << SEP(" S32:", -9, 4) << Ff(t, "S[3][2]") << SEP(" S33:", -9, 4) << Ff(t, "S[3][3]") << std::endl;
	}


	pdbFile << RM3("") << std::endl
			<< RM3(" NCS DETAILS") << std::endl;

	auto& ncs_dom = db["struct_ncs_dom"];
	if (ncs_dom.empty())
		pdbFile << RM3("  NUMBER OF NCS GROUPS : NULL") << std::endl;
	else
	{
		std::set<std::string> ncs_groups;
		for (auto i: ncs_dom)
			ncs_groups.insert(i["pdbx_ens_id"].as<std::string>());

		pdbFile << RM3("  NUMBER OF NCS GROUPS : ") << ncs_groups.size() << std::endl;
//			
//			for (auto ens_id: ncs_groups)
//			{
//				auto lim = db["struct_ncs_dom_lim"].find(cif::Key("pdbx_ens_id") == ens_id);
//				
//				set<std::string> chains;
//				set<int> component_ids;
//				
//				for (auto l: lim)
//				{
//					chains.insert(l["beg_auth_asym_id"]);
//					component_ids.insert(l["pdbx_component_id"].as<int>());
//				}
//				
//				pdbFile << RM3("") << std::endl
//						<< RM3(" NCS GROUP NUMBER               : ") << ens_id << std::endl
//						<< RM3("    CHAIN NAMES                    : ") << ba::join(chains, " ") << std::endl
//						<< RM3("    NUMBER OF COMPONENTS NCS GROUP : ") << component_ids.size() << std::endl
//						<< RM3("      COMPONENT C  SSSEQI  TO  C   SSSEQI   CODE") << std::endl;
//				
//				for (auto l: lim)
//				{
//					pdbFile << RM3("         ", -2)		<< Fi(l, "pdbx_component_id")
//							<< SEP(" ", -5)			<< Fs(l, "beg_auth_asym_id")
//							<< SEP("  ", -5)			<< Fi(l, "beg_auth_seq_id")
//							<< SEP("   ", -5)			<< Fs(l, "end_auth_asym_id")
//							<< SEP("   ", -5)			<< Fi(l, "end_auth_seq_id")
//							<< SEP("  ", -5)			<< Fs(l, "pdbx_refine_code")
//							<< std::endl;
//				}
//				
//				pdbFile << RM3("                  GROUP CHAIN        COUNT   RMS     WEIGHT") << std::endl;
//				for (auto l: db["refine_ls_restr_ncs"].find(cif::Key("pdbx_ens_id") == ens_id))
//				{
//					std::string type = l["pdbx_type"];
//					ba::to_upper(type);
//					
//					std::string unit;
//					if (ba::ends_with(type, "POSITIONAL"))
//						unit = "    (A): ";
//					else if (ba::ends_with(type, "THERMAL"))
//						unit = " (A**2): ";
//					else
//						unit = "       : ";
//					
//					pdbFile << RM3("  ", 18)			<< type
//							<< SEP("", -2)				<< Fi(l, "pdbx_ens_id")
//							<< SEP("    ", 1)			<< Fs(l, "pdbx_auth_asym_id")
//							<< SEP(unit.c_str(), -6)	<< Fi(l, "pdbx_number")
//							<< SEP(" ;", -6, 3)		<< Ff(l, "rms_dev_position")
//							<< SEP(" ;", -6, 3)		<< Ff(l, "weight_position")
//							<< std::endl;
//				}
//			}
	}

//	pdbFile << RM3("") << std::endl
//			<< RM3(" BULK SOLVENT MODELLING.") << std::endl
//			<< RM3("  METHOD USED : ") << Fs(refine, "solvent_model_details") << std::endl
//			<< RM3("  PARAMETERS FOR MASK CALCULATION") << std::endl
//			<< RM3("  VDW PROBE RADIUS   : ", 5, 2) << Ff(refine, "pdbx_solvent_vdw_probe_radii") << std::endl
//			<< RM3("  ION PROBE RADIUS   : ", 5, 2) << Ff(refine, "pdbx_solvent_ion_probe_radii") << std::endl
//			<< RM3("  SHRINKAGE RADIUS   : ", 5, 2) << Ff(refine, "pdbx_solvent_shrinkage_radii") << std::endl
//			
//			<< RM3("") << std::endl;

	pdbFile << RM3("") << std::endl;
}

void WriteRemark3XPlor(std::ostream& pdbFile, Datablock& db)
{
	auto refine = db["refine"].front();
	auto ls_shell = db["refine_ls_shell"].front();
	auto hist = db["refine_hist"].front();
	auto reflns = db["reflns"].front();
	auto analyze = db["refine_analyze"].front();
	auto& ls_restr = db["refine_ls_restr"];
	auto ls_restr_ncs = db["refine_ls_restr_ncs"].front();
	auto pdbx_xplor_file = db["pdbx_xplor_file"].front();
	
	pdbFile	<< RM3("") << std::endl
			<< RM3(" DATA USED IN REFINEMENT.") << std::endl
			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_high") << std::endl
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_low") << std::endl
			<< RM3("  DATA CUTOFF            (SIGMA(F)) : ", 6, 3)	<< Ff(refine, "pdbx_ls_sigma_F") << std::endl
			<< RM3("  DATA CUTOFF HIGH         (ABS(F)) : ", 6, 3)	<< Ff(refine, "pdbx_data_cutoff_high_absF") << std::endl
			<< RM3("  DATA CUTOFF LOW          (ABS(F)) : ", 6, 3)	<< Ff(refine, "pdbx_data_cutoff_low_absF") << std::endl
			<< RM3("  COMPLETENESS (WORKING+TEST)   (%) : ", 5, 2)	<< Ff(refine, "ls_percent_reflns_obs") << std::endl
			<< RM3("  NUMBER OF REFLECTIONS             : ", 12, 6)	<< Fi(refine, "ls_number_reflns_obs") << std::endl
	
			<< RM3("") << std::endl
			<< RM3(" FIT TO DATA USED IN REFINEMENT.") << std::endl
			<< RM3("  CROSS-VALIDATION METHOD          : ")	<< Fs(refine, "pdbx_ls_cross_valid_method") << std::endl
			<< RM3("  FREE R VALUE TEST SET SELECTION  : ")	<< Fs(refine, "pdbx_R_Free_selection_details") << std::endl
			<< RM3("  R VALUE            (WORKING SET) : ", 7, 3)	<< Ff(refine, "ls_R_factor_R_work") << std::endl
			<< RM3("  FREE R VALUE                     : ", 7, 3)	<< Ff(refine, "ls_R_factor_R_free") << std::endl
			<< RM3("  FREE R VALUE TEST SET SIZE   (%) : ", 7, 3)	<< Ff(refine, "ls_percent_reflns_R_free") << std::endl
			<< RM3("  FREE R VALUE TEST SET COUNT      : ", 12, 6)	<< Fi(refine, "ls_number_reflns_R_free") << std::endl
			<< RM3("  ESTIMATED ERROR OF FREE R VALUE  : ", 7, 3)	<< Ff(refine, "ls_R_factor_R_free_error") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" FIT IN THE HIGHEST RESOLUTION BIN.") << std::endl
			<< RM3("  TOTAL NUMBER OF BINS USED           : ", 12, 6)	<< Fi(ls_shell, "pdbx_total_number_of_bins_used") << std::endl
			<< RM3("  BIN RESOLUTION RANGE HIGH       (A) : ", 5, 2)	<< Ff(ls_shell, "d_res_high") << std::endl
			<< RM3("  BIN RESOLUTION RANGE LOW        (A) : ", 5, 2)	<< Ff(ls_shell, "d_res_low") << std::endl
			<< RM3("  BIN COMPLETENESS (WORKING+TEST) (%) : ", 5, 1)	<< Ff(ls_shell, "percent_reflns_obs") << std::endl
			<< RM3("  REFLECTIONS IN BIN    (WORKING SET) : ", 12, 6)	<< Fi(ls_shell, "number_reflns_R_work") << std::endl
			<< RM3("  BIN R VALUE           (WORKING SET) : ", 7, 3)	<< Ff(ls_shell, "R_factor_R_work") << std::endl
			<< RM3("  BIN FREE R VALUE                    : ", 7, 3)	<< Ff(ls_shell, "R_factor_R_free") << std::endl
			<< RM3("  BIN FREE R VALUE TEST SET SIZE  (%) : ", 5, 1)	<< Ff(ls_shell, "percent_reflns_R_free") << std::endl
			<< RM3("  BIN FREE R VALUE TEST SET COUNT     : ", 12, 6)	<< Fi(ls_shell, "number_reflns_R_free") << std::endl
			<< RM3("  ESTIMATED ERROR OF BIN FREE R VALUE : ", 7, 3)	<< Ff(ls_shell, "R_factor_R_free_error") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.") << std::endl
			<< RM3("  PROTEIN ATOMS            : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_protein") << std::endl
			<< RM3("  NUCLEIC ACID ATOMS       : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_nucleic_acid") << std::endl
			<< RM3("  HETEROGEN ATOMS          : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_ligand") << std::endl
			<< RM3("  SOLVENT ATOMS            : ", 12, 6)	<< Fi(hist, "number_atoms_solvent") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" B VALUES.") << std::endl
			<< RM3("  FROM WILSON PLOT           (A**2) : ", 7, 2)	<< Ff(reflns, "B_iso_Wilson_estimate") << std::endl
			<< RM3("  MEAN B VALUE      (OVERALL, A**2) : ", 7, 2)	<< Ff(refine, "B_iso_mean") << std::endl
			
			<< RM3("  OVERALL ANISOTROPIC B VALUE.") << std::endl
			<< RM3("   B11 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][1]") << std::endl
			<< RM3("   B22 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[2][2]") << std::endl
			<< RM3("   B33 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[3][3]") << std::endl
			<< RM3("   B12 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][2]") << std::endl
			<< RM3("   B13 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][3]") << std::endl
			<< RM3("   B23 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[2][3]") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" ESTIMATED COORDINATE ERROR.") << std::endl
			<< RM3("  ESD FROM LUZZATI PLOT        (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_coordinate_error_obs") << std::endl
			<< RM3("  ESD FROM SIGMAA              (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_sigma_a_obs") << std::endl
			<< RM3("  LOW RESOLUTION CUTOFF        (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_d_res_low_obs") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" CROSS-VALIDATED ESTIMATED COORDINATE ERROR.") << std::endl
			<< RM3("  ESD FROM C-V LUZZATI PLOT    (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_coordinate_error_free") << std::endl
			<< RM3("  ESD FROM C-V SIGMAA          (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_sigma_a_free") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" RMS DEVIATIONS FROM IDEAL VALUES.") << std::endl
			<< RM3("  BOND LENGTHS                 (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "x_bond_d", "dev_ideal") << std::endl
			<< RM3("  BOND ANGLES            (DEGREES) : ", 7, 2)	<< Ff(ls_restr, cif::Key("type") == "x_angle_deg", "dev_ideal") << std::endl
			<< RM3("  DIHEDRAL ANGLES        (DEGREES) : ", 7, 2)	<< Ff(ls_restr, cif::Key("type") == "x_dihedral_angle_d", "dev_ideal") << std::endl
			<< RM3("  IMPROPER ANGLES        (DEGREES) : ", 7, 2)	<< Ff(ls_restr, cif::Key("type") == "x_improper_angle_d", "dev_ideal") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" ISOTROPIC THERMAL MODEL : ") << Fs(refine, "pdbx_isotropic_thermal_model") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" ISOTROPIC THERMAL FACTOR RESTRAINTS.    RMS    SIGMA") << std::endl
			<< RM3("  MAIN-CHAIN BOND              (A**2) : ", 6, 2) << Ff(ls_restr, cif::Key("type") == "x_mcbond_it", "dev_ideal") << SEP("; ", 6, 2) 
																	 << Ff(ls_restr, cif::Key("type") == "x_mcbond_it", "dev_ideal_target") << std::endl
			<< RM3("  MAIN-CHAIN ANGLE             (A**2) : ", 6, 2) << Ff(ls_restr, cif::Key("type") == "x_mcangle_it", "dev_ideal") << SEP("; ", 6, 2) 
																	 << Ff(ls_restr, cif::Key("type") == "x_mcangle_it", "dev_ideal_target") << std::endl
			<< RM3("  SIDE-CHAIN BOND              (A**2) : ", 6, 2) << Ff(ls_restr, cif::Key("type") == "x_scbond_it", "dev_ideal") << SEP("; ", 6, 2) 
																	 << Ff(ls_restr, cif::Key("type") == "x_scbond_it", "dev_ideal_target") << std::endl
			<< RM3("  SIDE-CHAIN ANGLE             (A**2) : ", 6, 2) << Ff(ls_restr, cif::Key("type") == "x_scangle_it", "dev_ideal") << SEP("; ", 6, 2) 
																	 << Ff(ls_restr, cif::Key("type") == "x_scangle_it", "dev_ideal_target") << std::endl
			<< RM3("") << std::endl
			<< RM3(" NCS MODEL : ")	<< Fs(ls_restr_ncs, "ncs_model_details") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" NCS RESTRAINTS.                         RMS   SIGMA/WEIGHT") << std::endl
			
			// TODO: using only group 1 here, should this be fixed???
			<< RM3("  GROUP  1  POSITIONAL            (A) : ", 4, 2) << Ff(ls_restr_ncs, "rms_dev_position") << SEP("; ", 6, 2) 
																		<< Ff(ls_restr_ncs, "weight_position") << SEP("; ", 6, 2) << std::endl
			<< RM3("  GROUP  1  B-FACTOR           (A**2) : ", 4, 2) << Ff(ls_restr_ncs, "rms_dev_B_iso") << SEP("; ", 6, 2) 
																		<< Ff(ls_restr_ncs, "weight_B_iso") << SEP("; ", 6, 2) << std::endl
			
			// TODO: using only files from serial_no 1 here
			<< RM3("") << std::endl
			<< RM3(" PARAMETER FILE   1  : ") << Fs(pdbx_xplor_file, "param_file") << std::endl
			<< RM3(" TOPOLOGY FILE   1   : ") << Fs(pdbx_xplor_file, "topol_file") << std::endl
			
			<< RM3("") << std::endl;
}

void WriteRemark3NuclSQ(std::ostream& pdbFile, Datablock& db)
{
	auto refine = db["refine"].front();
	auto pdbx_refine = db["pdbx_refine"].front();
	auto hist = db["refine_hist"].front();
	auto reflns = db["reflns"].front();
	auto analyze = db["refine_analyze"].front();
	auto& ls_restr = db["refine_ls_restr"];
	
	pdbFile	<< RM3("") << std::endl
			<< RM3(" DATA USED IN REFINEMENT.") << std::endl

			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_high") << std::endl
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_low") << std::endl
			<< RM3("  DATA CUTOFF            (SIGMA(F)) : ", 6, 3)	<< Ff(refine, "pdbx_ls_sigma_F") << std::endl
			<< RM3("  COMPLETENESS FOR RANGE        (%) : ", 5, 2)	<< Ff(refine, "ls_percent_reflns_obs") << std::endl
			<< RM3("  NUMBER OF REFLECTIONS             : ", 12, 6)	<< Fi(refine, "ls_number_reflns_obs") << std::endl

			<< RM3("") << std::endl
			<< RM3(" FIT TO DATA USED IN REFINEMENT.") << std::endl
			<< RM3("  CROSS-VALIDATION METHOD          : ")	<< Fs(refine, "pdbx_ls_cross_valid_method") << std::endl
			<< RM3("  FREE R VALUE TEST SET SELECTION  : ")	<< Fs(refine, "pdbx_R_Free_selection_details") << std::endl
			<< RM3("  R VALUE     (WORKING + TEST SET) : ", 7, 3)	<< Ff(refine, "ls_R_factor_obs") << std::endl
			<< RM3("  R VALUE            (WORKING SET) : ", 7, 3)	<< Ff(refine, "ls_R_factor_R_work") << std::endl
			<< RM3("  FREE R VALUE                     : ", 7, 3)	<< Ff(refine, "ls_R_factor_R_free") << std::endl
			<< RM3("  FREE R VALUE TEST SET SIZE   (%) : ", 7, 3)	<< Ff(refine, "ls_percent_reflns_R_free") << std::endl
			<< RM3("  FREE R VALUE TEST SET COUNT      : ", 12, 6)	<< Fi(refine, "ls_number_reflns_R_free") << std::endl

			<< RM3("") << std::endl
			<< RM3(" FIT/AGREEMENT OF MODEL WITH ALL DATA.") << std::endl
			<< RM3("  R VALUE   (WORKING + TEST SET, NO CUTOFF) : ")	<< Fs(refine, "ls_R_factor_all") << std::endl
			<< RM3("  R VALUE          (WORKING SET, NO CUTOFF) : ")	<< Fs(pdbx_refine, "R_factor_obs_no_cutoff") << std::endl

			<< RM3("  FREE R VALUE                  (NO CUTOFF) : ")	<< Fs(pdbx_refine, "free_R_factor_no_cutoff") << std::endl
			<< RM3("  FREE R VALUE TEST SET SIZE (%, NO CUTOFF) : ")	<< Fs(pdbx_refine, "free_R_val_test_set_size_perc_no_cutoff") << std::endl
			<< RM3("  FREE R VALUE TEST SET COUNT   (NO CUTOFF) : ")	<< Fs(pdbx_refine, "free_R_val_test_set_ct_no_cutoff") << std::endl
			<< RM3("  TOTAL NUMBER OF REFLECTIONS   (NO CUTOFF) : ")	<< Fs(refine, "ls_number_reflns_all") << std::endl

			<< RM3("") << std::endl
			<< RM3(" NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.") << std::endl
			<< RM3("  PROTEIN ATOMS            : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_protein") << std::endl
			<< RM3("  NUCLEIC ACID ATOMS       : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_nucleic_acid") << std::endl
			<< RM3("  HETEROGEN ATOMS          : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_ligand") << std::endl
			<< RM3("  SOLVENT ATOMS            : ", 12, 6)	<< Fi(hist, "number_atoms_solvent") << std::endl
//			<< RM3("  ALL ATOMS                : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_protein") << std::endl

			<< RM3("") << std::endl
			<< RM3(" B VALUES.") << std::endl
//			<< RM3("  B VALUE TYPE                      : ", 7, 2)	<< Fs(refine, "pdbx_TLS_residual_ADP_flag") << std::endl
			<< RM3("  FROM WILSON PLOT           (A**2) : ", 7, 2)	<< Ff(reflns, "B_iso_Wilson_estimate") << std::endl
			<< RM3("  MEAN B VALUE      (OVERALL, A**2) : ", 7, 2)	<< Ff(refine, "B_iso_mean") << std::endl
			<< RM3("  OVERALL ANISOTROPIC B VALUE.") << std::endl
			<< RM3("   B11 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][1]") << std::endl
			<< RM3("   B22 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[2][2]") << std::endl
			<< RM3("   B33 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[3][3]") << std::endl
			<< RM3("   B12 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][2]") << std::endl
			<< RM3("   B13 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][3]") << std::endl
			<< RM3("   B23 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[2][3]") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" ESTIMATED COORDINATE ERROR.") << std::endl
			<< RM3("  ESD FROM LUZZATI PLOT        (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_coordinate_error_obs") << std::endl
			<< RM3("  ESD FROM SIGMAA              (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_sigma_a_obs") << std::endl
			<< RM3("  LOW RESOLUTION CUTOFF        (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_d_res_low_obs") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" RMS DEVIATIONS FROM IDEAL VALUES.") << std::endl
			<< RM3("  DISTANCE RESTRAINTS.                    RMS     SIGMA") << std::endl
			<< RM3("   SUGAR-BASE BOND DISTANCE        (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "n_sugar_bond_d", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "n_sugar_bond_d", "dev_ideal_target") << std::endl
			<< RM3("   SUGAR-BASE BOND ANGLE DISTANCE  (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "n_sugar_bond_angle_d", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "n_sugar_bond_angle_d", "dev_ideal_target") << std::endl
			<< RM3("   PHOSPHATE BONDS DISTANCE        (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "n_phos_bond_d", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "n_phos_bond_d", "dev_ideal_target") << std::endl
			<< RM3("   PHOSPHATE BOND ANGLE, H-BOND    (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "n_phos_bond_angle_d", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "n_phos_bond_angle_d", "dev_ideal_target") << std::endl

			<< RM3("") << std::endl
			<< RM3("  PLANE RESTRAINT                  (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "n_plane_restr", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "n_plane_restr", "dev_ideal_target") << std::endl
			<< RM3("  CHIRAL-CENTER RESTRAINT       (A**3) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "n_chiral_restr", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "n_chiral_restr", "dev_ideal_target") << std::endl

			<< RM3("") << std::endl
			<< RM3("  NON-BONDED CONTACT RESTRAINTS.") << std::endl
			<< RM3("   SINGLE TORSION CONTACT          (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "n_singtor_nbd", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "n_singtor_nbd", "dev_ideal_target") << std::endl
			<< RM3("   MULTIPLE TORSION CONTACT        (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "n_multtor_nbd", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "n_multtor_nbd", "dev_ideal_target") << std::endl

			<< RM3("") << std::endl
			<< RM3(" ISOTROPIC THERMAL FACTOR RESTRAINTS.    RMS     SIGMA") << std::endl
			<< RM3("  SUGAR-BASE BONDS             (A**2) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "n_sugar_bond_it", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "n_sugar_bond_it", "dev_ideal_target") << std::endl
			<< RM3("  SUGAR-BASE ANGLES            (A**2) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "n_sugar_angle_it", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "n_sugar_angle_it", "dev_ideal_target") << std::endl
			<< RM3("  PHOSPHATE BONDS              (A**2) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "n_phos_bond_it", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "n_phos_bond_it", "dev_ideal_target") << std::endl
			<< RM3("  PHOSPHATE BOND ANGLE, H-BOND (A**2) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "n_phos_angle_it", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "n_phos_angle_it", "dev_ideal_target") << std::endl

			<< RM3("") << std::endl;
}

void WriteRemark3ProlSQ(std::ostream& pdbFile, Datablock& db)
{
	auto refine = db["refine"].front();
	auto pdbx_refine = db["pdbx_refine"].front();
	auto hist = db["refine_hist"].front();
	auto reflns = db["reflns"].front();
	auto analyze = db["refine_analyze"].front();
	auto& ls_restr = db["refine_ls_restr"];
	
	pdbFile	<< RM3("") << std::endl
			<< RM3(" DATA USED IN REFINEMENT.") << std::endl

			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_high") << std::endl
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2)	<< Ff(refine, "ls_d_res_low") << std::endl
			<< RM3("  DATA CUTOFF            (SIGMA(F)) : ", 6, 3)	<< Ff(refine, "pdbx_ls_sigma_F") << std::endl
			<< RM3("  COMPLETENESS FOR RANGE        (%) : ", 5, 2)	<< Ff(refine, "ls_percent_reflns_obs") << std::endl
			<< RM3("  NUMBER OF REFLECTIONS             : ", 12, 6)	<< Fi(refine, "ls_number_reflns_obs") << std::endl

			<< RM3("") << std::endl
			<< RM3(" FIT TO DATA USED IN REFINEMENT.") << std::endl
			<< RM3("  CROSS-VALIDATION METHOD          : ")	<< Fs(refine, "pdbx_ls_cross_valid_method") << std::endl
			<< RM3("  FREE R VALUE TEST SET SELECTION  : ")	<< Fs(refine, "pdbx_R_Free_selection_details") << std::endl
			<< RM3("  R VALUE     (WORKING + TEST SET) : ", 7, 3)	<< Ff(refine, "ls_R_factor_obs") << std::endl
			<< RM3("  R VALUE            (WORKING SET) : ", 7, 3)	<< Ff(refine, "ls_R_factor_R_work") << std::endl
			<< RM3("  FREE R VALUE                     : ", 7, 3)	<< Ff(refine, "ls_R_factor_R_free") << std::endl
			<< RM3("  FREE R VALUE TEST SET SIZE   (%) : ", 7, 3)	<< Ff(refine, "ls_percent_reflns_R_free") << std::endl
			<< RM3("  FREE R VALUE TEST SET COUNT      : ", 12, 6)	<< Fi(refine, "ls_number_reflns_R_free") << std::endl

			<< RM3("") << std::endl
			<< RM3(" FIT/AGREEMENT OF MODEL WITH ALL DATA.") << std::endl
			<< RM3("  R VALUE   (WORKING + TEST SET, NO CUTOFF) : ")	<< Fs(refine, "ls_R_factor_all") << std::endl
			<< RM3("  R VALUE          (WORKING SET, NO CUTOFF) : ")	<< Fs(pdbx_refine, "R_factor_obs_no_cutoff") << std::endl

			<< RM3("  FREE R VALUE                  (NO CUTOFF) : ")	<< Fs(pdbx_refine, "free_R_factor_no_cutoff") << std::endl
			<< RM3("  FREE R VALUE TEST SET SIZE (%, NO CUTOFF) : ")	<< Fs(pdbx_refine, "free_R_val_test_set_size_perc_no_cutoff") << std::endl
			<< RM3("  FREE R VALUE TEST SET COUNT   (NO CUTOFF) : ")	<< Fs(pdbx_refine, "free_R_val_test_set_ct_no_cutoff") << std::endl
			<< RM3("  TOTAL NUMBER OF REFLECTIONS   (NO CUTOFF) : ")	<< Fs(refine, "ls_number_reflns_all") << std::endl

			<< RM3("") << std::endl
			<< RM3(" NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.") << std::endl
			<< RM3("  PROTEIN ATOMS            : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_protein") << std::endl
			<< RM3("  NUCLEIC ACID ATOMS       : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_nucleic_acid") << std::endl
			<< RM3("  HETEROGEN ATOMS          : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_ligand") << std::endl
			<< RM3("  SOLVENT ATOMS            : ", 12, 6)	<< Fi(hist, "number_atoms_solvent") << std::endl
//			<< RM3("  ALL ATOMS                : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_protein") << std::endl

			<< RM3("") << std::endl
			<< RM3(" B VALUES.") << std::endl
//			<< RM3("  B VALUE TYPE                      : ", 7, 2)	<< Fs(refine, "pdbx_TLS_residual_ADP_flag") << std::endl
			<< RM3("  FROM WILSON PLOT           (A**2) : ", 7, 2)	<< Ff(reflns, "B_iso_Wilson_estimate") << std::endl
			<< RM3("  MEAN B VALUE      (OVERALL, A**2) : ", 7, 2)	<< Ff(refine, "B_iso_mean") << std::endl
			<< RM3("  OVERALL ANISOTROPIC B VALUE.") << std::endl
			<< RM3("   B11 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][1]") << std::endl
			<< RM3("   B22 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[2][2]") << std::endl
			<< RM3("   B33 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[3][3]") << std::endl
			<< RM3("   B12 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][2]") << std::endl
			<< RM3("   B13 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[1][3]") << std::endl
			<< RM3("   B23 (A**2) : ", -7, 2)	<< Ff(refine, "aniso_B[2][3]") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" ESTIMATED COORDINATE ERROR.") << std::endl
			<< RM3("  ESD FROM LUZZATI PLOT        (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_coordinate_error_obs") << std::endl
			<< RM3("  ESD FROM SIGMAA              (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_sigma_a_obs") << std::endl
			<< RM3("  LOW RESOLUTION CUTOFF        (A) : ", 7, 2)	<< Ff(analyze, "Luzzati_d_res_low_obs") << std::endl
			
			<< RM3("") << std::endl
			<< RM3(" RMS DEVIATIONS FROM IDEAL VALUES.") << std::endl
			<< RM3("  DISTANCE RESTRAINTS.                    RMS    SIGMA") << std::endl
			<< RM3("   BOND LENGTH                     (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_bond_d", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_bond_d", "dev_ideal_target") << std::endl
			<< RM3("   ANGLE DISTANCE                  (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_angle_d", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_angle_d", "dev_ideal_target") << std::endl
			<< RM3("   INTRAPLANAR 1-4 DISTANCE        (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_planar_d", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_planar_d", "dev_ideal_target") << std::endl
			<< RM3("   H-BOND OR METAL COORDINATION    (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_hb_or_metal_coord", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_hb_or_metal_coord", "dev_ideal_target") << std::endl

			<< RM3("") << std::endl
			<< RM3("  PLANE RESTRAINT                 (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_plane_restr", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_plane_restr", "dev_ideal_target") << std::endl
			<< RM3("  CHIRAL-CENTER RESTRAINT      (A**3) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_chiral_restr", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_chiral_restr", "dev_ideal_target") << std::endl

			<< RM3("") << std::endl
			<< RM3("  NON-BONDED CONTACT RESTRAINTS.") << std::endl
			<< RM3("   SINGLE TORSION                  (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_singtor_nbd", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_singtor_nbd", "dev_ideal_target") << std::endl
			<< RM3("   MULTIPLE TORSION                (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_multtor_nbd", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_multtor_nbd", "dev_ideal_target") << std::endl
			<< RM3("   H-BOND (X...Y)                  (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_xyhbond_nbd", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_xyhbond_nbd", "dev_ideal_target") << std::endl
			<< RM3("   H-BOND (X-H...Y)                (A) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_xhyhbond_nbd", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_xhyhbond_nbd", "dev_ideal_target") << std::endl
			
			<< RM3("") << std::endl
			<< RM3("  CONFORMATIONAL TORSION ANGLE RESTRAINTS.") << std::endl
			<< RM3("   SPECIFIED                 (DEGREES) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_special_tor", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_special_tor", "dev_ideal_target") << std::endl
			<< RM3("   PLANAR                    (DEGREES) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_planar_tor", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_planar_tor", "dev_ideal_target") << std::endl
			<< RM3("   STAGGERED                 (DEGREES) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_staggered_tor", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_staggered_tor", "dev_ideal_target") << std::endl
			<< RM3("   TRANSVERSE                (DEGREES) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_transverse_tor", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_transverse_tor", "dev_ideal_target") << std::endl

			<< RM3("") << std::endl
			<< RM3("  ISOTROPIC THERMAL FACTOR RESTRAINTS. RMS SIGMA") << std::endl
			<< RM3("   MAIN-CHAIN BOND              (A**2) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_mcbond_it", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_mcbond_it", "dev_ideal_target") << std::endl
			<< RM3("   MAIN-CHAIN ANGLE             (A**2) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_mcangle_it", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_mcangle_it", "dev_ideal_target") << std::endl
			<< RM3("   SIDE-CHAIN BOND              (A**2) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_scbond_it", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_scbond_it", "dev_ideal_target") << std::endl
			<< RM3("   SIDE-CHAIN ANGLE             (A**2) : ", 7, 3)	<< Ff(ls_restr, cif::Key("type") == "p_scangle_it", "dev_ideal") << " ; "
																		<< Ff(ls_restr, cif::Key("type") == "p_scangle_it", "dev_ideal_target") << std::endl

			<< RM3("") << std::endl;
}

void WriteRemark3(std::ostream& pdbFile, Datablock& db)
{
	std::string program, authors;

	if (not db["pdbx_nmr_software"].empty())
	{
		auto software = db["pdbx_nmr_software"].find(cif::Key("classification") == "refinement");
		if (software.size() == 1)
			cif::tie(program, authors) = software.front().get("name", "authors");
		else if (software.size() > 1)
		{
			for (auto r: software)
			{
				if (program.empty() == false)
				{
					program += ", ";
					authors += ", ";
				}
				
				program += r["name"].as<std::string>();
				authors += r["authors"].as<std::string>() + " (" + r["name"].as<std::string>() + ")";
			}
		}
	}
	
	if (program.empty())
		program = cifSoftware(db, eRefinement);
	
	if (authors.empty())
		authors = "NULL";
	
	if (not program.empty())
	{
		pdbFile << RM3("") << std::endl
				<< RM3("REFINEMENT.") << std::endl;
		
		int l = 0;
		for (auto s: cif::wordWrap(program, 52))
			pdbFile << RM3(++l == 1 ? "  PROGRAM     : " : "                ") << s << std::endl;

		l = 0;
		for (auto s: cif::wordWrap(authors, 52))
			pdbFile << RM3(++l == 1 ? "  AUTHORS     : " : "                ") << s << std::endl;
	}
	
	if (not db["refine"].empty())
	{
		auto s = program.find(' ');
		if (s != std::string::npos)
			program.erase(s, std::string::npos);
		
		if (cif::iequals(program, "BUSTER") or cif::iequals(program, "BUSTER-TNT") or cif::iequals(program, "TNT"))
			WriteRemark3BusterTNT(pdbFile, db);
		else if (cif::iequals(program, "CNS") or cif::iequals(program, "CNX"))
			WriteRemark3CNS(pdbFile, db);
		else if (cif::iequals(program, "X-PLOR"))
			WriteRemark3XPlor(pdbFile, db);
		else if (cif::iequals(program, "REFMAC"))
			WriteRemark3Refmac(pdbFile, db);
		else if (cif::iequals(program, "SHELXL"))
			WriteRemark3Shelxl(pdbFile, db);
		else if (cif::iequals(program, "PHENIX"))
			WriteRemark3Phenix(pdbFile, db);
		else if (cif::iequals(program, "NUCLSQ"))
			WriteRemark3NuclSQ(pdbFile, db);
		else if (cif::iequals(program, "PROLSQ"))
			WriteRemark3ProlSQ(pdbFile, db);
	}

	for (auto r: db["refine"])
	{
		std::string remarks = r["details"].as<std::string>();
		if (remarks.empty())
			remarks = "NULL";
		
		WriteOneContinuedLine(pdbFile, "REMARK   3 ", 0, "OTHER REFINEMENT REMARKS: " + remarks);
		break;
	}
}

void WriteRemark200(std::ostream& pdbFile, Datablock& db)
{
	typedef RM<200> RM;
	
	try
	{
		for (auto diffrn: db["diffrn"])
		{
			std::string diffrn_id = diffrn["id"].as<std::string>();
			std::string crystal_id = diffrn["crystal_id"].as<std::string>();
			
			auto diffrn_radiation = db["diffrn_radiation"][cif::Key("diffrn_id") == diffrn_id];
			auto diffrn_radiation_wavelength = db["diffrn_radiation_wavelength"][cif::Key("id") == diffrn_radiation["wavelength_id"].as<std::string>()];
			auto diffrn_source = db["diffrn_source"][cif::Key("diffrn_id") == diffrn_id];
			auto diffrn_detector = db["diffrn_detector"][cif::Key("diffrn_id") == diffrn_id];
			auto exptl = db["exptl"][cif::Key("entry_id") == db.getName()];
			auto exptl_crystal = db["exptl_crystal"][cif::Key("id") == crystal_id];
			auto exptl_crystal_grow = db["exptl_crystal_grow"][cif::Key("crystal_id") == crystal_id];
			auto computing = db["computing"][cif::Key("entry_id") == db.getName()];
			auto reflns = db["reflns"][cif::Key("entry_id") == db.getName()];

			std::string pdbx_diffrn_id = reflns["pdbx_diffrn_id"].as<std::string>();
			
			auto reflns_shell = db["reflns_shell"][cif::Key("pdbx_diffrn_id") == pdbx_diffrn_id];
			auto refine = db["refine"][cif::Key("pdbx_diffrn_id") == pdbx_diffrn_id];
		
			std::string date = diffrn_detector["pdbx_collection_date"].as<std::string>();
			if (date.empty())
				date = "NULL";
			else
				date = cif2pdbDate(date);
				
			std::string iis = cifSoftware(db, eDataReduction);
			std::string dss = cifSoftware(db, eDataScaling);

			std::string source = diffrn_source["source"].as<std::string>();
			std::string synchrotron, type;
			
			if (source.empty())
				synchrotron = "NULL";
			else if (cif::iequals(source, "SYNCHROTRON"))
			{
				synchrotron = "Y";
				source = diffrn_source["pdbx_synchrotron_site"].as<std::string>();
				if (source.empty()) source = "NULL";
				type = "NULL";
			}
			else
			{
				synchrotron = "N";
				type = diffrn_source["type"].as<std::string>();
				if (type.empty()) type = "NULL";
			}
			
			if (source.empty()) source = "NULL";
			if (type.empty()) type = "NULL";
		
			pdbFile << RM("") << std::endl
				<< RM("EXPERIMENTAL DETAILS") << std::endl
				<< RM(" EXPERIMENT TYPE                : ")			<< Fs(exptl, "method") << std::endl
				<< RM(" DATE OF DATA COLLECTION        : ")			<< date << std::endl
				<< RM(" TEMPERATURE           (KELVIN) : ", 5, 1)	<< Ff(diffrn, "ambient_temp") << std::endl
				<< RM(" PH                             : ", 4, 1)	<< Ff(exptl_crystal_grow, "ph") << std::endl
				<< RM(" NUMBER OF CRYSTALS USED        : ")			<< Fi(exptl, "crystals_number") << std::endl
				<< RM("") << std::endl
				<< RM(" SYNCHROTRON              (Y/N) : ") << synchrotron << std::endl
				<< RM(" RADIATION SOURCE               : ") << source << std::endl
				<< RM(" BEAMLINE                       : ") << Fs(diffrn_source, "pdbx_synchrotron_beamline") << std::endl
				<< RM(" X-RAY GENERATOR MODEL          : ") << type << std::endl
				<< RM(" MONOCHROMATIC OR LAUE    (M/L) : ") << Fs(diffrn_radiation, "pdbx_monochromatic_or_laue_m_l") << std::endl
				<< RM(" WAVELENGTH OR RANGE        (A) : ", 7, 4) << Ff(diffrn_radiation_wavelength, "wavelength") << std::endl
				<< RM(" MONOCHROMATOR                  : ") << Fs(diffrn_radiation, "monochromator") << std::endl
				<< RM(" OPTICS                         : ") << Fs(diffrn_detector, "details") << std::endl
				<< RM("") << std::endl
				<< RM(" DETECTOR TYPE                  : ") << Fs(diffrn_detector, "detector") << std::endl
				<< RM(" DETECTOR MANUFACTURER          : ") << Fs(diffrn_detector, "type") << std::endl
				<< RM(" INTENSITY-INTEGRATION SOFTWARE : ") << iis << std::endl
				<< RM(" DATA SCALING SOFTWARE          : ") << dss << std::endl
				<< RM(" ") << std::endl
				<< RM(" NUMBER OF UNIQUE REFLECTIONS   : ") << Fi(reflns, "number_obs") << std::endl
				<< RM(" RESOLUTION RANGE HIGH      (A) : ", 7, 3) << Ff(reflns, "d_resolution_high") << std::endl
				<< RM(" RESOLUTION RANGE LOW       (A) : ", 7, 3) << Ff(reflns, "d_resolution_low") << std::endl
				<< RM(" REJECTION CRITERIA  (SIGMA(I)) : ", 7, 3) << Ff(reflns, "observed_criterion_sigma_I") << std::endl
				<< RM("") << std::endl
				<< RM("OVERALL.") << std::endl
				<< RM(" COMPLETENESS FOR RANGE     (%) : ", 7, 1) << Ff(reflns, "percent_possible_obs") << std::endl
				<< RM(" DATA REDUNDANCY                : ", 7, 3) << Ff(reflns, "pdbx_redundancy") << std::endl
				<< RM(" R MERGE                    (I) : ", 7, 5) << Ff(reflns, "pdbx_Rmerge_I_obs") << std::endl
				<< RM(" R SYM                      (I) : ", 7, 5) << Ff(reflns, "pdbx_Rsym_value") << std::endl
				<< RM(" <I/SIGMA(I)> FOR THE DATA SET  : ", 7, 4) << Ff(reflns, "pdbx_netI_over_sigmaI") << std::endl
				<< RM("") << std::endl
				<< RM("IN THE HIGHEST RESOLUTION SHELL.") << std::endl
				<< RM(" HIGHEST RESOLUTION SHELL, RANGE HIGH (A) : ", 7, 2) << Ff(reflns_shell, "d_res_high") << std::endl
				<< RM(" HIGHEST RESOLUTION SHELL, RANGE LOW  (A) : ", 7, 2) << Ff(reflns_shell, "d_res_low") << std::endl
				<< RM(" COMPLETENESS FOR SHELL     (%) : ", 7, 1) << Ff(reflns_shell, "percent_possible_all") << std::endl
				<< RM(" DATA REDUNDANCY IN SHELL       : ", 7, 2) << Ff(reflns_shell, "pdbx_redundancy") << std::endl
				<< RM(" R MERGE FOR SHELL          (I) : ", 7, 5) << Ff(reflns_shell, "Rmerge_I_obs") << std::endl
				<< RM(" R SYM FOR SHELL            (I) : ", 7, 5) << Ff(reflns_shell, "pdbx_Rsym_value") << std::endl
				<< RM(" <I/SIGMA(I)> FOR SHELL         : ", 7, 3) << Ff(reflns_shell, "meanI_over_sigI_obs") << std::endl
				<< RM("") << std::endl;

			struct { Row r; const char* field; const char* dst; }
			kTail[] = {
				{ diffrn_radiation, "pdbx_diffrn_protocol", "DIFFRACTION PROTOCOL: "},
				{ refine, "pdbx_method_to_determine_struct", "METHOD USED TO DETERMINE THE STRUCTURE: "},
				{ computing, "structure_solution", "SOFTWARE USED: "},
				{ refine, "pdbx_starting_model", "STARTING MODEL: "},
				{ exptl_crystal, "description", "\nREMARK: " }
			};

			for (auto& t: kTail)
			{
				std::string s = t.r[t.field].as<std::string>();
				
				if (s.empty())
				{
					if (strcmp(t.field, "structure_solution") == 0)
						s = cifSoftware(db, ePhasing);
					else
						s = "NULL";
				}
				
				WriteOneContinuedLine(pdbFile, "REMARK 200", 0, t.dst + s);
			}
		
			break;
		}
	}
	catch (const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}
}

void WriteRemark280(std::ostream& pdbFile, Datablock& db)
{
	typedef RM<280> RM;
	
	try
	{
		for (auto exptl_crystal: db["exptl_crystal"])
		{
			std::string crystal_id = exptl_crystal["id"].as<std::string>();
			auto exptl_crystal_grow = db["exptl_crystal_grow"][cif::Key("crystal_id") == crystal_id];

			pdbFile
				<< RM("") << std::endl
				<< RM("CRYSTAL") << std::endl
				<< RM("SOLVENT CONTENT, VS   (%): ", 6, 2) << Ff(exptl_crystal, "density_percent_sol") << std::endl
				<< RM("MATTHEWS COEFFICIENT, VM (ANGSTROMS**3/DA): ", 6, 2) << Ff(exptl_crystal, "density_Matthews") << std::endl
				<< RM("") << std::endl;

			std::vector<std::string> conditions;
			auto add = [&conditions](const std::string c)
			{
				if (find(conditions.begin(), conditions.end(), c) == conditions.end())
					conditions.push_back(c);
			};
			
			const char* keys[] = { "pdbx_details", "ph", "method", "temp" };
			
			for (size_t i = 0; i < (sizeof(keys) / sizeof(const char*)); ++i)
			{
				const char* c = keys[i];
				
				std::string v = exptl_crystal_grow[c].as<std::string>();
				if (not v.empty())
				{
					ba::to_upper(v);
					
					switch (i)
					{
						case 1:	add("PH " + v);							break;
						case 3: add("TEMPERATURE " + v + "K");			break;

						default:
							for (std::string::size_type b = 0, e = v.find(", "); b != std::string::npos; b = (e == std::string::npos ? e : e + 2), e = v.find(", ", b))
								add(v.substr(b, e - b));
							break;
					}

				}
			}

			WriteOneContinuedLine(pdbFile, "REMARK 280", 0, "CRYSTALLIZATION CONDITIONS: " + (conditions.empty() ? "NULL" : ba::join(conditions, ", ")));

			break;
		}
	}
	catch (const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}
}

void WriteRemark350(std::ostream& pdbFile, Datablock& db)
{
	auto& c1 = db["pdbx_struct_assembly"];
	if (c1.empty())
		return;

	std::vector<std::string> biomolecules, details;
	for (auto bm: c1)
	{
		std::string id = bm["id"].as<std::string>();
		biomolecules.push_back(id);
		
		for (auto r: db["struct_biol"].find(cif::Key("id") == id))
		{
			std::string s = r["details"].as<std::string>();
			if (not s.empty())
				details.push_back(s);
		}
	}
	
	// write out the mandatory REMARK 300 first
	
	pdbFile << RM<300>("") << std::endl
			<< RM<300>("BIOMOLECULE: ") << ba::join(biomolecules, ", ") << std::endl
			<< RM<300>("SEE REMARK 350 FOR THE AUTHOR PROVIDED AND/OR PROGRAM") << std::endl
			<< RM<300>("GENERATED ASSEMBLY INFORMATION FOR THE STRUCTURE IN") << std::endl
			<< RM<300>("THIS ENTRY. THE REMARK MAY ALSO PROVIDE INFORMATION ON") << std::endl
			<< RM<300>("BURIED SURFACE AREA.") << std::endl;

	if (not details.empty())
	{
		pdbFile << RM<300>("REMARK:") << std::endl;
		
		for (auto detail: details)
			WriteOneContinuedLine(pdbFile, "REMARK 300", 0, detail);
	}

	typedef RM<350> RM;
	
	pdbFile << RM("") << std::endl
			<< RM("COORDINATES FOR A COMPLETE MULTIMER REPRESENTING THE KNOWN") << std::endl
			<< RM("BIOLOGICALLY SIGNIFICANT OLIGOMERIZATION STATE OF THE") << std::endl
			<< RM("MOLECULE CAN BE GENERATED BY APPLYING BIOMT TRANSFORMATIONS") << std::endl
			<< RM("GIVEN BELOW.  BOTH NON-CRYSTALLOGRAPHIC AND") << std::endl
			<< RM("CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN.") << std::endl;

	for (auto bm: c1)
	{
		std::string id, details, method, oligomer;
		cif::tie(id, details, method, oligomer) = bm.get("id", "details", "method_details", "oligomeric_details");
		
		pdbFile << RM("") << std::endl
			 	<< RM("BIOMOLECULE: ") << id << std::endl;
		
		ba::to_upper(oligomer);
		
		if (details == "author_defined_assembly" or details == "author_and_software_defined_assembly")
			pdbFile << RM("AUTHOR DETERMINED BIOLOGICAL UNIT: ") << oligomer << std::endl;
		
		if (details == "software_defined_assembly" or details == "author_and_software_defined_assembly")
			pdbFile << RM("SOFTWARE DETERMINED QUATERNARY STRUCTURE: ") << oligomer << std::endl;
		
		if (not method.empty())
			pdbFile << RM("SOFTWARE USED: ") << method << std::endl;
		
		for (std::string type: { "ABSA (A^2)", "SSA (A^2)", "MORE" })
		{
			for (auto prop: db["pdbx_struct_assembly_prop"].find(cif::Key("biol_id") == id and cif::Key("type") == type))
			{
				std::string value = prop["value"].as<std::string>();
			
				if (cif::iequals(type, "ABSA (A^2)"))
					pdbFile << RM("TOTAL BURIED SURFACE AREA: ") << value << " ANGSTROM**2" << std::endl;
				else if (cif::iequals(type, "SSA (A^2)"))
					pdbFile << RM("SURFACE AREA OF THE COMPLEX: ") << value << " ANGSTROM**2" << std::endl;
				else if (cif::iequals(type, "MORE"))
					pdbFile << RM("CHANGE IN SOLVENT FREE ENERGY: ") << value << " KCAL/MOL" << std::endl;
			}
		}
		
		auto gen = db["pdbx_struct_assembly_gen"][cif::Key("assembly_id") == id];
		
		std::vector<std::string> asyms;
		std::string asym_id_list, oper_id_list;
		cif::tie(asym_id_list, oper_id_list) = gen.get("asym_id_list", "oper_expression");
		
		ba::split(asyms, asym_id_list, ba::is_any_of(","));
		
		std::vector<std::string> chains = MapAsymIDs2ChainIDs(asyms, db);
		pdbFile << RM("APPLY THE FOLLOWING TO CHAINS: ") << ba::join(chains, ", ") << std::endl; 
		

		for (auto i = make_split_iterator(oper_id_list, ba::token_finder(ba::is_any_of(","), ba::token_compress_on)); not i.eof(); ++i)
		{
			std::string oper_id{ i->begin(), i->end() };

			auto r = db["pdbx_struct_oper_list"][cif::Key("id") == oper_id];
			
			pdbFile << RM("  BIOMT1 ", -3) <<	Fs(r, "id")
					<< SEP(" ", -9, 6) <<		Ff(r, "matrix[1][1]")
					<< SEP(" ", -9, 6) <<		Ff(r, "matrix[1][2]")
					<< SEP(" ", -9, 6) <<		Ff(r, "matrix[1][3]")
					<< SEP(" ", -14, 5) <<	Ff(r, "vector[1]")
					<< std::endl
					<< RM("  BIOMT2 ", -3) <<	Fs(r, "id")
					<< SEP(" ", -9, 6) <<		Ff(r, "matrix[2][1]")
					<< SEP(" ", -9, 6) <<		Ff(r, "matrix[2][2]")
					<< SEP(" ", -9, 6) <<		Ff(r, "matrix[2][3]")
					<< SEP(" ", -14, 5) <<	Ff(r, "vector[2]")
					<< std::endl
					<< RM("  BIOMT3 ", -3) <<	Fs(r, "id")
					<< SEP(" ", -9, 6) <<		Ff(r, "matrix[3][1]")
					<< SEP(" ", -9, 6) <<		Ff(r, "matrix[3][2]")
					<< SEP(" ", -9, 6) <<		Ff(r, "matrix[3][3]")
					<< SEP(" ", -14, 5) <<	Ff(r, "vector[3]")
					<< std::endl;
		}
	}
}

void WriteRemark400(std::ostream& pdbFile, Datablock& db)
{
	for (auto& r: db["pdbx_entry_details"])
	{
		std::string compound_details = r["compound_details"].as<std::string>();
		if (not compound_details.empty())
			WriteOneContinuedLine(pdbFile, "REMARK 400", 0, "\nCOMPOUND\n" + compound_details);
	}
}

void WriteRemark450(std::ostream& pdbFile, Datablock& db)
{
	for (auto& r: db["pdbx_entry_details"])
	{
		std::string source_details = r["source_details"].as<std::string>();
		if (not source_details.empty())
			WriteOneContinuedLine(pdbFile, "REMARK 450", 0, "\nSOURCE\n" + source_details, 11);
		break;
	}
}

void WriteRemark465(std::ostream& pdbFile, Datablock& db)
{
	bool first = true;
	typedef RM<465> RM;
	boost::format fmt("REMARK 465 %3.3s %3.3s %1.1s %5.5d%1.1s");
	
	auto& c = db["pdbx_unobs_or_zero_occ_residues"];
	std::vector<Row> missing(c.begin(), c.end());
	stable_sort(missing.begin(), missing.end(), [](Row a, Row b) -> bool
	{
		int modelNrA, seqIDA, modelNrB, seqIDB;
		std::string asymIDA, asymIDB;
		
		cif::tie(modelNrA, asymIDA, seqIDA) = a.get("PDB_model_num", "auth_asym_id", "auth_seq_id");
		cif::tie(modelNrB, asymIDB, seqIDB) = b.get("PDB_model_num", "auth_asym_id", "auth_seq_id");
		
		int d = modelNrA - modelNrB;
		if (d == 0)
			d = asymIDA.compare(asymIDB);
		if (d == 0)
			d = seqIDA - seqIDB;
		
		return d < 0;
	});

	for (auto r: missing)
	{
		if (first)
		{
			pdbFile << RM("") << std::endl
					<< RM("MISSING RESIDUES") << std::endl
					<< RM("THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE") << std::endl
					<< RM("EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN") << std::endl
					<< RM("IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)") << std::endl
					<< RM("") << std::endl
					<< RM("  M RES C SSSEQI") << std::endl;
			first = false;
		}
		
		std::string modelNr, resName, chainID, iCode;
		int seqNr;
		
		cif::tie(modelNr, resName, chainID, iCode, seqNr) =
			r.get("PDB_model_num", "auth_comp_id", "auth_asym_id", "PDB_ins_code", "auth_seq_id");
		
		pdbFile << fmt % modelNr % resName % chainID % seqNr % iCode << std::endl;
	}
}

void WriteRemark470(std::ostream& pdbFile, Datablock& db)
{
	typedef RM<470> RM;
	boost::format fmt("REMARK 470 %3.3s %3.3s %1.1s%4.4d%1.1s  ");

	// wow...
	typedef std::tuple<std::string,std::string,int,std::string,std::string> key_type;
	std::map<key_type,std::deque<std::string>> data;

	for (auto r: db["pdbx_unobs_or_zero_occ_atoms"])
	{
		std::string modelNr, resName, chainID, iCode, atomID;
		int seqNr;
		
		cif::tie(modelNr, resName, chainID, iCode, seqNr, atomID) =
			r.get("PDB_model_num", "auth_comp_id", "auth_asym_id", "PDB_ins_code", "auth_seq_id", "auth_atom_id");

		key_type k{ modelNr, chainID, seqNr, iCode, resName };

		auto i = data.find(k);
		if (i == data.end())
			data[k] = std::deque<std::string>{ atomID };
		else
			i->second.push_back(atomID);
	}
	
	if (not data.empty())
	{
		pdbFile	<< RM("") << std::endl
				<< RM("MISSING ATOM") << std::endl
				<< RM("THE FOLLOWING RESIDUES HAVE MISSING ATOMS (M=MODEL NUMBER;") << std::endl
				<< RM("RES=RESIDUE NAME; C=CHAIN IDENTIFIER; SSEQ=SEQUENCE NUMBER;") << std::endl
				<< RM("I=INSERTION CODE):") << std::endl
				<< RM("  M RES CSSEQI  ATOMS") << std::endl;
		
		for (auto& a: data)
		{
			std::string modelNr, resName, chainID, iCode;
			int seqNr;

			tie(modelNr, chainID, seqNr, iCode, resName) = a.first;
			
			while (not a.second.empty())
			{
				pdbFile << fmt % modelNr % resName % chainID % seqNr % iCode << "  ";
				
				for (size_t i = 0; i < 6 and not a.second.empty(); ++i)
				{
					pdbFile << cif2pdbAtomName(a.second.front(), resName, db) << ' ';	
					a.second.pop_front();
				}

				pdbFile << std::endl;
			}
		}
		
	}
}

void WriteRemark610(std::ostream& pdbFile, Datablock& db)
{
	// #warning("unimplemented!");
}

void WriteRemark800(std::ostream& pdbFile, Datablock& db)
{
	int nr = 0;
	for (auto r: db["struct_site"])
	{
		pdbFile << "REMARK 800" << std::endl;
		if (++nr == 1)
		{
			pdbFile << "REMARK 800 SITE" << std::endl;
			++nr;
		}
		
		std::string ident, code, desc;
		cif::tie(ident, code, desc) = r.get("id", "pdbx_evidence_code", "details");

		ba::to_upper(code);

		for (auto l: { "SITE_IDENTIFIER: " + ident, "EVIDENCE_CODE: " + code, "SITE_DESCRIPTION: " + desc })
		{
			for (auto s: cif::wordWrap(l, 69))
				pdbFile << "REMARK 800 " << s << std::endl;
		};
	}
}

void WriteRemark999(std::ostream& pdbFile, Datablock& db)
{
	for (auto& r: db["pdbx_entry_details"])
	{
		std::string sequence_details = r["sequence_details"].as<std::string>();
		if (not sequence_details.empty())
			WriteOneContinuedLine(pdbFile, "REMARK 999", 0, "\nSEQUENCE\n" + sequence_details, 11);
		break;
	}
}

void WriteRemarks(std::ostream& pdbFile, Datablock& db)
{
	WriteRemark1(pdbFile, db);
	WriteRemark2(pdbFile, db);
	WriteRemark3(pdbFile, db);

	WriteRemark200(pdbFile, db);
	WriteRemark280(pdbFile, db);
	
	WriteRemark350(pdbFile, db);
	
	WriteRemark400(pdbFile, db);

	WriteRemark465(pdbFile, db);
	WriteRemark470(pdbFile, db);
	
	WriteRemark610(pdbFile, db);

	WriteRemark800(pdbFile, db);
	WriteRemark999(pdbFile, db);
}

int WritePrimaryStructure(std::ostream& pdbFile, Datablock& db)
{
	int numSeq = 0;
	
	// DBREF
	
	for (auto r: db["struct_ref"])
	{
		std::string id, db_name, db_code;
		cif::tie(id, db_name, db_code) = r.get("id", "db_name", "db_code");
		
		for (auto r1: db["struct_ref_seq"].find(cif::Key("ref_id") == id))
		{
			std::string idCode, chainID, insertBegin, insertEnd, dbAccession, dbinsBeg, dbinsEnd;
			std::string seqBegin, seqEnd, dbseqBegin, dbseqEnd;
			
			cif::tie(idCode, chainID, seqBegin, insertBegin, seqEnd, insertEnd, dbAccession, dbseqBegin, dbinsBeg, dbseqEnd, dbinsEnd)
				= r1.get("pdbx_PDB_id_code", "pdbx_strand_id", "pdbx_auth_seq_align_beg", "pdbx_seq_align_beg_ins_code", "pdbx_auth_seq_align_end",
					"pdbx_seq_align_end_ins_code", "pdbx_db_accession", "db_align_beg", "pdbx_db_align_beg_ins_code", "db_align_end", "pdbx_db_align_end_ins_code");

			if (dbAccession.length() > 8 or db_code.length() > 12 or atoi(dbseqEnd.c_str()) >= 100000)
				pdbFile << (boost::format(
								"DBREF1 %4.4s %1.1s %4.4s%1.1s %4.4s%1.1s %-6.6s               %-20.20s")
								% idCode
								% chainID
								% seqBegin
								% insertBegin
								% seqEnd
								% insertEnd
								% db_name
								% db_code
								).str()
						<< std::endl
						<< (boost::format(
								"DBREF2 %4.4s %1.1s     %-22.22s     %10.10s  %10.10s")
								% idCode
								% chainID
								% dbAccession
								% dbseqBegin
								% dbseqEnd).str()
						<< std::endl;
			else
				pdbFile << (boost::format(
							"DBREF  %4.4s %1.1s %4.4s%1.1s %4.4s%1.1s %-6.6s %-8.8s %-12.12s %5.5s%1.1s %5.5s%1.1s")
							% idCode
							% chainID
							% seqBegin
							% insertBegin
							% seqEnd
							% insertEnd
							% db_name
							% dbAccession
							% db_code
							% dbseqBegin
							% dbinsBeg
							% dbseqEnd
							% dbinsEnd).str() << std::endl;
		}
	}
	
	// SEQADV
	
	for (auto r: db["struct_ref_seq_dif"])
	{
		std::string idCode, resName, chainID, seqNum, iCode, database, dbAccession, dbRes, dbSeq, conflict;
		
		cif::tie(idCode, resName, chainID, seqNum, iCode, database, dbAccession, dbRes, dbSeq, conflict)
		 	= r.get("pdbx_PDB_id_code", "mon_id", "pdbx_pdb_strand_id", "pdbx_auth_seq_num", "pdbx_pdb_ins_code",
		 			"pdbx_seq_db_name", "pdbx_seq_db_accession_code", "db_mon_id", "pdbx_seq_db_seq_num",
		 			"details");
		
		ba::to_upper(conflict);
		
		pdbFile << (boost::format(
					"SEQADV %4.4s %3.3s %1.1s %4.4s%1.1s %-4.4s %-9.9s %3.3s %5.5s %-21.21s")
					% idCode
					% resName
					% chainID
					% seqNum
					% iCode
					% database
					% dbAccession
					% dbRes
					% dbSeq
					% conflict).str() << std::endl;
	}
	
	// SEQRES
	
	std::map<char,std::vector<std::string>> seqres;
	std::map<char,int> seqresl;
	for (auto r: db["pdbx_poly_seq_scheme"])
	{
		std::string chainID, res;
		cif::tie(chainID, res) = r.get("pdb_strand_id", "mon_id");
		if (chainID.empty() or res.length() > 3 or res.length() < 1)
			throw std::runtime_error("invalid pdbx_poly_seq_scheme record, chain: " + chainID + " res: " + res);
		seqres[chainID[0]].push_back(std::string(3 - res.length(), ' ') + res);
		++seqresl[chainID[0]];
	}
	
	for (auto s: seqres)
	{
		char chainID;
		std::vector<std::string> seq;
		tie(chainID, seq) = s;
		
		int n = 1;
		while (seq.empty() == false)
		{
			auto t = seq.size();
			if (t > 13)
				t = 13;
			
			auto r = boost::make_iterator_range(seq.begin(), seq.begin() + t);
			
			pdbFile << (boost::format(
						"SEQRES %3.3d %c %4.4d  %-51.51s          ")
						% n++
						% chainID
						% seqresl[chainID]
						% ba::join(r, " ")).str()
						<< std::endl;
			
			++numSeq;
			
			seq.erase(seq.begin(), seq.begin() + t);
		}
	}
	
	// MODRES
	
	for (auto r: db["pdbx_struct_mod_residue"])
	{
		std::string chainID, seqNum, resName, iCode, stdRes, comment;
		
		cif::tie(chainID, seqNum, resName, iCode, stdRes, comment) =
			r.get("auth_asym_id", "auth_seq_id", "auth_comp_id", "PDB_ins_code", "parent_comp_id", "details");
		
		pdbFile << (boost::format(
					"MODRES %4.4s %3.3s %1.1s %4.4s%1.1s %3.3s  %-41.41s")
					% db.getName()
					% resName
					% chainID
					% seqNum
					% iCode
					% stdRes
					% comment).str()
					<< std::endl;
	}
	
	return numSeq;
}

int WriteHeterogen(std::ostream& pdbFile, Datablock& db)
{
	int numHet = 0;
	
	std::string water_entity_id, water_comp_id;
	for (auto r: db["entity"].find(cif::Key("type") == std::string("water")))
	{
		water_entity_id = r["id"].as<std::string>();
		break;
	}
	
	std::map<std::string,std::string> het;
	
	for (auto r: db["chem_comp"])
	{
		std::string id, name, mon_nstd_flag;
		cif::tie(id, name, mon_nstd_flag) = r.get("id", "name", "mon_nstd_flag");
		
		if (mon_nstd_flag == "y")
			continue;
		
		het[id] = name;
	}

	for (auto r: db["pdbx_entity_nonpoly"])
	{
		std::string entity_id, name, comp_id;
		cif::tie(entity_id, name, comp_id) = r.get("entity_id", "name", "comp_id");
		
		if (entity_id == water_entity_id)
			water_comp_id = comp_id;

		if (het.count(comp_id) == 0)
			het[comp_id] = name;
	}
	
	struct HET {
		bool	water;
		std::string	hetID;
		char	chainID;
		int		seqNum;
		char	iCode;
		int		numHetAtoms;
		std::string	text;	// ignored
	};
	std::vector<HET> hets;
	
//	// construct component number map
//	map<int,int> component_nr;
//	std::string lChainID, lCompID, lICode;
//	int lSeqNum;
//	
//	for (auto r: db["atom_site"])
//	{
//		std::string chainID, compID, iCode;
//		int seqNum;
//		
//		cif::tie(seqNum, comp_id, chain_id, iCode) =
//			r.get("auth_seq_id", "auth_comp_id", "auth_asym_id", "pdbx_PDB_ins_code");
//		
//		if (chainID != lChainID or compID != lCompID or seqNum != lSeqNum or iCode != lICode)
//			
//	}
	
	// count the HETATM's
//	for (auto r: db["atom_site"].find(cif::Key("group_PDB") == std::string("HETATM")))
	std::set<std::string> missingHetNames;

	for (auto r: db["atom_site"])
	{
		int seqNum;
		std::string entity_id, comp_id, chain_id, iCode, modelNr;
		
		cif::tie(entity_id, seqNum, comp_id, chain_id, iCode, modelNr) =
			r.get("label_entity_id", "auth_seq_id", "auth_comp_id", "auth_asym_id", "pdbx_PDB_ins_code", "pdbx_PDB_model_num");

		if (mmcif::kAAMap.count(comp_id) or mmcif::kBaseMap.count(comp_id))
			continue;
		
		if (chain_id.length() != 1)
			throw std::runtime_error("Cannot produce PDB file, auth_asym_id not valid");
		
		if (entity_id != water_entity_id and het.count(comp_id) == 0)
			missingHetNames.insert(comp_id); 
		
		auto h = find_if(hets.begin(), hets.end(),
			[=](const HET& het) -> bool
			{
				return het.hetID == comp_id and het.chainID == chain_id[0] and het.seqNum == seqNum; 
			});
		
		if (h == hets.end())
		{
			hets.push_back({entity_id == water_entity_id, comp_id, chain_id[0], seqNum,
				(iCode.empty() ? ' ' : iCode[0]), 1
			});

		}
		else
			h->numHetAtoms += 1;
	}
	
	if (cif::VERBOSE > 1 and not missingHetNames.empty())
		std::cerr << "Missing het name(s) for " << ba::join(missingHetNames, ", ") << std::endl;
	
	boost::format kHET("HET    %3.3s  %1.1s%4.4d%1.1s  %5.5d");
	for (auto h: hets)
	{
		if (h.water)
			continue;
		pdbFile << (kHET % h.hetID % h.chainID % h.seqNum % h.iCode % h.numHetAtoms) << std::endl;
		++numHet;
	}
	
	for (auto n: het)
	{
		std::string id, name;
		tie(id, name) = n;

		if (id == water_comp_id)
			continue;

		ba::to_upper(name);
		
		int c = 1;
		
		boost::format kHETNAM("HETNAM  %2.2s %3.3s ");
		for (;;)
		{
			pdbFile << (kHETNAM % (c > 1 ? std::to_string(c) : std::string()) % id);
			++c;
			
			if (name.length() > 55)
			{
				bool done = false;
				for (auto e = name.begin() + 54; e != name.begin(); --e)
				{
					if (ispunct(*e))
					{
						pdbFile << std::string(name.begin(), e) << std::endl;
						name.erase(name.begin(), e);
						done = true;
						break;
					}
				}
				
				if (not done)
				{
					pdbFile << std::string(name.begin(), name.begin() + 55) << std::endl;
					name.erase(name.begin(), name.begin() + 55);
				}

				continue;
			}
			
			pdbFile << name << std::endl;
			break;
		}
	}
	
	for (auto n: het)
	{
		std::string id, name;
		tie(id, name) = n;

		if (id == water_comp_id)
			continue;

		std::string syn = db["chem_comp"][cif::Key("id") == id]["pdbx_synonyms"].as<std::string>();
		if (syn.empty())
			continue;
		
		WriteOneContinuedLine(pdbFile, "HETSYN", 4, id + ' ' + syn, 11);
	}
	
	// FORMUL

	std::vector<std::string> formulas;

	boost::format kFORMUL("FORMUL  %2.2d  %3.3s %2.2s%c");
	for (auto h: het)
	{
		std::string hetID = h.first;
		int componentNr = 0;
		
		std::string first_het_asym_id;
		for (auto p: db["pdbx_poly_seq_scheme"].find(cif::Key("mon_id") == hetID))
		{
			first_het_asym_id = p["asym_id"].as<std::string>();
			break;
		}
		
		if (first_het_asym_id.empty())
		{
			for (auto p: db["pdbx_nonpoly_scheme"].find(cif::Key("mon_id") == hetID))
			{
				first_het_asym_id = p["asym_id"].as<std::string>();
				break;
			}
		}
		
		if (not first_het_asym_id.empty())
		{
			for (auto a: db["struct_asym"])
			{
				++componentNr;
				if (a["id"] == first_het_asym_id)
					break;
			}
		}
		
		int nr = count_if(hets.begin(), hets.end(), [hetID](auto& h) -> bool { return h.hetID == hetID; });
		
		for (auto r: db["chem_comp"].find(cif::Key("id") == hetID))
		{
			std::string formula = r["formula"].as<std::string>();
			if (nr > 1)
				formula = std::to_string(nr) + '(' + formula + ')';
			
			int c = 1;
			for (;;)
			{
				std::stringstream fs;
				
				fs << (kFORMUL % componentNr % hetID % (c > 1 ? std::to_string(c) : std::string()) % (hetID == water_comp_id ? '*' : ' '));
				++c;
				
				if (formula.length() > 51)
				{
					bool done = false;
					for (auto e = formula.begin() + 50; e != formula.begin(); --e)
					{
						if (ispunct(*e))
						{
							pdbFile << std::string(formula.begin(), e) << std::endl;
							formula.erase(formula.begin(), e);
							done = true;
							break;
						}
					}
					
					if (not done)
					{
						pdbFile << std::string(formula.begin(), formula.begin() + 55) << std::endl;
						formula.erase(formula.begin(), formula.begin() + 55);
					}
	
					continue;
				}
				
				fs << formula << std::endl;
				
				formulas.push_back(fs.str());
				break;
			}

			break;
		}
	}
	
	sort(formulas.begin(), formulas.end(), [](const std::string& a, const std::string& b) -> bool
		{ return stoi(a.substr(8, 2)) < stoi(b.substr(8, 2)); });
	
	for (auto& f: formulas)
		pdbFile << f;
	
	return numHet;
}

std::tuple<int,int> WriteSecondaryStructure(std::ostream& pdbFile, Datablock& db)
{
	int numHelix = 0, numSheet = 0;
	
	// HELIX
	boost::format kHELIX("HELIX  %3.3d %3.3s %3.3s %c %4.4d%1.1s %3.3s %c %4.4d%1.1s%2.2d%-30.30s %5.5d");
	for (auto r: db["struct_conf"].find(cif::Key("conf_type_id") == "HELX_P"))
	{
		std::string pdbx_PDB_helix_id, beg_label_comp_id, pdbx_beg_PDB_ins_code,
			end_label_comp_id, pdbx_end_PDB_ins_code, beg_auth_comp_id,
			beg_auth_asym_id, end_auth_comp_id, end_auth_asym_id, details;
		int pdbx_PDB_helix_class, pdbx_PDB_helix_length, beg_auth_seq_id, end_auth_seq_id;
		
		cif::tie(pdbx_PDB_helix_id, beg_label_comp_id, pdbx_beg_PDB_ins_code,
			end_label_comp_id, pdbx_end_PDB_ins_code, beg_auth_comp_id,
			beg_auth_asym_id, end_auth_comp_id, end_auth_asym_id, details,
			pdbx_PDB_helix_class, pdbx_PDB_helix_length, beg_auth_seq_id, end_auth_seq_id) = 
			r.get("pdbx_PDB_helix_id", "beg_label_comp_id", "pdbx_beg_PDB_ins_code",
			"end_label_comp_id", "pdbx_end_PDB_ins_code", "beg_auth_comp_id",
			"beg_auth_asym_id", "end_auth_comp_id", "end_auth_asym_id", "details",
			"pdbx_PDB_helix_class", "pdbx_PDB_helix_length", "beg_auth_seq_id", "end_auth_seq_id");
		
		++numHelix;
		pdbFile << (kHELIX % numHelix % pdbx_PDB_helix_id
			% beg_label_comp_id % beg_auth_asym_id % beg_auth_seq_id % pdbx_beg_PDB_ins_code
			% end_label_comp_id % end_auth_asym_id % end_auth_seq_id % pdbx_end_PDB_ins_code
			% pdbx_PDB_helix_class % details % pdbx_PDB_helix_length) << std::endl;
	}
	
	// SHEET
	boost::format
				kSHEET0("SHEET  %3.3d %3.3s%2.2d %3.3s %1.1s%4.4d%1.1s %3.3s %1.1s%4.4d%1.1s%2.2d"),
				kSHEET1("SHEET  %3.3d %3.3s%2.2d %3.3s %1.1s%4.4d%1.1s %3.3s %1.1s%4.4d%1.1s%2.2d "
						"%-4.4s%3.3s %1.1s%4.4d%1.1s %-4.4s%3.3s %1.1s%4.4d%1.1s");

	for (auto r: db["struct_sheet"])
	{
		std::string sheetID;
		int numStrands = 0;
		
		cif::tie(sheetID, numStrands) = r.get("id", "number_strands");

		bool first = true;
		
		for (auto o: db["struct_sheet_order"].find(cif::Key("sheet_id") == sheetID))
		{
			int sense = 0;
			std::string s, rangeID1, rangeID2;

			cif::tie(s, rangeID1, rangeID2) = o.get("sense", "range_id_1", "range_id_2");
			if (s == "anti-parallel")
				sense = -1;
			else if (s == "parallel")
				sense = 1;
			
			if (first)
			{
				std::string initResName, initChainID, initICode, endResName, endChainID, endICode;
				int initSeqNum, endSeqNum;

				auto r1 = db["struct_sheet_range"][cif::Key("sheet_id") == sheetID and cif::Key("id") == rangeID1];
	
				cif::tie(initResName, initICode, endResName, endICode,
					initResName, initChainID, initSeqNum, endResName, endChainID, endSeqNum)
					 = r1.get("beg_label_comp_id", "pdbx_beg_PDB_ins_code", "end_label_comp_id",
				 	"pdbx_end_PDB_ins_code", "beg_auth_comp_id", "beg_auth_asym_id", "beg_auth_seq_id",
				 	"end_auth_comp_id", "end_auth_asym_id", "end_auth_seq_id");
			
				pdbFile << (kSHEET0 
					% rangeID1
					% sheetID
					% numStrands
					% initResName
					% initChainID
					% initSeqNum
					% initICode
					% endResName
					% endChainID
					% endSeqNum
					% endICode
					% 0) << std::endl;

				first = false;
			}
			
			std::string initResName, initChainID, initICode, endResName, endChainID, endICode, curAtom, curResName, curChainID, curICode, prevAtom, prevResName, prevChainID, prevICode;
			int initSeqNum, endSeqNum, curResSeq, prevResSeq;

			auto r2 = db["struct_sheet_range"][cif::Key("sheet_id") == sheetID and cif::Key("id") == rangeID2];

			cif::tie(initResName, initICode, endResName, endICode,
				initResName, initChainID, initSeqNum, endResName, endChainID, endSeqNum)
				 = r2.get("beg_label_comp_id", "pdbx_beg_PDB_ins_code", "end_label_comp_id",
			 	"pdbx_end_PDB_ins_code", "beg_auth_comp_id", "beg_auth_asym_id", "beg_auth_seq_id",
			 	"end_auth_comp_id", "end_auth_asym_id", "end_auth_seq_id");
			
			auto h = db["pdbx_struct_sheet_hbond"].find(cif::Key("sheet_id") == sheetID and cif::Key("range_id_1") == rangeID1 and cif::Key("range_id_2") == rangeID2);
			
			if (h.empty())
			{
				pdbFile << (kSHEET0
					% rangeID2
					% sheetID
					% numStrands
					% initResName
					% initChainID
					% initSeqNum
					% initICode
					% endResName
					% endChainID
					% endSeqNum
					% endICode
					% sense) << std::endl;
			}
			else
			{
				std::string compID[2];
				cif::tie(compID[0], compID[1]) = h.front().get("range_2_label_comp_id", "range_1_label_comp_id");
				
				cif::tie(curAtom, curResName, curResSeq, curChainID, curICode, prevAtom, prevResName, prevResSeq, prevChainID, prevICode)
					= h.front().get("range_2_auth_atom_id", "range_2_auth_comp_id", "range_2_auth_seq_id", "range_2_auth_asym_id", "range_2_PDB_ins_code",
						"range_1_auth_atom_id", "range_1_auth_comp_id", "range_1_auth_seq_id", "range_1_auth_asym_id", "range_1_PDB_ins_code");
				
				curAtom = cif2pdbAtomName(curAtom, compID[0], db);
				prevAtom = cif2pdbAtomName(prevAtom, compID[1], db);
				
				pdbFile << (kSHEET1
					% rangeID2
					% sheetID
					% numStrands
					% initResName
					% initChainID
					% initSeqNum
					% initICode
					% endResName
					% endChainID
					% endSeqNum
					% endICode
					% sense
					% curAtom
					% curResName
					% curChainID
					% curResSeq
					% curICode
					% prevAtom
					% prevResName
					% prevChainID
					% prevResSeq
					% prevICode) << std::endl;
			}
			
			++numSheet;
		}
	}
	
	return std::make_tuple(numHelix, numSheet);
}

void WriteConnectivity(std::ostream& pdbFile, cif::Datablock& db)
{
	// SSBOND
	// have to filter out alts
	std::set<std::tuple<char,int,char,char,int,char>> ssSeen;
	
	int nr = 1;
	boost::format kSSBOND("SSBOND %3.3d CYS %1.1s %4.4d%1.1s   CYS %1.1s %4.4d%1.1s                       %6.6s %6.6s %5.2f");
	for (auto r: db["struct_conn"].find(cif::Key("conn_type_id") == "disulf"))
	{
		std::string chainID1, icode1, chainID2, icode2, sym1, sym2;
		float Length;
		int seqNum1, seqNum2;
		
		cif::tie(
			chainID1, seqNum1, icode1, chainID2, seqNum2, icode2, sym1, sym2, Length) = 
			r.get("ptnr1_auth_asym_id", "ptnr1_auth_seq_id", "pdbx_ptnr1_PDB_ins_code",
				  "ptnr2_auth_asym_id", "ptnr2_auth_seq_id", "pdbx_ptnr2_PDB_ins_code",
				  "ptnr1_symmetry", "ptnr2_symmetry", "pdbx_dist_value");

		auto n = ssSeen.emplace(chainID1[0], seqNum1, icode1[0], chainID2[0], seqNum2, icode2[0]);
		if (n.second == false)
			continue;
		
		sym1 = cif2pdbSymmetry(sym1);
		sym2 = cif2pdbSymmetry(sym2);
		
		pdbFile << (kSSBOND
			% nr
			% chainID1
			% seqNum1
			% icode1
			% chainID2
			% seqNum2
			% icode2
			% sym1
			% sym2
			% Length) << std::endl;
		
		++nr;
	}
	
	// LINK
	
	boost::format kLINK("LINK        %-4.4s%1.1s%3.3s %1.1s%4.4d%1.1s               %-4.4s%1.1s%3.3s %1.1s%4.4d%1.1s  %6.6s %6.6s");
	for (auto r: db["struct_conn"].find(cif::Key("conn_type_id") == "metalc" or cif::Key("conn_type_id") == "covale"))
	{
		std::string name1, altLoc1, resName1, chainID1, iCode1, name2, altLoc2, resName2, chainID2, iCode2, sym1, sym2, Length;
		int resSeq1, resSeq2;
		
		cif::tie(name1, altLoc1, resName1, chainID1, resSeq1, iCode1, name2, altLoc2, resName2, chainID2, resSeq2, iCode2, sym1, sym2, Length) = 
			r.get("ptnr1_label_atom_id", "pdbx_ptnr1_label_alt_id", "ptnr1_label_comp_id", "ptnr1_auth_asym_id", "ptnr1_auth_seq_id", "pdbx_ptnr1_PDB_ins_code",
				  "ptnr2_label_atom_id", "pdbx_ptnr2_label_alt_id", "ptnr2_label_comp_id", "ptnr2_auth_asym_id", "ptnr2_auth_seq_id", "pdbx_ptnr2_PDB_ins_code",
				  "ptnr1_symmetry", "ptnr2_symmetry", "pdbx_dist_value");
		
		std::string compID[2];
		
		cif::tie(compID[0], compID[1]) = r.get("ptnr1_label_comp_id", "ptnr2_label_comp_id");

		name1 = cif2pdbAtomName(name1, compID[0], db);
		name2 = cif2pdbAtomName(name2, compID[1], db);

		sym1 = cif2pdbSymmetry(sym1);
		sym2 = cif2pdbSymmetry(sym2);
		
		pdbFile << (kLINK
			% name1
			% altLoc1
			% resName1
			% chainID1
			% resSeq1
			% iCode1
			% name2
			% altLoc2
			% resName2
			% chainID2
			% resSeq2
			% iCode2
			% sym1
			% sym2);
		
		if (not Length.empty())
			pdbFile << boost::format(" %5.2f") % stod(Length);
			
		pdbFile << std::endl;
	}

	// CISPEP
	
	boost::format kCISPEP("CISPEP %3.3d %3.3s %1.1s %4.4d%1.1s   %3.3s %1.1s %4.4d%1.1s       %3.3d       %6.2f");
	for (auto r: db["struct_mon_prot_cis"])
	{
		std::string serNum, pep1, chainID1, icode1, pep2, chainID2, icode2, modNum;
		int seqNum1, seqNum2;
		float measure;
		
		cif::tie(serNum, pep1, chainID1, seqNum1, icode1, pep2, chainID2, seqNum2, icode2, modNum, measure) =
			r.get("pdbx_id", "label_comp_id", "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code",
				  "pdbx_label_comp_id_2", "pdbx_auth_asym_id_2", "pdbx_auth_seq_id_2", "pdbx_PDB_ins_code_2",
				  "pdbx_PDB_model_num", "pdbx_omega_angle");
		
		pdbFile << (kCISPEP
			% serNum
			% pep1
			% chainID1
			% seqNum1
			% icode1
			% pep2
			% chainID2
			% seqNum2
			% icode2
			% modNum
			% measure) << std::endl;
	}
}

int WriteMiscellaneousFeatures(std::ostream& pdbFile, Datablock& db)
{
	int numSite = 0;
	
	// SITE
	
	std::map<std::string,std::deque<std::string>> sites;
	
	boost::format kSITERES("%3.3s %1.1s%4.4d%1.1s ");
	for (auto r: db["struct_site_gen"])
	{
		std::string siteID, resName, chainID, iCode;
		int seq;
		
		cif::tie(siteID, resName, chainID, seq, iCode) =
			r.get("site_id", "auth_comp_id", "auth_asym_id", "auth_seq_id", "pdbx_auth_ins_code");

		sites[siteID].push_back((kSITERES % resName % chainID % seq % iCode).str());
	}
	
	boost::format kSITE("SITE   %3.3d %3.3s %2.2d ");
	for (auto s: sites)
	{
		std::string siteID = std::get<0>(s);
		std::deque<std::string>& res = std::get<1>(s);
		
		int numRes = res.size();
		
		int nr = 1;
		while (res.empty() == false)
		{
			pdbFile << (kSITE % nr % siteID % numRes);
			
			for (int i = 0; i < 4; ++i)
			{
				if (not res.empty())
				{
					pdbFile << res.front();
					res.pop_front();
				}
				else
					pdbFile << std::string(11, ' ');
			}

			pdbFile << std::endl;
			++nr;
			++numSite;
		}
	}
	
	return numSite;
}

void WriteCrystallographic(std::ostream& pdbFile, Datablock& db)
{
	auto r = db["symmetry"][cif::Key("entry_id") == db.getName()];
	std::string symmetry = r["space_group_name_H-M"].as<std::string>();

	r = db["cell"][cif::Key("entry_id") == db.getName()];
	
	boost::format kCRYST1("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11.11s%4.4d");
	
	pdbFile << (kCRYST1
		% r["length_a"].as<double>()
		% r["length_b"].as<double>()
		% r["length_c"].as<double>()
		% r["angle_alpha"].as<double>()
		% r["angle_beta"].as<double>()
		% r["angle_gamma"].as<double>()
		% symmetry
		% r["Z_PDB"].as<double>()) << std::endl;
}

int WriteCoordinateTransformation(std::ostream& pdbFile, Datablock& db)
{
	int result = 0;
	
	for (auto r: db["database_PDB_matrix"])
	{
		boost::format kORIGX("ORIGX%1.1d    %10.6f%10.6f%10.6f     %10.5f");
		pdbFile << (kORIGX % 1 % r["origx[1][1]"].as<float>() % r["origx[1][2]"].as<float>() % r["origx[1][3]"].as<float>() % r["origx_vector[1]"].as<float>()) << std::endl;
		pdbFile << (kORIGX % 2 % r["origx[2][1]"].as<float>() % r["origx[2][2]"].as<float>() % r["origx[2][3]"].as<float>() % r["origx_vector[2]"].as<float>()) << std::endl;
		pdbFile << (kORIGX % 3 % r["origx[3][1]"].as<float>() % r["origx[3][2]"].as<float>() % r["origx[3][3]"].as<float>() % r["origx_vector[3]"].as<float>()) << std::endl;
		result += 3;
		break;
	}

	for (auto r: db["atom_sites"])
	{
		boost::format kSCALE("SCALE%1.1d    %10.6f%10.6f%10.6f     %10.5f");
		pdbFile << (kSCALE % 1 % r["fract_transf_matrix[1][1]"].as<float>() % r["fract_transf_matrix[1][2]"].as<float>() % r["fract_transf_matrix[1][3]"].as<float>() % r["fract_transf_vector[1]"].as<float>()) << std::endl;
		pdbFile << (kSCALE % 2 % r["fract_transf_matrix[2][1]"].as<float>() % r["fract_transf_matrix[2][2]"].as<float>() % r["fract_transf_matrix[2][3]"].as<float>() % r["fract_transf_vector[2]"].as<float>()) << std::endl;
		pdbFile << (kSCALE % 3 % r["fract_transf_matrix[3][1]"].as<float>() % r["fract_transf_matrix[3][2]"].as<float>() % r["fract_transf_matrix[3][3]"].as<float>() % r["fract_transf_vector[3]"].as<float>()) << std::endl;
		result += 3;
		break;
	}

	int nr = 1;
	boost::format kMTRIX("MTRIX%1.1d %3.3d%10.6f%10.6f%10.6f     %10.5f    %1.1s");
	for (auto r: db["struct_ncs_oper"])
	{
		std::string given = r["code"] == "given" ? "1" : "";
		
		pdbFile << (kMTRIX % 1 % nr % r["matrix[1][1]"].as<float>() % r["matrix[1][2]"].as<float>() % r["matrix[1][3]"].as<float>() % r["vector[1]"].as<float>() % given) << std::endl;
		pdbFile << (kMTRIX % 2 % nr % r["matrix[2][1]"].as<float>() % r["matrix[2][2]"].as<float>() % r["matrix[2][3]"].as<float>() % r["vector[2]"].as<float>() % given) << std::endl;
		pdbFile << (kMTRIX % 3 % nr % r["matrix[3][1]"].as<float>() % r["matrix[3][2]"].as<float>() % r["matrix[3][3]"].as<float>() % r["vector[3]"].as<float>() % given) << std::endl;
		
		++nr;
		result += 3;
	}
	
	return result;
}

std::tuple<int,int> WriteCoordinatesForModel(std::ostream& pdbFile, Datablock& db,
	const std::map<std::string,std::tuple<std::string,int,std::string>>& last_resseq_for_chain_map,
	std::set<std::string>& TERminatedChains, int model_nr)
{
	int numCoord = 0, numTer = 0;
	
	boost::format kATOM(	"%-6.6s%5.5d %-4.4s%1.1s%3.3s %1.1s%4.4d%1.1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s%2.2s"),
		kANISOU(			"ANISOU%5.5d %-4.4s%1.1s%3.3s %1.1s%4.4d%1.1s %7.7d%7.7d%7.7d%7.7d%7.7d%7.7d      %2.2s%2.2s"),
		kTER(				"TER   %5.5d      %3.3s %1.1s%4.4d%1.1s");
	
	auto& atom_site = db["atom_site"];
	auto& atom_site_anisotrop = db["atom_site_anisotrop"];
	
	int serial = 1;
	auto ri = atom_site.begin();
	
	std::string id, group, name, altLoc, resName, chainID, iCode, element;
	int resSeq, charge;
	
	for (;;)
	{
		std::string nextResName, nextChainID, nextICode, modelNum;
		int nextResSeq = 0;
		
		if (ri != atom_site.end())
			cif::tie(nextResName, nextChainID, nextICode, nextResSeq, modelNum) =
				ri->get("label_comp_id", "auth_asym_id", "pdbx_PDB_ins_code", "auth_seq_id", "pdbx_PDB_model_num");
		
		if (modelNum.empty() == false and stol(modelNum) != model_nr)
		{
			++ri;
			continue;
		}

		if (chainID.empty() == false and TERminatedChains.count(chainID) == 0)
		{
			bool terminate = nextChainID != chainID;

			if (not terminate)
				terminate =
					(nextResSeq != resSeq or iCode != nextICode) and
					(last_resseq_for_chain_map.count(chainID) == false or last_resseq_for_chain_map.at(chainID) == make_tuple(resName, resSeq, iCode));

			if (terminate)
			{
				pdbFile << (kTER
					% serial
					% resName
					% chainID
					% resSeq
					% iCode) << std::endl;
	
				++serial;
				TERminatedChains.insert(chainID);
				
				++numTer;
			}
		}
		
		if (ri == atom_site.end())
			break;
		
		auto r = *ri++;

		try
		{
			if (r["pdbx_PDB_model_num"].as<int>() != model_nr)
				continue;
		}
		catch (...) { /* perhaps no model number here */ }

		float x, y, z, occupancy, tempFactor;
		
		cif::tie(id, group, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, element, charge) =
			r.get("id", "group_PDB", "label_atom_id", "label_alt_id", "auth_comp_id", "auth_asym_id", "auth_seq_id",
				"pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "type_symbol", "pdbx_formal_charge");
		
		if (name.length() < 4 and (element.length() == 1 or std::toupper(name[0]) != std::toupper(element[0]) or std::toupper(name[1]) != std::toupper(element[1])))
			name.insert(name.begin(), ' ');
		
		std::string sCharge;
		if (charge != 0)
			sCharge = std::to_string(charge) + (charge > 0 ? '+' : '-');

		pdbFile << (kATOM
			% group
			% serial
			% name
			% altLoc
			% resName
			% chainID
			% resSeq
			% iCode
			% x
			% y
			% z
			% occupancy
			% tempFactor
			% element
			% sCharge) << std::endl;
		
		++numCoord;
		
		auto ai = atom_site_anisotrop[cif::Key("id") == id];
		if (not ai.empty())
//		
//		auto ai = find_if(atom_site_anisotrop.begin(), atom_site_anisotrop.end(), [id](Row r) -> bool { return r["id"] == id; });
//		if (ai != atom_site_anisotrop.end())
		{
			float u11, u22, u33, u12, u13, u23;
			
			cif::tie(u11, u22, u33, u12, u13, u23) =
				ai.get("U[1][1]", "U[2][2]", "U[3][3]", "U[1][2]", "U[1][3]", "U[2][3]");
			
			pdbFile << (kANISOU
				% serial
				% name
				% altLoc
				% resName
				% chainID
				% resSeq
				% iCode
				% lrintf(u11 * 10000)
				% lrintf(u22 * 10000)
				% lrintf(u33 * 10000)
				% lrintf(u12 * 10000)
				% lrintf(u13 * 10000)
				% lrintf(u23 * 10000)
				% element
				% sCharge) << std::endl;
		}
		
		++serial;
	}
	
	return std::make_tuple(numCoord, numTer);
}

std::tuple<int,int> WriteCoordinate(std::ostream& pdbFile, Datablock& db)
{
	// residues known from seqres
//	map<tuple<std::string,int,std::string>,std::string> res2chain_map;
	std::map<std::string,std::tuple<std::string,int,std::string>> last_resseq_for_chain_map;
	
	for (auto r: db["pdbx_poly_seq_scheme"])
	{
		std::string chainID, resName, iCode;
		int resSeq;
		
		if (r["auth_seq_num"].empty())
			continue;
		
		cif::tie(chainID, resName, resSeq, iCode) = r.get("pdb_strand_id", "pdb_mon_id", "auth_seq_num", "pdb_ins_code");
		
		last_resseq_for_chain_map[chainID] = make_tuple(resName, resSeq, iCode);
//		res2chain_map[make_tuple(resName, resSeq, iCode)] = chainID;
	}
	
	// collect known model numbers
	std::set<int> models;
	try
	{
		for (auto r: db["atom_site"])
			models.insert(r["pdbx_PDB_model_num"].as<int>());
	}
	catch (...)
	{
	}

	std::tuple<int,int> result;

	if (models.empty() or models == std::set<int>{ 0 })
	{
		std::set<std::string> TERminatedChains;
		result = WriteCoordinatesForModel(pdbFile, db, last_resseq_for_chain_map, TERminatedChains, 0);
	}
	else
	{
		boost::format kModel("MODEL     %4.4d");
		
		for (int model_nr: models)
		{
			if (models.size() > 1)
				pdbFile << (kModel % model_nr) << std::endl;
			
			std::set<std::string> TERminatedChains;
			auto n = WriteCoordinatesForModel(pdbFile, db, last_resseq_for_chain_map, TERminatedChains, model_nr);
			if (model_nr == 1)
				result = n;
			
			if (models.size() > 1)
				pdbFile << "ENDMDL" << std::endl;
		}
	}
	
	return result;
}

void WritePDBFile(std::ostream& pdbFile, cif::File& cifFile)
{
	io::filtering_ostream out;
	out.push(FillOutLineFilter());
	out.push(pdbFile);

	auto filter = out.component<FillOutLineFilter>(0);
	assert(filter);

	auto& db = cifFile.firstDatablock();
	
	int numRemark = 0, numHet = 0, numHelix = 0, numSheet = 0, numTurn = 0, numSite = 0, numXform = 0, numCoord = 0, numTer = 0, numConect = 0, numSeq = 0;
	
								WriteTitle(out, db);
	
	int savedLineCount = filter->GetLineCount();
//	numRemark = 				WriteRemarks(out, db);
								WriteRemarks(out, db);
	numRemark = filter->GetLineCount() - savedLineCount;

	numSeq = 					WritePrimaryStructure(out, db);
	numHet = 					WriteHeterogen(out, db);
	std::tie(numHelix, numSheet) =	WriteSecondaryStructure(out, db);
								WriteConnectivity(out, db);
	numSite =					WriteMiscellaneousFeatures(out, db);
								WriteCrystallographic(out, db);
	numXform =					WriteCoordinateTransformation(out, db);
	std::tie(numCoord, numTer) =		WriteCoordinate(out, db);

	boost::format kMASTER("MASTER    %5.5d    0%5.5d%5.5d%5.5d%5.5d%5.5d%5.5d%5.5d%5.5d%5.5d%5.5d");
	
	out	<< (kMASTER % numRemark % numHet % numHelix % numSheet % numTurn % numSite % numXform % numCoord % numTer % numConect % numSeq) << std::endl
		<< "END" << std::endl;
}

void WritePDBHeaderLines(std::ostream& os, cif::File& cifFile)
{
	// io::filtering_ostream out;
	// out.push(FillOutLineFilter());
	// out.push(os);

	// auto filter = out.component<FillOutLineFilter>(0);
	// assert(filter);

	auto& db = cifFile.firstDatablock();

	WriteHeaderLines(os, db);
}

std::string FixStringLength(const std::string& s, std::string::size_type l)
{
	auto result = s;

	if (result.length() > l)
		result = result.substr(0, l - 4) + "... ";
	else if (result.length() < l)
		result.append(l - result.length(), ' ');

	return result;
}

std::string GetPDBHEADERLine(cif::File& cifFile, std::string::size_type truncate_at)
{
	auto& db = cifFile.firstDatablock();

	//    0         1         2         3         4         5         6         7         8
	//    HEADER    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxDDDDDDDDD   IIII
	const char kHeader[] =
		 "HEADER    %1$-40.40s"                           "%2$-9.9s""   %3$-4.4s";

	// HEADER
	
	std::string keywords;
	auto& cat1 = db["struct_keywords"];
	
	for (auto r: cat1)
	{
		keywords = r["pdbx_keywords"].as<std::string>();
		if (keywords.length() > truncate_at - 40)
			keywords = keywords.substr(0, truncate_at - 44) + " ...";
	}
	
	std::string date;
	for (auto r: db["pdbx_database_status"])
	{
		date = r["recvd_initial_deposition_date"].as<std::string>();
		if (date.empty())
			continue;
		date = cif2pdbDate(date);
		break;
	}
	
	if (date.empty())
	{
		for (auto r: db["database_PDB_rev"])
		{
			date = r["date_original"].as<std::string>();
			if (date.empty())
				continue;
			date = cif2pdbDate(date);
			break;
		}
	}
	
	return FixStringLength((boost::format(kHeader) % keywords % date % db.getName()).str(), truncate_at);
}

std::string GetPDBCOMPNDLine(cif::File& cifFile, std::string::size_type truncate_at)
{
	auto& db = cifFile.firstDatablock();

	// COMPND
	using namespace std::placeholders;
	
	int molID = 0;
	std::vector<std::string> cmpnd;
	
	for (auto r: db["entity"])
	{
		if (r["type"] != "polymer")
			continue;
		
		std::string entityID = r["id"].as<std::string>();
		
		++molID;
		cmpnd.push_back("MOL_ID: " + std::to_string(molID));
		
		std::string molecule = r["pdbx_description"].as<std::string>();
		cmpnd.push_back("MOLECULE: " + molecule);

		auto poly = db["entity_poly"].find(cif::Key("entity_id") == entityID);
		if (not poly.empty())
		{
			std::string chains = poly.front()["pdbx_strand_id"].as<std::string>();
			ba::replace_all(chains, ",", ", ");
			cmpnd.push_back("CHAIN: " + chains);
		}

		std::string fragment = r["pdbx_fragment"].as<std::string>();
		if (not fragment.empty())
			cmpnd.push_back("FRAGMENT: " + fragment);
		
		for (auto sr: db["entity_name_com"].find(cif::Key("entity_id") == entityID))
		{
			std::string syn = sr["name"].as<std::string>();
			if (not syn.empty())
				cmpnd.push_back("SYNONYM: " + syn);
		}
		
		std::string mutation = r["pdbx_mutation"].as<std::string>();
		if (not mutation.empty())
			cmpnd.push_back("MUTATION: " + mutation);
		
		std::string ec = r["pdbx_ec"].as<std::string>();
		if (not ec.empty())
			cmpnd.push_back("EC: " + ec);
		
		if (r["src_method"] == "man" or r["src_method"] == "syn")
			cmpnd.push_back("ENGINEERED: YES");

		std::string details = r["details"].as<std::string>();
		if (not details.empty())
			cmpnd.push_back("OTHER_DETAILS: " + details);
	}

	return FixStringLength("COMPND    " + ba::join(cmpnd, "; "), truncate_at);
}

std::string GetPDBSOURCELine(cif::File& cifFile, std::string::size_type truncate_at)
{
	auto& db = cifFile.firstDatablock();

	// SOURCE
	
	int molID = 0;
	std::vector<std::string> source;
	
	for (auto r: db["entity"])
	{
		if (r["type"] != "polymer")
			continue;
		
		std::string entityID = r["id"].as<std::string>();
		
		++molID;
		source.push_back("MOL_ID: " + std::to_string(molID));
		
		if (r["src_method"] == "syn")
			source.push_back("SYNTHETIC: YES");
		
		auto& gen = db["entity_src_gen"];
		const std::pair<const char*,const char*> kGenSourceMapping[] = {
			{ "gene_src_common_name",			"ORGANISM_COMMON" },
			{ "pdbx_gene_src_gene",				"GENE" },
			{ "gene_src_strain",				"STRAIN" },
			{ "pdbx_gene_src_cell_line",		"CELL_LINE" },
			{ "pdbx_gene_src_organelle",		"ORGANELLE" },
			{ "pdbx_gene_src_cellular_location","CELLULAR_LOCATION" },
			{ "pdbx_gene_src_scientific_name",	"ORGANISM_SCIENTIFIC" },
			{ "pdbx_gene_src_ncbi_taxonomy_id",	"ORGANISM_TAXID" },
			{ "pdbx_host_org_scientific_name",	"EXPRESSION_SYSTEM" },
			{ "pdbx_host_org_ncbi_taxonomy_id",	"EXPRESSION_SYSTEM_TAXID" },
			{ "pdbx_host_org_strain",			"EXPRESSION_SYSTEM_STRAIN" },
			{ "pdbx_host_org_variant",			"EXPRESSION_SYSTEM_VARIANT" },
			{ "pdbx_host_org_cellular_location","EXPRESSION_SYSTEM_CELLULAR_LOCATION" },
			{ "pdbx_host_org_vector_type",		"EXPRESSION_SYSTEM_VECTOR_TYPE" },
			{ "pdbx_host_org_vector",			"EXPRESSION_SYSTEM_VECTOR" },
			{ "pdbx_host_org_gene",				"EXPRESSION_SYSTEM_GENE" },
			{ "plasmid_name",					"EXPRESSION_SYSTEM_PLASMID" },
			{ "details",						"OTHER_DETAILS" }
		};

		for (auto gr: gen.find(cif::Key("entity_id") == entityID))
		{
			for (auto m: kGenSourceMapping)
			{
				std::string cname, sname;
				tie(cname, sname) = m;
				
				std::string s = gr[cname].as<std::string>();
				if (not s.empty())
					source.push_back(sname + ": " + s);
			}
		}
		
		auto& nat = db["entity_src_nat"];
		const std::pair<const char*, const char*> kNatSourceMapping[] = {
			{ "common_name",				"ORGANISM_COMMON" },
			{ "strain",						"STRAIN" },
			{ "pdbx_organism_scientific",	"ORGANISM_SCIENTIFIC" },
			{ "pdbx_ncbi_taxonomy_id",		"ORGANISM_TAXID" },
			{ "pdbx_cellular_location",		"CELLULAR_LOCATION" },
			{ "pdbx_plasmid_name",			"PLASMID" },
			{ "pdbx_organ",					"ORGAN" },
			{ "details",					"OTHER_DETAILS" }
		};
		
		for (auto nr: nat.find(cif::Key("entity_id") == entityID))
		{
			for (auto m: kNatSourceMapping)
			{
				std::string cname, sname;
				tie(cname, sname) = m;
				
				std::string s = nr[cname].as<std::string>();
				if (not s.empty())
					source.push_back(sname + ": " + s);
			}
		}
	}

	return FixStringLength("SOURCE    " + ba::join(source, "; "), truncate_at);
}

std::string GetPDBAUTHORLine(cif::File& cifFile, std::string::size_type truncate_at)
{
	auto& db = cifFile.firstDatablock();

	// AUTHOR
	std::vector<std::string> author;
	for (auto r: db["audit_author"])
		author.push_back(cif2pdbAuth(r["name"].as<std::string>()));

	return FixStringLength("AUTHOR    " + ba::join(author, "; "), truncate_at);
}
