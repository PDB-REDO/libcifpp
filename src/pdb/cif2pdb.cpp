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

#include "cif++.hpp"

#include <cmath>
#include <deque>
#include <iomanip>
#include <map>
#include <regex>
#include <set>


namespace cif::pdb
{

using namespace std::literals;

// --------------------------------------------------------------------
// conversion routines between cif and pdb format

std::string cif2pdbDate(const std::string &d)
{
	const std::regex rx(R"((\d{4})-(\d{2})(?:-(\d{2}))?)");
	const char *kMonths[12] = {
		"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
	};

	std::smatch m;
	std::string result;

	if (std::regex_match(d, m, rx))
	{
		int year = std::stoi(m[1].str());
		int month = std::stoi(m[2].str());

		if (m[3].matched)
			result = cif::format("%02.2d-%3.3s-%02.2d", stoi(m[3].str()), kMonths[month - 1], (year % 100)).str();
		else
			result = cif::format("%3.3s-%02.2d", kMonths[month - 1], (year % 100)).str();
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

std::string cif2pdbAtomName(std::string name, std::string resName, const datablock &db)
{
	if (name.length() < 4)
	{
		for (auto r : db["atom_site"].find(key("label_atom_id") == name and key("label_comp_id") == resName))
		{
			std::string element = r["type_symbol"].as<std::string>();

			if (element.length() == 1 or not iequals(name, element))
				name.insert(name.begin(), ' ');

			break;
		}
	}

	return name;
}

enum SoftwareType
{
	eRefinement,
	eDataScaling,
	eDataExtraction,
	eDataReduction,
	ePhasing
};

std::string cifSoftware(const datablock &db, SoftwareType sw)
{
	std::string result = "NULL";

	try
	{
		switch (sw)
		{
			case eRefinement: result = db["computing"].find_first<std::string>(key("entry_id") == db.name(), "structure_refinement"); break;
			case eDataScaling: result = db["computing"].find_first<std::string>(key("entry_id") == db.name(), "pdbx_data_reduction_ds"); break;
			case eDataReduction: result = db["computing"].find_first<std::string>(key("entry_id") == db.name(), "pdbx_data_reduction_ii"); break;
			default: break;
		}

		if (result.empty() or result == "NULL")
		{
			auto &software = db["software"];

			row_handle r;

			switch (sw)
			{
				case eRefinement: r = software.find_first(key("classification") == "refinement"); break;
				case eDataScaling: r = software.find_first(key("classification") == "data scaling"); break;
				case eDataExtraction: r = software.find_first(key("classification") == "data extraction"); break;
				case eDataReduction: r = software.find_first(key("classification") == "data reduction"); break;
				case ePhasing: r = software.find_first(key("classification") == "phasing"); break;
			}

			if (not r.empty())
				result = r["name"].as<std::string>() + " " + r["version"].as<std::string>();
		}

		trim(result);
		to_upper(result);

		if (result.empty())
			result = "NULL";
	}
	catch (...)
	{
	}

	return result;
}

// Map asym ID's back to PDB Chain ID's
std::vector<std::string> MapAsymIDs2ChainIDs(const std::vector<std::string> &asymIDs, const datablock &db)
{
	std::set<std::string> result;

	for (auto asym : asymIDs)
	{
		for (auto r : db["pdbx_poly_seq_scheme"].find(key("asym_id") == asym))
		{
			result.insert(r["pdb_strand_id"].as<std::string>());
			break;
		}

		for (auto r : db["pdbx_nonpoly_scheme"].find(key("asym_id") == asym))
		{
			result.insert(r["pdb_strand_id"].as<std::string>());
			break;
		}
	}

	return { result.begin(), result.end() };
}

// support for wrapping text using a 'continuation marker'
std::size_t WriteContinuedLine(std::ostream &pdbFile, std::string header, int &count, int cLen, std::string text, std::string::size_type lStart = 0)
{
	if (lStart == 0)
	{
		if (cLen == 0)
			lStart = header.length() + 1;
		else
			lStart = header.length() + cLen;
	}

	std::string::size_type maxLength = 80 - lStart - 1;

	std::vector<std::string> lines = word_wrap(text, maxLength);

	for (auto &line : lines)
	{
		// to_upper(line);

		pdbFile << header;

		if (++count <= 1 or cLen == 0)
		{
			pdbFile << std::string(lStart - header.length(), ' ');
			if (count == 1)
				lStart = header.length() + cLen + 1;
		}
		else
			pdbFile << std::fixed << std::setw(cLen) << std::right << count << ' ';

		pdbFile << line << '\n';
	}

	return lines.size();
}

std::size_t WriteOneContinuedLine(std::ostream &pdbFile, std::string header, int cLen, std::string line, int lStart = 0)
{
	int count = 0;
	return WriteContinuedLine(pdbFile, header, count, cLen, line, lStart);
}

std::size_t WriteCitation(std::ostream &pdbFile, const datablock &db, row_handle r, int reference)
{
	std::size_t result = 0;

	std::string s1;

	if (reference > 0)
	{
		pdbFile << "REMARK   1 REFERENCE " << std::to_string(reference) << '\n';
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
	for (auto r1 : db["citation_author"].find(key("citation_id") == id))
		authors.push_back(cif2pdbAuth(r1["name"].as<std::string>()));

	if (not authors.empty())
		result += WriteOneContinuedLine(pdbFile, s1 + "AUTH", 2, join(authors, ","), 19);

	result += WriteOneContinuedLine(pdbFile, s1 + "TITL", 2, title, 19);

	if (not pubname.empty())
	{
		to_upper(pubname);

		const std::string kRefHeader = s1 + "REF %2.2s %-28.28s  %2.2s%4.4s %5.5s %4.4s";
		pdbFile << cif::format(kRefHeader, "" /* continuation */, pubname, (volume.empty() ? "" : "V."), volume, pageFirst, year)
				<< '\n';
		++result;
	}

	if (not issn.empty())
	{
		const std::string kRefHeader = s1 + "REFN                   ISSN %-25.25s";
		pdbFile << cif::format(kRefHeader, issn) << '\n';
		++result;
	}

	//		if (not issn.empty() or astm.empty())
	//		{
	////    0         1         2         3         4         5         6         7         8
	////    HEADER    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxDDDDDDDDD   IIII
	// const char kRefHeader[] =
	//          "REMARK   1  REFN    %4.4s %-6.6s  %2.2s %-25.25s";
	//
	//			pdbFile << (boost::cif::format(kRefHeader)
	//						% (astm.empty() ? "" : "ASTN")
	//						% astm
	//						% country
	//						% issn).str()
	//					<< '\n';
	//		}

	if (not pmid.empty())
	{
		const std::string kPMID = s1 + "PMID   %-60.60s ";
		pdbFile << cif::format(kPMID, pmid) << '\n';
		++result;
	}

	if (not doi.empty())
	{
		const std::string kDOI = s1 + "DOI    %-60.60s ";
		pdbFile << cif::format(kDOI, doi) << '\n';
		++result;
	}

	return result;
}

void write_header_lines(std::ostream &pdbFile, const datablock &db)
{
	//    0         1         2         3         4         5         6         7         8
	//    HEADER    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxDDDDDDDDD   IIII
	const char kHeader[] =
		"HEADER    %-40.40s"
		"%-9.9s"
		"   %-4.4s";

	// HEADER

	std::string keywords;
	auto &cat1 = db["struct_keywords"];

	for (auto r : cat1)
	{
		keywords = r["pdbx_keywords"].as<std::string>();
		break;
	}

	std::string date;
	for (auto r : db["pdbx_database_status"])
	{
		date = r["recvd_initial_deposition_date"].as<std::string>();
		if (date.empty())
			continue;
		date = cif2pdbDate(date);
		break;
	}

	if (date.empty())
	{
		for (auto r : db["database_PDB_rev"])
		{
			date = r["date_original"].as<std::string>();
			if (date.empty())
				continue;
			date = cif2pdbDate(date);
			break;
		}
	}

	pdbFile << cif::format(kHeader, keywords, date, db.name()) << '\n';

	// TODO: implement
	// OBSLTE (skip for now)

	// TITLE
	for (auto r : db["struct"])
	{
		std::string title = r["title"].as<std::string>();
		trim(title);
		WriteOneContinuedLine(pdbFile, "TITLE   ", 2, title);
		break;
	}

	// COMPND
	using namespace std::placeholders;

	int molID = 0;
	std::vector<std::string> cmpnd;

	for (auto r : db["entity"])
	{
		if (r["type"] != "polymer")
			continue;

		std::string entityID = r["id"].as<std::string>();

		++molID;
		cmpnd.push_back("MOL_ID: " + std::to_string(molID));

		std::string molecule = r["pdbx_description"].as<std::string>();
		cmpnd.push_back("MOLECULE: " + molecule);

		auto poly = db["entity_poly"].find(key("entity_id") == entityID);
		if (not poly.empty())
		{
			std::string chains = poly.front()["pdbx_strand_id"].as<std::string>();
			replace_all(chains, ",", ", ");
			cmpnd.push_back("CHAIN: " + chains);
		}

		std::string fragment = r["pdbx_fragment"].as<std::string>();
		if (not fragment.empty())
			cmpnd.push_back("FRAGMENT: " + fragment);

		for (auto sr : db["entity_name_com"].find(key("entity_id") == entityID))
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

	WriteOneContinuedLine(pdbFile, "COMPND ", 3, join(cmpnd, ";\n"));

	// SOURCE

	molID = 0;
	std::vector<std::string> source;

	for (auto r : db["entity"])
	{
		if (r["type"] != "polymer")
			continue;

		std::string entityID = r["id"].as<std::string>();

		++molID;
		source.push_back("MOL_ID: " + std::to_string(molID));

		if (r["src_method"] == "syn")
			source.push_back("SYNTHETIC: YES");

		auto &gen = db["entity_src_gen"];
		const std::pair<const char *, const char *> kGenSourceMapping[] = {
			{ "gene_src_common_name", "ORGANISM_COMMON" },
			{ "pdbx_gene_src_gene", "GENE" },
			{ "gene_src_strain", "STRAIN" },
			{ "pdbx_gene_src_cell_line", "CELL_LINE" },
			{ "pdbx_gene_src_organelle", "ORGANELLE" },
			{ "pdbx_gene_src_cellular_location", "CELLULAR_LOCATION" },
			{ "pdbx_gene_src_scientific_name", "ORGANISM_SCIENTIFIC" },
			{ "pdbx_gene_src_ncbi_taxonomy_id", "ORGANISM_TAXID" },
			{ "pdbx_host_org_scientific_name", "EXPRESSION_SYSTEM" },
			{ "pdbx_host_org_ncbi_taxonomy_id", "EXPRESSION_SYSTEM_TAXID" },
			{ "pdbx_host_org_strain", "EXPRESSION_SYSTEM_STRAIN" },
			{ "pdbx_host_org_variant", "EXPRESSION_SYSTEM_VARIANT" },
			{ "pdbx_host_org_cellular_location", "EXPRESSION_SYSTEM_CELLULAR_LOCATION" },
			{ "pdbx_host_org_vector_type", "EXPRESSION_SYSTEM_VECTOR_TYPE" },
			{ "pdbx_host_org_vector", "EXPRESSION_SYSTEM_VECTOR" },
			{ "pdbx_host_org_gene", "EXPRESSION_SYSTEM_GENE" },
			{ "plasmid_name", "EXPRESSION_SYSTEM_PLASMID" },
			{ "details", "OTHER_DETAILS" }
		};

		for (auto gr : gen.find(key("entity_id") == entityID))
		{
			for (const auto &[cname, sname] : kGenSourceMapping)
			{
				std::string s = gr[cname].as<std::string>();
				if (not s.empty())
					source.push_back(sname + ": "s + s);
			}
		}

		auto &nat = db["entity_src_nat"];
		const std::pair<const char *, const char *> kNatSourceMapping[] = {
			{ "common_name", "ORGANISM_COMMON" },
			{ "strain", "STRAIN" },
			{ "pdbx_organism_scientific", "ORGANISM_SCIENTIFIC" },
			{ "pdbx_ncbi_taxonomy_id", "ORGANISM_TAXID" },
			{ "pdbx_cellular_location", "CELLULAR_LOCATION" },
			{ "pdbx_plasmid_name", "PLASMID" },
			{ "pdbx_organ", "ORGAN" },
			{ "details", "OTHER_DETAILS" }
		};

		for (auto nr : nat.find(key("entity_id") == entityID))
		{
			for (const auto &[cname, sname] : kNatSourceMapping)
			{
				std::string s = nr[cname].as<std::string>();
				if (not s.empty())
					source.push_back(sname + ": "s + s);
			}
		}
	}

	WriteOneContinuedLine(pdbFile, "SOURCE ", 3, join(source, ";\n"));

	// KEYWDS

	keywords.clear();
	for (auto r : cat1)
	{
		if (not r["text"].empty())
			keywords += r["text"].as<std::string>();
		else
			keywords += r["pdbx_keywords"].as<std::string>();
	}

	if (not keywords.empty())
		WriteOneContinuedLine(pdbFile, "KEYWDS  ", 2, keywords);

	// EXPDTA

	auto &dbexpt = db["exptl"];
	if (not dbexpt.empty())
	{
		std::vector<std::string> method;
		for (auto r : dbexpt)
			method.push_back(r["method"].as<std::string>());
		if (not method.empty())
			WriteOneContinuedLine(pdbFile, "EXPDTA  ", 2, join(method, "; "));
	}

	// NUMMDL
	// TODO...

	// MDLTYP
	// TODO...

	// AUTHOR
	std::vector<std::string> authors;
	for (auto r : db["audit_author"])
		authors.push_back(cif2pdbAuth(r["name"].as<std::string>()));
	if (not authors.empty())
		WriteOneContinuedLine(pdbFile, "AUTHOR  ", 2, join(authors, ","));
}

void WriteTitle(std::ostream &pdbFile, const datablock &db)
{
	write_header_lines(pdbFile, db);

	// REVDAT
	const char kRevDatFmt[] = "REVDAT %3d%2.2s %9.9s %4.4s    %1d      ";
	auto &cat2 = db["database_PDB_rev"];
	std::vector<row_handle> rev(cat2.begin(), cat2.end());
	sort(rev.begin(), rev.end(), [](row_handle a, row_handle b) -> bool
		{ return a["num"].as<int>() > b["num"].as<int>(); });
	for (auto r : rev)
	{
		int revNum, modType;
		std::string date, replaces;

		cif::tie(revNum, modType, date, replaces) = r.get("num", "mod_type", "date", "replaces");

		date = cif2pdbDate(date);

		std::vector<std::string> types;

		for (auto r1 : db["database_PDB_rev_record"].find(key("rev_num") == revNum))
			types.push_back(r1["type"].as<std::string>());

		int continuation = 0;
		do
		{
			std::string cs = ++continuation > 1 ? std::to_string(continuation) : std::string();

			pdbFile << cif::format(kRevDatFmt, revNum, cs, date, db.name(), modType);
			for (std::size_t i = 0; i < 4; ++i)
				pdbFile << cif::format(" %-6.6s", (i < types.size() ? types[i] : std::string()));
			pdbFile << '\n';

			if (types.size() > 4)
				types.erase(types.begin(), types.begin() + 4);
			else
				types.clear();
		} while (types.empty() == false);
	}

	// SPRSDE
	// TODO...

	// JRNL
	for (auto r : db["citation"])
	{
		WriteCitation(pdbFile, db, r, 0);
		break;
	}
}

void WriteRemark1(std::ostream &pdbFile, const datablock &db)
{
	int reference = 0;

	for (auto r : db["citation"])
	{
		if (reference > 0)
		{
			if (reference == 1)
				pdbFile << "REMARK   1\n";

			WriteCitation(pdbFile, db, r, reference);
		}

		++reference;
	}
}

void WriteRemark2(std::ostream &pdbFile, const datablock &db)
{
	auto &refine = db["refine"];
	if (refine.empty())
	{
		pdbFile << "REMARK   2\n"
				<< "REMARK   2 RESOLUTION. NOT APPLICABLE.\n";
	}
	else
	{
		try
		{
			float resHigh = refine.front()["ls_d_res_high"].as<float>();
			pdbFile << "REMARK   2\n"
					<< cif::format("REMARK   2 RESOLUTION. %7.2f ANGSTROMS.", resHigh) << '\n';
		}
		catch (...)
		{ /* skip it */
		}
	}
}

// --------------------------------------------------------------------
// Code to help format RERMARK 3 data

class FBase
{
  public:
	virtual ~FBase() {}

	virtual void out(std::ostream &os) = 0;

  protected:
	FBase(row_handle r, const char *f)
		: mRow(r)
		, mField(f)
	{
	}
	FBase(const category &cat, condition &&cond, const char *f)
		: mField(f)
	{
		auto r = cat.find(std::move(cond));
		if (not r.empty())
			mRow = r.front();
	}

	std::string_view text() const
	{
		return mRow.empty() or mRow[mField].empty() ? "" : mRow[mField].text();
	}

	row_handle mRow;
	const char *mField;
};

class Fi : public FBase
{
  public:
	Fi(row_handle r, const char *f)
		: FBase(r, f)
	{
	}
	Fi(const category &cat, condition &&cond, const char *f)
		: FBase(cat, std::move(cond), f)
	{
	}

	virtual void out(std::ostream &os)
	{
		std::string s{ text() };

		if (s.empty())
		{
			os << "NULL";
			if (os.width() > 4)
				os << std::string(os.width() - 4, ' ');
		}
		else
		{
			long l = 0;
			auto r = std::from_chars(s.data(), s.data() + s.length(), l);
			if ((bool)r.ec)
			{
				if (VERBOSE > 0)
					std::cerr << "Failed to write '" << s << "' as a long from field " << mField << ", this indicates an error in the code for writing PDB files\n";
				os << s;
			}
			else
				os << l;
		}
	}
};

class Ff : public FBase
{
  public:
	Ff(row_handle r, const char *f)
		: FBase(r, f)
	{
	}
	Ff(const category &cat, condition &&cond, const char *f)
		: FBase(cat, std::move(cond), f)
	{
	}

	virtual void out(std::ostream &os)
	{
		if (mRow.empty() or mRow[mField].empty())
		{
			os << "NULL";
			if (os.width() > 4)
				os << std::string(os.width() - 4, ' ');
		}
		else
		{
			std::string s{ text() };

			double d = 0;
			auto r = cif::from_chars(s.data(), s.data() + s.length(), d);
			if ((bool)r.ec)
			{
				if (VERBOSE > 0)
					std::cerr << "Failed to write '" << s << "' as a double from field " << mField << ", this indicates an error in the code for writing PDB files\n";
				os << s;
			}
			else
				os << d;
		}
	}
};

class Fs : public FBase
{
  public:
	Fs(row_handle r, const char *f, int remarkNr = 3)
		: FBase(r, f)
		, mNr(remarkNr)
	{
	}
	Fs(const category &cat, condition &&cond, const char *f, int remarkNr = 3)
		: FBase(cat, std::move(cond), f)
		, mNr(remarkNr)
	{
	}

	virtual void out(std::ostream &os)
	{
		std::string s{ text() };
		std::size_t width = os.width();

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
			os << '\n';

			std::stringstream ss;
			ss << "REMARK " << std::setw(3) << std::right << mNr << ' ';
			WriteOneContinuedLine(os, ss.str(), 0, s);
		}
	}

	int mNr = 3;
};

std::ostream &operator<<(std::ostream &os, FBase &&fld)
{
	fld.out(os);
	return os;
}

template <int N>
struct RM
{
	RM(const char *desc, int width = 0, int precision = 6)
		: mDesc(desc)
		, mWidth(width)
		, mPrecision(precision)
	{
	}
	const char *mDesc;
	int mWidth, mPrecision;
};

typedef RM<3> RM3;

template <int N>
std::ostream &operator<<(std::ostream &os, RM<N> &&rm)
{
	os << "REMARK " << std::setw(3) << std::right << N << " " << rm.mDesc << (rm.mWidth > 0 ? std::left : std::right) << std::fixed << std::setw(std::abs(rm.mWidth)) << std::setprecision(rm.mPrecision);
	return os;
}

struct SEP
{
	SEP(const char *txt, int width, int precision = 6)
		: mText(txt)
		, mWidth(width)
		, mPrecision(precision)
	{
	}
	const char *mText;
	int mWidth, mPrecision;
};

std::ostream &operator<<(std::ostream &os, SEP &&sep)
{
	os << sep.mText << (sep.mWidth > 0 ? std::left : std::right) << std::fixed << std::setw(std::abs(sep.mWidth)) << std::setprecision(sep.mPrecision);
	return os;
}

// --------------------------------------------------------------------

void WriteRemark3BusterTNT(std::ostream &pdbFile, const datablock &db)
{
	auto refine = db["refine"].front();
	auto ls_shell = db["refine_ls_shell"].front();
	auto hist = db["refine_hist"].front();
	auto reflns = db["reflns"].front();
	auto analyze = db["refine_analyze"].front();
	auto &ls_restr = db["refine_ls_restr"];
	//	auto ls_restr_ncs = db["refine_ls_restr_ncs"].front();
	//	auto pdbx_xplor_file = db["pdbx_xplor_file"].front();
	//	auto pdbx_refine = db["pdbx_refine"].front();

	pdbFile << RM3("") << '\n'
			<< RM3(" DATA USED IN REFINEMENT.") << '\n'
			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_high") << '\n'
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_low") << '\n'
			<< RM3("  DATA CUTOFF            (SIGMA(F)) : ", 6, 3) << Ff(refine, "pdbx_ls_sigma_F") << '\n'
			<< RM3("  COMPLETENESS FOR RANGE        (%) : ", 6, 1) << Ff(refine, "ls_percent_reflns_obs") << '\n'
			<< RM3("  NUMBER OF REFLECTIONS             : ", 12, 6) << Fi(refine, "ls_number_reflns_obs") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT TO DATA USED IN REFINEMENT.") << '\n'
			<< RM3("  CROSS-VALIDATION METHOD          : ") << Fs(refine, "pdbx_ls_cross_valid_method") << '\n'
			<< RM3("  FREE R VALUE TEST SET SELECTION  : ") << Fs(refine, "pdbx_R_Free_selection_details") << '\n'
			<< RM3("  R VALUE     (WORKING + TEST SET) : ", 7, 3) << Ff(refine, "ls_R_factor_obs") << '\n'
			<< RM3("  R VALUE            (WORKING SET) : ", 7, 3) << Ff(refine, "ls_R_factor_R_work") << '\n'
			<< RM3("  FREE R VALUE                     : ", 7, 3) << Ff(refine, "ls_R_factor_R_free") << '\n'
			<< RM3("  FREE R VALUE TEST SET SIZE   (%) : ", 7, 3) << Ff(refine, "ls_percent_reflns_R_free") << '\n'
			<< RM3("  FREE R VALUE TEST SET COUNT      : ", 12, 6) << Fi(refine, "ls_number_reflns_R_free") << '\n'
			<< RM3("  ESTIMATED ERROR OF FREE R VALUE  : ", 7, 3) << Ff(refine, "ls_R_factor_R_free_error") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT IN THE HIGHEST RESOLUTION BIN.") << '\n'
			<< RM3("  TOTAL NUMBER OF BINS USED               : ", 12, 6) << Fi(ls_shell, "pdbx_total_number_of_bins_used") << '\n'
			<< RM3("  BIN RESOLUTION RANGE HIGH   (ANGSTROMS) : ", 5, 2) << Ff(ls_shell, "d_res_high") << '\n'
			<< RM3("  BIN RESOLUTION RANGE LOW    (ANGSTROMS) : ", 5, 2) << Ff(ls_shell, "d_res_low") << '\n'
			<< RM3("  BIN COMPLETENESS     (WORKING+TEST) (%) : ", 6, 2) << Ff(ls_shell, "percent_reflns_obs") << '\n'
			<< RM3("  REFLECTIONS IN BIN (WORKING + TEST SET) : ", 12, 6) << Fi(ls_shell, "number_reflns_all") << '\n'
			<< RM3("  BIN R VALUE        (WORKING + TEST SET) : ", 8, 4) << Ff(ls_shell, "R_factor_all") << '\n'
			<< RM3("  REFLECTIONS IN BIN        (WORKING SET) : ", 12, 6) << Fi(ls_shell, "number_reflns_R_work") << '\n'
			<< RM3("  BIN R VALUE               (WORKING SET) : ", 8, 4) << Ff(ls_shell, "R_factor_R_work") << '\n'
			<< RM3("  BIN FREE R VALUE                        : ", 8, 4) << Ff(ls_shell, "R_factor_R_free") << '\n'
			<< RM3("  BIN FREE R VALUE TEST SET SIZE      (%) : ", 6, 2) << Ff(ls_shell, "percent_reflns_R_free") << '\n'
			<< RM3("  BIN FREE R VALUE TEST SET COUNT         : ", 12, 7) << Fi(ls_shell, "number_reflns_R_free") << '\n'
			<< RM3("  ESTIMATED ERROR OF BIN FREE R VALUE     : ", 7, 3) << Ff(ls_shell, "R_factor_R_free_error") << '\n'

			<< RM3("") << '\n'
			<< RM3(" NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.") << '\n'
			<< RM3("  PROTEIN ATOMS            : ", 12, 6) << Fi(hist, "pdbx_number_atoms_protein") << '\n'
			<< RM3("  NUCLEIC ACID ATOMS       : ", 12, 6) << Fi(hist, "pdbx_number_atoms_nucleic_acid") << '\n'
			<< RM3("  HETEROGEN ATOMS          : ", 12, 6) << Fi(hist, "pdbx_number_atoms_ligand") << '\n'
			<< RM3("  SOLVENT ATOMS            : ", 12, 6) << Fi(hist, "number_atoms_solvent") << '\n'

			<< RM3("") << '\n'
			<< RM3(" B VALUES.") << '\n'
			//			<< RM3("  B VALUE TYPE                      : ")		<< Fs(refine, "pdbx_TLS_residual_ADP_flag") << '\n'
			<< RM3("  FROM WILSON PLOT           (A**2) : ", 7, 2) << Ff(reflns, "B_iso_Wilson_estimate") << '\n'
			<< RM3("  MEAN B VALUE      (OVERALL, A**2) : ", 7, 2) << Ff(refine, "B_iso_mean") << '\n'

			<< RM3("  OVERALL ANISOTROPIC B VALUE.") << '\n'
			<< RM3("   B11 (A**2) : ", -8, 5) << Ff(refine, "aniso_B[1][1]") << '\n'
			<< RM3("   B22 (A**2) : ", -8, 5) << Ff(refine, "aniso_B[2][2]") << '\n'
			<< RM3("   B33 (A**2) : ", -8, 5) << Ff(refine, "aniso_B[3][3]") << '\n'
			<< RM3("   B12 (A**2) : ", -8, 5) << Ff(refine, "aniso_B[1][2]") << '\n'
			<< RM3("   B13 (A**2) : ", -8, 5) << Ff(refine, "aniso_B[1][3]") << '\n'
			<< RM3("   B23 (A**2) : ", -8, 5) << Ff(refine, "aniso_B[2][3]") << '\n'

			<< RM3("") << '\n'
			<< RM3(" ESTIMATED COORDINATE ERROR.") << '\n'
			<< RM3("  ESD FROM LUZZATI PLOT                    (A) : ", 7, 3) << Ff(analyze, "Luzzati_coordinate_error_obs") << '\n'
			<< RM3("  DPI (BLOW EQ-10) BASED ON R VALUE        (A) : ", 5, 3) << Ff(refine, "pdbx_overall_SU_R_Blow_DPI") << '\n'
			<< RM3("  DPI (BLOW EQ-9) BASED ON FREE R VALUE    (A) : ", 5, 3) << Ff(refine, "pdbx_overall_SU_R_free_Blow_DPI") << '\n'
			<< RM3("  DPI (CRUICKSHANK) BASED ON R VALUE       (A) : ", 5, 3) << Ff(refine, "overall_SU_R_Cruickshank_DPI") << '\n'
			<< RM3("  DPI (CRUICKSHANK) BASED ON FREE R VALUE  (A) : ", 5, 3) << Ff(refine, "pdbx_overall_SU_R_free_Cruickshank_DPI") << '\n'

			<< RM3("") << '\n'
			<< RM3("  REFERENCES: BLOW, D. (2002) ACTA CRYST D58, 792-797") << '\n'
			<< RM3("              CRUICKSHANK, D.W.J. (1999) ACTA CRYST D55, 583-601") << '\n'

			<< RM3("") << '\n'
			<< RM3("  CORRELATION COEFFICIENTS.") << '\n'
			<< RM3("  CORRELATION COEFFICIENT FO-FC      : ", 5, 3) << Ff(refine, "correlation_coeff_Fo_to_Fc") << '\n'
			<< RM3("  CORRELATION COEFFICIENT FO-FC FREE : ", 5, 3) << Ff(refine, "correlation_coeff_Fo_to_Fc_free") << '\n'

			<< RM3("") << '\n'
			<< RM3("  NUMBER OF GEOMETRIC FUNCTION TERMS DEFINED : 15") << '\n'
			<< RM3("  TERM                          COUNT    WEIGHT   FUNCTION.") << '\n'
			<< RM3("   BOND LENGTHS              : ", 7, 0) << Ff(ls_restr, key("type") == "t_bond_d", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_bond_d", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_bond_d", "pdbx_restraint_function") << '\n'
			<< RM3("   BOND ANGLES               : ", 7, 0) << Ff(ls_restr, key("type") == "t_angle_deg", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_angle_deg", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_angle_deg", "pdbx_restraint_function") << '\n'
			<< RM3("   TORSION ANGLES            : ", 7, 0) << Ff(ls_restr, key("type") == "t_dihedral_angle_d", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_dihedral_angle_d", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_dihedral_angle_d", "pdbx_restraint_function") << '\n'
			<< RM3("   TRIGONAL CARBON PLANES    : ", 7, 0) << Ff(ls_restr, key("type") == "t_trig_c_planes", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_trig_c_planes", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_trig_c_planes", "pdbx_restraint_function") << '\n'
			<< RM3("   GENERAL PLANES            : ", 7, 0) << Ff(ls_restr, key("type") == "t_gen_planes", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_gen_planes", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_gen_planes", "pdbx_restraint_function") << '\n'
			<< RM3("   ISOTROPIC THERMAL FACTORS : ", 7, 0) << Ff(ls_restr, key("type") == "t_it", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_it", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_it", "pdbx_restraint_function") << '\n'
			<< RM3("   BAD NON-BONDED CONTACTS   : ", 7, 0) << Ff(ls_restr, key("type") == "t_nbd", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_nbd", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_nbd", "pdbx_restraint_function") << '\n'
			<< RM3("   IMPROPER TORSIONS         : ", 7, 0) << Ff(ls_restr, key("type") == "t_improper_torsion", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_improper_torsion", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_improper_torsion", "pdbx_restraint_function") << '\n'
			<< RM3("   PSEUDOROTATION ANGLES     : ", 7, 0) << Ff(ls_restr, key("type") == "t_pseud_angle", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_pseud_angle", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_pseud_angle", "pdbx_restraint_function") << '\n'
			<< RM3("   CHIRAL IMPROPER TORSION   : ", 7, 0) << Ff(ls_restr, key("type") == "t_chiral_improper_torsion", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_chiral_improper_torsion", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_chiral_improper_torsion", "pdbx_restraint_function") << '\n'
			<< RM3("   SUM OF OCCUPANCIES        : ", 7, 0) << Ff(ls_restr, key("type") == "t_sum_occupancies", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_sum_occupancies", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_sum_occupancies", "pdbx_restraint_function") << '\n'
			<< RM3("   UTILITY DISTANCES         : ", 7, 0) << Ff(ls_restr, key("type") == "t_utility_distance", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_utility_distance", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_utility_distance", "pdbx_restraint_function") << '\n'
			<< RM3("   UTILITY ANGLES            : ", 7, 0) << Ff(ls_restr, key("type") == "t_utility_angle", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_utility_angle", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_utility_angle", "pdbx_restraint_function") << '\n'
			<< RM3("   UTILITY TORSION           : ", 7, 0) << Ff(ls_restr, key("type") == "t_utility_torsion", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_utility_torsion", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_utility_torsion", "pdbx_restraint_function") << '\n'
			<< RM3("   IDEAL-DIST CONTACT TERM   : ", 7, 0) << Ff(ls_restr, key("type") == "t_ideal_dist_contact", "number")
			<< SEP("; ", 7, 3) << Ff(ls_restr, key("type") == "t_ideal_dist_contact", "weight")
			<< SEP("; ", 12) << Fs(ls_restr, key("type") == "t_ideal_dist_contact", "pdbx_restraint_function") << '\n'

			<< RM3("") << '\n'
			<< RM3(" RMS DEVIATIONS FROM IDEAL VALUES.") << '\n'
			<< RM3("  BOND LENGTHS                       (A) : ", 7, 3) << Ff(ls_restr, key("type") == "t_bond_d", "dev_ideal") << '\n'
			<< RM3("  BOND ANGLES                  (DEGREES) : ", 7, 2) << Ff(ls_restr, key("type") == "t_angle_deg", "dev_ideal") << '\n'
			<< RM3("  PEPTIDE OMEGA TORSION ANGLES (DEGREES) : ", 7, 2) << Ff(ls_restr, key("type") == "t_omega_torsion", "dev_ideal") << '\n'
			<< RM3("  OTHER TORSION ANGLES         (DEGREES) : ", 7, 2) << Ff(ls_restr, key("type") == "t_other_torsion", "dev_ideal") << '\n';

	auto &tls = db["pdbx_refine_tls"];

	pdbFile << RM3("") << '\n'
			<< RM3(" TLS DETAILS") << '\n'
			<< RM3("  NUMBER OF TLS GROUPS  : ") << (tls.size() ? std::to_string(tls.size()) : "NULL") << '\n';

	for (auto t : tls)
	{
		std::string id = t["id"].as<std::string>();
		auto g = db["pdbx_refine_tls_group"].find_first(key("refine_tls_id") == id);

		pdbFile << RM3("") << '\n'
				<< RM3("  TLS GROUP : ") << id << '\n'
				<< RM3("   SELECTION: ") << Fs(g, "selection_details") << '\n';

		pdbFile << RM3("   ORIGIN FOR THE GROUP (A):", -9, 4) << Ff(t, "origin_x")
				<< SEP("", -9, 4) << Ff(t, "origin_y")
				<< SEP("", -9, 4) << Ff(t, "origin_z") << '\n'
				<< RM3("   T TENSOR") << '\n'
				<< RM3("     T11:", -9, 4) << Ff(t, "T[1][1]") << SEP(" T22:", -9, 4) << Ff(t, "T[2][2]") << '\n'
				<< RM3("     T33:", -9, 4) << Ff(t, "T[3][3]") << SEP(" T12:", -9, 4) << Ff(t, "T[1][2]") << '\n'
				<< RM3("     T13:", -9, 4) << Ff(t, "T[1][3]") << SEP(" T23:", -9, 4) << Ff(t, "T[2][3]") << '\n'
				<< RM3("   L TENSOR") << '\n'
				<< RM3("     L11:", -9, 4) << Ff(t, "L[1][1]") << SEP(" L22:", -9, 4) << Ff(t, "L[2][2]") << '\n'
				<< RM3("     L33:", -9, 4) << Ff(t, "L[3][3]") << SEP(" L12:", -9, 4) << Ff(t, "L[1][2]") << '\n'
				<< RM3("     L13:", -9, 4) << Ff(t, "L[1][3]") << SEP(" L23:", -9, 4) << Ff(t, "L[2][3]") << '\n'
				<< RM3("   S TENSOR") << '\n'
				<< RM3("     S11:", -9, 4) << Ff(t, "S[1][1]") << SEP(" S12:", -9, 4) << Ff(t, "S[1][2]") << SEP(" S13:", -9, 4) << Ff(t, "S[1][3]") << '\n'
				<< RM3("     S21:", -9, 4) << Ff(t, "S[2][1]") << SEP(" S22:", -9, 4) << Ff(t, "S[2][2]") << SEP(" S23:", -9, 4) << Ff(t, "S[2][3]") << '\n'
				<< RM3("     S31:", -9, 4) << Ff(t, "S[3][1]") << SEP(" S32:", -9, 4) << Ff(t, "S[3][2]") << SEP(" S33:", -9, 4) << Ff(t, "S[3][3]") << '\n';
	}

	pdbFile << RM3("") << '\n';
}

// --------------------------------------------------------------------

void WriteRemark3CNS(std::ostream &pdbFile, const datablock &db)
{
	auto refine = db["refine"].front();
	auto ls_shell = db["refine_ls_shell"].front();
	auto hist = db["refine_hist"].front();
	auto reflns = db["reflns"].front();
	auto analyze = db["refine_analyze"].front();
	auto &ls_restr = db["refine_ls_restr"];
	auto ls_restr_ncs = db["refine_ls_restr_ncs"].front();
	//	auto pdbx_xplor_file = db["pdbx_xplor_file"].front();
	//	auto pdbx_refine = db["pdbx_refine"].front();

	pdbFile << RM3("") << '\n'
			<< RM3("REFINEMENT TARGET : ") << Fs(refine, "pdbx_stereochemistry_target_values") << '\n'
			<< RM3("") << '\n'
			<< RM3(" DATA USED IN REFINEMENT.") << '\n'
			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_high") << '\n'
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_low") << '\n'
			<< RM3("  DATA CUTOFF            (SIGMA(F)) : ", 6, 3) << Ff(refine, "pdbx_ls_sigma_F") << '\n'
			<< RM3("  DATA CUTOFF HIGH         (ABS(F)) : ", 6, 3) << Ff(refine, "pdbx_data_cutoff_high_absF") << '\n'
			<< RM3("  DATA CUTOFF LOW          (ABS(F)) : ", 7, 4) << Ff(refine, "pdbx_data_cutoff_low_absF") << '\n'
			<< RM3("  COMPLETENESS (WORKING+TEST)   (%) : ", 4, 1) << Ff(refine, "ls_percent_reflns_obs") << '\n'
			<< RM3("  NUMBER OF REFLECTIONS             : ", 12, 6) << Fi(refine, "ls_number_reflns_obs") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT TO DATA USED IN REFINEMENT.") << '\n'
			<< RM3("  CROSS-VALIDATION METHOD          : ") << Fs(refine, "pdbx_ls_cross_valid_method") << '\n'
			<< RM3("  FREE R VALUE TEST SET SELECTION  : ") << Fs(refine, "pdbx_R_Free_selection_details") << '\n'
			//			<< RM3("  R VALUE     (WORKING + TEST SET) : ", 7, 5)	<< Ff(refine, "ls_R_factor_obs") << '\n'
			<< RM3("  R VALUE            (WORKING SET) : ", 7, 3) << Ff(refine, "ls_R_factor_R_work") << '\n'
			<< RM3("  FREE R VALUE                     : ", 7, 3) << Ff(refine, "ls_R_factor_R_free") << '\n'
			<< RM3("  FREE R VALUE TEST SET SIZE   (%) : ", 7, 3) << Ff(refine, "ls_percent_reflns_R_free") << '\n'
			<< RM3("  FREE R VALUE TEST SET COUNT      : ", 12, 6) << Fi(refine, "ls_number_reflns_R_free") << '\n'
			<< RM3("  ESTIMATED ERROR OF FREE R VALUE  : ", 7, 3) << Ff(refine, "ls_R_factor_R_free_error") << '\n'

			//			<< RM3("") << '\n'
	        //			<< RM3(" FIT/AGREEMENT OF MODEL WITH ALL DATA.") << '\n'
	        //			<< RM3("  R VALUE   (WORKING + TEST SET, NO CUTOFF) : ", 7, 3)	<< Ff(pdbx_refine, "R_factor_all_no_cutoff") << '\n'
	        //			<< RM3("  R VALUE          (WORKING SET, NO CUTOFF) : ", 7, 3)	<< Ff(pdbx_refine, "R_factor_obs_no_cutoff") << '\n'
	        //			<< RM3("  FREE R VALUE                  (NO CUTOFF) : ", 7, 3)	<< Ff(pdbx_refine, "free_R_factor_no_cutoff") << '\n'
	        //			<< RM3("  FREE R VALUE TEST SET SIZE (%, NO CUTOFF) : ", 7, 3)	<< Ff(pdbx_refine, "free_R_val_test_set_size_perc_no_cutoff") << '\n'
	        //			<< RM3("  FREE R VALUE TEST SET COUNT   (NO CUTOFF) : ", 12, 6)	<< Fi(pdbx_refine, "free_R_val_test_set_ct_no_cutoff") << '\n'
	        //			<< RM3("  TOTAL NUMBER OF REFLECTIONS   (NO CUTOFF) : ", 12, 6)	<< Fi(refine, "ls_number_reflns_all") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT IN THE HIGHEST RESOLUTION BIN.") << '\n'
			<< RM3("  TOTAL NUMBER OF BINS USED           : ", 12, 6) << Fi(ls_shell, "pdbx_total_number_of_bins_used") << '\n'
			<< RM3("  BIN RESOLUTION RANGE HIGH       (A) : ", 5, 2) << Ff(ls_shell, "d_res_high") << '\n'
			<< RM3("  BIN RESOLUTION RANGE LOW        (A) : ", 5, 2) << Ff(ls_shell, "d_res_low") << '\n'
			<< RM3("  BIN COMPLETENESS (WORKING+TEST) (%) : ", 6, 2) << Ff(ls_shell, "percent_reflns_obs") << '\n'
			<< RM3("  REFLECTIONS IN BIN    (WORKING SET) : ", 12, 6) << Fi(ls_shell, "number_reflns_R_work") << '\n'
			<< RM3("  BIN R VALUE           (WORKING SET) : ", 8, 4) << Ff(ls_shell, "R_factor_R_work") << '\n'
			<< RM3("  BIN FREE R VALUE                    : ", 8, 4) << Ff(ls_shell, "R_factor_R_free") << '\n'
			<< RM3("  BIN FREE R VALUE TEST SET SIZE  (%) : ", 6, 2) << Ff(ls_shell, "percent_reflns_R_free") << '\n'
			<< RM3("  BIN FREE R VALUE TEST SET COUNT     : ", 12, 7) << Fi(ls_shell, "number_reflns_R_free") << '\n'
			<< RM3("  ESTIMATED ERROR OF BIN FREE R VALUE : ", 7, 3) << Ff(ls_shell, "R_factor_R_free_error") << '\n'

			<< RM3("") << '\n'
			<< RM3(" NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.") << '\n'
			<< RM3("  PROTEIN ATOMS            : ", 12, 6) << Fi(hist, "pdbx_number_atoms_protein") << '\n'
			<< RM3("  NUCLEIC ACID ATOMS       : ", 12, 6) << Fi(hist, "pdbx_number_atoms_nucleic_acid") << '\n'
			<< RM3("  HETEROGEN ATOMS          : ", 12, 6) << Fi(hist, "pdbx_number_atoms_ligand") << '\n'
			<< RM3("  SOLVENT ATOMS            : ", 12, 6) << Fi(hist, "number_atoms_solvent") << '\n'

			<< RM3("") << '\n'
			<< RM3(" B VALUES.") << '\n'
			<< RM3("  B VALUE TYPE                      : ") << Fs(refine, "pdbx_TLS_residual_ADP_flag") << '\n'
			<< RM3("  FROM WILSON PLOT           (A**2) : ", 7, 2) << Ff(reflns, "B_iso_Wilson_estimate") << '\n'
			<< RM3("  MEAN B VALUE      (OVERALL, A**2) : ", 7, 2) << Ff(refine, "B_iso_mean") << '\n'

			<< RM3("  OVERALL ANISOTROPIC B VALUE.") << '\n'
			<< RM3("   B11 (A**2) : ", -8, 5) << Ff(refine, "aniso_B[1][1]") << '\n'
			<< RM3("   B22 (A**2) : ", -8, 5) << Ff(refine, "aniso_B[2][2]") << '\n'
			<< RM3("   B33 (A**2) : ", -8, 5) << Ff(refine, "aniso_B[3][3]") << '\n'
			<< RM3("   B12 (A**2) : ", -8, 5) << Ff(refine, "aniso_B[1][2]") << '\n'
			<< RM3("   B13 (A**2) : ", -8, 5) << Ff(refine, "aniso_B[1][3]") << '\n'
			<< RM3("   B23 (A**2) : ", -8, 5) << Ff(refine, "aniso_B[2][3]") << '\n'

			<< RM3("") << '\n'
			<< RM3(" ESTIMATED COORDINATE ERROR.") << '\n'
			<< RM3("  ESD FROM LUZZATI PLOT        (A) : ", 7, 2) << Ff(analyze, "Luzzati_coordinate_error_obs") << '\n'
			<< RM3("  ESD FROM SIGMAA              (A) : ", 7, 2) << Ff(analyze, "Luzzati_sigma_a_obs") << '\n'
			<< RM3("  LOW RESOLUTION CUTOFF        (A) : ", 7, 2) << Ff(analyze, "Luzzati_d_res_low_obs") << '\n'

			<< RM3("") << '\n'
			<< RM3(" CROSS-VALIDATED ESTIMATED COORDINATE ERROR.") << '\n'
			<< RM3("  ESD FROM C-V LUZZATI PLOT    (A) : ", 7, 2) << Ff(analyze, "Luzzati_coordinate_error_free") << '\n'
			<< RM3("  ESD FROM C-V SIGMAA          (A) : ", 7, 2) << Ff(analyze, "Luzzati_sigma_a_free") << '\n'

			<< RM3("") << '\n'
			<< RM3(" RMS DEVIATIONS FROM IDEAL VALUES.") << '\n'
			<< RM3("  BOND LENGTHS                 (A) : ", 7, 3) << Ff(ls_restr, key("type") == "c_bond_d", "dev_ideal") << '\n'
			<< RM3("  BOND ANGLES            (DEGREES) : ", 7, 2) << Ff(ls_restr, key("type") == "c_angle_deg", "dev_ideal") << '\n'
			<< RM3("  DIHEDRAL ANGLES        (DEGREES) : ", 7, 2) << Ff(ls_restr, key("type") == "c_dihedral_angle_d", "dev_ideal") << '\n'
			<< RM3("  IMPROPER ANGLES        (DEGREES) : ", 7, 2) << Ff(ls_restr, key("type") == "c_improper_angle_d", "dev_ideal") << '\n'

			<< RM3("") << '\n'
			<< RM3(" ISOTROPIC THERMAL MODEL : ") << Fs(refine, "pdbx_isotropic_thermal_model") << '\n'

			<< RM3("") << '\n'
			<< RM3(" ISOTROPIC THERMAL FACTOR RESTRAINTS.    RMS    SIGMA") << '\n'
			<< RM3("  MAIN-CHAIN BOND              (A**2) : ", 7, 3) << Ff(ls_restr, key("type") == "c_mcbond_it", "dev_ideal") << SEP("; ", 7, 3)
			<< Ff(ls_restr, key("type") == "c_mcbond_it", "dev_ideal_target") << '\n'
			<< RM3("  MAIN-CHAIN ANGLE             (A**2) : ", 7, 3) << Ff(ls_restr, key("type") == "c_mcangle_it", "dev_ideal") << SEP("; ", 7, 3)
			<< Ff(ls_restr, key("type") == "c_mcangle_it", "dev_ideal_target") << '\n'
			<< RM3("  SIDE-CHAIN BOND              (A**2) : ", 7, 3) << Ff(ls_restr, key("type") == "c_scbond_it", "dev_ideal") << SEP("; ", 7, 3)
			<< Ff(ls_restr, key("type") == "c_scbond_it", "dev_ideal_target") << '\n'
			<< RM3("  SIDE-CHAIN ANGLE             (A**2) : ", 7, 3) << Ff(ls_restr, key("type") == "c_scangle_it", "dev_ideal") << SEP("; ", 7, 3)
			<< Ff(ls_restr, key("type") == "c_scangle_it", "dev_ideal_target") << '\n'

			<< RM3("") << '\n'
			<< RM3(" BULK SOLVENT MODELING.") << '\n'
			<< RM3("  METHOD USED        : ") << Fs(refine, "solvent_model_details") << '\n'
			<< RM3("  KSOL               : ", 5, 2) << Ff(refine, "solvent_model_param_ksol") << '\n'
			<< RM3("  BSOL               : ", 5, 2) << Ff(refine, "solvent_model_param_bsol") << '\n'

			<< RM3("") << '\n'
			<< RM3(" NCS MODEL : ") << Fs(ls_restr_ncs, "ncs_model_details") << '\n'

			<< RM3("") << '\n'
			<< RM3(" NCS RESTRAINTS.                         RMS   SIGMA/WEIGHT") << '\n'

			// TODO: using only group 1 here, should this be fixed???
			<< RM3("  GROUP  1  POSITIONAL            (A) : ", 4, 2) << Ff(ls_restr_ncs, "rms_dev_position") << SEP("; ", 6, 2)
			<< Ff(ls_restr_ncs, "weight_position") << SEP("; ", 6, 2) << '\n'
			<< RM3("  GROUP  1  B-FACTOR           (A**2) : ", 4, 2) << Ff(ls_restr_ncs, "rms_dev_B_iso") << SEP("; ", 6, 2)
			<< Ff(ls_restr_ncs, "weight_B_iso") << SEP("; ", 6, 2) << '\n'

			// TODO: using only files from serial_no 1 here
	        //			<< RM3("") << '\n'
	        //			<< RM3(" PARAMETER FILE   1  : ") << Fs(pdbx_xplor_file, "param_file") << '\n'
	        //			<< RM3(" TOPOLOGY FILE   1   : ") << Fs(pdbx_xplor_file, "topol_file") << '\n'

			<< RM3("") << '\n';
}

// --------------------------------------------------------------------

void WriteRemark3Refmac(std::ostream &pdbFile, const datablock &db)
{
	auto refine = db["refine"].front();
	auto ls_shell = db["refine_ls_shell"].front();
	auto hist = db["refine_hist"].front();
	auto reflns = db["reflns"].front();
	//	auto analyze = db["refine_analyze"].front();
	auto &ls_restr = db["refine_ls_restr"];
	//	auto pdbx_xplor_file = db["pdbx_xplor_file"].front();

	auto c = [](const char *t) -> condition
	{ return key("type") == t; };

	pdbFile << RM3("") << '\n'
			<< RM3("REFINEMENT TARGET : ") << Fs(refine, "pdbx_stereochemistry_target_values") << '\n'
			<< RM3("") << '\n'
			<< RM3(" DATA USED IN REFINEMENT.") << '\n'
			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_high") << '\n'
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_low") << '\n'
			<< RM3("  DATA CUTOFF            (SIGMA(F)) : ", 6, 3) << Ff(refine, "pdbx_ls_sigma_F") << '\n'
			<< RM3("  COMPLETENESS FOR RANGE        (%) : ", 5, 2) << Ff(refine, "ls_percent_reflns_obs") << '\n'
			<< RM3("  NUMBER OF REFLECTIONS             : ", 12, 6) << Fi(refine, "ls_number_reflns_obs") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT TO DATA USED IN REFINEMENT.") << '\n'
			<< RM3("  CROSS-VALIDATION METHOD          : ") << Fs(refine, "pdbx_ls_cross_valid_method") << '\n'
			<< RM3("  FREE R VALUE TEST SET SELECTION  : ") << Fs(refine, "pdbx_R_Free_selection_details") << '\n'
			<< RM3("  R VALUE     (WORKING + TEST SET) : ", 7, 5) << Ff(refine, "ls_R_factor_obs") << '\n'
			<< RM3("  R VALUE            (WORKING SET) : ", 7, 5) << Ff(refine, "ls_R_factor_R_work") << '\n'
			<< RM3("  FREE R VALUE                     : ", 7, 5) << Ff(refine, "ls_R_factor_R_free") << '\n'
			<< RM3("  FREE R VALUE TEST SET SIZE   (%) : ", 7, 1) << Ff(refine, "ls_percent_reflns_R_free") << '\n'
			<< RM3("  FREE R VALUE TEST SET COUNT      : ", 12, 6) << Fi(refine, "ls_number_reflns_R_free") << '\n'
			<< RM3("  ESTIMATED ERROR OF FREE R VALUE  : ", 7, 3) << Ff(refine, "ls_R_factor_R_free_error") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT IN THE HIGHEST RESOLUTION BIN.") << '\n'
			<< RM3("  TOTAL NUMBER OF BINS USED           : ") << Fi(ls_shell, "pdbx_total_number_of_bins_used") << '\n'
			<< RM3("  BIN RESOLUTION RANGE HIGH       (A) : ", 5, 3) << Ff(ls_shell, "d_res_high") << '\n'
			<< RM3("  BIN RESOLUTION RANGE LOW        (A) : ", 5, 3) << Ff(ls_shell, "d_res_low") << '\n'
			<< RM3("  REFLECTION IN BIN     (WORKING SET) : ") << Fi(ls_shell, "number_reflns_R_work") << '\n'
			<< RM3("  BIN COMPLETENESS (WORKING+TEST) (%) : ", 5, 2) << Ff(ls_shell, "percent_reflns_obs") << '\n'
			<< RM3("  BIN R VALUE           (WORKING SET) : ", 7, 3) << Ff(ls_shell, "R_factor_R_work") << '\n'
			<< RM3("  BIN FREE R VALUE SET COUNT          : ") << Fi(ls_shell, "number_reflns_R_free") << '\n'
			<< RM3("  BIN FREE R VALUE                    : ", 7, 3) << Ff(ls_shell, "R_factor_R_free") << '\n'

			<< RM3("") << '\n'
			<< RM3(" NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.") << '\n'
			<< RM3("  PROTEIN ATOMS            : ") << Fi(hist, "pdbx_number_atoms_protein") << '\n'
			<< RM3("  NUCLEIC ACID ATOMS       : ") << Fi(hist, "pdbx_number_atoms_nucleic_acid") << '\n'
			<< RM3("  HETEROGEN ATOMS          : ") << Fi(hist, "pdbx_number_atoms_ligand") << '\n'
			<< RM3("  SOLVENT ATOMS            : ") << Fi(hist, "number_atoms_solvent") << '\n'

			<< RM3("") << '\n'
			<< RM3(" B VALUES.") << '\n'
			<< RM3("  B VALUE TYPE                      : ") << Fs(refine, "pdbx_TLS_residual_ADP_flag") << '\n'
			<< RM3("  FROM WILSON PLOT           (A**2) : ", 8, 3) << Ff(reflns, "B_iso_Wilson_estimate") << '\n'
			<< RM3("  MEAN B VALUE      (OVERALL, A**2) : ", 8, 3) << Ff(refine, "B_iso_mean") << '\n'

			<< RM3("  OVERALL ANISOTROPIC B VALUE.") << '\n'
			<< RM3("   B11 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][1]") << '\n'
			<< RM3("   B22 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[2][2]") << '\n'
			<< RM3("   B33 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[3][3]") << '\n'
			<< RM3("   B12 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][2]") << '\n'
			<< RM3("   B13 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][3]") << '\n'
			<< RM3("   B23 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[2][3]") << '\n'

			<< RM3("") << '\n'
			<< RM3(" ESTIMATED OVERALL COORDINATE ERROR.") << '\n'
			<< RM3("  ESU BASED ON R VALUE                            (A): ", 6, 3) << Ff(refine, "pdbx_overall_ESU_R") << '\n'
			<< RM3("  ESU BASED ON FREE R VALUE                       (A): ", 6, 3) << Ff(refine, "pdbx_overall_ESU_R_Free") << '\n'
			<< RM3("  ESU BASED ON MAXIMUM LIKELIHOOD                 (A): ", 6, 3) << Ff(refine, "overall_SU_ML") << '\n'
			<< RM3("  ESU FOR B VALUES BASED ON MAXIMUM LIKELIHOOD (A**2): ", 6, 3) << Ff(refine, "overall_SU_B") << '\n'

			<< RM3("") << '\n'
			<< RM3(" CORRELATION COEFFICIENTS.") << '\n'
			<< RM3("  CORRELATION COEFFICIENT FO-FC      : ", 6, 3) << Ff(refine, "correlation_coeff_Fo_to_Fc") << '\n'
			<< RM3("  CORRELATION COEFFICIENT FO-FC FREE : ", 6, 3) << Ff(refine, "correlation_coeff_Fo_to_Fc_free") << '\n'

			<< RM3("") << '\n'
			<< RM3(" RMS DEVIATIONS FROM IDEAL VALUES        COUNT    RMS    WEIGHT") << '\n'
			<< RM3("  BOND LENGTHS REFINED ATOMS        (A): ", -5) << Fi(ls_restr, c("r_bond_refined_d"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_bond_refined_d"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_bond_refined_d"), "dev_ideal_target") << '\n'
			<< RM3("  BOND LENGTHS OTHERS               (A): ", -5) << Fi(ls_restr, c("r_bond_other_d"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_bond_other_d"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_bond_other_d"), "dev_ideal_target") << '\n'
			<< RM3("  BOND ANGLES REFINED ATOMS   (DEGREES): ", -5) << Fi(ls_restr, c("r_angle_refined_deg"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_angle_refined_deg"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_angle_refined_deg"), "dev_ideal_target") << '\n'
			<< RM3("  BOND ANGLES OTHERS          (DEGREES): ", -5) << Fi(ls_restr, c("r_angle_other_deg"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_angle_other_deg"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_angle_other_deg"), "dev_ideal_target") << '\n'
			<< RM3("  TORSION ANGLES, PERIOD 1    (DEGREES): ", -5) << Fi(ls_restr, c("r_dihedral_angle_1_deg"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_dihedral_angle_1_deg"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_dihedral_angle_1_deg"), "dev_ideal_target") << '\n'
			<< RM3("  TORSION ANGLES, PERIOD 2    (DEGREES): ", -5) << Fi(ls_restr, c("r_dihedral_angle_2_deg"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_dihedral_angle_2_deg"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_dihedral_angle_2_deg"), "dev_ideal_target") << '\n'
			<< RM3("  TORSION ANGLES, PERIOD 3    (DEGREES): ", -5) << Fi(ls_restr, c("r_dihedral_angle_3_deg"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_dihedral_angle_3_deg"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_dihedral_angle_3_deg"), "dev_ideal_target") << '\n'
			<< RM3("  TORSION ANGLES, PERIOD 4    (DEGREES): ", -5) << Fi(ls_restr, c("r_dihedral_angle_4_deg"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_dihedral_angle_4_deg"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_dihedral_angle_4_deg"), "dev_ideal_target") << '\n'
			<< RM3("  CHIRAL-CENTER RESTRAINTS       (A**3): ", -5) << Fi(ls_restr, c("r_chiral_restr"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_chiral_restr"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_chiral_restr"), "dev_ideal_target") << '\n'
			<< RM3("  GENERAL PLANES REFINED ATOMS      (A): ", -5) << Fi(ls_restr, c("r_gen_planes_refined"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_gen_planes_refined"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_gen_planes_refined"), "dev_ideal_target") << '\n'
			<< RM3("  GENERAL PLANES OTHERS             (A): ", -5) << Fi(ls_restr, c("r_gen_planes_other"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_gen_planes_other"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_gen_planes_other"), "dev_ideal_target") << '\n'
			<< RM3("  NON-BONDED CONTACTS REFINED ATOMS (A): ", -5) << Fi(ls_restr, c("r_nbd_refined"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_nbd_refined"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_nbd_refined"), "dev_ideal_target") << '\n'
			<< RM3("  NON-BONDED CONTACTS OTHERS        (A): ", -5) << Fi(ls_restr, c("r_nbd_other"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_nbd_other"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_nbd_other"), "dev_ideal_target") << '\n'
			<< RM3("  NON-BONDED TORSION REFINED ATOMS  (A): ", -5) << Fi(ls_restr, c("r_nbtor_refined"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_nbtor_refined"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_nbtor_refined"), "dev_ideal_target") << '\n'
			<< RM3("  NON-BONDED TORSION OTHERS         (A): ", -5) << Fi(ls_restr, c("r_nbtor_other"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_nbtor_other"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_nbtor_other"), "dev_ideal_target") << '\n'
			<< RM3("  H-BOND (X...Y) REFINED ATOMS      (A): ", -5) << Fi(ls_restr, c("r_xyhbond_nbd_refined"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_xyhbond_nbd_refined"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_xyhbond_nbd_refined"), "dev_ideal_target") << '\n'
			<< RM3("  H-BOND (X...Y) OTHERS             (A): ", -5) << Fi(ls_restr, c("r_xyhbond_nbd_other"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_xyhbond_nbd_other"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_xyhbond_nbd_other"), "dev_ideal_target") << '\n'
			<< RM3("  POTENTIAL METAL-ION REFINED ATOMS (A): ", -5) << Fi(ls_restr, c("r_metal_ion_refined"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_metal_ion_refined"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_metal_ion_refined"), "dev_ideal_target") << '\n'
			<< RM3("  POTENTIAL METAL-ION OTHERS        (A): ", -5) << Fi(ls_restr, c("r_metal_ion_other"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_metal_ion_other"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_metal_ion_other"), "dev_ideal_target") << '\n'
			<< RM3("  SYMMETRY VDW REFINED ATOMS        (A): ", -5) << Fi(ls_restr, c("r_symmetry_vdw_refined"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_symmetry_vdw_refined"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_symmetry_vdw_refined"), "dev_ideal_target") << '\n'
			<< RM3("  SYMMETRY VDW OTHERS               (A): ", -5) << Fi(ls_restr, c("r_symmetry_vdw_other"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_symmetry_vdw_other"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_symmetry_vdw_other"), "dev_ideal_target") << '\n'
			<< RM3("  SYMMETRY H-BOND REFINED ATOMS     (A): ", -5) << Fi(ls_restr, c("r_symmetry_hbond_refined"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_symmetry_hbond_refined"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_symmetry_hbond_refined"), "dev_ideal_target") << '\n'
			<< RM3("  SYMMETRY H-BOND OTHERS            (A): ", -5) << Fi(ls_restr, c("r_symmetry_hbond_other"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_symmetry_hbond_other"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_symmetry_hbond_other"), "dev_ideal_target") << '\n'
			<< RM3("  SYMMETRY METAL-ION REFINED ATOMS  (A): ", -5) << Fi(ls_restr, c("r_symmetry_metal_ion_refined"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_symmetry_metal_ion_refined"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_symmetry_metal_ion_refined"), "dev_ideal_target") << '\n'
			<< RM3("  SYMMETRY METAL-ION OTHERS         (A): ", -5) << Fi(ls_restr, c("r_symmetry_metal_ion_other"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_symmetry_metal_ion_other"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_symmetry_metal_ion_other"), "dev_ideal_target") << '\n'

			<< RM3("") << '\n'
			<< RM3(" ISOTROPIC THERMAL FACTOR RESTRAINTS.     COUNT   RMS    WEIGHT") << '\n'
			<< RM3("  MAIN-CHAIN BOND REFINED ATOMS  (A**2): ", -5) << Fi(ls_restr, c("r_mcbond_it"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_mcbond_it"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_mcbond_it"), "dev_ideal_target") << '\n'
			<< RM3("  MAIN-CHAIN BOND OTHER ATOMS    (A**2): ", -5) << Fi(ls_restr, c("r_mcbond_other"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_mcbond_other"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_mcbond_other"), "dev_ideal_target") << '\n'
			<< RM3("  MAIN-CHAIN ANGLE REFINED ATOMS (A**2): ", -5) << Fi(ls_restr, c("r_mcangle_it"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_mcangle_it"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_mcangle_it"), "dev_ideal_target") << '\n'
			<< RM3("  MAIN-CHAIN ANGLE OTHER ATOMS   (A**2): ", -5) << Fi(ls_restr, c("r_mcangle_other"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_mcangle_other"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_mcangle_other"), "dev_ideal_target") << '\n'
			<< RM3("  SIDE-CHAIN BOND REFINED ATOMS  (A**2): ", -5) << Fi(ls_restr, c("r_scbond_it"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_scbond_it"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_scbond_it"), "dev_ideal_target") << '\n'
			<< RM3("  SIDE-CHAIN BOND OTHER ATOMS    (A**2): ", -5) << Fi(ls_restr, c("r_scbond_other"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_scbond_other"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_scbond_other"), "dev_ideal_target") << '\n'
			<< RM3("  SIDE-CHAIN ANGLE REFINED ATOMS (A**2): ", -5) << Fi(ls_restr, c("r_scangle_it"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_scangle_it"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_scangle_it"), "dev_ideal_target") << '\n'
			<< RM3("  SIDE-CHAIN ANGLE OTHER ATOMS   (A**2): ", -5) << Fi(ls_restr, c("r_scangle_other"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_scangle_other"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_scangle_other"), "dev_ideal_target") << '\n'
			<< RM3("  LONG RANGE B REFINED ATOMS     (A**2): ", -5) << Fi(ls_restr, c("r_long_range_B_refined"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_long_range_B_refined"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_long_range_B_refined"), "dev_ideal_target") << '\n'
			<< RM3("  LONG RANGE B OTHER ATOMS       (A**2): ", -5) << Fi(ls_restr, c("r_long_range_B_other"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_long_range_B_other"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_long_range_B_other"), "dev_ideal_target") << '\n'

			<< RM3("") << '\n'
			<< RM3(" ANISOTROPIC THERMAL FACTOR RESTRAINTS.   COUNT   RMS    WEIGHT") << '\n'
			<< RM3("  RIGID-BOND RESTRAINTS          (A**2): ", -5) << Fi(ls_restr, c("r_rigid_bond_restr"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_rigid_bond_restr"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_rigid_bond_restr"), "dev_ideal_target") << '\n'
			<< RM3("  SPHERICITY; FREE ATOMS         (A**2): ", -5) << Fi(ls_restr, c("r_sphericity_free"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_sphericity_free"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_sphericity_free"), "dev_ideal_target") << '\n'
			<< RM3("  SPHERICITY; BONDED ATOMS       (A**2): ", -5) << Fi(ls_restr, c("r_sphericity_bonded"), "number") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_sphericity_bonded"), "dev_ideal") << SEP(" ;", -6, 3)
			<< Ff(ls_restr, c("r_sphericity_bonded"), "dev_ideal_target") << '\n'

			<< RM3("") << '\n'
			<< RM3(" NCS RESTRAINTS STATISTICS") << '\n';

	auto &ncs_dom = db["struct_ncs_dom"];
	if (ncs_dom.empty())
		pdbFile << RM3("  NUMBER OF DIFFERENT NCS GROUPS : NULL") << '\n';
	else
	{
		std::set<std::string> ncs_groups;
		for (auto i : ncs_dom)
			ncs_groups.insert(i["pdbx_ens_id"].as<std::string>());

		pdbFile << RM3("  NUMBER OF DIFFERENT NCS GROUPS : ") << ncs_groups.size() << '\n';

		for (auto ens_id : ncs_groups)
		{
			auto lim = db["struct_ncs_dom_lim"].find(key("pdbx_ens_id") == ens_id);

			std::set<std::string> chains;
			std::set<int> component_ids;

			for (auto l : lim)
			{
				chains.insert(l["beg_auth_asym_id"].as<std::string>());
				component_ids.insert(l["pdbx_component_id"].as<int>());
			}

			pdbFile << RM3("") << '\n'
					<< RM3(" NCS GROUP NUMBER               : ") << ens_id << '\n'
					<< RM3("    CHAIN NAMES                    : ") << join(chains, " ") << '\n'
					<< RM3("    NUMBER OF COMPONENTS NCS GROUP : ") << component_ids.size() << '\n'
					<< RM3("      COMPONENT C  SSSEQI  TO  C   SSSEQI   CODE") << '\n';

			for (auto l : lim)
			{
				pdbFile << RM3("         ", -2) << Fi(l, "pdbx_component_id")
						<< SEP(" ", -5) << Fs(l, "beg_auth_asym_id")
						<< SEP("  ", -5) << Fi(l, "beg_auth_seq_id")
						<< SEP("   ", -5) << Fs(l, "end_auth_asym_id")
						<< SEP("   ", -5) << Fi(l, "end_auth_seq_id")
						<< SEP("  ", -5) << Fs(l, "pdbx_refine_code")
						<< '\n';
			}

			pdbFile << RM3("                  GROUP CHAIN        COUNT   RMS     WEIGHT") << '\n';
			for (auto l : db["refine_ls_restr_ncs"].find(key("pdbx_ens_id") == ens_id))
			{
				std::string type = l["pdbx_type"].as<std::string>();
				to_upper(type);

				std::string unit;
				if (ends_with(type, "POSITIONAL"))
					unit = "    (A): ";
				else if (ends_with(type, "THERMAL"))
					unit = " (A**2): ";
				else
					unit = "       : ";

				pdbFile << RM3("  ", 18) << type
						<< SEP("", -2) << Fs(l, "pdbx_ens_id")
						<< SEP("    ", 1) << Fs(l, "pdbx_auth_asym_id")
						<< SEP(unit.c_str(), -6) << Fi(l, "pdbx_number")
						<< SEP(" ;", -6, 3) << Ff(l, "rms_dev_position")
						<< SEP(" ;", -6, 3) << Ff(l, "weight_position")
						<< '\n';
			}
		}
	}

	pdbFile << RM3("") << '\n'
			<< RM3(" TWIN DETAILS") << '\n';

	auto &twins = db["pdbx_reflns_twin"];
	if (twins.empty())
		pdbFile << RM3("  NUMBER OF TWIN DOMAINS  : NULL") << '\n';
	else
	{
		pdbFile << RM3("  NUMBER OF TWIN DOMAINS  :    ") << twins.size() << '\n';

		int nr = 1;
		for (auto twin : twins)
		{
			pdbFile << RM3("     TWIN DOMAIN   : ") << nr++ << '\n'
					<< RM3("     TWIN OPERATOR : ") << Fs(twin, "operator") << '\n'
					<< RM3("     TWIN FRACTION : ") << SEP("", -6, 3) << Ff(twin, "fraction") << '\n';
		}
	}

	auto &tls = db["pdbx_refine_tls"];

	pdbFile << RM3("") << '\n'
			<< RM3(" TLS DETAILS") << '\n'
			<< RM3("  NUMBER OF TLS GROUPS  : ") << (tls.size() ? std::to_string(tls.size()) : "NULL") << '\n';

	for (auto t : tls)
	{
		std::string id = t["id"].as<std::string>();
		auto g = db["pdbx_refine_tls_group"].find(key("refine_tls_id") == id);

		pdbFile << RM3("") << '\n'
				<< RM3("  TLS GROUP : ") << id << '\n'
				<< RM3("   NUMBER OF COMPONENTS GROUP : ") << g.size() << '\n'
				<< RM3("   COMPONENTS        C SSSEQI   TO  C SSSEQI") << '\n';

		for (auto gi : g)
		{
			pdbFile << RM3("   RESIDUE RANGE :   ") << Fs(gi, "beg_auth_asym_id")
					<< SEP("", -6) << Fi(gi, "beg_auth_seq_id")
					<< SEP("", -9) << Fs(gi, "end_auth_asym_id")
					<< SEP("", -6) << Fi(gi, "end_auth_seq_id")
					<< '\n';
		}

		pdbFile << RM3("   ORIGIN FOR THE GROUP (A):", -9, 4) << Ff(t, "origin_x")
				<< SEP("", -9, 4) << Ff(t, "origin_y")
				<< SEP("", -9, 4) << Ff(t, "origin_z") << '\n'
				<< RM3("   T TENSOR") << '\n'
				<< RM3("     T11:", -9, 4) << Ff(t, "T[1][1]") << SEP(" T22:", -9, 4) << Ff(t, "T[2][2]") << '\n'
				<< RM3("     T33:", -9, 4) << Ff(t, "T[3][3]") << SEP(" T12:", -9, 4) << Ff(t, "T[1][2]") << '\n'
				<< RM3("     T13:", -9, 4) << Ff(t, "T[1][3]") << SEP(" T23:", -9, 4) << Ff(t, "T[2][3]") << '\n'
				<< RM3("   L TENSOR") << '\n'
				<< RM3("     L11:", -9, 4) << Ff(t, "L[1][1]") << SEP(" L22:", -9, 4) << Ff(t, "L[2][2]") << '\n'
				<< RM3("     L33:", -9, 4) << Ff(t, "L[3][3]") << SEP(" L12:", -9, 4) << Ff(t, "L[1][2]") << '\n'
				<< RM3("     L13:", -9, 4) << Ff(t, "L[1][3]") << SEP(" L23:", -9, 4) << Ff(t, "L[2][3]") << '\n'
				<< RM3("   S TENSOR") << '\n'
				<< RM3("     S11:", -9, 4) << Ff(t, "S[1][1]") << SEP(" S12:", -9, 4) << Ff(t, "S[1][2]") << SEP(" S13:", -9, 4) << Ff(t, "S[1][3]") << '\n'
				<< RM3("     S21:", -9, 4) << Ff(t, "S[2][1]") << SEP(" S22:", -9, 4) << Ff(t, "S[2][2]") << SEP(" S23:", -9, 4) << Ff(t, "S[2][3]") << '\n'
				<< RM3("     S31:", -9, 4) << Ff(t, "S[3][1]") << SEP(" S32:", -9, 4) << Ff(t, "S[3][2]") << SEP(" S33:", -9, 4) << Ff(t, "S[3][3]") << '\n';
	}

	pdbFile << RM3("") << '\n'
			<< RM3(" BULK SOLVENT MODELLING.") << '\n'
			<< RM3("  METHOD USED : ") << Fs(refine, "solvent_model_details") << '\n'
			<< RM3("  PARAMETERS FOR MASK CALCULATION") << '\n'
			<< RM3("  VDW PROBE RADIUS   : ", 5, 2) << Ff(refine, "pdbx_solvent_vdw_probe_radii") << '\n'
			<< RM3("  ION PROBE RADIUS   : ", 5, 2) << Ff(refine, "pdbx_solvent_ion_probe_radii") << '\n'
			<< RM3("  SHRINKAGE RADIUS   : ", 5, 2) << Ff(refine, "pdbx_solvent_shrinkage_radii") << '\n'

			<< RM3("") << '\n';
}

void WriteRemark3Shelxl(std::ostream &pdbFile, const datablock &db)
{
	auto refine = db["refine"].front();
	//	auto ls_shell = db["refine_ls_shell"].front();
	auto refine_hist = db["refine_hist"].front();
	//	auto reflns = db["reflns"].front();
	auto refine_analyze = db["refine_analyze"].front();
	auto &ls_restr = db["refine_ls_restr"];
	//	auto pdbx_xplor_file = db["pdbx_xplor_file"].front();
	auto pdbx_refine = db["pdbx_refine"].front();

	auto c = [](const char *t) -> condition
	{ return key("type") == t; };

	pdbFile << RM3("") << '\n'
			<< RM3(" DATA USED IN REFINEMENT.") << '\n'
			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_high") << '\n'
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_low") << '\n'
			<< RM3("  DATA CUTOFF            (SIGMA(F)) : ", 6, 3) << Ff(refine, "pdbx_ls_sigma_F") << '\n'
			<< RM3("  COMPLETENESS FOR RANGE        (%) : ", 5, 2) << Ff(refine, "ls_percent_reflns_obs") << '\n'
			<< RM3("  CROSS-VALIDATION METHOD           : ") << Fs(refine, "pdbx_ls_cross_valid_method") << '\n'
			<< RM3("  FREE R VALUE TEST SET SELECTION   : ") << Fs(refine, "pdbx_R_Free_selection_details") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT TO DATA USED IN REFINEMENT (NO CUTOFF).") << '\n'
			<< RM3("  R VALUE   (WORKING + TEST SET, NO CUTOFF) : ", 7, 3) << Ff(pdbx_refine, "R_factor_all_no_cutoff") << '\n'
			<< RM3("  R VALUE          (WORKING SET, NO CUTOFF) : ", 7, 3) << Ff(pdbx_refine, "R_factor_obs_no_cutoff") << '\n'
			<< RM3("  FREE R VALUE                  (NO CUTOFF) : ", 7, 3) << Ff(pdbx_refine, "free_R_factor_no_cutoff") << '\n'
			<< RM3("  FREE R VALUE TEST SET SIZE (%, NO CUTOFF) : ", 7, 3) << Ff(pdbx_refine, "free_R_val_test_set_size_perc_no_cutoff") << '\n'
			<< RM3("  FREE R VALUE TEST SET COUNT   (NO CUTOFF) : ", 12, 6) << Fi(pdbx_refine, "free_R_val_test_set_ct_no_cutoff") << '\n'
			<< RM3("  TOTAL NUMBER OF REFLECTIONS   (NO CUTOFF) : ", 12, 6) << Fi(refine, "ls_number_reflns_all") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT/AGREEMENT OF MODEL FOR DATA WITH F>4SIG(F).") << '\n'
			<< RM3("  R VALUE   (WORKING + TEST SET, F>4SIG(F)) : ", 7, 3) << Ff(pdbx_refine, "R_factor_all_4sig_cutoff") << '\n'
			<< RM3("  R VALUE          (WORKING SET, F>4SIG(F)) : ", 7, 3) << Ff(pdbx_refine, "R_factor_obs_4sig_cutoff") << '\n'
			<< RM3("  FREE R VALUE                  (F>4SIG(F)) : ", 7, 3) << Ff(pdbx_refine, "free_R_factor_4sig_cutoff") << '\n'
			<< RM3("  FREE R VALUE TEST SET SIZE (%, F>4SIG(F)) : ", 7, 3) << Ff(pdbx_refine, "free_R_val_test_set_size_perc_4sig_cutoff") << '\n'
			<< RM3("  FREE R VALUE TEST SET COUNT   (F>4SIG(F)) : ") << Fi(pdbx_refine, "free_R_val_test_set_ct_4sig_cutoff") << '\n'
			<< RM3("  TOTAL NUMBER OF REFLECTIONS   (F>4SIG(F)) : ") << Fi(pdbx_refine, "number_reflns_obs_4sig_cutoff") << '\n'

			<< RM3("") << '\n'
			<< RM3(" NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.") << '\n'
			<< RM3("  PROTEIN ATOMS      : ") << Fi(refine_hist, "pdbx_number_atoms_protein") << '\n'
			<< RM3("  NUCLEIC ACID ATOMS : ") << Fi(refine_hist, "pdbx_number_atoms_nucleic_acid") << '\n'
			<< RM3("  HETEROGEN ATOMS    : ") << Fi(refine_hist, "pdbx_number_atoms_ligand") << '\n'
			<< RM3("  SOLVENT ATOMS      : ") << Fi(refine_hist, "number_atoms_solvent") << '\n'

			<< RM3("") << '\n'
			<< RM3(" MODEL REFINEMENT.") << '\n'
			<< RM3("  OCCUPANCY SUM OF NON-HYDROGEN ATOMS      : ", 7, 3) << Ff(refine_analyze, "occupancy_sum_non_hydrogen") << '\n'
			<< RM3("  OCCUPANCY SUM OF HYDROGEN ATOMS          : ", 7, 3) << Ff(refine_analyze, "occupancy_sum_hydrogen") << '\n'
			<< RM3("  NUMBER OF DISCRETELY DISORDERED RESIDUES : ") << Fi(refine_analyze, "number_disordered_residues") << '\n'
			<< RM3("  NUMBER OF LEAST-SQUARES PARAMETERS       : ") << Fi(refine, "ls_number_parameters") << '\n'
			<< RM3("  NUMBER OF RESTRAINTS                     : ") << Fi(refine, "ls_number_restraints") << '\n'

			<< RM3("") << '\n'
			<< RM3(" RMS DEVIATIONS FROM RESTRAINT TARGET VALUES.") << '\n'
			<< RM3("  BOND LENGTHS                         (A) : ", 7, 3) << Ff(ls_restr, c("s_bond_d"), "dev_ideal") << '\n'
			<< RM3("  ANGLE DISTANCES                      (A) : ", 7, 3) << Ff(ls_restr, c("s_angle_d"), "dev_ideal") << '\n'
			<< RM3("  SIMILAR DISTANCES (NO TARGET VALUES) (A) : ", 7, 3) << Ff(ls_restr, c("s_similar_dist"), "dev_ideal") << '\n'
			<< RM3("  DISTANCES FROM RESTRAINT PLANES      (A) : ", 7, 3) << Ff(ls_restr, c("s_from_restr_planes"), "dev_ideal") << '\n'
			<< RM3("  ZERO CHIRAL VOLUMES               (A**3) : ", 7, 3) << Ff(ls_restr, c("s_zero_chiral_vol"), "dev_ideal") << '\n'
			<< RM3("  NON-ZERO CHIRAL VOLUMES           (A**3) : ", 7, 3) << Ff(ls_restr, c("s_non_zero_chiral_vol"), "dev_ideal") << '\n'
			<< RM3("  ANTI-BUMPING DISTANCE RESTRAINTS     (A) : ", 7, 3) << Ff(ls_restr, c("s_anti_bump_dis_restr"), "dev_ideal") << '\n'
			<< RM3("  RIGID-BOND ADP COMPONENTS         (A**2) : ", 7, 3) << Ff(ls_restr, c("s_rigid_bond_adp_cmpnt"), "dev_ideal") << '\n'
			<< RM3("  SIMILAR ADP COMPONENTS            (A**2) : ", 7, 3) << Ff(ls_restr, c("s_similar_adp_cmpnt"), "dev_ideal") << '\n'
			<< RM3("  APPROXIMATELY ISOTROPIC ADPS      (A**2) : ", 7, 3) << Ff(ls_restr, c("s_approx_iso_adps"), "dev_ideal") << '\n'

			<< RM3("") << '\n'
			<< RM3(" BULK SOLVENT MODELING.") << '\n'
			<< RM3("  METHOD USED: ") << Fs(refine, "solvent_model_details") << '\n'

			<< RM3("") << '\n'
			<< RM3(" STEREOCHEMISTRY TARGET VALUES : ") << Fs(refine, "pdbx_stereochemistry_target_values") << '\n'
			<< RM3("  SPECIAL CASE: ") << Fs(refine, "pdbx_stereochem_target_val_spec_case") << '\n'

			<< RM3("") << '\n';
}

void WriteRemark3Phenix(std::ostream &pdbFile, const datablock &db)
{
	auto refine = db["refine"].front();
	//	auto ls_shell = db["refine_ls_shell"].front();
	//	auto hist = db["refine_hist"].front();
	auto reflns = db["reflns"].front();
	//	auto analyze = db["refine_analyze"].front();
	auto &ls_restr = db["refine_ls_restr"];
	//	auto pdbx_xplor_file = db["pdbx_xplor_file"].front();
	auto pdbx_reflns_twin = db["pdbx_reflns_twin"].front();

	auto c = [](const char *t) -> condition
	{ return key("type") == t; };

	pdbFile << RM3("") << '\n'
			<< RM3("   REFINEMENT TARGET : ") << Fs(refine, "pdbx_stereochemistry_target_values") << '\n'
			<< RM3("") << '\n'
			<< RM3(" DATA USED IN REFINEMENT.") << '\n'
			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_high") << '\n'
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_low") << '\n'
			<< RM3("  MIN(FOBS/SIGMA_FOBS)              : ", 6, 3) << Ff(refine, "pdbx_ls_sigma_F") << '\n'
			<< RM3("  COMPLETENESS FOR RANGE        (%) : ", 5, 2) << Ff(refine, "ls_percent_reflns_obs") << '\n'
			<< RM3("  NUMBER OF REFLECTIONS             : ", 12, 6) << Fi(refine, "ls_number_reflns_obs") << '\n'
			<< RM3("") << '\n'
			<< RM3(" FIT TO DATA USED IN REFINEMENT.") << '\n'
			<< RM3("  R VALUE     (WORKING + TEST SET) : ", 7, 5) << Ff(refine, "ls_R_factor_obs") << '\n'
			<< RM3("  R VALUE            (WORKING SET) : ", 7, 5) << Ff(refine, "ls_R_factor_R_work") << '\n'
			<< RM3("  FREE R VALUE                     : ", 7, 5) << Ff(refine, "ls_R_factor_R_free") << '\n'
			<< RM3("  FREE R VALUE TEST SET SIZE   (%) : ", 7, 3) << Ff(refine, "ls_percent_reflns_R_free") << '\n'
			<< RM3("  FREE R VALUE TEST SET COUNT      : ", 12, 6) << Fi(refine, "ls_number_reflns_R_free") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT TO DATA USED IN REFINEMENT (IN BINS).") << '\n'
			<< RM3("  BIN  RESOLUTION RANGE  COMPL.    NWORK NFREE   RWORK  RFREE") << '\n';

	int bin = 1;
	std::vector<row_handle> bins;
	for (auto r : db["refine_ls_shell"])
		bins.push_back(r);
	//	reverse(bins.begin(), bins.end());
	try
	{
		sort(bins.begin(), bins.end(), [](row_handle a, row_handle b) -> bool
			{ return a["d_res_high"].as<float>() > b["d_res_high"].as<float>(); });
	}
	catch (...)
	{
	}

	for (auto r : bins)
	{
		float d_res_low, d_res_high, percent_reflns_obs, R_factor_R_work, R_factor_R_free;
		int number_reflns_R_work, number_reflns_R_free;

		tie(d_res_low, d_res_high, percent_reflns_obs, number_reflns_R_work,
			number_reflns_R_free, R_factor_R_work, R_factor_R_free) =
			r.get("d_res_low", "d_res_high", "percent_reflns_obs", "number_reflns_R_work",
				"number_reflns_R_free", "R_factor_R_work", "R_factor_R_free");

		percent_reflns_obs /= 100;

		pdbFile << RM3("  ") << cif::format("%3d %7.4f - %7.4f    %4.2f %8d %5d  %6.4f %6.4f", bin++, d_res_low, d_res_high, percent_reflns_obs, number_reflns_R_work, number_reflns_R_free, R_factor_R_work, R_factor_R_free) << '\n';
	}

	pdbFile << RM3("") << '\n'
			<< RM3(" BULK SOLVENT MODELLING.") << '\n'
			<< RM3("  METHOD USED        : ") << Fs(refine, "solvent_model_details") << '\n'
			<< RM3("  SOLVENT RADIUS     : ", 5, 2) << Ff(refine, "pdbx_solvent_vdw_probe_radii") << '\n'
			<< RM3("  SHRINKAGE RADIUS   : ", 5, 2) << Ff(refine, "pdbx_solvent_shrinkage_radii") << '\n'
			<< RM3("  K_SOL              : ", 5, 2) << Ff(refine, "solvent_model_param_ksol") << '\n'
			<< RM3("  B_SOL              : ", 5, 2) << Ff(refine, "solvent_model_param_bsol") << '\n'

			<< RM3("") << '\n'
			<< RM3(" ERROR ESTIMATES.") << '\n'
			<< RM3("  COORDINATE ERROR (MAXIMUM-LIKELIHOOD BASED)     : ", 6, 3) << Ff(refine, "overall_SU_ML") << '\n'
			<< RM3("  PHASE ERROR (DEGREES, MAXIMUM-LIKELIHOOD BASED) : ", 6, 3) << Ff(refine, "pdbx_overall_phase_error") << '\n'

			<< RM3("") << '\n'
			<< RM3(" B VALUES.") << '\n'
			<< RM3("  B VALUE TYPE                      : ") << Fs(refine, "pdbx_TLS_residual_ADP_flag") << '\n'
			<< RM3("  FROM WILSON PLOT           (A**2) : ", 7, 4) << Ff(reflns, "B_iso_Wilson_estimate") << '\n'
			<< RM3("  MEAN B VALUE      (OVERALL, A**2) : ", 7, 4) << Ff(refine, "B_iso_mean") << '\n'
			<< RM3("  OVERALL ANISOTROPIC B VALUE.") << '\n'
			<< RM3("   B11 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][1]") << '\n'
			<< RM3("   B22 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[2][2]") << '\n'
			<< RM3("   B33 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[3][3]") << '\n'
			<< RM3("   B12 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][2]") << '\n'
			<< RM3("   B13 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][3]") << '\n'
			<< RM3("   B23 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[2][3]") << '\n'

			<< RM3("") << '\n'
			<< RM3(" TWINNING INFORMATION.") << '\n'
			<< RM3("  FRACTION: ") << Fs(pdbx_reflns_twin, "fraction") << '\n'
			<< RM3("  OPERATOR: ") << Fs(pdbx_reflns_twin, "operator") << '\n'

			<< RM3("") << '\n'
			<< RM3(" DEVIATIONS FROM IDEAL VALUES.") << '\n'
			<< RM3("                RMSD          COUNT") << '\n'
			<< RM3("  BOND      : ", -6, 3) << Ff(ls_restr, c("f_bond_d"), "dev_ideal") << SEP("        ", -7)
			<< Fi(ls_restr, c("f_bond_d"), "number")
			<< '\n'
			<< RM3("  ANGLE     : ", -6, 3) << Ff(ls_restr, c("f_angle_d"), "dev_ideal") << SEP("        ", -7)
			<< Fi(ls_restr, c("f_angle_d"), "number")
			<< '\n'
			<< RM3("  CHIRALITY : ", -6, 3) << Ff(ls_restr, c("f_chiral_restr"), "dev_ideal") << SEP("        ", -7)
			<< Fi(ls_restr, c("f_chiral_restr"), "number")
			<< '\n'
			<< RM3("  PLANARITY : ", -6, 3) << Ff(ls_restr, c("f_plane_restr"), "dev_ideal") << SEP("        ", -7)
			<< Fi(ls_restr, c("f_plane_restr"), "number")
			<< '\n'
			<< RM3("  DIHEDRAL  : ", -6, 3) << Ff(ls_restr, c("f_dihedral_angle_d"), "dev_ideal") << SEP("        ", -7)
			<< Fi(ls_restr, c("f_dihedral_angle_d"), "number")
			<< '\n';

	auto &tls = db["pdbx_refine_tls"];

	pdbFile << RM3("") << '\n'
			<< RM3(" TLS DETAILS") << '\n'
			<< RM3("  NUMBER OF TLS GROUPS  : ") << (tls.size() ? std::to_string(tls.size()) : "NULL") << '\n';

	for (auto t : tls)
	{
		std::string id = t["id"].as<std::string>();

		auto pdbx_refine_tls_group = db["pdbx_refine_tls_group"].find_first(key("refine_tls_id") == id);

		pdbFile << RM3("  TLS GROUP : ") << id << '\n'
				<< RM3("   SELECTION: ") << Fs(pdbx_refine_tls_group, "selection_details") << '\n'
				<< RM3("   ORIGIN FOR THE GROUP (A):", -9, 4) << Ff(t, "origin_x")
				<< SEP("", -9, 4) << Ff(t, "origin_y")
				<< SEP("", -9, 4) << Ff(t, "origin_z") << '\n'
				<< RM3("   T TENSOR") << '\n'
				<< RM3("     T11:", -9, 4) << Ff(t, "T[1][1]") << SEP(" T22:", -9, 4) << Ff(t, "T[2][2]") << '\n'
				<< RM3("     T33:", -9, 4) << Ff(t, "T[3][3]") << SEP(" T12:", -9, 4) << Ff(t, "T[1][2]") << '\n'
				<< RM3("     T13:", -9, 4) << Ff(t, "T[1][3]") << SEP(" T23:", -9, 4) << Ff(t, "T[2][3]") << '\n'
				<< RM3("   L TENSOR") << '\n'
				<< RM3("     L11:", -9, 4) << Ff(t, "L[1][1]") << SEP(" L22:", -9, 4) << Ff(t, "L[2][2]") << '\n'
				<< RM3("     L33:", -9, 4) << Ff(t, "L[3][3]") << SEP(" L12:", -9, 4) << Ff(t, "L[1][2]") << '\n'
				<< RM3("     L13:", -9, 4) << Ff(t, "L[1][3]") << SEP(" L23:", -9, 4) << Ff(t, "L[2][3]") << '\n'
				<< RM3("   S TENSOR") << '\n'
				<< RM3("     S11:", -9, 4) << Ff(t, "S[1][1]") << SEP(" S12:", -9, 4) << Ff(t, "S[1][2]") << SEP(" S13:", -9, 4) << Ff(t, "S[1][3]") << '\n'
				<< RM3("     S21:", -9, 4) << Ff(t, "S[2][1]") << SEP(" S22:", -9, 4) << Ff(t, "S[2][2]") << SEP(" S23:", -9, 4) << Ff(t, "S[2][3]") << '\n'
				<< RM3("     S31:", -9, 4) << Ff(t, "S[3][1]") << SEP(" S32:", -9, 4) << Ff(t, "S[3][2]") << SEP(" S33:", -9, 4) << Ff(t, "S[3][3]") << '\n';
	}

	pdbFile << RM3("") << '\n'
			<< RM3(" NCS DETAILS") << '\n';

	auto &ncs_dom = db["struct_ncs_dom"];
	if (ncs_dom.empty())
		pdbFile << RM3("  NUMBER OF NCS GROUPS : NULL") << '\n';
	else
	{
		std::set<std::string> ncs_groups;
		for (auto i : ncs_dom)
			ncs_groups.insert(i["pdbx_ens_id"].as<std::string>());

		pdbFile << RM3("  NUMBER OF NCS GROUPS : ") << ncs_groups.size() << '\n';
		//
		//			for (auto ens_id: ncs_groups)
		//			{
		//				auto lim = db["struct_ncs_dom_lim"].find(key("pdbx_ens_id") == ens_id);
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
		//				pdbFile << RM3("") << '\n'
		//						<< RM3(" NCS GROUP NUMBER               : ") << ens_id << '\n'
		//						<< RM3("    CHAIN NAMES                    : ") << join(chains, " ") << '\n'
		//						<< RM3("    NUMBER OF COMPONENTS NCS GROUP : ") << component_ids.size() << '\n'
		//						<< RM3("      COMPONENT C  SSSEQI  TO  C   SSSEQI   CODE") << '\n';
		//
		//				for (auto l: lim)
		//				{
		//					pdbFile << RM3("         ", -2)		<< Fi(l, "pdbx_component_id")
		//							<< SEP(" ", -5)			<< Fs(l, "beg_auth_asym_id")
		//							<< SEP("  ", -5)			<< Fi(l, "beg_auth_seq_id")
		//							<< SEP("   ", -5)			<< Fs(l, "end_auth_asym_id")
		//							<< SEP("   ", -5)			<< Fi(l, "end_auth_seq_id")
		//							<< SEP("  ", -5)			<< Fs(l, "pdbx_refine_code")
		//							<< '\n';
		//				}
		//
		//				pdbFile << RM3("                  GROUP CHAIN        COUNT   RMS     WEIGHT") << '\n';
		//				for (auto l: db["refine_ls_restr_ncs"].find(key("pdbx_ens_id") == ens_id))
		//				{
		//					std::string type = l["pdbx_type"];
		//					to_upper(type);
		//
		//					std::string unit;
		//					if (ends_with(type, "POSITIONAL"))
		//						unit = "    (A): ";
		//					else if (ends_with(type, "THERMAL"))
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
		//							<< '\n';
		//				}
		//			}
	}

	//	pdbFile << RM3("") << '\n'
	//			<< RM3(" BULK SOLVENT MODELLING.") << '\n'
	//			<< RM3("  METHOD USED : ") << Fs(refine, "solvent_model_details") << '\n'
	//			<< RM3("  PARAMETERS FOR MASK CALCULATION") << '\n'
	//			<< RM3("  VDW PROBE RADIUS   : ", 5, 2) << Ff(refine, "pdbx_solvent_vdw_probe_radii") << '\n'
	//			<< RM3("  ION PROBE RADIUS   : ", 5, 2) << Ff(refine, "pdbx_solvent_ion_probe_radii") << '\n'
	//			<< RM3("  SHRINKAGE RADIUS   : ", 5, 2) << Ff(refine, "pdbx_solvent_shrinkage_radii") << '\n'
	//
	//			<< RM3("") << '\n';

	pdbFile << RM3("") << '\n';
}

void WriteRemark3XPlor(std::ostream &pdbFile, const datablock &db)
{
	auto refine = db["refine"].front();
	auto ls_shell = db["refine_ls_shell"].front();
	auto hist = db["refine_hist"].front();
	auto reflns = db["reflns"].front();
	auto analyze = db["refine_analyze"].front();
	auto &ls_restr = db["refine_ls_restr"];
	auto ls_restr_ncs = db["refine_ls_restr_ncs"].front();
	auto pdbx_xplor_file = db["pdbx_xplor_file"].front();

	pdbFile << RM3("") << '\n'
			<< RM3(" DATA USED IN REFINEMENT.") << '\n'
			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_high") << '\n'
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_low") << '\n'
			<< RM3("  DATA CUTOFF            (SIGMA(F)) : ", 6, 3) << Ff(refine, "pdbx_ls_sigma_F") << '\n'
			<< RM3("  DATA CUTOFF HIGH         (ABS(F)) : ", 6, 3) << Ff(refine, "pdbx_data_cutoff_high_absF") << '\n'
			<< RM3("  DATA CUTOFF LOW          (ABS(F)) : ", 6, 3) << Ff(refine, "pdbx_data_cutoff_low_absF") << '\n'
			<< RM3("  COMPLETENESS (WORKING+TEST)   (%) : ", 5, 2) << Ff(refine, "ls_percent_reflns_obs") << '\n'
			<< RM3("  NUMBER OF REFLECTIONS             : ", 12, 6) << Fi(refine, "ls_number_reflns_obs") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT TO DATA USED IN REFINEMENT.") << '\n'
			<< RM3("  CROSS-VALIDATION METHOD          : ") << Fs(refine, "pdbx_ls_cross_valid_method") << '\n'
			<< RM3("  FREE R VALUE TEST SET SELECTION  : ") << Fs(refine, "pdbx_R_Free_selection_details") << '\n'
			<< RM3("  R VALUE            (WORKING SET) : ", 7, 3) << Ff(refine, "ls_R_factor_R_work") << '\n'
			<< RM3("  FREE R VALUE                     : ", 7, 3) << Ff(refine, "ls_R_factor_R_free") << '\n'
			<< RM3("  FREE R VALUE TEST SET SIZE   (%) : ", 7, 3) << Ff(refine, "ls_percent_reflns_R_free") << '\n'
			<< RM3("  FREE R VALUE TEST SET COUNT      : ", 12, 6) << Fi(refine, "ls_number_reflns_R_free") << '\n'
			<< RM3("  ESTIMATED ERROR OF FREE R VALUE  : ", 7, 3) << Ff(refine, "ls_R_factor_R_free_error") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT IN THE HIGHEST RESOLUTION BIN.") << '\n'
			<< RM3("  TOTAL NUMBER OF BINS USED           : ", 12, 6) << Fi(ls_shell, "pdbx_total_number_of_bins_used") << '\n'
			<< RM3("  BIN RESOLUTION RANGE HIGH       (A) : ", 5, 2) << Ff(ls_shell, "d_res_high") << '\n'
			<< RM3("  BIN RESOLUTION RANGE LOW        (A) : ", 5, 2) << Ff(ls_shell, "d_res_low") << '\n'
			<< RM3("  BIN COMPLETENESS (WORKING+TEST) (%) : ", 5, 1) << Ff(ls_shell, "percent_reflns_obs") << '\n'
			<< RM3("  REFLECTIONS IN BIN    (WORKING SET) : ", 12, 6) << Fi(ls_shell, "number_reflns_R_work") << '\n'
			<< RM3("  BIN R VALUE           (WORKING SET) : ", 7, 3) << Ff(ls_shell, "R_factor_R_work") << '\n'
			<< RM3("  BIN FREE R VALUE                    : ", 7, 3) << Ff(ls_shell, "R_factor_R_free") << '\n'
			<< RM3("  BIN FREE R VALUE TEST SET SIZE  (%) : ", 5, 1) << Ff(ls_shell, "percent_reflns_R_free") << '\n'
			<< RM3("  BIN FREE R VALUE TEST SET COUNT     : ", 12, 6) << Fi(ls_shell, "number_reflns_R_free") << '\n'
			<< RM3("  ESTIMATED ERROR OF BIN FREE R VALUE : ", 7, 3) << Ff(ls_shell, "R_factor_R_free_error") << '\n'

			<< RM3("") << '\n'
			<< RM3(" NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.") << '\n'
			<< RM3("  PROTEIN ATOMS            : ", 12, 6) << Fi(hist, "pdbx_number_atoms_protein") << '\n'
			<< RM3("  NUCLEIC ACID ATOMS       : ", 12, 6) << Fi(hist, "pdbx_number_atoms_nucleic_acid") << '\n'
			<< RM3("  HETEROGEN ATOMS          : ", 12, 6) << Fi(hist, "pdbx_number_atoms_ligand") << '\n'
			<< RM3("  SOLVENT ATOMS            : ", 12, 6) << Fi(hist, "number_atoms_solvent") << '\n'

			<< RM3("") << '\n'
			<< RM3(" B VALUES.") << '\n'
			<< RM3("  FROM WILSON PLOT           (A**2) : ", 7, 2) << Ff(reflns, "B_iso_Wilson_estimate") << '\n'
			<< RM3("  MEAN B VALUE      (OVERALL, A**2) : ", 7, 2) << Ff(refine, "B_iso_mean") << '\n'

			<< RM3("  OVERALL ANISOTROPIC B VALUE.") << '\n'
			<< RM3("   B11 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][1]") << '\n'
			<< RM3("   B22 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[2][2]") << '\n'
			<< RM3("   B33 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[3][3]") << '\n'
			<< RM3("   B12 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][2]") << '\n'
			<< RM3("   B13 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][3]") << '\n'
			<< RM3("   B23 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[2][3]") << '\n'

			<< RM3("") << '\n'
			<< RM3(" ESTIMATED COORDINATE ERROR.") << '\n'
			<< RM3("  ESD FROM LUZZATI PLOT        (A) : ", 7, 2) << Ff(analyze, "Luzzati_coordinate_error_obs") << '\n'
			<< RM3("  ESD FROM SIGMAA              (A) : ", 7, 2) << Ff(analyze, "Luzzati_sigma_a_obs") << '\n'
			<< RM3("  LOW RESOLUTION CUTOFF        (A) : ", 7, 2) << Ff(analyze, "Luzzati_d_res_low_obs") << '\n'

			<< RM3("") << '\n'
			<< RM3(" CROSS-VALIDATED ESTIMATED COORDINATE ERROR.") << '\n'
			<< RM3("  ESD FROM C-V LUZZATI PLOT    (A) : ", 7, 2) << Ff(analyze, "Luzzati_coordinate_error_free") << '\n'
			<< RM3("  ESD FROM C-V SIGMAA          (A) : ", 7, 2) << Ff(analyze, "Luzzati_sigma_a_free") << '\n'

			<< RM3("") << '\n'
			<< RM3(" RMS DEVIATIONS FROM IDEAL VALUES.") << '\n'
			<< RM3("  BOND LENGTHS                 (A) : ", 7, 3) << Ff(ls_restr, key("type") == "x_bond_d", "dev_ideal") << '\n'
			<< RM3("  BOND ANGLES            (DEGREES) : ", 7, 2) << Ff(ls_restr, key("type") == "x_angle_deg", "dev_ideal") << '\n'
			<< RM3("  DIHEDRAL ANGLES        (DEGREES) : ", 7, 2) << Ff(ls_restr, key("type") == "x_dihedral_angle_d", "dev_ideal") << '\n'
			<< RM3("  IMPROPER ANGLES        (DEGREES) : ", 7, 2) << Ff(ls_restr, key("type") == "x_improper_angle_d", "dev_ideal") << '\n'

			<< RM3("") << '\n'
			<< RM3(" ISOTROPIC THERMAL MODEL : ") << Fs(refine, "pdbx_isotropic_thermal_model") << '\n'

			<< RM3("") << '\n'
			<< RM3(" ISOTROPIC THERMAL FACTOR RESTRAINTS.    RMS    SIGMA") << '\n'
			<< RM3("  MAIN-CHAIN BOND              (A**2) : ", 6, 2) << Ff(ls_restr, key("type") == "x_mcbond_it", "dev_ideal") << SEP("; ", 6, 2)
			<< Ff(ls_restr, key("type") == "x_mcbond_it", "dev_ideal_target") << '\n'
			<< RM3("  MAIN-CHAIN ANGLE             (A**2) : ", 6, 2) << Ff(ls_restr, key("type") == "x_mcangle_it", "dev_ideal") << SEP("; ", 6, 2)
			<< Ff(ls_restr, key("type") == "x_mcangle_it", "dev_ideal_target") << '\n'
			<< RM3("  SIDE-CHAIN BOND              (A**2) : ", 6, 2) << Ff(ls_restr, key("type") == "x_scbond_it", "dev_ideal") << SEP("; ", 6, 2)
			<< Ff(ls_restr, key("type") == "x_scbond_it", "dev_ideal_target") << '\n'
			<< RM3("  SIDE-CHAIN ANGLE             (A**2) : ", 6, 2) << Ff(ls_restr, key("type") == "x_scangle_it", "dev_ideal") << SEP("; ", 6, 2)
			<< Ff(ls_restr, key("type") == "x_scangle_it", "dev_ideal_target") << '\n'
			<< RM3("") << '\n'
			<< RM3(" NCS MODEL : ") << Fs(ls_restr_ncs, "ncs_model_details") << '\n'

			<< RM3("") << '\n'
			<< RM3(" NCS RESTRAINTS.                         RMS   SIGMA/WEIGHT") << '\n'

			// TODO: using only group 1 here, should this be fixed???
			<< RM3("  GROUP  1  POSITIONAL            (A) : ", 4, 2) << Ff(ls_restr_ncs, "rms_dev_position") << SEP("; ", 6, 2)
			<< Ff(ls_restr_ncs, "weight_position") << SEP("; ", 6, 2) << '\n'
			<< RM3("  GROUP  1  B-FACTOR           (A**2) : ", 4, 2) << Ff(ls_restr_ncs, "rms_dev_B_iso") << SEP("; ", 6, 2)
			<< Ff(ls_restr_ncs, "weight_B_iso") << SEP("; ", 6, 2) << '\n'

			// TODO: using only files from serial_no 1 here
			<< RM3("") << '\n'
			<< RM3(" PARAMETER FILE   1  : ") << Fs(pdbx_xplor_file, "param_file") << '\n'
			<< RM3(" TOPOLOGY FILE   1   : ") << Fs(pdbx_xplor_file, "topol_file") << '\n'

			<< RM3("") << '\n';
}

void WriteRemark3NuclSQ(std::ostream &pdbFile, const datablock &db)
{
	auto refine = db["refine"].front();
	auto pdbx_refine = db["pdbx_refine"].front();
	auto hist = db["refine_hist"].front();
	auto reflns = db["reflns"].front();
	auto analyze = db["refine_analyze"].front();
	auto &ls_restr = db["refine_ls_restr"];

	pdbFile << RM3("") << '\n'
			<< RM3(" DATA USED IN REFINEMENT.") << '\n'

			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_high") << '\n'
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_low") << '\n'
			<< RM3("  DATA CUTOFF            (SIGMA(F)) : ", 6, 3) << Ff(refine, "pdbx_ls_sigma_F") << '\n'
			<< RM3("  COMPLETENESS FOR RANGE        (%) : ", 5, 2) << Ff(refine, "ls_percent_reflns_obs") << '\n'
			<< RM3("  NUMBER OF REFLECTIONS             : ", 12, 6) << Fi(refine, "ls_number_reflns_obs") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT TO DATA USED IN REFINEMENT.") << '\n'
			<< RM3("  CROSS-VALIDATION METHOD          : ") << Fs(refine, "pdbx_ls_cross_valid_method") << '\n'
			<< RM3("  FREE R VALUE TEST SET SELECTION  : ") << Fs(refine, "pdbx_R_Free_selection_details") << '\n'
			<< RM3("  R VALUE     (WORKING + TEST SET) : ", 7, 3) << Ff(refine, "ls_R_factor_obs") << '\n'
			<< RM3("  R VALUE            (WORKING SET) : ", 7, 3) << Ff(refine, "ls_R_factor_R_work") << '\n'
			<< RM3("  FREE R VALUE                     : ", 7, 3) << Ff(refine, "ls_R_factor_R_free") << '\n'
			<< RM3("  FREE R VALUE TEST SET SIZE   (%) : ", 7, 3) << Ff(refine, "ls_percent_reflns_R_free") << '\n'
			<< RM3("  FREE R VALUE TEST SET COUNT      : ", 12, 6) << Fi(refine, "ls_number_reflns_R_free") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT/AGREEMENT OF MODEL WITH ALL DATA.") << '\n'
			<< RM3("  R VALUE   (WORKING + TEST SET, NO CUTOFF) : ") << Fs(refine, "ls_R_factor_all") << '\n'
			<< RM3("  R VALUE          (WORKING SET, NO CUTOFF) : ") << Fs(pdbx_refine, "R_factor_obs_no_cutoff") << '\n'

			<< RM3("  FREE R VALUE                  (NO CUTOFF) : ") << Fs(pdbx_refine, "free_R_factor_no_cutoff") << '\n'
			<< RM3("  FREE R VALUE TEST SET SIZE (%, NO CUTOFF) : ") << Fs(pdbx_refine, "free_R_val_test_set_size_perc_no_cutoff") << '\n'
			<< RM3("  FREE R VALUE TEST SET COUNT   (NO CUTOFF) : ") << Fs(pdbx_refine, "free_R_val_test_set_ct_no_cutoff") << '\n'
			<< RM3("  TOTAL NUMBER OF REFLECTIONS   (NO CUTOFF) : ") << Fs(refine, "ls_number_reflns_all") << '\n'

			<< RM3("") << '\n'
			<< RM3(" NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.") << '\n'
			<< RM3("  PROTEIN ATOMS            : ", 12, 6) << Fi(hist, "pdbx_number_atoms_protein") << '\n'
			<< RM3("  NUCLEIC ACID ATOMS       : ", 12, 6) << Fi(hist, "pdbx_number_atoms_nucleic_acid") << '\n'
			<< RM3("  HETEROGEN ATOMS          : ", 12, 6) << Fi(hist, "pdbx_number_atoms_ligand") << '\n'
			<< RM3("  SOLVENT ATOMS            : ", 12, 6) << Fi(hist, "number_atoms_solvent") << '\n'
			//			<< RM3("  ALL ATOMS                : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_protein") << '\n'

			<< RM3("") << '\n'
			<< RM3(" B VALUES.") << '\n'
			//			<< RM3("  B VALUE TYPE                      : ", 7, 2)	<< Fs(refine, "pdbx_TLS_residual_ADP_flag") << '\n'
			<< RM3("  FROM WILSON PLOT           (A**2) : ", 7, 2) << Ff(reflns, "B_iso_Wilson_estimate") << '\n'
			<< RM3("  MEAN B VALUE      (OVERALL, A**2) : ", 7, 2) << Ff(refine, "B_iso_mean") << '\n'
			<< RM3("  OVERALL ANISOTROPIC B VALUE.") << '\n'
			<< RM3("   B11 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][1]") << '\n'
			<< RM3("   B22 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[2][2]") << '\n'
			<< RM3("   B33 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[3][3]") << '\n'
			<< RM3("   B12 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][2]") << '\n'
			<< RM3("   B13 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][3]") << '\n'
			<< RM3("   B23 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[2][3]") << '\n'

			<< RM3("") << '\n'
			<< RM3(" ESTIMATED COORDINATE ERROR.") << '\n'
			<< RM3("  ESD FROM LUZZATI PLOT        (A) : ", 7, 2) << Ff(analyze, "Luzzati_coordinate_error_obs") << '\n'
			<< RM3("  ESD FROM SIGMAA              (A) : ", 7, 2) << Ff(analyze, "Luzzati_sigma_a_obs") << '\n'
			<< RM3("  LOW RESOLUTION CUTOFF        (A) : ", 7, 2) << Ff(analyze, "Luzzati_d_res_low_obs") << '\n'

			<< RM3("") << '\n'
			<< RM3(" RMS DEVIATIONS FROM IDEAL VALUES.") << '\n'
			<< RM3("  DISTANCE RESTRAINTS.                    RMS     SIGMA") << '\n'
			<< RM3("   SUGAR-BASE BOND DISTANCE        (A) : ", 7, 3) << Ff(ls_restr, key("type") == "n_sugar_bond_d", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "n_sugar_bond_d", "dev_ideal_target") << '\n'
			<< RM3("   SUGAR-BASE BOND ANGLE DISTANCE  (A) : ", 7, 3) << Ff(ls_restr, key("type") == "n_sugar_bond_angle_d", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "n_sugar_bond_angle_d", "dev_ideal_target") << '\n'
			<< RM3("   PHOSPHATE BONDS DISTANCE        (A) : ", 7, 3) << Ff(ls_restr, key("type") == "n_phos_bond_d", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "n_phos_bond_d", "dev_ideal_target") << '\n'
			<< RM3("   PHOSPHATE BOND ANGLE, H-BOND    (A) : ", 7, 3) << Ff(ls_restr, key("type") == "n_phos_bond_angle_d", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "n_phos_bond_angle_d", "dev_ideal_target") << '\n'

			<< RM3("") << '\n'
			<< RM3("  PLANE RESTRAINT                  (A) : ", 7, 3) << Ff(ls_restr, key("type") == "n_plane_restr", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "n_plane_restr", "dev_ideal_target") << '\n'
			<< RM3("  CHIRAL-CENTER RESTRAINT       (A**3) : ", 7, 3) << Ff(ls_restr, key("type") == "n_chiral_restr", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "n_chiral_restr", "dev_ideal_target") << '\n'

			<< RM3("") << '\n'
			<< RM3("  NON-BONDED CONTACT RESTRAINTS.") << '\n'
			<< RM3("   SINGLE TORSION CONTACT          (A) : ", 7, 3) << Ff(ls_restr, key("type") == "n_singtor_nbd", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "n_singtor_nbd", "dev_ideal_target") << '\n'
			<< RM3("   MULTIPLE TORSION CONTACT        (A) : ", 7, 3) << Ff(ls_restr, key("type") == "n_multtor_nbd", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "n_multtor_nbd", "dev_ideal_target") << '\n'

			<< RM3("") << '\n'
			<< RM3(" ISOTROPIC THERMAL FACTOR RESTRAINTS.    RMS     SIGMA") << '\n'
			<< RM3("  SUGAR-BASE BONDS             (A**2) : ", 7, 3) << Ff(ls_restr, key("type") == "n_sugar_bond_it", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "n_sugar_bond_it", "dev_ideal_target") << '\n'
			<< RM3("  SUGAR-BASE ANGLES            (A**2) : ", 7, 3) << Ff(ls_restr, key("type") == "n_sugar_angle_it", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "n_sugar_angle_it", "dev_ideal_target") << '\n'
			<< RM3("  PHOSPHATE BONDS              (A**2) : ", 7, 3) << Ff(ls_restr, key("type") == "n_phos_bond_it", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "n_phos_bond_it", "dev_ideal_target") << '\n'
			<< RM3("  PHOSPHATE BOND ANGLE, H-BOND (A**2) : ", 7, 3) << Ff(ls_restr, key("type") == "n_phos_angle_it", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "n_phos_angle_it", "dev_ideal_target") << '\n'

			<< RM3("") << '\n';
}

void WriteRemark3ProlSQ(std::ostream &pdbFile, const datablock &db)
{
	auto refine = db["refine"].front();
	auto pdbx_refine = db["pdbx_refine"].front();
	auto hist = db["refine_hist"].front();
	auto reflns = db["reflns"].front();
	auto analyze = db["refine_analyze"].front();
	auto &ls_restr = db["refine_ls_restr"];

	pdbFile << RM3("") << '\n'
			<< RM3(" DATA USED IN REFINEMENT.") << '\n'

			<< RM3("  RESOLUTION RANGE HIGH (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_high") << '\n'
			<< RM3("  RESOLUTION RANGE LOW  (ANGSTROMS) : ", 5, 2) << Ff(refine, "ls_d_res_low") << '\n'
			<< RM3("  DATA CUTOFF            (SIGMA(F)) : ", 6, 3) << Ff(refine, "pdbx_ls_sigma_F") << '\n'
			<< RM3("  COMPLETENESS FOR RANGE        (%) : ", 5, 2) << Ff(refine, "ls_percent_reflns_obs") << '\n'
			<< RM3("  NUMBER OF REFLECTIONS             : ", 12, 6) << Fi(refine, "ls_number_reflns_obs") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT TO DATA USED IN REFINEMENT.") << '\n'
			<< RM3("  CROSS-VALIDATION METHOD          : ") << Fs(refine, "pdbx_ls_cross_valid_method") << '\n'
			<< RM3("  FREE R VALUE TEST SET SELECTION  : ") << Fs(refine, "pdbx_R_Free_selection_details") << '\n'
			<< RM3("  R VALUE     (WORKING + TEST SET) : ", 7, 3) << Ff(refine, "ls_R_factor_obs") << '\n'
			<< RM3("  R VALUE            (WORKING SET) : ", 7, 3) << Ff(refine, "ls_R_factor_R_work") << '\n'
			<< RM3("  FREE R VALUE                     : ", 7, 3) << Ff(refine, "ls_R_factor_R_free") << '\n'
			<< RM3("  FREE R VALUE TEST SET SIZE   (%) : ", 7, 3) << Ff(refine, "ls_percent_reflns_R_free") << '\n'
			<< RM3("  FREE R VALUE TEST SET COUNT      : ", 12, 6) << Fi(refine, "ls_number_reflns_R_free") << '\n'

			<< RM3("") << '\n'
			<< RM3(" FIT/AGREEMENT OF MODEL WITH ALL DATA.") << '\n'
			<< RM3("  R VALUE   (WORKING + TEST SET, NO CUTOFF) : ") << Fs(refine, "ls_R_factor_all") << '\n'
			<< RM3("  R VALUE          (WORKING SET, NO CUTOFF) : ") << Fs(pdbx_refine, "R_factor_obs_no_cutoff") << '\n'

			<< RM3("  FREE R VALUE                  (NO CUTOFF) : ") << Fs(pdbx_refine, "free_R_factor_no_cutoff") << '\n'
			<< RM3("  FREE R VALUE TEST SET SIZE (%, NO CUTOFF) : ") << Fs(pdbx_refine, "free_R_val_test_set_size_perc_no_cutoff") << '\n'
			<< RM3("  FREE R VALUE TEST SET COUNT   (NO CUTOFF) : ") << Fs(pdbx_refine, "free_R_val_test_set_ct_no_cutoff") << '\n'
			<< RM3("  TOTAL NUMBER OF REFLECTIONS   (NO CUTOFF) : ") << Fs(refine, "ls_number_reflns_all") << '\n'

			<< RM3("") << '\n'
			<< RM3(" NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.") << '\n'
			<< RM3("  PROTEIN ATOMS            : ", 12, 6) << Fi(hist, "pdbx_number_atoms_protein") << '\n'
			<< RM3("  NUCLEIC ACID ATOMS       : ", 12, 6) << Fi(hist, "pdbx_number_atoms_nucleic_acid") << '\n'
			<< RM3("  HETEROGEN ATOMS          : ", 12, 6) << Fi(hist, "pdbx_number_atoms_ligand") << '\n'
			<< RM3("  SOLVENT ATOMS            : ", 12, 6) << Fi(hist, "number_atoms_solvent") << '\n'
			//			<< RM3("  ALL ATOMS                : ", 12, 6)	<< Fi(hist, "pdbx_number_atoms_protein") << '\n'

			<< RM3("") << '\n'
			<< RM3(" B VALUES.") << '\n'
			//			<< RM3("  B VALUE TYPE                      : ", 7, 2)	<< Fs(refine, "pdbx_TLS_residual_ADP_flag") << '\n'
			<< RM3("  FROM WILSON PLOT           (A**2) : ", 7, 2) << Ff(reflns, "B_iso_Wilson_estimate") << '\n'
			<< RM3("  MEAN B VALUE      (OVERALL, A**2) : ", 7, 2) << Ff(refine, "B_iso_mean") << '\n'
			<< RM3("  OVERALL ANISOTROPIC B VALUE.") << '\n'
			<< RM3("   B11 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][1]") << '\n'
			<< RM3("   B22 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[2][2]") << '\n'
			<< RM3("   B33 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[3][3]") << '\n'
			<< RM3("   B12 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][2]") << '\n'
			<< RM3("   B13 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[1][3]") << '\n'
			<< RM3("   B23 (A**2) : ", -7, 2) << Ff(refine, "aniso_B[2][3]") << '\n'

			<< RM3("") << '\n'
			<< RM3(" ESTIMATED COORDINATE ERROR.") << '\n'
			<< RM3("  ESD FROM LUZZATI PLOT        (A) : ", 7, 2) << Ff(analyze, "Luzzati_coordinate_error_obs") << '\n'
			<< RM3("  ESD FROM SIGMAA              (A) : ", 7, 2) << Ff(analyze, "Luzzati_sigma_a_obs") << '\n'
			<< RM3("  LOW RESOLUTION CUTOFF        (A) : ", 7, 2) << Ff(analyze, "Luzzati_d_res_low_obs") << '\n'

			<< RM3("") << '\n'
			<< RM3(" RMS DEVIATIONS FROM IDEAL VALUES.") << '\n'
			<< RM3("  DISTANCE RESTRAINTS.                    RMS    SIGMA") << '\n'
			<< RM3("   BOND LENGTH                     (A) : ", 7, 3) << Ff(ls_restr, key("type") == "p_bond_d", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_bond_d", "dev_ideal_target") << '\n'
			<< RM3("   ANGLE DISTANCE                  (A) : ", 7, 3) << Ff(ls_restr, key("type") == "p_angle_d", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_angle_d", "dev_ideal_target") << '\n'
			<< RM3("   INTRAPLANAR 1-4 DISTANCE        (A) : ", 7, 3) << Ff(ls_restr, key("type") == "p_planar_d", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_planar_d", "dev_ideal_target") << '\n'
			<< RM3("   H-BOND OR METAL COORDINATION    (A) : ", 7, 3) << Ff(ls_restr, key("type") == "p_hb_or_metal_coord", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_hb_or_metal_coord", "dev_ideal_target") << '\n'

			<< RM3("") << '\n'
			<< RM3("  PLANE RESTRAINT                 (A) : ", 7, 3) << Ff(ls_restr, key("type") == "p_plane_restr", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_plane_restr", "dev_ideal_target") << '\n'
			<< RM3("  CHIRAL-CENTER RESTRAINT      (A**3) : ", 7, 3) << Ff(ls_restr, key("type") == "p_chiral_restr", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_chiral_restr", "dev_ideal_target") << '\n'

			<< RM3("") << '\n'
			<< RM3("  NON-BONDED CONTACT RESTRAINTS.") << '\n'
			<< RM3("   SINGLE TORSION                  (A) : ", 7, 3) << Ff(ls_restr, key("type") == "p_singtor_nbd", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_singtor_nbd", "dev_ideal_target") << '\n'
			<< RM3("   MULTIPLE TORSION                (A) : ", 7, 3) << Ff(ls_restr, key("type") == "p_multtor_nbd", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_multtor_nbd", "dev_ideal_target") << '\n'
			<< RM3("   H-BOND (X...Y)                  (A) : ", 7, 3) << Ff(ls_restr, key("type") == "p_xyhbond_nbd", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_xyhbond_nbd", "dev_ideal_target") << '\n'
			<< RM3("   H-BOND (X-H...Y)                (A) : ", 7, 3) << Ff(ls_restr, key("type") == "p_xhyhbond_nbd", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_xhyhbond_nbd", "dev_ideal_target") << '\n'

			<< RM3("") << '\n'
			<< RM3("  CONFORMATIONAL TORSION ANGLE RESTRAINTS.") << '\n'
			<< RM3("   SPECIFIED                 (DEGREES) : ", 7, 3) << Ff(ls_restr, key("type") == "p_special_tor", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_special_tor", "dev_ideal_target") << '\n'
			<< RM3("   PLANAR                    (DEGREES) : ", 7, 3) << Ff(ls_restr, key("type") == "p_planar_tor", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_planar_tor", "dev_ideal_target") << '\n'
			<< RM3("   STAGGERED                 (DEGREES) : ", 7, 3) << Ff(ls_restr, key("type") == "p_staggered_tor", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_staggered_tor", "dev_ideal_target") << '\n'
			<< RM3("   TRANSVERSE                (DEGREES) : ", 7, 3) << Ff(ls_restr, key("type") == "p_transverse_tor", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_transverse_tor", "dev_ideal_target") << '\n'

			<< RM3("") << '\n'
			<< RM3("  ISOTROPIC THERMAL FACTOR RESTRAINTS. RMS SIGMA") << '\n'
			<< RM3("   MAIN-CHAIN BOND              (A**2) : ", 7, 3) << Ff(ls_restr, key("type") == "p_mcbond_it", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_mcbond_it", "dev_ideal_target") << '\n'
			<< RM3("   MAIN-CHAIN ANGLE             (A**2) : ", 7, 3) << Ff(ls_restr, key("type") == "p_mcangle_it", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_mcangle_it", "dev_ideal_target") << '\n'
			<< RM3("   SIDE-CHAIN BOND              (A**2) : ", 7, 3) << Ff(ls_restr, key("type") == "p_scbond_it", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_scbond_it", "dev_ideal_target") << '\n'
			<< RM3("   SIDE-CHAIN ANGLE             (A**2) : ", 7, 3) << Ff(ls_restr, key("type") == "p_scangle_it", "dev_ideal") << " ; "
			<< Ff(ls_restr, key("type") == "p_scangle_it", "dev_ideal_target") << '\n'

			<< RM3("") << '\n';
}

void WriteRemark3(std::ostream &pdbFile, const datablock &db)
{
	std::string program, authors;

	if (not db["pdbx_nmr_software"].empty())
	{
		auto software = db["pdbx_nmr_software"].find(key("classification") == "refinement");
		if (software.size() == 1)
			cif::tie(program, authors) = software.front().get("name", "authors");
		else if (software.size() > 1)
		{
			for (auto r : software)
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
		pdbFile << RM3("") << '\n'
				<< RM3("REFINEMENT.") << '\n';

		int l = 0;
		for (auto s : word_wrap(program, 52))
			pdbFile << RM3(++l == 1 ? "  PROGRAM     : " : "                ") << s << '\n';

		l = 0;
		for (auto s : word_wrap(authors, 52))
			pdbFile << RM3(++l == 1 ? "  AUTHORS     : " : "                ") << s << '\n';
	}

	if (not db["refine"].empty())
	{
		auto s = program.find(' ');
		if (s != std::string::npos)
			program.erase(s, std::string::npos);

		if (iequals(program, "BUSTER") or iequals(program, "BUSTER-TNT") or iequals(program, "TNT"))
			WriteRemark3BusterTNT(pdbFile, db);
		else if (iequals(program, "CNS") or iequals(program, "CNX"))
			WriteRemark3CNS(pdbFile, db);
		else if (iequals(program, "X-PLOR"))
			WriteRemark3XPlor(pdbFile, db);
		else if (iequals(program, "REFMAC"))
			WriteRemark3Refmac(pdbFile, db);
		else if (iequals(program, "SHELXL"))
			WriteRemark3Shelxl(pdbFile, db);
		else if (iequals(program, "PHENIX"))
			WriteRemark3Phenix(pdbFile, db);
		else if (iequals(program, "NUCLSQ"))
			WriteRemark3NuclSQ(pdbFile, db);
		else if (iequals(program, "PROLSQ"))
			WriteRemark3ProlSQ(pdbFile, db);
	}

	for (auto r : db["refine"])
	{
		std::string remarks = r["details"].as<std::string>();
		if (remarks.empty())
			remarks = "NULL";

		WriteOneContinuedLine(pdbFile, "REMARK   3 ", 0, "OTHER REFINEMENT REMARKS: " + remarks);
		break;
	}
}

void WriteRemark200(std::ostream &pdbFile, const datablock &db)
{
	typedef RM<200> RM;

	try
	{
		for (auto diffrn : db["diffrn"])
		{
			std::string diffrn_id = diffrn["id"].as<std::string>();
			std::string crystal_id = diffrn["crystal_id"].as<std::string>();

			auto diffrn_radiation = db["diffrn_radiation"].find_first(key("diffrn_id") == diffrn_id);
			auto diffrn_radiation_wavelength = db["diffrn_radiation_wavelength"].find_first(key("id") == diffrn_radiation["wavelength_id"].as<std::string>());
			auto diffrn_source = db["diffrn_source"].find_first(key("diffrn_id") == diffrn_id);
			auto diffrn_detector = db["diffrn_detector"].find_first(key("diffrn_id") == diffrn_id);
			auto exptl = db["exptl"].find_first(key("entry_id") == db.name());
			auto exptl_crystal = db["exptl_crystal"].find_first(key("id") == crystal_id);
			auto exptl_crystal_grow = db["exptl_crystal_grow"].find_first(key("crystal_id") == crystal_id);
			auto computing = db["computing"].find_first(key("entry_id") == db.name());
			auto reflns = db["reflns"].find_first(key("entry_id") == db.name());

			std::string pdbx_diffrn_id = reflns["pdbx_diffrn_id"].as<std::string>();

			auto reflns_shell = db["reflns_shell"].find_first(key("pdbx_diffrn_id") == pdbx_diffrn_id);
			auto refine = db["refine"].find_first(key("pdbx_diffrn_id") == pdbx_diffrn_id);

			std::string date =
				diffrn_detector.empty() ? "NULL" : cif2pdbDate(diffrn_detector["pdbx_collection_date"].as<std::string>());

			std::string iis = cifSoftware(db, eDataReduction);
			std::string dss = cifSoftware(db, eDataScaling);

			std::string source = diffrn_source["source"].as<std::string>();
			std::string synchrotron, type;

			if (source.empty())
				synchrotron = "NULL";
			else if (iequals(source, "SYNCHROTRON"))
			{
				synchrotron = "Y";
				source = diffrn_source["pdbx_synchrotron_site"].as<std::string>();
				if (source.empty())
					source = "NULL";
				type = "NULL";
			}
			else
			{
				synchrotron = "N";
				type = diffrn_source["type"].as<std::string>();
				if (type.empty())
					type = "NULL";
			}

			if (source.empty())
				source = "NULL";
			if (type.empty())
				type = "NULL";

			pdbFile << RM("") << '\n'
					<< RM("EXPERIMENTAL DETAILS") << '\n'
					<< RM(" EXPERIMENT TYPE                : ") << Fs(exptl, "method") << '\n'
					<< RM(" DATE OF DATA COLLECTION        : ") << date << '\n'
					<< RM(" TEMPERATURE           (KELVIN) : ", 5, 1) << Ff(diffrn, "ambient_temp") << '\n'
					<< RM(" PH                             : ", 4, 1) << Ff(exptl_crystal_grow, "ph") << '\n'
					<< RM(" NUMBER OF CRYSTALS USED        : ") << Fi(exptl, "crystals_number") << '\n'
					<< RM("") << '\n'
					<< RM(" SYNCHROTRON              (Y/N) : ") << synchrotron << '\n'
					<< RM(" RADIATION SOURCE               : ") << source << '\n'
					<< RM(" BEAMLINE                       : ") << Fs(diffrn_source, "pdbx_synchrotron_beamline") << '\n'
					<< RM(" X-RAY GENERATOR MODEL          : ") << type << '\n'
					<< RM(" MONOCHROMATIC OR LAUE    (M/L) : ") << Fs(diffrn_radiation, "pdbx_monochromatic_or_laue_m_l") << '\n'
					<< RM(" WAVELENGTH OR RANGE        (A) : ", 7, 4) << Ff(diffrn_radiation_wavelength, "wavelength") << '\n'
					<< RM(" MONOCHROMATOR                  : ") << Fs(diffrn_radiation, "monochromator") << '\n'
					<< RM(" OPTICS                         : ") << Fs(diffrn_detector, "details") << '\n'
					<< RM("") << '\n'
					<< RM(" DETECTOR TYPE                  : ") << Fs(diffrn_detector, "detector") << '\n'
					<< RM(" DETECTOR MANUFACTURER          : ") << Fs(diffrn_detector, "type") << '\n'
					<< RM(" INTENSITY-INTEGRATION SOFTWARE : ") << iis << '\n'
					<< RM(" DATA SCALING SOFTWARE          : ") << dss << '\n'
					<< RM(" ") << '\n'
					<< RM(" NUMBER OF UNIQUE REFLECTIONS   : ") << Fi(reflns, "number_obs") << '\n'
					<< RM(" RESOLUTION RANGE HIGH      (A) : ", 7, 3) << Ff(reflns, "d_resolution_high") << '\n'
					<< RM(" RESOLUTION RANGE LOW       (A) : ", 7, 3) << Ff(reflns, "d_resolution_low") << '\n'
					<< RM(" REJECTION CRITERIA  (SIGMA(I)) : ", 7, 3) << Ff(reflns, "observed_criterion_sigma_I") << '\n'
					<< RM("") << '\n'
					<< RM("OVERALL.") << '\n'
					<< RM(" COMPLETENESS FOR RANGE     (%) : ", 7, 1) << Ff(reflns, "percent_possible_obs") << '\n'
					<< RM(" DATA REDUNDANCY                : ", 7, 3) << Ff(reflns, "pdbx_redundancy") << '\n'
					<< RM(" R MERGE                    (I) : ", 7, 5) << Ff(reflns, "pdbx_Rmerge_I_obs") << '\n'
					<< RM(" R SYM                      (I) : ", 7, 5) << Ff(reflns, "pdbx_Rsym_value") << '\n'
					<< RM(" <I/SIGMA(I)> FOR THE DATA SET  : ", 7, 4) << Ff(reflns, "pdbx_netI_over_sigmaI") << '\n'
					<< RM("") << '\n'
					<< RM("IN THE HIGHEST RESOLUTION SHELL.") << '\n'
					<< RM(" HIGHEST RESOLUTION SHELL, RANGE HIGH (A) : ", 7, 2) << Ff(reflns_shell, "d_res_high") << '\n'
					<< RM(" HIGHEST RESOLUTION SHELL, RANGE LOW  (A) : ", 7, 2) << Ff(reflns_shell, "d_res_low") << '\n'
					<< RM(" COMPLETENESS FOR SHELL     (%) : ", 7, 1) << Ff(reflns_shell, "percent_possible_all") << '\n'
					<< RM(" DATA REDUNDANCY IN SHELL       : ", 7, 2) << Ff(reflns_shell, "pdbx_redundancy") << '\n'
					<< RM(" R MERGE FOR SHELL          (I) : ", 7, 5) << Ff(reflns_shell, "Rmerge_I_obs") << '\n'
					<< RM(" R SYM FOR SHELL            (I) : ", 7, 5) << Ff(reflns_shell, "pdbx_Rsym_value") << '\n'
					<< RM(" <I/SIGMA(I)> FOR SHELL         : ", 7, 3) << Ff(reflns_shell, "meanI_over_sigI_obs") << '\n'
					<< RM("") << '\n';

			struct
			{
				row_handle r;
				const char *field;
				const char *dst;
			} kTail[] = {
				{ diffrn_radiation, "pdbx_diffrn_protocol", "DIFFRACTION PROTOCOL: " },
				{ refine, "pdbx_method_to_determine_struct", "METHOD USED TO DETERMINE THE STRUCTURE: " },
				{ computing, "structure_solution", "SOFTWARE USED: " },
				{ refine, "pdbx_starting_model", "STARTING MODEL: " },
				{ exptl_crystal, "description", "\nREMARK: " }
			};

			for (auto &t : kTail)
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
	catch (const std::exception &ex)
	{
		if (VERBOSE >= 0)
			std::cerr << ex.what() << '\n';
	}
}

void WriteRemark280(std::ostream &pdbFile, const datablock &db)
{
	typedef RM<280> RM;

	try
	{
		for (auto exptl_crystal : db["exptl_crystal"])
		{
			std::string crystal_id = exptl_crystal["id"].as<std::string>();
			auto exptl_crystal_grow = db["exptl_crystal_grow"].find_first(key("crystal_id") == crystal_id);

			pdbFile
				<< RM("") << '\n'
				<< RM("CRYSTAL") << '\n'
				<< RM("SOLVENT CONTENT, VS   (%): ", 6, 2) << Ff(exptl_crystal, "density_percent_sol") << '\n'
				<< RM("MATTHEWS COEFFICIENT, VM (ANGSTROMS**3/DA): ", 6, 2) << Ff(exptl_crystal, "density_Matthews") << '\n'
				<< RM("") << '\n';

			std::vector<std::string> conditions;
			auto add = [&conditions](const std::string c)
			{
				if (find(conditions.begin(), conditions.end(), c) == conditions.end())
					conditions.push_back(c);
			};

			const char *keys[] = { "pdbx_details", "ph", "method", "temp" };

			for (std::size_t i = 0; i < (sizeof(keys) / sizeof(const char *)); ++i)
			{
				const char *c = keys[i];

				std::string v = exptl_crystal_grow[c].as<std::string>();
				if (not v.empty())
				{
					to_upper(v);

					switch (i)
					{
						case 1: add("PH " + v); break;
						case 3: add("TEMPERATURE " + v + "K"); break;

						default:
							for (std::string::size_type b = 0, e = v.find(", "); b != std::string::npos; b = (e == std::string::npos ? e : e + 2), e = v.find(", ", b))
								add(v.substr(b, e - b));
							break;
					}
				}
			}

			WriteOneContinuedLine(pdbFile, "REMARK 280", 0, "CRYSTALLIZATION CONDITIONS: " + (conditions.empty() ? "NULL" : join(conditions, ", ")));

			break;
		}
	}
	catch (const std::exception &ex)
	{
		if (VERBOSE >= 0)
			std::cerr << ex.what() << '\n';
	}
}

void WriteRemark350(std::ostream &pdbFile, const datablock &db)
{
	auto &c1 = db["pdbx_struct_assembly"];
	if (c1.empty())
		return;

	std::vector<std::string> biomolecules, details;
	for (auto bm : c1)
	{
		std::string id = bm["id"].as<std::string>();
		biomolecules.push_back(id);

		for (auto r : db["struct_biol"].find(key("id") == id))
		{
			std::string s = r["details"].as<std::string>();
			if (not s.empty())
				details.push_back(s);
		}
	}

	// write out the mandatory REMARK 300 first

	pdbFile << RM<300>("") << '\n'
			<< RM<300>("BIOMOLECULE: ") << join(biomolecules, ", ") << '\n'
			<< RM<300>("SEE REMARK 350 FOR THE AUTHOR PROVIDED AND/OR PROGRAM") << '\n'
			<< RM<300>("GENERATED ASSEMBLY INFORMATION FOR THE STRUCTURE IN") << '\n'
			<< RM<300>("THIS ENTRY. THE REMARK MAY ALSO PROVIDE INFORMATION ON") << '\n'
			<< RM<300>("BURIED SURFACE AREA.") << '\n';

	if (not details.empty())
	{
		pdbFile << RM<300>("REMARK:") << '\n';

		for (auto detail : details)
			WriteOneContinuedLine(pdbFile, "REMARK 300", 0, detail);
	}

	typedef RM<350> RM;

	pdbFile << RM("") << '\n'
			<< RM("COORDINATES FOR A COMPLETE MULTIMER REPRESENTING THE KNOWN") << '\n'
			<< RM("BIOLOGICALLY SIGNIFICANT OLIGOMERIZATION STATE OF THE") << '\n'
			<< RM("MOLECULE CAN BE GENERATED BY APPLYING BIOMT TRANSFORMATIONS") << '\n'
			<< RM("GIVEN BELOW.  BOTH NON-CRYSTALLOGRAPHIC AND") << '\n'
			<< RM("CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN.") << '\n';

	for (auto bm : c1)
	{
		std::string id, detail, method, oligomer;
		cif::tie(id, detail, method, oligomer) = bm.get("id", "details", "method_details", "oligomeric_details");

		pdbFile << RM("") << '\n'
				<< RM("BIOMOLECULE: ") << id << '\n';

		to_upper(oligomer);

		if (detail == "author_defined_assembly" or detail == "author_and_software_defined_assembly")
			pdbFile << RM("AUTHOR DETERMINED BIOLOGICAL UNIT: ") << oligomer << '\n';

		if (detail == "software_defined_assembly" or detail == "author_and_software_defined_assembly")
			pdbFile << RM("SOFTWARE DETERMINED QUATERNARY STRUCTURE: ") << oligomer << '\n';

		if (not method.empty())
			pdbFile << RM("SOFTWARE USED: ") << method << '\n';

		for (std::string type : { "ABSA (A^2)", "SSA (A^2)", "MORE" })
		{
			for (auto prop : db["pdbx_struct_assembly_prop"].find(key("biol_id") == id and key("type") == type))
			{
				std::string value = prop["value"].as<std::string>();

				if (iequals(type, "ABSA (A^2)"))
					pdbFile << RM("TOTAL BURIED SURFACE AREA: ") << value << " ANGSTROM**2\n";
				else if (iequals(type, "SSA (A^2)"))
					pdbFile << RM("SURFACE AREA OF THE COMPLEX: ") << value << " ANGSTROM**2\n";
				else if (iequals(type, "MORE"))
					pdbFile << RM("CHANGE IN SOLVENT FREE ENERGY: ") << value << " KCAL/MOL\n";
			}
		}

		auto gen = db["pdbx_struct_assembly_gen"].find_first(key("assembly_id") == id);

		if (gen)
		{
			std::string asym_id_list, oper_id_list;
			cif::tie(asym_id_list, oper_id_list) = gen.get("asym_id_list", "oper_expression");

			auto asyms = split<std::string>(asym_id_list, ",");

			std::vector<std::string> chains = MapAsymIDs2ChainIDs(asyms, db);
			pdbFile << RM("APPLY THE FOLLOWING TO CHAINS: ") << join(chains, ", ") << '\n';

			for (auto oper_id : split<std::string>(oper_id_list, ",", true))
			{
				auto r = db["pdbx_struct_oper_list"].find_first(key("id") == oper_id);

				pdbFile << RM("  BIOMT1 ", -3) << Fs(r, "id")
						<< SEP(" ", -9, 6) << Ff(r, "matrix[1][1]")
						<< SEP(" ", -9, 6) << Ff(r, "matrix[1][2]")
						<< SEP(" ", -9, 6) << Ff(r, "matrix[1][3]")
						<< SEP(" ", -14, 5) << Ff(r, "vector[1]")
						<< '\n'
						<< RM("  BIOMT2 ", -3) << Fs(r, "id")
						<< SEP(" ", -9, 6) << Ff(r, "matrix[2][1]")
						<< SEP(" ", -9, 6) << Ff(r, "matrix[2][2]")
						<< SEP(" ", -9, 6) << Ff(r, "matrix[2][3]")
						<< SEP(" ", -14, 5) << Ff(r, "vector[2]")
						<< '\n'
						<< RM("  BIOMT3 ", -3) << Fs(r, "id")
						<< SEP(" ", -9, 6) << Ff(r, "matrix[3][1]")
						<< SEP(" ", -9, 6) << Ff(r, "matrix[3][2]")
						<< SEP(" ", -9, 6) << Ff(r, "matrix[3][3]")
						<< SEP(" ", -14, 5) << Ff(r, "vector[3]")
						<< '\n';
			}
		}
	}
}

void WriteRemark400(std::ostream &pdbFile, const datablock &db)
{
	for (auto r : db["pdbx_entry_details"])
	{
		std::string compound_details = r["compound_details"].as<std::string>();
		if (not compound_details.empty())
			WriteOneContinuedLine(pdbFile, "REMARK 400", 0, "\nCOMPOUND\n" + compound_details);
	}
}

void WriteRemark450(std::ostream &pdbFile, const datablock &db)
{
	for (auto r : db["pdbx_entry_details"])
	{
		std::string source_details = r["source_details"].as<std::string>();
		if (not source_details.empty())
			WriteOneContinuedLine(pdbFile, "REMARK 450", 0, "\nSOURCE\n" + source_details, 11);
		break;
	}
}

void WriteRemark465(std::ostream &pdbFile, const datablock &db)
{
	bool first = true;
	typedef RM<465> RM;

	auto &c = db["pdbx_unobs_or_zero_occ_residues"];
	std::vector<row_handle> missing(c.begin(), c.end());
	stable_sort(missing.begin(), missing.end(), [](row_handle a, row_handle b) -> bool
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
		
		return d < 0; });

	for (auto r : missing)
	{
		if (first)
		{
			pdbFile << RM("") << '\n'
					<< RM("MISSING RESIDUES") << '\n'
					<< RM("THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE") << '\n'
					<< RM("EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN") << '\n'
					<< RM("IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)") << '\n'
					<< RM("") << '\n'
					<< RM("  M RES C SSSEQI") << '\n';
			first = false;
		}

		std::string modelNr, resName, chainID, iCode;
		int seqNr;

		cif::tie(modelNr, resName, chainID, iCode, seqNr) =
			r.get("PDB_model_num", "auth_comp_id", "auth_asym_id", "PDB_ins_code", "auth_seq_id");

		pdbFile << cif::format("REMARK 465 %3.3s %3.3s %1.1s %5d%1.1s", modelNr, resName, chainID, seqNr, iCode) << '\n';
	}
}

void WriteRemark470(std::ostream &pdbFile, const datablock &db)
{
	typedef RM<470> RM;

	// wow...
	typedef std::tuple<std::string, std::string, int, std::string, std::string> key_type;
	std::map<key_type, std::deque<std::string>> data;

	for (auto r : db["pdbx_unobs_or_zero_occ_atoms"])
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
		pdbFile << RM("") << '\n'
				<< RM("MISSING ATOM") << '\n'
				<< RM("THE FOLLOWING RESIDUES HAVE MISSING ATOMS (M=MODEL NUMBER;") << '\n'
				<< RM("RES=RESIDUE NAME; C=CHAIN IDENTIFIER; SSEQ=SEQUENCE NUMBER;") << '\n'
				<< RM("I=INSERTION CODE):") << '\n'
				<< RM("  M RES CSSEQI  ATOMS") << '\n';

		for (auto &a : data)
		{
			std::string modelNr, resName, chainID, iCode;
			int seqNr;

			std::tie(modelNr, chainID, seqNr, iCode, resName) = a.first;

			while (not a.second.empty())
			{
				pdbFile << cif::format("REMARK 470 %3.3s %3.3s %1.1s%4d%1.1s  ", modelNr, resName, chainID, seqNr, iCode) << "  ";

				for (std::size_t i = 0; i < 6 and not a.second.empty(); ++i)
				{
					pdbFile << cif2pdbAtomName(a.second.front(), resName, db) << ' ';
					a.second.pop_front();
				}

				pdbFile << '\n';
			}
		}
	}
}

void WriteRemark610(std::ostream &pdbFile, const datablock &db)
{
	// #warning("unimplemented!");
}

void WriteRemark800(std::ostream &pdbFile, const datablock &db)
{
	int nr = 0;
	for (auto r : db["struct_site"])
	{
		pdbFile << "REMARK 800\n";
		if (++nr == 1)
		{
			pdbFile << "REMARK 800 SITE\n";
			++nr;
		}

		std::string ident, code, desc;
		cif::tie(ident, code, desc) = r.get("id", "pdbx_evidence_code", "details");

		to_upper(code);

		for (auto l : { "SITE_IDENTIFIER: " + ident, "EVIDENCE_CODE: " + code, "SITE_DESCRIPTION: " + desc })
		{
			for (auto s : word_wrap(l, 69))
				pdbFile << "REMARK 800 " << s << '\n';
		};
	}
}

void WriteRemark999(std::ostream &pdbFile, const datablock &db)
{
	for (auto r : db["pdbx_entry_details"])
	{
		std::string sequence_details = r["sequence_details"].as<std::string>();
		if (not sequence_details.empty())
			WriteOneContinuedLine(pdbFile, "REMARK 999", 0, "\nSEQUENCE\n" + sequence_details, 11);
		break;
	}
}

void WriteRemarks(std::ostream &pdbFile, const datablock &db)
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

int WritePrimaryStructure(std::ostream &pdbFile, const datablock &db)
{
	int numSeq = 0;

	// DBREF

	for (auto r : db["struct_ref"])
	{
		std::string id, db_name, db_code;
		cif::tie(id, db_name, db_code) = r.get("id", "db_name", "db_code");

		for (auto r1 : db["struct_ref_seq"].find(key("ref_id") == id))
		{
			std::string idCode, chainID, insertBegin, insertEnd, dbAccession, dbinsBeg, dbinsEnd;
			std::string seqBegin, seqEnd, dbseqBegin, dbseqEnd;

			cif::tie(idCode, chainID, seqBegin, insertBegin, seqEnd, insertEnd, dbAccession, dbseqBegin, dbinsBeg, dbseqEnd, dbinsEnd) = r1.get("pdbx_PDB_id_code", "pdbx_strand_id", "pdbx_auth_seq_align_beg", "pdbx_seq_align_beg_ins_code", "pdbx_auth_seq_align_end",
				"pdbx_seq_align_end_ins_code", "pdbx_db_accession", "db_align_beg", "pdbx_db_align_beg_ins_code", "db_align_end", "pdbx_db_align_end_ins_code");

			if (dbAccession.length() > 8 or db_code.length() > 12 or atoi(dbseqEnd.c_str()) >= 100000)
				pdbFile << cif::format(
							   "DBREF1 %4.4s %1.1s %4.4s%1.1s %4.4s%1.1s %-6.6s               %-20.20s",
							   idCode, chainID, seqBegin, insertBegin, seqEnd, insertEnd, db_name, db_code)
						<< '\n'
						<< cif::format(
							   "DBREF2 %4.4s %1.1s     %-22.22s     %10.10s  %10.10s",
							   idCode, chainID, dbAccession, dbseqBegin, dbseqEnd)
						<< '\n';
			else
				pdbFile << cif::format(
							   "DBREF  %4.4s %1.1s %4.4s%1.1s %4.4s%1.1s %-6.6s %-8.8s %-12.12s %5.5s%1.1s %5.5s%1.1s",
							   idCode, chainID, seqBegin, insertBegin, seqEnd, insertEnd, db_name, dbAccession, db_code, dbseqBegin, dbinsBeg, dbseqEnd, dbinsEnd)
						<< '\n';
		}
	}

	// SEQADV

	for (auto r : db["struct_ref_seq_dif"])
	{
		std::string idCode, resName, chainID, seqNum, iCode, database, dbAccession, dbRes, dbSeq, conflict;

		cif::tie(idCode, resName, chainID, seqNum, iCode, database, dbAccession, dbRes, dbSeq, conflict) = r.get("pdbx_PDB_id_code", "mon_id", "pdbx_pdb_strand_id", "pdbx_auth_seq_num", "pdbx_pdb_ins_code",
			"pdbx_seq_db_name", "pdbx_seq_db_accession_code", "db_mon_id", "pdbx_seq_db_seq_num",
			"details");

		to_upper(conflict);

		pdbFile << cif::format(
					   "SEQADV %4.4s %3.3s %1.1s %4.4s%1.1s %-4.4s %-9.9s %3.3s %5.5s %-21.21s",
					   idCode, resName, chainID, seqNum, iCode, database, dbAccession, dbRes, dbSeq, conflict)
					   .str()
				<< '\n';
	}

	// SEQRES

	std::map<char, std::vector<std::string>> seqres;
	std::map<char, int> seqresl;
	for (auto r : db["pdbx_poly_seq_scheme"])
	{
		std::string chainID, res;
		cif::tie(chainID, res) = r.get("pdb_strand_id", "mon_id");
		if (chainID.empty() or res.length() > 3 or res.length() < 1)
			throw std::runtime_error("invalid pdbx_poly_seq_scheme record, chain: " + chainID + " res: " + res);
		seqres[chainID[0]].push_back(std::string(3 - res.length(), ' ') + res);
		++seqresl[chainID[0]];
	}

	for (auto &&[chainID, seq] : seqres)
	{
		int n = 1;
		while (seq.empty() == false)
		{
			auto t = seq.size();
			if (t > 13)
				t = 13;

			pdbFile << cif::format(
						   "SEQRES %3d %1.1s %4d  %-51.51s          ",
						   n++, std::string{ chainID }, seqresl[chainID], join(seq.begin(), seq.begin() + t, " "))
					<< '\n';

			++numSeq;

			seq.erase(seq.begin(), seq.begin() + t);
		}
	}

	// MODRES

	for (auto r : db["pdbx_struct_mod_residue"])
	{
		std::string chainID, seqNum, resName, iCode, stdRes, comment;

		cif::tie(chainID, seqNum, resName, iCode, stdRes, comment) =
			r.get("auth_asym_id", "auth_seq_id", "auth_comp_id", "PDB_ins_code", "parent_comp_id", "details");

		pdbFile << cif::format(
					   "MODRES %4.4s %3.3s %1.1s %4.4s%1.1s %3.3s  %-41.41s",
					   db.name(), resName, chainID, seqNum, iCode, stdRes, comment)
					   .str()
				<< '\n';
	}

	return numSeq;
}

int WriteHeterogen(std::ostream &pdbFile, const datablock &db)
{
	int numHet = 0;

	std::string water_entity_id, water_comp_id;
	for (auto r : db["entity"].find(key("type") == std::string("water")))
	{
		water_entity_id = r["id"].as<std::string>();
		break;
	}

	std::map<std::string, std::string> het;

	for (auto r : db["chem_comp"])
	{
		std::string id, name, mon_nstd_flag;
		cif::tie(id, name, mon_nstd_flag) = r.get("id", "name", "mon_nstd_flag");

		if (mon_nstd_flag == "y")
			continue;

		het[id] = name;
	}

	for (auto r : db["pdbx_entity_nonpoly"])
	{
		std::string entity_id, name, comp_id;
		cif::tie(entity_id, name, comp_id) = r.get("entity_id", "name", "comp_id");

		if (entity_id == water_entity_id)
			water_comp_id = comp_id;

		if (het.count(comp_id) == 0)
			het[comp_id] = name;
	}

	struct HET
	{
		bool water;
		std::string hetID;
		char chainID;
		int seqNum;
		char iCode;
		int numHetAtoms;
		std::string text; // ignored
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
	//		tie(seqNum, comp_id, chain_id, iCode) =
	//			r.get("auth_seq_id", "auth_comp_id", "auth_asym_id", "pdbx_PDB_ins_code");
	//
	//		if (chainID != lChainID or compID != lCompID or seqNum != lSeqNum or iCode != lICode)
	//
	//	}

	// count the HETATM's
	//	for (auto r: db["atom_site"].find(key("group_PDB") == std::string("HETATM")))
	std::set<std::string> missingHetNames;

	for (auto r : db["atom_site"])
	{
		int seqNum;
		std::string entity_id, comp_id, chain_id, iCode, modelNr;

		cif::tie(entity_id, seqNum, comp_id, chain_id, iCode, modelNr) =
			r.get("label_entity_id", "auth_seq_id", "auth_comp_id", "auth_asym_id", "pdbx_PDB_ins_code", "pdbx_PDB_model_num");

		if (compound_factory::kAAMap.count(comp_id) or compound_factory::kBaseMap.count(comp_id))
			continue;

		if (chain_id.length() != 1)
			throw std::runtime_error("Cannot produce PDB file, auth_asym_id not valid");

		if (entity_id != water_entity_id and het.count(comp_id) == 0)
			missingHetNames.insert(comp_id);

		auto h = find_if(hets.begin(), hets.end(),
			[=](const HET &het) -> bool
			{
				return het.hetID == comp_id and het.chainID == chain_id[0] and het.seqNum == seqNum;
			});

		if (h == hets.end())
		{
			hets.push_back({ entity_id == water_entity_id, comp_id, chain_id[0], seqNum,
				(iCode.empty() ? ' ' : iCode[0]), 1 });
		}
		else
			h->numHetAtoms += 1;
	}

	if (VERBOSE > 1 and not missingHetNames.empty())
		std::cerr << "Missing het name(s) for " << join(missingHetNames, ", ") << '\n';

	for (auto h : hets)
	{
		if (h.water)
			continue;
		pdbFile << cif::format("HET    %3.3s  %c%4d%c  %5d", h.hetID, h.chainID, h.seqNum, h.iCode, h.numHetAtoms) << '\n';
		++numHet;
	}

	for (auto &&[id, name] : het)
	{
		if (id == water_comp_id)
			continue;

		to_upper(name);

		int c = 1;

		for (;;)
		{
			pdbFile << cif::format("HETNAM  %2.2s %3.3s ", (c > 1 ? std::to_string(c) : std::string()), id);
			++c;

			if (name.length() > 55)
			{
				bool done = false;
				for (auto e = name.begin() + 54; e != name.begin(); --e)
				{
					if (ispunct(*e))
					{
						pdbFile << std::string(name.begin(), e) << '\n';
						name.erase(name.begin(), e);
						done = true;
						break;
					}
				}

				if (not done)
				{
					pdbFile << std::string(name.begin(), name.begin() + 55) << '\n';
					name.erase(name.begin(), name.begin() + 55);
				}

				continue;
			}

			pdbFile << name << '\n';
			break;
		}
	}

	for (auto &&[id, name] : het)
	{
		if (id == water_comp_id)
			continue;

		std::string syn = db["chem_comp"].find_first<std::string>(key("id") == id, "pdbx_synonyms");
		if (syn.empty())
			continue;

		WriteOneContinuedLine(pdbFile, "HETSYN", 4, id + ' ' + syn, 11);
	}

	// FORMUL

	std::vector<std::string> formulas;

	for (auto h : het)
	{
		std::string hetID = h.first;
		int componentNr = 0;

		std::string first_het_asym_id;
		for (auto p : db["pdbx_poly_seq_scheme"].find(key("mon_id") == hetID))
		{
			first_het_asym_id = p["asym_id"].as<std::string>();
			break;
		}

		if (first_het_asym_id.empty())
		{
			for (auto p : db["pdbx_nonpoly_scheme"].find(key("mon_id") == hetID))
			{
				first_het_asym_id = p["asym_id"].as<std::string>();
				break;
			}
		}

		if (not first_het_asym_id.empty())
		{
			for (auto a : db["struct_asym"])
			{
				++componentNr;
				if (a["id"] == first_het_asym_id)
					break;
			}
		}

		auto nr = count_if(hets.begin(), hets.end(), [hetID](auto &h) -> bool
			{ return h.hetID == hetID; });

		for (auto r : db["chem_comp"].find(key("id") == hetID))
		{
			std::string formula = r["formula"].as<std::string>();
			if (nr > 1)
				formula = std::to_string(nr) + '(' + formula + ')';

			int c = 1;
			for (;;)
			{
				std::stringstream fs;

				fs << cif::format("FORMUL  %2d  %3.3s %2.2s%c", componentNr, hetID, (c > 1 ? std::to_string(c) : std::string()), (hetID == water_comp_id ? '*' : ' '));
				++c;

				if (formula.length() > 51)
				{
					bool done = false;
					for (auto e = formula.begin() + 50; e != formula.begin(); --e)
					{
						if (ispunct(*e))
						{
							pdbFile << std::string(formula.begin(), e) << '\n';
							formula.erase(formula.begin(), e);
							done = true;
							break;
						}
					}

					if (not done)
					{
						pdbFile << std::string(formula.begin(), formula.begin() + 55) << '\n';
						formula.erase(formula.begin(), formula.begin() + 55);
					}

					continue;
				}

				fs << formula << '\n';

				formulas.push_back(fs.str());
				break;
			}

			break;
		}
	}

	sort(formulas.begin(), formulas.end(), [](const std::string &a, const std::string &b) -> bool
		{ return stoi(a.substr(8, 2)) < stoi(b.substr(8, 2)); });

	for (auto &f : formulas)
		pdbFile << f;

	return numHet;
}

std::tuple<int, int> WriteSecondaryStructure(std::ostream &pdbFile, const datablock &db)
{
	int numHelix = 0, numSheet = 0;

	// HELIX
	for (auto r : db["struct_conf"].find(key("conf_type_id") == "HELX_P"))
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
		pdbFile << cif::format("HELIX  %3d %3.3s %3.3s %1.1s %4d%1.1s %3.3s %1.1s %4d%1.1s%2d%-30.30s %5d",
					   numHelix, pdbx_PDB_helix_id, beg_label_comp_id, beg_auth_asym_id, beg_auth_seq_id, pdbx_beg_PDB_ins_code, end_label_comp_id, end_auth_asym_id, end_auth_seq_id, pdbx_end_PDB_ins_code, pdbx_PDB_helix_class, details, pdbx_PDB_helix_length)
				<< '\n';
	}

	for (auto r : db["struct_sheet"])
	{
		std::string sheetID;
		int numStrands = 0;

		cif::tie(sheetID, numStrands) = r.get("id", "number_strands");

		bool first = true;

		for (auto o : db["struct_sheet_order"].find(key("sheet_id") == sheetID))
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

				auto r1 = db["struct_sheet_range"].find_first(key("sheet_id") == sheetID and key("id") == rangeID1);

				cif::tie(initResName, initICode, endResName, endICode,
					initResName, initChainID, initSeqNum, endResName, endChainID, endSeqNum) = r1.get("beg_label_comp_id", "pdbx_beg_PDB_ins_code", "end_label_comp_id",
					"pdbx_end_PDB_ins_code", "beg_auth_comp_id", "beg_auth_asym_id", "beg_auth_seq_id",
					"end_auth_comp_id", "end_auth_asym_id", "end_auth_seq_id");

				pdbFile << cif::format("SHEET  %3.3s %3.3s%2d %3.3s %1.1s%4d%1.1s %3.3s %1.1s%4d%1.1s%2d", rangeID1, sheetID, numStrands, initResName, initChainID, initSeqNum, initICode, endResName, endChainID, endSeqNum, endICode, 0) << '\n';

				first = false;
			}

			std::string initResName, initChainID, initICode, endResName, endChainID, endICode, curAtom, curResName, curChainID, curICode, prevAtom, prevResName, prevChainID, prevICode;
			int initSeqNum, endSeqNum, curResSeq, prevResSeq;

			auto r2 = db["struct_sheet_range"].find_first(key("sheet_id") == sheetID and key("id") == rangeID2);

			cif::tie(initResName, initICode, endResName, endICode,
				initResName, initChainID, initSeqNum, endResName, endChainID, endSeqNum) = r2.get("beg_label_comp_id", "pdbx_beg_PDB_ins_code", "end_label_comp_id",
				"pdbx_end_PDB_ins_code", "beg_auth_comp_id", "beg_auth_asym_id", "beg_auth_seq_id",
				"end_auth_comp_id", "end_auth_asym_id", "end_auth_seq_id");

			auto h = db["pdbx_struct_sheet_hbond"].find(key("sheet_id") == sheetID and key("range_id_1") == rangeID1 and key("range_id_2") == rangeID2);

			if (h.empty())
			{
				pdbFile << cif::format("SHEET  %3.3s %3.3s%2d %3.3s %1.1s%4d%1.1s %3.3s %1.1s%4d%1.1s%2d", rangeID2, sheetID, numStrands, initResName, initChainID, initSeqNum, initICode, endResName, endChainID, endSeqNum, endICode, sense) << '\n';
			}
			else
			{
				std::string compID[2];
				cif::tie(compID[0], compID[1]) = h.front().get("range_2_label_comp_id", "range_1_label_comp_id");

				cif::tie(curAtom, curResName, curResSeq, curChainID, curICode, prevAtom, prevResName, prevResSeq, prevChainID, prevICode) = h.front().get("range_2_auth_atom_id", "range_2_auth_comp_id", "range_2_auth_seq_id", "range_2_auth_asym_id", "range_2_PDB_ins_code",
					"range_1_auth_atom_id", "range_1_auth_comp_id", "range_1_auth_seq_id", "range_1_auth_asym_id", "range_1_PDB_ins_code");

				curAtom = cif2pdbAtomName(curAtom, compID[0], db);
				prevAtom = cif2pdbAtomName(prevAtom, compID[1], db);

				pdbFile << cif::format("SHEET  %3.3s %3.3s%2d %3.3s %1.1s%4d%1.1s %3.3s %1.1s%4d%1.1s%2d "
										"%-4.4s%3.3s %1.1s%4d%1.1s %-4.4s%3.3s %1.1s%4d%1.1s",
							   rangeID2, sheetID, numStrands, initResName, initChainID, initSeqNum, initICode, endResName, endChainID, endSeqNum, endICode, sense, curAtom, curResName, curChainID, curResSeq, curICode, prevAtom, prevResName, prevChainID, prevResSeq, prevICode)
						<< '\n';
			}

			++numSheet;
		}
	}

	return std::make_tuple(numHelix, numSheet);
}

void WriteConnectivity(std::ostream &pdbFile, const datablock &db)
{
	// SSBOND
	// have to filter out alts
	std::set<std::tuple<char, int, char, char, int, char>> ssSeen;

	int nr = 1;
	for (auto r : db["struct_conn"].find(key("conn_type_id") == "disulf"))
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

		pdbFile << cif::format("SSBOND %3d CYS %1.1s %4d%1.1s   CYS %1.1s %4d%1.1s                       %6.6s %6.6s %5.2f", nr, chainID1, seqNum1, icode1, chainID2, seqNum2, icode2, sym1, sym2, Length) << '\n';

		++nr;
	}

	// LINK

	for (auto r : db["struct_conn"].find(key("conn_type_id") == "metalc" or key("conn_type_id") == "covale"))
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

		pdbFile << cif::format("LINK        %-4.4s%1.1s%3.3s %1.1s%4d%1.1s               %-4.4s%1.1s%3.3s %1.1s%4d%1.1s  %6.6s %6.6s", name1, altLoc1, resName1, chainID1, resSeq1, iCode1, name2, altLoc2, resName2, chainID2, resSeq2, iCode2, sym1, sym2);

		if (not Length.empty())
			pdbFile << cif::format(" %5.2f", stod(Length));

		pdbFile << '\n';
	}

	// CISPEP

	for (auto r : db["struct_mon_prot_cis"])
	{
		std::string serNum, pep1, chainID1, icode1, pep2, chainID2, icode2, modNum;
		int seqNum1, seqNum2;
		float measure;

		cif::tie(serNum, pep1, chainID1, seqNum1, icode1, pep2, chainID2, seqNum2, icode2, modNum, measure) =
			r.get("pdbx_id", "label_comp_id", "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code",
				"pdbx_label_comp_id_2", "pdbx_auth_asym_id_2", "pdbx_auth_seq_id_2", "pdbx_PDB_ins_code_2",
				"pdbx_PDB_model_num", "pdbx_omega_angle");

		pdbFile << cif::format("CISPEP %3.3s %3.3s %1.1s %4d%1.1s   %3.3s %1.1s %4d%1.1s       %3.3s       %6.2f",
			serNum, pep1, chainID1, seqNum1, icode1, pep2, chainID2, seqNum2, icode2, modNum, measure) << '\n';
	}
}

int WriteMiscellaneousFeatures(std::ostream &pdbFile, const datablock &db)
{
	int numSite = 0;

	// SITE

	std::map<std::string, std::deque<std::string>> sites;

	for (auto r : db["struct_site_gen"])
	{
		std::string siteID, resName, chainID, iCode;
		int seq;

		cif::tie(siteID, resName, chainID, seq, iCode) =
			r.get("site_id", "auth_comp_id", "auth_asym_id", "auth_seq_id", "pdbx_auth_ins_code");

		sites[siteID].push_back(cif::format("%3.3s %1.1s%4d%1.1s ", resName, chainID, seq, iCode).str());
	}

	for (auto s : sites)
	{
		std::string siteID = std::get<0>(s);
		std::deque<std::string> &res = std::get<1>(s);

		std::size_t numRes = res.size();

		int nr = 1;
		while (res.empty() == false)
		{
			pdbFile << cif::format("SITE   %3d %3.3s %2d ", nr, siteID, numRes);

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

			pdbFile << '\n';
			++nr;
			++numSite;
		}
	}

	return numSite;
}

void WriteCrystallographic(std::ostream &pdbFile, const datablock &db)
{
	auto r = db["symmetry"].find_first(key("entry_id") == db.name());
	std::string symmetry = r["space_group_name_H-M"].as<std::string>();

	r = db["cell"].find_first(key("entry_id") == db.name());

	pdbFile << cif::format("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11.11s%4d", r["length_a"].as<double>(), r["length_b"].as<double>(), r["length_c"].as<double>(), r["angle_alpha"].as<double>(), r["angle_beta"].as<double>(), r["angle_gamma"].as<double>(), symmetry, r["Z_PDB"].as<int>()) << '\n';
}

int WriteCoordinateTransformation(std::ostream &pdbFile, const datablock &db)
{
	int result = 0;

	for (auto r : db["database_PDB_matrix"])
	{
		pdbFile << cif::format("ORIGX%1d    %10.6f%10.6f%10.6f     %10.5f", 1, r["origx[1][1]"].as<float>(), r["origx[1][2]"].as<float>(), r["origx[1][3]"].as<float>(), r["origx_vector[1]"].as<float>()) << '\n';
		pdbFile << cif::format("ORIGX%1d    %10.6f%10.6f%10.6f     %10.5f", 2, r["origx[2][1]"].as<float>(), r["origx[2][2]"].as<float>(), r["origx[2][3]"].as<float>(), r["origx_vector[2]"].as<float>()) << '\n';
		pdbFile << cif::format("ORIGX%1d    %10.6f%10.6f%10.6f     %10.5f", 3, r["origx[3][1]"].as<float>(), r["origx[3][2]"].as<float>(), r["origx[3][3]"].as<float>(), r["origx_vector[3]"].as<float>()) << '\n';
		result += 3;
		break;
	}

	for (auto r : db["atom_sites"])
	{
		pdbFile << cif::format("SCALE%1d    %10.6f%10.6f%10.6f     %10.5f", 1, r["fract_transf_matrix[1][1]"].as<float>(), r["fract_transf_matrix[1][2]"].as<float>(), r["fract_transf_matrix[1][3]"].as<float>(), r["fract_transf_vector[1]"].as<float>()) << '\n';
		pdbFile << cif::format("SCALE%1d    %10.6f%10.6f%10.6f     %10.5f", 2, r["fract_transf_matrix[2][1]"].as<float>(), r["fract_transf_matrix[2][2]"].as<float>(), r["fract_transf_matrix[2][3]"].as<float>(), r["fract_transf_vector[2]"].as<float>()) << '\n';
		pdbFile << cif::format("SCALE%1d    %10.6f%10.6f%10.6f     %10.5f", 3, r["fract_transf_matrix[3][1]"].as<float>(), r["fract_transf_matrix[3][2]"].as<float>(), r["fract_transf_matrix[3][3]"].as<float>(), r["fract_transf_vector[3]"].as<float>()) << '\n';
		result += 3;
		break;
	}

	int nr = 1;
	for (auto r : db["struct_ncs_oper"])
	{
		std::string given = r["code"] == "given" ? "1" : "";

		pdbFile << cif::format("MTRIX%1d %3d%10.6f%10.6f%10.6f     %10.5f    %1.1s", 1, nr, r["matrix[1][1]"].as<float>(), r["matrix[1][2]"].as<float>(), r["matrix[1][3]"].as<float>(), r["vector[1]"].as<float>(), given) << '\n';
		pdbFile << cif::format("MTRIX%1d %3d%10.6f%10.6f%10.6f     %10.5f    %1.1s", 2, nr, r["matrix[2][1]"].as<float>(), r["matrix[2][2]"].as<float>(), r["matrix[2][3]"].as<float>(), r["vector[2]"].as<float>(), given) << '\n';
		pdbFile << cif::format("MTRIX%1d %3d%10.6f%10.6f%10.6f     %10.5f    %1.1s", 3, nr, r["matrix[3][1]"].as<float>(), r["matrix[3][2]"].as<float>(), r["matrix[3][3]"].as<float>(), r["vector[3]"].as<float>(), given) << '\n';

		++nr;
		result += 3;
	}

	return result;
}

std::tuple<int, int> WriteCoordinatesForModel(std::ostream &pdbFile, const datablock &db,
	const std::map<std::string, std::tuple<std::string, int, std::string>> &last_resseq_for_chain_map,
	std::set<std::string> &terminatedChains, int model_nr)
{
	using namespace cif::literals;

	int numCoord = 0, numTer = 0;

	auto &atom_site = db["atom_site"];
	auto &atom_site_anisotrop = db["atom_site_anisotrop"];
	auto &entity = db["entity"];
	// auto &pdbx_poly_seq_scheme = db["pdbx_poly_seq_scheme"];
	// auto &pdbx_nonpoly_scheme = db["pdbx_nonpoly_scheme"];
	auto &pdbx_branch_scheme = db["pdbx_branch_scheme"];

	int serial = 1;
	auto ri = atom_site.begin();

	std::string id, group, name, altLoc, resName, chainID, iCode, element;
	int resSeq = 0, charge;

	for (;;)
	{
		std::string nextResName, nextChainID, nextICode, modelNum;
		int nextResSeq = 0;

		if (ri != atom_site.end())
			cif::tie(nextResName, nextChainID, nextICode, nextResSeq, modelNum) =
				(*ri).get("label_comp_id", "auth_asym_id", "pdbx_PDB_ins_code", "auth_seq_id", "pdbx_PDB_model_num");

		if (modelNum.empty() == false)
		{
			int nr = 0;
			auto r = std::from_chars(modelNum.data(), modelNum.data() + modelNum.length(), nr);
			if ((bool)r.ec)
			{
				if (VERBOSE > 0)
					std::cerr << "Model number '" << modelNum << "' is not a valid integer\n";
			}

			if (nr != model_nr)
			{
				++ri;
				continue;
			}
		}

		if (chainID.empty() == false and terminatedChains.count(chainID) == 0)
		{
			bool terminate = nextChainID != chainID;

			if (not terminate)
				terminate =
					(nextResSeq != resSeq or iCode != nextICode) and
					(last_resseq_for_chain_map.count(chainID) == false or last_resseq_for_chain_map.at(chainID) == make_tuple(resName, resSeq, iCode));

			if (terminate)
			{
				pdbFile << cif::format("TER   %5d      %3.3s %1.1s%4d%1.1s",  serial,  resName,  chainID,  resSeq,  iCode) << '\n';

				++serial;
				terminatedChains.insert(chainID);

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
		catch (...)
		{ /* perhaps no model number here */
		}

		float x, y, z, occupancy, tempFactor;

		cif::tie(id, group, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, element, charge) =
			r.get("id", "group_PDB", "label_atom_id", "label_alt_id", "auth_comp_id", "auth_asym_id", "auth_seq_id",
				"pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "type_symbol", "pdbx_formal_charge");

		if (resName != "HOH")
		{
			int entity_id = r.get<int>("label_entity_id");
			try
			{
				auto type = entity.find1<std::string>("id"_key == entity_id, "type");

				if (type == "branched")	// find the real auth_seq_num, since sugars have their auth_seq_num reused as sugar number... sigh.
					resSeq = pdbx_branch_scheme.find1<int>("asym_id"_key == r.get<std::string>("label_asym_id") and "pdb_seq_num"_key == resSeq, "auth_seq_num");
				// else if (type == "non-polymer")	// same for non-polymers
				// 	resSeq = pdbx_nonpoly_scheme.find1<int>("asym_id"_key == r.get<std::string>("label_asym_id") and "pdb_seq_num"_key == resSeq, "auth_seq_num");
				// else if (type == "polymer")
				// 	resSeq = pdbx_poly_seq_scheme.find1<int>("asym_id"_key == r.get<std::string>("label_asym_id") and "pdb_seq_num"_key == resSeq, "auth_seq_num");
			}
			catch (const std::exception &ex)
			{
				std::cerr << "Oops, there was not exactly one entity with id " << entity_id << '\n';
			}
		}
		
		if (chainID.length() > 1)
			throw std::runtime_error("Chain ID " + chainID + " won't fit into a PDB file");

		if (name.length() < 4 and (element.length() == 1 or std::toupper(name[0]) != std::toupper(element[0]) or std::toupper(name[1]) != std::toupper(element[1])))
			name.insert(name.begin(), ' ');

		std::string sCharge;
		if (charge != 0)
			sCharge = std::to_string(charge) + (charge > 0 ? '+' : '-');

		pdbFile << cif::format("%-6.6s%5d %-4.4s%1.1s%3.3s %1.1s%4d%1.1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s%2.2s", group, serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, element, sCharge) << '\n';

		++numCoord;

		auto ai = atom_site_anisotrop.find_first(key("id") == id);
		if (not ai.empty())
		//
		//		auto ai = find_if(atom_site_anisotrop.begin(), atom_site_anisotrop.end(), [id](row_handle r) -> bool { return r["id"] == id; });
		//		if (ai != atom_site_anisotrop.end())
		{
			float u11, u22, u33, u12, u13, u23;

			tie(u11, u22, u33, u12, u13, u23) =
				ai.get("U[1][1]", "U[2][2]", "U[3][3]", "U[1][2]", "U[1][3]", "U[2][3]");

			pdbFile << cif::format("ANISOU%5d %-4.4s%1.1s%3.3s %1.1s%4d%1.1s %7d%7d%7d%7d%7d%7d      %2.2s%2.2s", serial, name, altLoc, resName, chainID, resSeq, iCode, std::lrintf(u11 * 10000), std::lrintf(u22 * 10000), std::lrintf(u33 * 10000), std::lrintf(u12 * 10000), std::lrintf(u13 * 10000), std::lrintf(u23 * 10000), element, sCharge) << '\n';
		}

		++serial;
	}

	return std::make_tuple(numCoord, numTer);
}

std::tuple<int, int> WriteCoordinate(std::ostream &pdbFile, const datablock &db)
{
	// residues known from seqres
	//	map<tuple<std::string,int,std::string>,std::string> res2chain_map;
	std::map<std::string, std::tuple<std::string, int, std::string>> last_resseq_for_chain_map;

	for (auto r : db["pdbx_poly_seq_scheme"])
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
		for (auto r : db["atom_site"])
			models.insert(r["pdbx_PDB_model_num"].as<int>());
	}
	catch (...)
	{
	}

	std::tuple<int, int> result;

	if (models.empty() or models == std::set<int>{ 0 })
	{
		std::set<std::string> TERminatedChains;
		result = WriteCoordinatesForModel(pdbFile, db, last_resseq_for_chain_map, TERminatedChains, 0);
	}
	else
	{
		for (int model_nr : models)
		{
			if (models.size() > 1)
				pdbFile << cif::format("MODEL     %4d",  model_nr) << '\n';

			std::set<std::string> TERminatedChains;
			auto n = WriteCoordinatesForModel(pdbFile, db, last_resseq_for_chain_map, TERminatedChains, model_nr);
			if (model_nr == 1)
				result = n;

			if (models.size() > 1)
				pdbFile << "ENDMDL\n";
		}
	}

	return result;
}

void WritePDBHeaderLines(std::ostream &os, const datablock &db)
{
	fill_out_streambuf fb(os);
	write_header_lines(os, db);
}

std::string FixStringLength(const std::string &s, std::string::size_type l)
{
	auto result = s;

	if (result.length() > l)
		result = result.substr(0, l - 4) + "... ";
	else if (result.length() < l)
		result.append(l - result.length(), ' ');

	return result;
}

std::string get_HEADER_line(const datablock &db, std::string::size_type truncate_at)
{
	//    0         1         2         3         4         5         6         7         8
	//    HEADER    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxDDDDDDDDD   IIII

	// HEADER

	std::string keywords;
	auto &cat1 = db["struct_keywords"];

	for (auto r : cat1)
	{
		keywords = r["pdbx_keywords"].as<std::string>();
		if (keywords.length() > truncate_at - 40)
			keywords = keywords.substr(0, truncate_at - 44) + " ...";
	}

	std::string date;
	for (auto r : db["pdbx_database_status"])
	{
		date = r["recvd_initial_deposition_date"].as<std::string>();
		if (date.empty())
			continue;
		date = cif2pdbDate(date);
		break;
	}

	if (date.empty())
	{
		for (auto r : db["database_PDB_rev"])
		{
			date = r["date_original"].as<std::string>();
			if (date.empty())
				continue;
			date = cif2pdbDate(date);
			break;
		}
	}

	return FixStringLength(cif::format("HEADER    %-40.40s%-9.9s   %-4.4s", keywords, date, db.name()).str(), truncate_at);
}

std::string get_COMPND_line(const datablock &db, std::string::size_type truncate_at)
{
	// COMPND
	using namespace std::placeholders;

	int molID = 0;
	std::vector<std::string> cmpnd;

	for (auto r : db["entity"])
	{
		if (r["type"] != "polymer")
			continue;

		std::string entityID = r["id"].as<std::string>();

		++molID;
		cmpnd.push_back("MOL_ID: " + std::to_string(molID));

		std::string molecule = r["pdbx_description"].as<std::string>();
		cmpnd.push_back("MOLECULE: " + molecule);

		auto poly = db["entity_poly"].find(key("entity_id") == entityID);
		if (not poly.empty())
		{
			std::string chains = poly.front()["pdbx_strand_id"].as<std::string>();
			replace_all(chains, ",", ", ");
			cmpnd.push_back("CHAIN: " + chains);
		}

		std::string fragment = r["pdbx_fragment"].as<std::string>();
		if (not fragment.empty())
			cmpnd.push_back("FRAGMENT: " + fragment);

		for (auto sr : db["entity_name_com"].find(key("entity_id") == entityID))
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

	return FixStringLength("COMPND    " + join(cmpnd, "; "), truncate_at);
}

std::string get_SOURCE_line(const datablock &db, std::string::size_type truncate_at)
{
	// SOURCE

	int molID = 0;
	std::vector<std::string> source;

	for (auto r : db["entity"])
	{
		if (r["type"] != "polymer")
			continue;

		std::string entityID = r["id"].as<std::string>();

		++molID;
		source.push_back("MOL_ID: " + std::to_string(molID));

		if (r["src_method"] == "syn")
			source.push_back("SYNTHETIC: YES");

		auto &gen = db["entity_src_gen"];
		const std::pair<const char *, const char *> kGenSourceMapping[] = {
			{ "gene_src_common_name", "ORGANISM_COMMON" },
			{ "pdbx_gene_src_gene", "GENE" },
			{ "gene_src_strain", "STRAIN" },
			{ "pdbx_gene_src_cell_line", "CELL_LINE" },
			{ "pdbx_gene_src_organelle", "ORGANELLE" },
			{ "pdbx_gene_src_cellular_location", "CELLULAR_LOCATION" },
			{ "pdbx_gene_src_scientific_name", "ORGANISM_SCIENTIFIC" },
			{ "pdbx_gene_src_ncbi_taxonomy_id", "ORGANISM_TAXID" },
			{ "pdbx_host_org_scientific_name", "EXPRESSION_SYSTEM" },
			{ "pdbx_host_org_ncbi_taxonomy_id", "EXPRESSION_SYSTEM_TAXID" },
			{ "pdbx_host_org_strain", "EXPRESSION_SYSTEM_STRAIN" },
			{ "pdbx_host_org_variant", "EXPRESSION_SYSTEM_VARIANT" },
			{ "pdbx_host_org_cellular_location", "EXPRESSION_SYSTEM_CELLULAR_LOCATION" },
			{ "pdbx_host_org_vector_type", "EXPRESSION_SYSTEM_VECTOR_TYPE" },
			{ "pdbx_host_org_vector", "EXPRESSION_SYSTEM_VECTOR" },
			{ "pdbx_host_org_gene", "EXPRESSION_SYSTEM_GENE" },
			{ "plasmid_name", "EXPRESSION_SYSTEM_PLASMID" }
		};

		for (auto gr : gen.find(key("entity_id") == entityID))
		{
			for (const auto &[cname, sname] : kGenSourceMapping)
			{
				std::string s = gr[cname].as<std::string>();
				if (not s.empty())
					source.push_back(sname + ": "s + s);
			}
		}

		auto &nat = db["entity_src_nat"];
		const std::pair<const char *, const char *> kNatSourceMapping[] = {
			{ "common_name", "ORGANISM_COMMON" },
			{ "strain", "STRAIN" },
			{ "pdbx_organism_scientific", "ORGANISM_SCIENTIFIC" },
			{ "pdbx_ncbi_taxonomy_id", "ORGANISM_TAXID" },
			{ "pdbx_cellular_location", "CELLULAR_LOCATION" },
			{ "pdbx_plasmid_name", "PLASMID" },
			{ "pdbx_organ", "ORGAN" },
			{ "details", "OTHER_DETAILS" }
		};

		for (auto nr : nat.find(key("entity_id") == entityID))
		{
			for (const auto &[cname, sname] : kNatSourceMapping)
			{
				std::string s = nr[cname].as<std::string>();
				if (not s.empty())
					source.push_back(sname + ": "s + s);
			}
		}
	}

	return FixStringLength("SOURCE    " + join(source, "; "), truncate_at);
}

std::string get_AUTHOR_line(const datablock &db, std::string::size_type truncate_at)
{
	// AUTHOR
	std::vector<std::string> author;
	for (auto r : db["audit_author"])
		author.push_back(cif2pdbAuth(r["name"].as<std::string>()));

	return FixStringLength("AUTHOR    " + join(author, "; "), truncate_at);
}

// --------------------------------------------------------------------

void write(std::ostream &os, const datablock &db)
{
	fill_out_streambuf fb(os);

	int numRemark = 0, numHet = 0, numHelix = 0, numSheet = 0, numTurn = 0, numSite = 0, numXform = 0, numCoord = 0, numTer = 0, numConect = 0, numSeq = 0;

	WriteTitle(os, db);

	int savedLineCount = fb.get_line_count();
	//	numRemark = 				WriteRemarks(pdbFile, db);
	WriteRemarks(os, db);
	numRemark = fb.get_line_count() - savedLineCount;

	numSeq = WritePrimaryStructure(os, db);
	numHet = WriteHeterogen(os, db);
	std::tie(numHelix, numSheet) = WriteSecondaryStructure(os, db);
	WriteConnectivity(os, db);
	numSite = WriteMiscellaneousFeatures(os, db);
	WriteCrystallographic(os, db);
	numXform = WriteCoordinateTransformation(os, db);
	std::tie(numCoord, numTer) = WriteCoordinate(os, db);

	os << cif::format("MASTER    %5d    0%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d",  numRemark,  numHet,  numHelix,  numSheet,  numTurn,  numSite,  numXform,  numCoord,  numTer,  numConect,  numSeq) << '\n'
			<< "END\n";
}

void write(const std::filesystem::path &p, const datablock &db)
{
	gzio::ofstream out(p);

	bool writePDB = false;
	if (p.extension() == ".gz")
		writePDB = iequals(p.stem().extension().string(), ".pdb");
	else
		writePDB = iequals(p.extension().string(), ".pdb");

	if (writePDB)
		write(out, db);
	else
		db.write(out);
}


} // namespace cif::pdb
