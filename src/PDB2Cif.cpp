#include "libpr.h"

#include <map>
#include <set>
#include <boost/regex.hpp>

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "peptidedb.h"
#include "pdb2cif.h"
//#include "remark3/templates.h"
#include "libcif/atom_type.h"
#include "libcif/compound.h"

#include "pdb2cif-remark3.h"

using namespace std;
namespace ba = boost::algorithm;

using cif::datablock;
using cif::category;
using cif::row;
using cif::key;
using cif::iequals;

// --------------------------------------------------------------------
// attempt to come up with better error handling

namespace error
{
	enum pdb_errors
	{
		residue_not_found	= 1000
	};
	
	namespace detail
	{
		class pdb_category : public boost::system::error_category
		{
		  public:

			const char* name() const BOOST_SYSTEM_NOEXCEPT
			{
				return "pdb";
			}
			
			string message(int value) const
			{
				switch (value)
				{
					case residue_not_found:
						return "Residue not found";
					
					default:
						return "Error in PDB format";
				}
			}
		};
	}

	boost::system::error_category& pdb_category()
	{
		static detail::pdb_category impl;
		return impl;
	}
	
	inline boost::system::error_code make_error_code(pdb_errors e)
	{
		return boost::system::error_code(static_cast<int>(e), pdb_category());
	}
}

namespace boost {
namespace system {

template<> struct is_error_code_enum<error::pdb_errors>
{
	static const bool value = true;
};

}
}

// --------------------------------------------------------------------

const map<string,int> kMonths{
	{ "JAN",  1 }, { "FEB",  2 }, { "MAR",  3 },
	{ "APR",  4 }, { "MAY",  5 }, { "JUN",  6 },
	{ "JUL",  7 }, { "AUG",  8 }, { "SEP",  9 },
	{ "OCT", 10 }, { "NOV", 11 }, { "DEC", 12 },
};

const set<string> kSupportedRecords{
	"HEADER", "OBSLTE", "TITLE ", "SPLIT ", "CAVEAT", "COMPND", "SOURCE",
	"KEYWDS", "EXPDTA", "NUMMDL", "MDLTYP", "AUTHOR", "REVDAT", "SPRSDE",
	"JRNL  ", "REMARK", "DBREF ", "DBREF1", "DBREF2", "SEQADV", "SEQRES",
	"MODRES", "HET   ", "HETNAM", "HETSYN", "FORMUL", "HELIX ", "SHEET ",
	"SSBOND", "LINK  ", "CISPEP", "SITE  ", "CRYST1", "ORIGX1", "SCALE1",
	"MTRIX1", "ORIGX2", "SCALE2", "MTRIX2", "ORIGX3", "SCALE3", "MTRIX3",
	"MODEL ", "ATOM  ", "ANISOU", "TER   ", "HETATM", "ENDMDL", "CONECT",
	"MASTER", "END   "
};

bool is_water(const string& resname)
{
	return resname == "HOH" or resname == "H2O" or resname == "OH2" or resname == "WAT" or resname == "DOD";
}

// --------------------------------------------------------------------
//	Unfortunately, parsing a PDB file requires several passes over the
//	data. Therefore we first obtain all records where a record has the
//	value flattened out for continuation.

PDBRecord::PDBRecord(uint32 line_nr, const string& name, const string& value)
	: m_next(nullptr), m_line_nr(line_nr), m_vlen(value.length())
{
	assert(name.length() <= 10);
	
	strcpy(m_name, name.c_str());
	strcpy(m_value, value.c_str());
}

PDBRecord::~PDBRecord()
{
	delete m_next;
}

void* PDBRecord::operator new(size_t size, size_t v_len)
{
	return malloc(size + v_len + 1);
}

void PDBRecord::operator delete(void* p)
{
	free(p);
}

bool PDBRecord::is(const char* name) const
{
	return iequals(m_name, name);
}

char PDBRecord::v_c(size_t column)
{
	char result = ' ';
	if (column - 7 < m_vlen)
		result = m_value[column - 7];
	return result;
}

string PDBRecord::v_s(size_t column_first, size_t column_last)
{
	string result;
	
	if (column_last > m_vlen + 6)
		column_last = m_vlen + 6;
	
	if (column_first < m_vlen + 7)
	{
		result = string{m_value + column_first - 7, m_value + column_last - 7 + 1};
		ba::trim(result);
	}

	return result;
}

int PDBRecord::v_i(int column_first, int column_last)
{
	int result = 0;

	const char* e = m_value + m_vlen;
	if (e > m_value + column_last - 7 + 1)
		e = m_value + column_last - 7 + 1;
	
	enum { start, digit, tail } state = start;
	bool negate = false;

	for (const char* p = m_value + column_first - 7; p < e; ++p)
	{
		switch (state)
		{
			case start:
				if (*p == '+')
					state = digit;
				else if (*p == '-')
				{
					negate = true;
					state = digit;
				}
				else if (isdigit(*p))
				{
					result = *p - '0';
					state = digit;
				}
				else if (not isspace(*p))
					throw runtime_error("Not a valid integer in PDB record");
				break;
				
			case digit:
				if (isspace(*p))
					state = tail;
				else if (not isdigit(*p))
					throw runtime_error("Not a valid integer in PDB record");
				else
					result = result * 10 + *p - '0';
				break;
			
			case tail:
				if (not isspace(*p))
					throw runtime_error("Not a valid integer in PDB record");
				break;
		}
	}
	
	if (negate)
		result = -result;
	
	return result;
}

string PDBRecord::v_f(size_t column_first, size_t column_last)
{
	// for now... TODO: check format?
	return v_s(column_first, column_last);
}

// --------------------------------------------------------------------

class SpecificationListParser
{
  public:
	SpecificationListParser(const string& text)
		: m_text(text), m_p(m_text.begin()) {}
	
	tuple<string,string> GetNextSpecification();
	
  private:
	string				m_text;
	string::iterator	m_p;
};

tuple<string,string> SpecificationListParser::GetNextSpecification()
{
	string id, value;
	
	string::iterator start = m_p, backup;

	enum { eStart, eID, eColon, eValue, eNL, eNL_ID, eSemiColon, eError, eDone } state = eStart;
	
	while (m_p != m_text.end() and state != eDone)
	{
		char ch = *m_p++;
		
		switch (state)
		{
			case eStart:
				if (isalnum(ch) or ch == '_')
				{
					id = { ch };
					value.clear();
					state = eID;
					start = m_p;
				}
				else if (not isspace(ch))
				{
					if (VERBOSE)
						cerr << "skipping invalid character in SOURCE ID: " << ch << endl;
				}
				break;
			
			case eID:
				if (isalnum(ch) or ch == '_')
					id += ch;
				else if (ch == ':')
					state = eColon;
				else
					state = eError;
				break;
			
			case eColon:
				if (ch == ';')
				{
					if (VERBOSE)
						cerr << "Empty value for SOURCE: " << id << endl;
					state = eStart;
				}
				else if (not isspace(ch))
				{
					value = { ch };
					state = eValue;
				}
				break;
			
			case eValue:
				if (ch == '\n')
				{
					backup = m_p;
					state = eNL;
				}
				else if (ch == ';')
				{
					backup = m_p;
					state = eSemiColon;
				}
				else
					value += ch;
				break;
			
			case eSemiColon:
				if (ch == '\n')
					state = eDone;
				else if (ch != ' ')
				{
					value.insert(value.end(), backup, m_p);
					state = eValue;
				}
				break;
			
			case eNL:
				if (isalnum(ch))
				{
					value += ' ';
					state = eNL_ID;
				}
				else if (isspace(ch))
					state = eValue;
				break;
			
			case eNL_ID:
				if (ch == ':')
				{
					m_p = backup;
					state = eDone;
				}
				else if (ch == ';')
					state = eSemiColon;
				else if (not (isalnum(ch) or ch == '_'))
				{
					value.insert(value.end(), backup, m_p);
					state = eValue;
				}
				break;
			
			case eError:
				if (ch == ';')
				{
					if (VERBOSE)
						cerr << "Skipping invalid header line: '" << string(start, m_p) << endl;
					state = eStart;
				}
				break;
			
			case eDone: break;	// keep compiler happy
		}
	}

	ba::trim(value);

	return make_tuple(id, value);
}

// --------------------------------------------------------------------

class PDBFileParser
{
  public:
	
	PDBFileParser()
		: m_data(nullptr), m_rec(nullptr)
	{
	}
	
	~PDBFileParser()
	{
		delete m_data;
	}
	
	void Parse(istream& is, cif::file& result);

  private:

	// ----------------------------------------------------------------
	
	struct DBREF
	{
		string	PDBIDCode;
		char	chainID;
		int		seqBegin;
		char	insertBegin;
		int		seqEnd;
		char	insertEnd;
		string	database;
		string	dbAccession;
		string	dbIdCode;
		int		dbSeqBegin;
		char	dbinsBeg;
		int		dbSeqEnd;
		char	dbinsEnd;
	};

	struct HET
	{
		string		hetID;
		char		chainID;
		int			seqNum;
		char		iCode;
		int			numHetAtoms;
		string		text;
		string		asym_id;
		set<int>	atoms;
	};
	
	struct UNOBS
	{
		int model_nr;
		string res;
		char chain;
		int seq;
		char iCode;
		vector<string> atoms;
	};
	
	// ----------------------------------------------------------------

	/*
		To get from PDB chains to CIF entity and poly records we take the following steps:
		
		First check if there is a Primary Structure Section. If there is, it should contain
		a valid DBREF/SEQRES pair that allows the reconstruction of numbering of residues.
		
		If that fails, we fall back to:
	
		1. Collect the chains from the PDB file.
		2. For each chain, split out the residues and waters, assign those to new entities
		3. If there are multiple chains containing residues, align those to find unique polymers
		4. Annotate the entity records with available information in the PDB file (COMPND e.g.)
		5. Create the mapping structures from PDB numbering to CIF numbering.
	*/

	struct PDBCompound
	{
		int						m_mol_id;
		string					m_title;
		set<char>				m_chains;
		map<string,string>		m_info;
		map<string,string>		m_source;
		int						m_count = 0;
	};
	
	struct PDBSeqRes
	{
		string					m_mon_id;
		int						m_seq_num;
		char					m_icode;

		int						m_db_seq_num;
		bool					m_seen = false;
		set<string>				m_alts;
		
		bool operator==(const PDBSeqRes& rhs) const
		{
			return m_seq_num == rhs.m_seq_num and m_mon_id == rhs.m_mon_id and m_icode == rhs.m_icode;
		}
	};
	
	struct PDBChain
	{
		PDBChain(const string& structure_id, char chain_id, int mol_id)
			: m_dbref{structure_id, chain_id}, m_waters(0), m_ter_index(0), m_mol_id(mol_id), m_next_seq_num(1), m_next_db_seq_num(1)
		{
		}
		
		DBREF					m_dbref;
		vector<PDBSeqRes>		m_seqres, m_het;
		int						m_waters;
		size_t					m_ter_index;
		
		int						m_mol_id;

		// scratch values for reading SEQRES records
		int						m_next_seq_num;
		int						m_next_db_seq_num;
		
		// scratch value for aligning
		struct AtomRes
		{
			string				m_mon_id;
			int					m_seq_num;
			char				m_icode;
			
			bool operator==(const AtomRes& rhs) const		{ return m_mon_id == rhs.m_mon_id and m_seq_num == rhs.m_seq_num and m_icode == rhs.m_icode; }
			bool operator!=(const AtomRes& rhs) const		{ return m_mon_id != rhs.m_mon_id or m_seq_num != rhs.m_seq_num or m_icode != rhs.m_icode; }
		};
		vector<AtomRes>			m_residues_seen;
		
		void AlignResToSeqRes();
		bool SameSequence(const PDBChain& rhs) const;
	};

	// ----------------------------------------------------------------

	PDBCompound& GetOrCreateCompound(int mol_id)
	{
		auto i = find_if(m_compounds.begin(), m_compounds.end(), [mol_id](PDBCompound& comp) -> bool { return comp.m_mol_id == mol_id; });
		if (i == m_compounds.end())
		{
			m_compounds.push_back(PDBCompound{ mol_id });
			
			m_molID2EntityID[mol_id] = to_string(m_next_entity_nr++);
			
			i = prev(m_compounds.end());
		}
		
		return *i;
	}

	// locate the PDBChain record for a chain ID, or create it with dummy data if missing
	PDBChain& GetChainForID(char chainID, int numRes = 0)
	{
		auto i = find_if(m_chains.begin(), m_chains.end(), [chainID](PDBChain& ch) -> bool { return ch.m_dbref.chainID == chainID; });
		
		if (i == m_chains.end())
		{
			// locate the compound for this chain, if any (does that happen?)
			int mol_id = 0;
			for (auto& cmp: m_compounds)
			{
				if (cmp.m_chains.count(chainID) > 0)
				{
					mol_id = cmp.m_mol_id;
					break;
				}
			}
			
			m_chains.emplace_back(m_structure_id, chainID, mol_id);

			i = prev(m_chains.end());
		}
		
		return *i;
	};
	
	void InsertChemComp(const string& chem_comp)
	{
		if (find(m_chem_comp.begin(), m_chem_comp.end(), chem_comp) == m_chem_comp.end())
			m_chem_comp.push_back(chem_comp);
	}

	void InsertAtomType(const string& atom_type)
	{
		if (find(m_atom_types.begin(), m_atom_types.end(), atom_type) == m_atom_types.end())
			m_atom_types.push_back(atom_type);
	}

	// ----------------------------------------------------------------

	template<typename Predicate>
	PDBRecord* FindRecord(Predicate&& pred)
	{
		PDBRecord* result;
		
		for (result = m_data; result != nullptr; result = result->m_next)
		{
			if (pred(*result))
				break;
		}

		return result;
	}

	PDBRecord* FindRecord(const char* name)
	{
		return FindRecord([name](PDBRecord& rec) -> bool { return rec.is(name); });
	}

	// ----------------------------------------------------------------

	char v_c(size_t column) const
	{
		return m_rec->v_c(column);
	}
	
	string v_s(size_t column_first, size_t column_last = numeric_limits<size_t>::max()) const
	{
		return m_rec->v_s(column_first, column_last);
	}

	string v_f(size_t column_first, size_t column_last) const
	{
		return m_rec->v_f(column_first, column_last);
	}

	int v_i(int column_first, int column_last) const
	{
		int result = 0;
		try
		{
			result = m_rec->v_i(column_first, column_last);
		}
		catch (const exception& ex)
		{
			Error(ex.what());
		}
		return result;
	}

	// ----------------------------------------------------------------
	
	// Map a PDB residue location to a seqnum in a struct_asym 
	tuple<string,int,bool> MapResidue(char chainID, int resSeq, char iCode) const
	{
		auto key = make_tuple(chainID, resSeq, iCode);
		if (not m_chainSeq2AsymSeq.count(key))
			throw runtime_error(string("Residue ") + chainID + to_string(resSeq) + iCode + " could not be mapped"); 
		
		return m_chainSeq2AsymSeq.at(key);
	}

	tuple<string,int,bool> MapResidue(char chainID, int resSeq, char iCode, boost::system::error_code& ec) const
	{
		auto key = make_tuple(chainID, resSeq, iCode);
		
		tuple<string,int,bool> result;
		
		if (not m_chainSeq2AsymSeq.count(key))
		{
			ec = error::make_error_code(error::pdb_errors::residue_not_found);
			if (VERBOSE)
				cerr << "Residue " << chainID << resSeq << iCode << " could not be mapped" << endl;
		}
		else
			result = m_chainSeq2AsymSeq.at(key);
		
		return result;
	}
	
	tuple<string,int,bool> MapResidue(char chainID, int resSeq, char iCode, const string& resName);

	// ----------------------------------------------------------------
	
	void PreParseInput(istream& is);

	void GetNextRecord();
	void Match(const string& expected);
	void Error(const string& msg) const
	{
		string line_nr;
		if (m_rec != nullptr)
			line_nr = " (at line " + to_string(m_rec->m_line_nr) + ')';
		
		throw runtime_error("Error parsing PDB file" + line_nr + ": " + msg);
	}

	void ParseTitle();
	void ParseCitation(const string& id);
	void ParseRemarks();
	
//	void ParseRemark3();
//	size_t ParseRemark3(const string& program, const Remark3Template templ[], size_t N);
//	string NextRemark3Line();

	void ParseRemark200();
	void ParseRemark350();

	void ParsePrimaryStructure();
	void ParseHeterogen();
	void ConstructEntities();
	void ParseSecondaryStructure();
	void ParseConnectivtyAnnotation();
	void ParseMiscellaneousFeatures();
	void ParseCrystallographic();
	void ParseCoordinateTransformation();
	void ParseCoordinate(int model_nr);
	void ParseConnectivty();
	void ParseBookkeeping();

	// ----------------------------------------------------------------

	category* get_category(string name)
	{
		datablock::iterator i;
		std::tie(i, ignore) = m_datablock->emplace(name);
		return &*i;
	}
	
	vector<string> SplitCSV(const string& value);

	string pdb2cif_date(string s)
	{
		smatch m;
		const regex
			rx1(R"((\d{2})-(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)-(\d{2}))"),
			rx2(R"((JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)-(\d{2}))");

		if (regex_match(s, m, rx1))
		{
			using namespace boost::gregorian;
			
			int day = stoi(m[1].str());
			auto mi = kMonths.find(m[2].str());
			if (mi == kMonths.end())
				Error("Invalid month");
			int month = mi->second;
			int year = 1900 + stoi(m[3].str());
			if (year < 1950)
				year += 100;
		
			date date_original(year, month, day);
			
			s = to_iso_extended_string(date_original);
		}
		else if (regex_match(s, m, rx2))
		{
			auto mi = kMonths.find(m[1].str());
			if (mi == kMonths.end())
				Error("Invalid month");
			int month = mi->second;
			int year = 1900 + stoi(m[2].str());
			if (year < 1950)
				year += 100;
		
			s = (boost::format("%04d-%02d") % year % month).str();
		}
		
		return s;
	}
	
	string pdb2cif_auth(string author)
	{
		ba::trim(author);
		
		const regex rx(R"(((?:[A-Z]+\.)+)(.+))");
		smatch m;
		if (regex_match(author, m, rx))
			author = m[2].str() + ", " + m[1].str();
		
		bool upper = true;
		for (auto& c: author)
		{
			if (ispunct(c) or isspace(c))
				upper = true;
			else if (upper)
				upper = false;
			else
				c = tolower(c);
		} 

		return author;
	}

	string pdb2cif_symmetry(string s)
	{
		static const regex sg_rx(R"((\d{1,3})(\d{3}))");
		
		if (not s.empty())
		{
			smatch m;
			if (not regex_match(s, m, sg_rx))
				throw runtime_error("invalid symmetry value");
	
			s = m[1].str() + "_" + m[2].str();
		}
		
		return s;
	}
	
	string pdb2cif_charge(string c)
	{
		regex rx(R"((\d+)(\+|-))");
		smatch m;
		
		if (regex_match(c, m, rx))
		{
			if (m[2].str() == "-")
				c = '-' + m[1].str();
			else
				c = m[1].str();
		}

		return c;
	}
	
	string cif_id_for_int(int nr) const
	{
		string result;
		
		if (nr >= 26 * 26 * 26)
			result = 'L' + to_string(nr);
		else
		{
			if (nr >= 26 * 26)
			{
				int v = nr / (26 * 26);
				result += 'A' - 1 + v;
				nr %= (26 * 26);
			}
			
			if (nr >= 26)
			{
				int v = nr / 26;
				result += 'A' - 1 + v;
				nr %= 26;
			}
			
			result += 'A' + nr;
		}
		
		assert(not result.empty());
		return result;
	}
	
	vector<char> alt_locs_for_atom(char chainID, int seqNum, char iCode, string atomName);
	void MapChainID2AsymIDS(char chainID, vector<string>& asym_ids);

	// ----------------------------------------------------------------
	
	PDBRecord*	m_data;
	PDBRecord*	m_rec;
	cif::datablock*	m_datablock = nullptr;

	string		m_structure_id;
	string		m_original_date;
	string		m_exp_method = "X-RAY DIFFRACTION";
	int			m_citation_author_nr = 1, m_citation_editor_nr = 1;
	int			m_next_mol_id = 1, m_next_entity_nr = 1;
	int			m_next_software_ord = 1;

	struct SEQADV
	{
		string	resName;
		char	chainID;
		int		seqNum;
		char	iCode;
		string	database;
		string	dbAccession;
		string	dbRes;
		int		dbSeq;
		string	conflict;
	};

	vector<SEQADV>						m_seqadvs;
	
	list<PDBCompound>					m_compounds;
	list<PDBChain>						m_chains;
	vector<HET> 						m_hets;
	map<string,string>					m_hetnams;
	map<string,string>					m_hetsyns;
	map<string,string>					m_formuls;
	string								m_water_het_id;
	vector<string>						m_chem_comp, m_atom_types;
	
	map<string,string>					m_remark200;
	string								m_refinement_software;
	int									m_atom_id = 0;
	int									m_pdbx_dif_ordinal = 0;

	vector<UNOBS>						m_unobs;

	// various maps between numbering schemes
	map<tuple<char,int,char>,tuple<string,int,bool>>	m_chainSeq2AsymSeq;

	map<int,string>							m_molID2EntityID;
	map<string,string>						m_het2EntityID;
	map<string,string>						m_asymID2EntityID;
	map<string,string>						m_mod2parent;
};

// --------------------------------------------------------------------

vector<char> PDBFileParser::alt_locs_for_atom(char inChainID, int inResSeq, char inICode, string inAtomName)
{
	// well, maybe this could be optimized...
	set<char> result;
	
	for (auto r = m_data; r != nullptr; r = r->m_next)
	{
		if (r->is("ATOM  ") or r->is("HETATM"))		//	 1 -  6        Record name   "ATOM  "
		{											//	 ...
			string name			= r->v_s(13, 16);	//	13 - 16        Atom          name         Atom name.
			char altLoc			= r->v_c(17);		//	17             Character     altLoc       Alternate location indicator.
			char chainID		= r->v_c(22);		//	22             Character     chainID      Chain identifier.
			int resSeq			= r->v_i(23, 26);	//	23 - 26        Integer       resSeq       Residue sequence number.
			char iCode			= r->v_c(27);		//	27             AChar         iCode        Code for insertion of residues.

			if (chainID == inChainID and resSeq == inResSeq and iCode == inICode and name == inAtomName and altLoc != ' ')
				result.insert(altLoc);
		}
	}
	
	return { result.begin(), result.end() };
}

void PDBFileParser::MapChainID2AsymIDS(char chainID, vector<string>& asym_ids)
{
	struct l : binary_function<string,string,bool>
	{
		bool operator()(const string& a, const string& b) const
		{
			int d = a.length() - b.length();
			if (d == 0)
				d = a.compare(b);
			return d < 0;
		}
	};
	
	set<string,l> asym(asym_ids.begin(), asym_ids.end());

	for (auto& m: m_chainSeq2AsymSeq)
	{
		if (get<0>(m.first) == chainID)
			asym.insert(get<0>(m.second));
	}
	
	asym_ids.assign(asym.begin(), asym.end());
}

// --------------------------------------------------------------------

void PDBFileParser::PreParseInput(istream& is)
{
	string lookahead;
	uint32 line_nr = 1;
	getline(is, lookahead);

	if (ba::starts_with(lookahead, "HEADER") == false)
		Error("This does not look like a PDB file, should start with a HEADER line");

	PDBRecord* last = nullptr;
	set<string> dropped;

	for (;;)
	{
		if (lookahead.empty())
			break;
	
		string type = lookahead.substr(0, 6);
		string value;
		if (lookahead.length() > 6)
			value = ba::trim_right_copy(lookahead.substr(6));
	
		getline(is, lookahead);

		uint32 cur_line_nr = line_nr;
		++line_nr;
		
		if (kSupportedRecords.count(type) == 0)
		{
			ba::trim(type);
			
			if (type != "END")	// special case
				dropped.insert(type);
			continue;
		}
		
		// see if we need to append continuation values
		if (type == "AUTHOR" or
			type == "EXPDTA" or
			type == "MDLTYP" or
			type == "KEYWDS" or
			type == "SPLIT " or
			type == "SPRSDE" or
			type == "TITLE ")
		{
			int n = 2;
			while (lookahead.substr(0, 6) == type and stoi(lookahead.substr(7, 3)) == n)
			{
				value += ba::trim_right_copy(lookahead.substr(10));
				getline(is, lookahead);
				++line_nr;
				++n;
			}
		}
		else if (type == "COMPND")
		{
			int n = 2;
			value += '\n';
			while (lookahead.substr(0, 6) == type and stoi(lookahead.substr(7, 3)) == n)
			{
				value += ba::trim_right_copy(lookahead.substr(10));
				value += '\n';
				getline(is, lookahead);
				++line_nr;
				++n;
			}
		}
		else if (type == "REVDAT")
		{
			int rev_nr = stoi(value.substr(1, 3));
			int n = 2;
			while (lookahead.substr(0, 6) == type and
				stoi(lookahead.substr(7, 3)) == rev_nr and
				stoi(lookahead.substr(10, 2)) == n)
			{
				value += lookahead.substr(38);
				getline(is, lookahead);
				++line_nr;
				++n;
			}
		}
		else if (type == "CAVEAT")
		{
			int n = 2;
			while (lookahead.substr(0, 6) == type and stoi(lookahead.substr(7, 3)) == n)
			{
				value += ba::trim_right_copy(lookahead.substr(13));
				getline(is, lookahead);
				++line_nr;
				++n;
			}
		}
		else if (type == "OBSLTE")
		{
			while (lookahead.substr(0, 6) == type)
			{
				value += lookahead.substr(31);
				getline(is, lookahead);
				++line_nr;
			}
		}
		else if (type == "SOURCE")
		{
			value += '\n';
			int n = 2;
			while (lookahead.substr(0, 6) == type and stoi(lookahead.substr(7, 3)) == n)
			{
				value += ba::trim_copy(lookahead.substr(10));
				value += '\n';
				getline(is, lookahead);
				++line_nr;
				++n;
			}
		}
		else if (type == "FORMUL")
		{
			int comp_nr = stoi(value.substr(1, 3));
			int n = 2;
			while (lookahead.substr(0, 6) == type and
				stoi(lookahead.substr(7, 3)) == comp_nr and
				lookahead.substr(16, 2) != "  " and
				stoi(lookahead.substr(16, 2)) == n)
			{
				value += lookahead.substr(19);
				getline(is, lookahead);
				++line_nr;
				++n;
			}
		}
		else if (type == "HETNAM" or
				 type == "HETSYN")
		{
			int n = 2;
			while (lookahead.substr(0, 6) == type and
				lookahead.substr(8, 2) != "  " and
				stoi(lookahead.substr(8, 2)) == n)
			{
				value += lookahead.substr(16);
				getline(is, lookahead);
				++line_nr;
				++n;
			}
		}
		else if (type == "SITE  ")
		{
			string site_name = value.substr(5, 3);
			ba::trim_right(value);
			size_t n = value.length() - 12;
			value += string(11 - (n % 11), ' ');
			
			while (lookahead.substr(0, 6) == type and lookahead.substr(11, 3) == site_name)
			{
				string s = lookahead.substr(18);
				ba::trim_right(s);
				s += string(11 - (s.length() % 11), ' ');
				value += s;
				
// TODO: improve this... either use numRes or don't lump together all text
//				value += " " + ba::trim_right_copy();
				getline(is, lookahead);
				++line_nr;
			}
		}
		else if (type == "REMARK")
		{
			type += value.substr(0, 4);

			// parse it now, makes life easier later on
			if (type == "REMARK 200" or type == "REMARK 240")
			{
				auto i = value.find(":");
		
				if (i != string::npos)
				{
					string k = value.substr(4, i - 4);
					string v = value.substr(i + 1);
					
					ba::trim(k);
					while (k.find("  ") != string::npos)
						ba::replace_all(k, "  ", " ");
					ba::trim(v);
					
					if (iequals(v, "NONE"))
						m_remark200[k] = ".";
					else if (not iequals(v, "NULL"))
						m_remark200[k] = v;
				}
			}
		}
		
		PDBRecord* cur = new(value.length()) PDBRecord(cur_line_nr, type, value);
		
		if (last == nullptr)
			last = m_data = cur;
		else
			last->m_next = cur;

		last = cur;
	}
	
	if (not dropped.empty())
	{
		cerr << "Dropped unsupported records: " << ba::join(dropped, ", ") << endl;
	}
	
	m_rec = m_data;
}

void PDBFileParser::GetNextRecord()
{
	if (m_rec != nullptr)
		m_rec = m_rec->m_next;

	if (m_rec == nullptr)
	{
		static PDBRecord* end = new(0)PDBRecord({ 0, "END   ", ""});
		m_rec = end;
	}
}

void PDBFileParser::Match(const string& expected)
{
	if (m_rec->m_name != expected)
		Error("Expected record " + expected + " but found " + m_rec->m_name);
}

vector<string> PDBFileParser::SplitCSV(const string& value)
{
	vector<string> vs;
	ba::split(vs, value, ba::is_any_of(","));
	for (auto& v: vs)
		ba::trim(v);
	return vs;
}

void PDBFileParser::ParseTitle()
{
	// strict ordering required
	
	// HEADER
	//	 1 -  6       Record name    "HEADER"
	//	11 - 50       String(40)     classification    Classifies the molecule(s).
	//	51 - 59       Date           depDate           Deposition date. This is the date the
	//	                                               coordinates  were received at the PDB.
	//	63 - 66       IDcode         idCode            This identifier is unique within the PDB.

	Match("HEADER");
	
	m_structure_id		= v_s(63, 66);
	string keywords		= v_s(11, 50);
	m_original_date		= pdb2cif_date(v_s(51, 59));

	m_datablock = new cif::datablock(m_structure_id);
	
	ba::trim(keywords);
	
	auto cat = get_category("entry");
//	cat->add_column("id");
	cat->emplace({ {"id", m_structure_id} });
	
	GetNextRecord();
	
	// OBSLTE
	if (m_rec->is("OBSLTE"))
	{
		//	 1 -  6       Record name   "OBSLTE"
		//	 9 - 10       Continuation  continuation  Allows concatenation of multiple records
		//	12 - 20       Date          repDate       Date that this datablock was replaced.
		//	22 - 25       IDcode        idCode        ID code of this datablock.
		//	32 - 35       IDcode        rIdCode       ID code of datablock that replaced this one.
		//	37 - 40       ...


		string old		= v_s(22, 25);
		string date		= pdb2cif_date(v_s(12, 20));
		cat = get_category("pdbx_database_PDB_obs");
		
		string value = m_rec->v_s(32);
		for (auto i = make_split_iterator(value, ba::token_finder(ba::is_any_of(" "), ba::token_compress_on)); not i.eof(); ++i)
		{
			cat->emplace({
				{ "id", "OBSLTE" },
				{ "date", date },
				{ "replace_pdb_id", old },
				{ "pdb_id", string(i->begin(), i->end()) }
			});
		}

		GetNextRecord();
	}
	
	// TITLE
//	Match("TITLE ");
	string title;
	if (m_rec->is("TITLE "))	//	 1 -  6       Record name    "TITLE "
	{							//	 9 - 10       Continuation   continuation  Allows concatenation of multiple records.
		title = v_s(11);		//	11 - 80       String         title         Title of the  experiment.
		GetNextRecord();
	}
	
	// SPLIT
	if (m_rec->is("SPLIT"))
	{
		//	 1 -  6        Record  name  "SPLIT "
		//	 9 - 10        Continuation  continuation  Allows concatenation of multiple records.
		//	12 - 15        IDcode        idCode        ID code of related datablock.
		if (VERBOSE)	
			Error("skipping unimplemented SPLIT record");
		GetNextRecord();
	}
	
	// CAVEAT
	if (m_rec->is("CAVEAT"))							//	  1 - 6       Record name   "CAVEAT"
	{													//	 9 - 10       Continuation  continuation   Allows concatenation of multiple records.
		get_category("database_PDB_caveat")->emplace({
			{ "id", v_s(12, 15) },						//	12 - 15       IDcode        idCode         PDB ID code of this datablock.                  
			{ "text", string{m_rec->v_s(20) } }    		//	20 - 79       String        comment        Free text giving the reason for the  CAVEAT.
		});
		
		GetNextRecord();
	}
	
	// COMPND
	Match("COMPND");
	//	 1 -  6       Record name     "COMPND"   
	//	 8 - 10       Continuation    continuation  Allows concatenation of multiple records.
	//	11 - 80       Specification   compound      Description of the molecular components.
	//	              list   

	string value{m_rec->v_s(11)};
	if (value.find(':') == string::npos)
	{
			// special case for dumb, stripped files		
		auto& comp = GetOrCreateCompound(1);
		comp.m_info["MOLECULE"] = value;
	}
	else
	{
		SpecificationListParser p(value);
		
		for (;;)
		{
			string key, value;
			std::tie(key, value) = p.GetNextSpecification();

			if (key.empty())
				break;
		
			if (not iequals(key, "MOL_ID") and m_compounds.empty())
			{
				cerr << "Ignoring invalid COMPND record" << endl;
				break;
			}
			
			if (key == "MOL_ID")
			{
				auto& comp = GetOrCreateCompound(stoi(value));
				comp.m_title = title;
			}
			else if (key == "CHAIN")
			{
				vector<string> chains;
				
				ba::split(chains, value, ba::is_any_of(","));
				for (auto& c: chains)
				{
					ba::trim(c);
					m_compounds.back().m_chains.insert(c[0]);
				}
			}
			else
				m_compounds.back().m_info[key] = value;
		}	
	}

	GetNextRecord();

	// SOURCE
//	Match("SOURCE");

	if (m_rec->is("SOURCE"))
	{
		//	 1 -  6      Record name    "SOURCE"       
		//	 8 - 10      Continuation   continuation   Allows concatenation of multiple records.
		//	11 - 79      Specification  srcName        Identifies the source of the
		//	             List                          macromolecule in a  token: value format.
		
		map<string,string>* source = nullptr;
	
//		value = { m_rec->v_s(11) };
//		for (auto si = ba::make_split_iterator(value, ba::token_finder(ba::is_any_of(";"), ba::token_compress_on)); not si.eof(); ++si)
//		{
//			string s(si->begin(), si->end());
//			if (s.empty())
//				continue;
//			
//			auto colon = s.find(": ");
//			if (colon == string::npos)
//			{
//				if (VERBOSE)
//					cerr << "invalid source field, missing colon (" << s << ')' << endl;
//				continue;
//			}
		SpecificationListParser p(v_s(11));
		
		for (;;)
		{
			string key, value;
			std::tie(key, value) = p.GetNextSpecification();

			if (key.empty())
				break;
			
			if (key == "MOL_ID")
			{
				for (auto& c: m_compounds)
				{
					if (c.m_mol_id == stoi(value))
					{
						source = &c.m_source;
						break;
					}
				}
				
				continue;
			}
			
			if (source == nullptr)
				Error("MOL_ID missing");
			
			(*source)[key] = value;
		}
		
		GetNextRecord();
	}

	// KEYWDS
//	Match("KEYWDS");
	string pdbx_keywords;

	if (m_rec->is("KEYWDS"))		//	 1 -  6       Record name    "KEYWDS"  
	{								//	 9 - 10       Continuation   continuation  Allows concatenation of records if necessary.
		pdbx_keywords = v_s(11);	//	11 - 79       List           keywds        Comma-separated list of keywords relevant
									//	                                           to the datablock.       
		GetNextRecord();
	}
	
	if (not (keywords.empty() and pdbx_keywords.empty()))
	{
		get_category("struct_keywords")->emplace({
			{ "entry_id",  m_structure_id },
			{ "pdbx_keywords", keywords },
			{ "text", pdbx_keywords }
		});
	}

	// EXPDTA
//	Match("EXPDTA");
	if (m_rec->is("EXPDTA"))
	{
		m_exp_method = v_s(11);
		
		cat = get_category("exptl");
		cat->emplace({
			{ "entry_id", m_structure_id },
			{ "method", m_exp_method },
			{ "crystals_number", m_remark200["NUMBER OF CRYSTALS USED"] }
		});
		
		GetNextRecord();
	}

	// NUMMDL
	if (m_rec->is("NUMMDL"))
	{
		if (VERBOSE)	
			Error("skipping unimplemented NUMMDL record");
		GetNextRecord();
	}
	
	// MDLTYP
	if (m_rec->is("MDLTYP"))
	{
		if (VERBOSE)	
			Error("skipping unimplemented MDLTYP record");
		GetNextRecord();
	}

	// AUTHOR
//	Match("AUTHOR");
	if (m_rec->is("AUTHOR"))
	{
		int n = 1;
		cat = get_category("audit_author");
		
		value = { m_rec->v_s(11) };
		for (auto si = ba::make_split_iterator(value, ba::token_finder(ba::is_any_of(","), ba::token_compress_on)); not si.eof(); ++si)
		{
			string author(si->begin(), si->end());
			
			cat->emplace({
				{ "name", pdb2cif_auth(author) },
				{ "pdbx_ordinal", n }
			});
			++n;
		}
	
		GetNextRecord();
	}

	// REVDAT
	bool firstRevDat = true;
	struct RevDat {
		int rev_num;
		string date, date_original, replaces;
		int mod_type;
		vector<string> types;

		bool operator<(const RevDat& rhs) const { return rev_num < rhs.rev_num; }
	};
	vector<RevDat> revdats;
	
	while (m_rec->is("REVDAT"))
	{
													//	 1 -  6       Record name    "REVDAT"                                             
		int rev_num = v_i(8, 10);					//	 8 - 10       Integer        modNum        Modification number.                   
													//	11 - 12       Continuation   continuation  Allows concatenation of multiple records.
		string date = pdb2cif_date(v_s(14, 22));	//	14 - 22       Date           modDate       Date of modification (or release  for   
													//	                                           new entries)  in DD-MMM-YY format. This is
													//	                                           not repeated on continued lines.
		string modId = v_s(24, 27);					//	24 - 27       IDCode         modId         ID code of this datablock. This is not repeated on 
													//	                                           continuation lines.    
		int mod_type = v_i(32, 32);					//	32            Integer        modType       An integer identifying the type of    
													//	                                           modification. For all  revisions, the
													//	                                           modification type is listed as 1 
		string value = v_s(40);						//	40 - 45       LString(6)     record        Modification detail. 
													//	47 - 52       LString(6)     record        Modification detail. 
													//	54 - 59       LString(6)     record        Modification detail. 
													//	61 - 66       LString(6)     record        Modification detail.
		
		revdats.push_back({rev_num, date, mod_type == 0 ? m_original_date : "", modId, mod_type});
		
		ba::split(revdats.back().types, value, ba::is_any_of(" "));

		if (firstRevDat)
		{
			cat = get_category("database_2");
			cat->emplace({
				{ "database_id", "PDB" },
				{ "database_code", modId }
			});
		}

		GetNextRecord();
		firstRevDat = false;
	}
	
/*
	This is internal stuff for PDB, don't write it

	sort(revdats.begin(), revdats.end());
	for (auto& revdat: revdats)
	{
		get_category("database_PDB_rev")->emplace({
			{ "num",			revdat.rev_num },
			{ "date",			revdat.date },
			{ "date_original",	revdat.date_original },
			{ "replaces", 		revdat.replaces },
			{ "mod_type",		revdat.mod_type }
		});
		
		for (auto& type: revdat.types)
		{
			if (type.empty())
				continue;
			get_category("database_PDB_rev_record")->emplace({
				{ "rev_num",	revdat.rev_num  },
				{ "type",		type }
			});
		}
	}
*/

	// SPRSDE
	if (m_rec->is("SPRSDE"))
	{
		if (VERBOSE)	
			cerr << "skipping unimplemented SPRSDE record" << endl;
		GetNextRecord();
	}
	
	// JRNL
	if (m_rec->is("JRNL  "))
		ParseCitation("primary");
}

void PDBFileParser::ParseCitation(const string& id)
{
	const char* rec = m_rec->m_name;
	
	string auth, titl, edit, publ, refn, pmid, doi;
	string pubname, volume, astm, country, issn, csd;
	int page_first = 0, page_last = 0, year = 0;

	auto extend = [](string& s, const string& p)
	{
		if (not s.empty())	
			s += ' ';
		s += ba::trim_copy(p);
	};

	while (m_rec->is(rec) and (id == "primary" or v_c(12) == ' '))
	{
		string k = v_s(13, 16);
		if (k == "AUTH")				extend(auth, v_s(20, 79));
		else if (k == "TITL")			extend(titl, v_s(20, 79));
		else if (k == "EDIT")			extend(edit, v_s(20, 79));
		else if (k == "REF")
		{
			if (pubname.empty())
			{
				extend(pubname, v_s(20, 47));
				if (v_s(50, 51) == "V.")
					volume = ba::trim_copy(v_s(52, 55));
				page_first = v_i(57, 61);
				year = v_i(63, 66);
			}
			else
				extend(pubname, v_s(20, 47));
		}
		else if (k == "PUBL")			extend(publ, v_s(20, 70));
		else if (k == "REFN")
		{
			if (v_s(20, 23) == "ASTN")
				astm = v_s(25, 30);
			country = v_s(33, 34);
			if (v_s(36, 39) == "ISSN")
				issn = v_s(41, 65);
		}
		else if (k == "PMID")			pmid = v_s(20, 79);
		else if (k == "DOI")			doi = v_s(20, 79);
		
		GetNextRecord();
	}
	
	auto cat = get_category("citation");
	cat->emplace({
		{ "id", id },
		{ "title", titl },
		{ "journal_abbrev", pubname },
		{ "journal_volume", volume },
		{ "page_first", page_first > 0 ? to_string(page_first) : "" },
		{ "page_last", page_last > 0 ? to_string(page_last) : "" },
		{ "year", year > 0 ? to_string(year) : "" },
		{ "journal_id_ASTM", astm },
		{ "country", country },
		{ "journal_id_ISSN", issn },
		{ "journal_id_CSD", csd },
		{ "book_publisher", publ },
		{ "pdbx_database_id_PubMed", pmid },
		{ "pdbx_database_id_DOI", doi }
	});

	if (not auth.empty())
	{
		cat = get_category("citation_author");
		for (auto si = ba::make_split_iterator(auth, ba::token_finder(ba::is_any_of(","), ba::token_compress_on)); not si.eof(); ++si)
		{
			string author(si->begin(), si->end());
			
			cat->emplace({
				{ "citation_id", id },
				{ "name", pdb2cif_auth(author) },
				{ "ordinal", m_citation_author_nr }
			});

			++m_citation_author_nr;
		}
	}

	if (not edit.empty())
	{
		cat = get_category("citation_editor");
		for (auto si = ba::make_split_iterator(edit, ba::token_finder(ba::is_any_of(","), ba::token_compress_on)); not si.eof(); ++si)
		{
			string editor(si->begin(), si->end());
			
			cat->emplace({
				{ "citation_id", id },
				{ "name", pdb2cif_auth(editor) },
				{ "ordinal", m_citation_editor_nr }
			});

			++m_citation_editor_nr;
		}
	}
}

void PDBFileParser::ParseRemarks()
{
	string sequence_details, compound_details, source_details;
	
	while (ba::starts_with(m_rec->m_name, "REMARK"))
	{
		int remarkNr = v_i(8, 10);
		
		switch (remarkNr)
		{
			case 1:
				while (m_rec->is("REMARK   1"))
				{
					if (m_rec->m_vlen > 15 and v_s(12, 20) == "REFERENCE")
					{
						string id = v_s(22, 70);
						GetNextRecord();
						
						ParseCitation(id);
					}
					else
						GetNextRecord();
				}
				break;
			
			case 3:
				// we skip REMARK 3 until we know the mapping
				while (m_rec->is("REMARK   3"))
					GetNextRecord();
				break;

			case 4:
				// who cares...
				while (m_rec->is("REMARK   4"))
					GetNextRecord();
				break;

			case 100:
			{
				const regex rx(R"(THE (\S+) ID CODE IS (\S+?)\.?\s*)");
				smatch m;
				string r = v_s(12);
				
				if (regex_match(r, m, rx))
				{
					auto cat = get_category("database_2");
					cat->emplace({
						{ "database_id", m[1].str() },
						{ "database_code", m[2].str() }
					});
				}
				
				GetNextRecord();
				break;
			}
			
			case 200:
			{
				// we already parsed most of this remark, but the "REMARK:" part might have been split.

				bool remark = false;
				
				do
				{
					string r = m_rec->v_s(12);

					if (ba::starts_with(r, "REMARK: "))
					{
						m_remark200["REMARK"] = r.substr(8);
						remark = true;
					}
					else if (remark)
					{
						if (r.empty())
							remark = false;
						else
							m_remark200["REMARK"] += r;
					}
					
					GetNextRecord();
				}
				while (m_rec->is("REMARK 200"));
				break;
			}

			case 280:
			{
				string density_Matthews, density_percent_sol, conditions;
				
				const regex rx1(R"(SOLVENT CONTENT, VS +\(%\): *(.+))"),
					rx2(R"(MATTHEWS COEFFICIENT, VM \(ANGSTROMS\*\*3/DA\): *(.+))");
				
				smatch m;

				do
				{
					string r = v_s(12);

					if (conditions.empty())
					{
						if (regex_match(r, m, rx1))
							density_percent_sol = m[1].str();
						else if (regex_match(r, m, rx2))
							density_Matthews = m[1].str();
						else if (ba::starts_with(r, "CRYSTALLIZATION CONDITIONS: "))
							conditions = r.substr(28);
					}
					else
						conditions = conditions + ' ' + r;
					
					GetNextRecord();
				}
				while (m_rec->is("REMARK 280"));
				
				string desc = m_remark200["REMARK"];
				if (desc == "NULL")
					desc.clear();
				
				get_category("exptl_crystal")->emplace({
					{ "id", 1 },
					{ "density_Matthews", iequals(density_Matthews, "NULL") ? "" : density_Matthews  },
					{ "density_percent_sol", iequals(density_percent_sol, "NULL") ? "" : density_percent_sol },
					{ "description", desc }
				});
				
				// now try to parse the conditions
				const regex rx3(R"(TEMPERATURE +(\d+)K)"), rx4(R"(PH *(?:: *)?(\d+(?:\.\d+)?))")/*, rx5(R"(\b(\d+)C\b)")*/;
				
				string temp, ph, method;

				for (auto i = make_split_iterator(conditions, ba::token_finder(ba::is_any_of(","), ba::token_compress_on)); not i.eof(); ++i)
				{
					string s(i->begin(), i->end());

					ba::trim(s);

					if (regex_search(s, m, rx3))
						temp = m[1].str();
					if (regex_search(s, m, rx4))
						ph = m[1].str();
					if (s.length() < 60 and
						(ba::icontains(s, "drop") or ba::icontains(s, "vapor") or ba::icontains(s, "batch")))
					{
						if (not method.empty())
							method = method + ", " + s;
						else
							method = s;
					}
				}
				
				if (not (method.empty() and temp.empty() and ph.empty() and (conditions.empty() or conditions == "NULL")))
				{
					get_category("exptl_crystal_grow")->emplace({
						{ "crystal_id", 1 },
						{ "method", method },
						{ "temp", temp },
						{ "pH", ph },
						{ "pdbx_details", conditions }
					});
				}

				break;
			}
			
//			case 290:
//				
//				break;
			
			case 350:
				// postponed since we don't have the required information yet
				for (; m_rec->is("REMARK 350"); GetNextRecord())
					;
				break;
			
			case 400:
			{
				stringstream s;
				GetNextRecord();
				if (v_s(12) == "COMPOUND")
					GetNextRecord();
				
				while (m_rec->is("REMARK 400"))
				{
					s << v_s(12) << endl;
					GetNextRecord();
				}
				
				compound_details = s.str();
				break;
			}
			
			case 450:
			{
				stringstream s;
				GetNextRecord();
				if (v_s(12) == "SOURCE")
					GetNextRecord();
				
				while (m_rec->is("REMARK 450"))
				{
					s << v_s(12) << endl;
					GetNextRecord();
				}
				
				source_details = s.str();
				break;
			}
			
			case 465:
			{
				bool headerSeen = false;
				regex rx(R"( *MODELS *(\d+)-(\d+))");
				int models[2] = { -1, -1 };
				
				for (; m_rec->is("REMARK 465"); GetNextRecord())
				{
					if (not headerSeen)
					{
						string line = v_s(12);
						smatch m;

						if (regex_match(line, m, rx))
						{
							models[0] = stoi(m[1].str());
							models[1] = stoi(m[2].str());
						}
						else
							headerSeen = ba::contains(line, "RES C SSSEQI");
						continue;
					}
					
					if (models[0] == models[1])
						models[0] = models[1] = v_i(12, 14);
					
					string res = v_s(16, 18);
					char chain = v_c(20);
					int seq = v_i(22, 26);
					char iCode = v_c(27);
					
					for (int model_nr = models[0]; model_nr <= models[1]; ++model_nr)
						m_unobs.push_back({model_nr, res, chain, seq, iCode});
				}
				
				break;
			}
			
			case 470:
			{
				bool headerSeen = false;
				regex rx(R"( *MODELS *(\d+)-(\d+))");
				int models[2] = { -1, -1 };
				
				for (; m_rec->is("REMARK 470"); GetNextRecord())
				{
					if (not headerSeen)
					{
						string line = v_s(12);
						smatch m;

						if (regex_match(line, m, rx))
						{
							models[0] = stoi(m[1].str());
							models[1] = stoi(m[2].str());
						}
						else
							headerSeen = ba::contains(line, "RES CSSEQI  ATOMS");
						continue;
					}
					
					if (models[0] == models[1])
						models[0] = models[1] = v_i(12, 14);
					
					string res = v_s(16, 18);
					char chain = v_c(20);
					int seq = v_i(21, 25);
					char iCode = v_c(26);
					
					vector<string> atoms;
					string atomStr = m_rec->v_s(29);
					for (auto i = make_split_iterator(atomStr, ba::token_finder(ba::is_any_of(" "), ba::token_compress_on)); not i.eof(); ++i)
						atoms.push_back({ i-> begin(), i->end() });
					
					for (int model_nr = models[0]; model_nr <= models[1]; ++model_nr)
						m_unobs.push_back({model_nr, res, chain, seq, iCode, atoms});
				}
				
				break;
			}
			
			case 500:
			{
				GetNextRecord();
				
				enum State { e_start, e_CCinSAU, e_CC, e_CBL, e_CBA, e_TA, e_CTg, e_PG, e_MCP, e_ChC } state = e_start;
				bool headerSeen = false;
				int id = 0;
				
				for (; m_rec->is("REMARK 500"); GetNextRecord())
				{
					string line = v_s(12);
					
					if (line == "GEOMETRY AND STEREOCHEMISTRY")
						continue;
					
					switch (state)
					{
						case e_start:
						{
							if (line.empty() or not ba::starts_with(line, "SUBTOPIC: "))
								continue;

							string subtopic = line.substr(10);

							if (subtopic == "CLOSE CONTACTS IN SAME ASYMMETRIC UNIT")	state = e_CCinSAU;
							else if (subtopic == "CLOSE CONTACTS")						state = e_CC;
							else if (subtopic == "COVALENT BOND LENGTHS")				state = e_CBL;
							else if (subtopic == "COVALENT BOND ANGLES")				state = e_CBA;
							else if (subtopic == "TORSION ANGLES")						state = e_TA;
							else if (subtopic == "NON-CIS, NON-TRANS")					state = e_CTg;
							else if (subtopic == "PLANAR GROUPS")						state = e_PG;
							else if (subtopic == "MAIN CHAIN PLANARITY")				state = e_MCP;
							else if (subtopic == "CHIRAL CENTERS")						state = e_ChC;
							else if (VERBOSE)
								Error("Unknown subtopic in REMARK 500: " + subtopic);
							
							headerSeen = false;
							id = 0;
							break;
						}
						
						case e_CCinSAU:
						{
							if (not headerSeen)
								headerSeen =
									line == "ATM1  RES C  SSEQI   ATM2  RES C  SSEQI           DISTANCE";
							else if (line.empty())
								state = e_start;
							else
							{
								string atom_1 = v_s(13, 16);
								string res_1 = v_s(19, 21);
								string alt_1 = v_s(17, 17);
								char chain_1 = v_c(23);
								int seq_1 = v_i(25, 29);
								string iCode_1 = v_s(30, 30);
								
								string atom_2 = v_s(34, 37);
								string alt_2 = v_s(38, 38);
								string res_2 = v_s(40, 42);
								char chain_2 = v_c(44);
								int seq_2 = v_i(46, 50);
								string iCode_2 = v_s(51, 51);
								
								string distance = v_f(63, 71);
								
								get_category("pdbx_validate_close_contact")->emplace({
									{ "id",				to_string(++id) },
									{ "PDB_model_num",	1 },
									{ "auth_atom_id_1",	atom_1 },
									{ "auth_asym_id_1", string{ chain_1 } },
									{ "auth_comp_id_1", res_1 },
									{ "auth_seq_id_1",	seq_1 },
									{ "PDB_ins_code_1", iCode_1 },
									{ "label_alt_id_1", alt_1 },
									{ "auth_atom_id_2",	atom_2 },
									{ "auth_asym_id_2", string { chain_2 } },
									{ "auth_comp_id_2", res_2 },
									{ "auth_seq_id_2",	seq_2 },
									{ "PDB_ins_code_2", iCode_2 },
									{ "label_alt_id_2", alt_2 },
									{ "dist", distance }
								});
							}
							break;
						}
						
						case e_CC:
						{
							if (not headerSeen)
								headerSeen = line == "ATM1  RES C  SSEQI   ATM2  RES C  SSEQI  SSYMOP   DISTANCE";
							else if (line.empty())
								state = e_start;
							else
							{
								string atom_1 = v_s(13, 16);
								string res_1 = v_s(19, 21);
								char chain_1 = v_c(23);
								int seq_1 = v_i(25, 29);
								
								string atom_2 = v_s(34, 37);
								string res_2 = v_s(40, 42);
								char chain_2 = v_c(44);
								int seq_2 = v_i(46, 50);
								
								string symop = pdb2cif_symmetry(v_s(54, 59));
								
								string distance = v_f(63, 71);
								
								get_category("pdbx_validate_symm_contact")->emplace({
									{ "id",				to_string(++id) },
									{ "PDB_model_num",	1 },
									{ "auth_atom_id_1",	atom_1 },
									{ "auth_asym_id_1", string{ chain_1 } },
									{ "auth_comp_id_1", res_1 },
									{ "auth_seq_id_1",	seq_1 },
//									{ "PDB_ins_code_1", "" },
//									{ "label_alt_id_1", "" },
									{ "site_symmetry_1", "1_555" },
									{ "auth_atom_id_2",	atom_2 },
									{ "auth_asym_id_2", string { chain_2 } },
									{ "auth_comp_id_2", res_2 },
									{ "auth_seq_id_2",	seq_2 },
//									{ "PDB_ins_code_2", "" },
//									{ "label_alt_id_2", "" },
									{ "site_symmetry_2", symop },
									{ "dist", distance }
								});
							}
							break;
						}
						
						case e_CBL:
						{
							if (not headerSeen)
							{
								if (ba::starts_with(line, "FORMAT: ") and line != "FORMAT: (10X,I3,1X,2(A3,1X,A1,I4,A1,1X,A4,3X),1X,F6.3)")
									Error("Unexpected format in REMARK 500");
								
								headerSeen = line == "M RES CSSEQI ATM1   RES CSSEQI ATM2   DEVIATION";
							}
							else if (line.empty())
								state = e_start;
							else
							{
								int model = v_i(11, 13);
								string resNam1 = v_s(15, 17);
								string chainID1 { v_c(19) };
								int seqNum1 = v_i(20, 24);
								string iCode1 { v_c(25) };
								string alt1 = v_s(30, 30);
								string atm1 = v_s(26, 29);

								string resNam2 = v_s(33, 35);
								string chainID2 { v_c(37) };
								int seqNum2 = v_i(38, 41);
								string iCode2 { v_c(42) };
								string alt2 = v_s(48, 48);
								string atm2 = v_s(44, 47);
								
								string deviation = v_f(51, 57);
								
								if (iCode1 == " ") iCode1.clear();
								if (iCode2 == " ") iCode2.clear();
								
								get_category("pdbx_validate_rmsd_bond")->emplace({
									{ "id",				to_string(++id) },
									{ "PDB_model_num",	model ? model : 1 },
									{ "auth_atom_id_1",	atm1 },
									{ "auth_asym_id_1", chainID1 },
									{ "auth_comp_id_1", resNam1 },
									{ "auth_seq_id_1",	seqNum1 },
									{ "PDB_ins_code_1",	iCode1 },
									{ "label_alt_id_1", alt1 },
									{ "auth_atom_id_2",	atm2 },
									{ "auth_asym_id_2", chainID2 },
									{ "auth_comp_id_2", resNam2 },
									{ "auth_seq_id_2",	seqNum2 },
									{ "PDB_ins_code_2",	iCode2 },
									{ "label_alt_id_2", alt2 },
									{ "bond_deviation",deviation }
								});
							}
							
							break;
						}
						
						case e_CBA:
							if (not headerSeen)
							{
								if (ba::starts_with(line, "FORMAT: ") and line != "FORMAT: (10X,I3,1X,A3,1X,A1,I4,A1,3(1X,A4,2X),12X,F5.1)")
									Error("Unexpected format in REMARK 500");
								
								headerSeen = line == "M RES CSSEQI ATM1   ATM2   ATM3";
							}
							else if (line.empty())
								state = e_start;
							else if (v_s(64) == "DEGREES")
							{
								int model = v_i(11, 13);
								string resNam = v_s(15, 17);
								string chainID { v_c(19) };
								int seqNum = v_i(20, 24);
								string iCode { v_c(25) };
								
								if (iCode == " ")
									iCode.clear();
								
								string atoms[3] = { v_s(27, 30), v_s(34, 37), v_s(41, 44) };
								string deviation = v_f(57, 62);
								
								get_category("pdbx_validate_rmsd_angle")->emplace({
									{ "id",				to_string(++id) },
									{ "PDB_model_num",	model ? model : 1 },
									{ "auth_atom_id_1",	atoms[0] },
									{ "auth_asym_id_1", chainID },
									{ "auth_comp_id_1", resNam },
									{ "auth_seq_id_1",	seqNum },
									{ "PDB_ins_code_1",	iCode },
									{ "auth_atom_id_2",	atoms[1] },
									{ "auth_asym_id_2", chainID },
									{ "auth_comp_id_2", resNam },
									{ "auth_seq_id_2",	seqNum },
									{ "PDB_ins_code_2",	iCode },
									{ "auth_atom_id_3",	atoms[2] },
									{ "auth_asym_id_3", chainID },
									{ "auth_comp_id_3", resNam },
									{ "auth_seq_id_3",	seqNum },
									{ "PDB_ins_code_3",	iCode },
									{ "angle_deviation",deviation }
								});
							}
							
							break;
						
						case e_TA:
							if (not headerSeen)
							{
								if (ba::starts_with(line, "FORMAT: ") and line != "FORMAT:(10X,I3,1X,A3,1X,A1,I4,A1,4X,F7.2,3X,F7.2)")
									Error("Unexpected format in REMARK 500");
								
								headerSeen = line == "M RES CSSEQI        PSI       PHI";
							}
							else if (line.empty())
								state = e_start;
							else
							{
								int model = v_i(11, 13);
								string resNam = v_s(15, 17);
								string chainID { v_c(19) };
								int seqNum = v_i(20, 24);
								string iCode { v_c(25) };
								
								if (iCode == " ")
									iCode.clear();
								
								string psi = v_f(27, 35);
								string phi = v_f(37, 45);
								
								get_category("pdbx_validate_torsion")->emplace({
									{ "id",				to_string(++id) },
									{ "PDB_model_num",	model ? model : 1 },
									{ "auth_comp_id", 	resNam },
									{ "auth_asym_id", 	chainID },
									{ "auth_seq_id",	seqNum },
									{ "PDB_ins_code",	iCode },
									{ "phi",			phi },
									{ "psi",			psi }
								});
							}
							break;
						
						case e_CTg:
							if (not headerSeen)
								headerSeen = line == "MODEL     OMEGA";
							else if (line.empty())
								state = e_start;
							else
							{
								int model = v_i(45, 48);

								string resNam1 = v_s(12, 14);
								string chainID1 { v_c(16) };
								int seqNum1 = v_i(17, 21);
								string iCode1 { v_c(22) };
								
								if (iCode1 == " ")
									iCode1.clear();
								
								string resNam2 = v_s(27, 29);
								string chainID2 { v_c(31) };
								int seqNum2 = v_i(32, 36);
								string iCode2 { v_c(37) };
								
								if (iCode2 == " ")
									iCode2.clear();
								
								string omega = v_f(54, 60);
								
								get_category("pdbx_validate_peptide_omega")->emplace({
									{ "id",				to_string(++id) },
									{ "PDB_model_num",	model ? model : 1 },
									{ "auth_comp_id_1",	resNam1 },
									{ "auth_asym_id_1",	chainID1 },
									{ "auth_seq_id_1",	seqNum1 },
									{ "PDB_ins_code_1",	iCode1 },
									{ "auth_comp_id_2",	resNam2 },
									{ "auth_asym_id_2",	chainID2 },
									{ "auth_seq_id_2",	seqNum2 },
									{ "PDB_ins_code_2",	iCode2 },
									{ "omega",			omega }
								});
							}
							break;
						
						case e_PG:
							if (not headerSeen)
								headerSeen = line == "M RES CSSEQI        RMS     TYPE";
							else if (line.empty())
								state = e_start;
							else
							{
								int model = v_i(11, 13);
								string resNam = v_s(15, 17);
								string chainID { v_c(19) };
								int seqNum = v_i(20, 24);
								string iCode { v_c(25) };
								
								if (iCode == " ")
									iCode.clear();
								
								string rmsd = v_f(32, 36);
								string type = v_s(41);
								
								get_category("pdbx_validate_planes")->emplace({
									{ "id",				to_string(++id) },
									{ "PDB_model_num",	model ? model : 1 },
									{ "auth_comp_id",	resNam },
									{ "auth_asym_id",	chainID },
									{ "auth_seq_id",	seqNum },
									{ "PDB_ins_code",	iCode },
									{ "rmsd",			rmsd },
									{ "type",			type }
								});
							}
							break;
						
						
						default:
							state = e_start;
							break;
					}
				}
				
				break;
			}
			
			case 610:
			{
				bool headerSeen = false;
				
				for (; m_rec->is("REMARK 610"); GetNextRecord())
				{
					if (not headerSeen)
					{
						string line = v_s(12);
						headerSeen = ba::contains(line, "RES C SSEQI");
						continue;
					}
					
					int model_nr = v_i(12, 14);
					if (model_nr == 0)
						model_nr = 1;
					string res = v_s(16, 18);
					char chain = v_c(20);
					int seq = v_i(22, 25);
					char iCode = v_c(26);
					
					auto compound = libcif::compound::create(res);
					if (compound == nullptr)
						continue;
					
					vector<string> atoms;
					for (auto atom: compound->atoms())
					{
						if (atom.type_symbol != libcif::H)
							atoms.push_back(atom.id);
					}
					
					m_unobs.push_back({model_nr, res, chain, seq, iCode, { atoms }});
				}
				
				break;
			}
			
			case 800:
			{
				const regex rx1(R"(SITE_IDENTIFIER: (.+))"),
					rx2(R"(EVIDENCE_CODE: (.+))"),
					rx3(R"(SITE_DESCRIPTION: (binding site for residue ([[:alnum:]]{1,3}) ([[:alnum:]]) (\d+)|.+))", regex_constants::icase);
				
				string id, evidence, desc;
				string pdbx_auth_asym_id, pdbx_auth_comp_id, pdbx_auth_seq_id, pdbx_auth_ins_code;
				smatch m;
				
				enum State { s_start, s_id, s_evidence, s_desc, s_desc_2 } state = s_start;
				
				auto store = [&]()
				{
					// Find the matching SITE record
					auto site = FindRecord([id](PDBRecord& r) -> bool
					{
						return r.is("SITE  ") and r.v_s(12, 14) == id;
					});
					
					if (site == nullptr)
						throw runtime_error("Invalid REMARK 800, no SITE record for id " + id);
					
					// next record, store what we have	
					get_category("struct_site")->emplace({
						{ "id", id },
						{ "details", desc },
						{ "pdbx_auth_asym_id", pdbx_auth_asym_id },
						{ "pdbx_auth_comp_id", pdbx_auth_comp_id },
						{ "pdbx_auth_seq_id", pdbx_auth_seq_id },
						{ "pdbx_num_residues", site->v_i(16, 17) },
						{ "pdbx_evidence_code", evidence }
					});
				};

				for ( ; m_rec->is("REMARK 800"); GetNextRecord())
				{
					string s = m_rec->v_s(12);
					if (s.empty())
						continue;

					switch (state)
					{
						case s_start:
 							if (s == "SITE")
 								state = s_id;
							else if (VERBOSE)
								Error("Invalid REMARK 800 record, expected SITE");
 							break;
						
						case s_id:
							if (regex_match(s, m, rx1))
							{
								id = m[1].str();
								state = s_evidence;
							}
							else if (VERBOSE)
								Error("Invalid REMARK 800 record, expected SITE_IDENTIFIER");
							break;
						
						case s_evidence:
							if (regex_match(s, m, rx2))
							{
								evidence = m[1].str();
								state = s_desc;
							}
							else if (VERBOSE)
								Error("Invalid REMARK 800 record, expected SITE_IDENTIFIER");
							break;
						
						case s_desc:
							if (regex_match(s, m, rx3))
							{
								desc = m[1].str();
								pdbx_auth_comp_id = m[2].str();
								pdbx_auth_asym_id = m[3].str();
								pdbx_auth_seq_id = m[4].str();

								state = s_desc_2;
							}
							break;
						
						case s_desc_2:
							if (regex_match(s, m, rx1))
							{
								store();
								
								id = m[1].str();
								state = s_evidence;
								evidence.clear();
								desc.clear();
							}
							else
								desc = desc + ' ' + s;
							break;
					}
				}
				
				if (not id.empty())
					store();

				break;
			}
			
			case 999:
			{
				stringstream s;
				GetNextRecord();
				if (v_s(12) == "SEQUENCE")
					GetNextRecord();
				
				while (m_rec->is("REMARK 999"))
				{
					s << v_s(12) << endl;
					GetNextRecord();
				}
				
				sequence_details = s.str();
				break;
			}

			// these are skipped 

			case 2:
			case 290:
			case 300:
				GetNextRecord();
				break;

			default:
			{
				string skipped = m_rec->m_name;
				
				stringstream s;
				
				if (not m_rec->v_s(11).empty())
					s << m_rec->v_s(11) << endl;
				GetNextRecord();
					
				while (m_rec->is(skipped.c_str()))
				{
					s << m_rec->v_s(11) << endl;
					GetNextRecord();
				}
				
				get_category("pdbx_database_remark")->emplace({
					{ "id", remarkNr },
					{ "text", s.str() }
				});
				
				break;
			}
		}
	}

	if (not (compound_details.empty() and sequence_details.empty() and source_details.empty()))
	{
		get_category("pdbx_entry_details")->emplace({
			{ "entry_id",			m_structure_id },
			{ "compound_details",	compound_details },
			{ "sequence_details",	sequence_details },
			{ "source_details",		source_details }
		});
	}

	// store remark 200 info (special case)
	if (not m_remark200.empty())
		ParseRemark200();
}

void PDBFileParser::ParseRemark200()
{
	auto rm200 = [&](const char* name, int diffrnNr) -> string
	{
		int nr = 0;
		string result;
		
		for (auto i = make_split_iterator(m_remark200[name],
			ba::token_finder(ba::is_any_of(";"), ba::token_compress_off)); not i.eof(); ++i)
		{
			if (++nr != diffrnNr)
				continue;
			
			result.assign(i->begin(), i->end());;
			ba::trim(result);
			
			if (result == "NULL")
				result.clear();

			break;
		}
		
		return result;
	};

	auto inRM200 = [this](initializer_list<const char*> s) -> bool
	{
		bool result = false;
		
		for (auto* n: s)
		{
			if (not this->m_remark200[n].empty())
			{
				result = true;
				break;
			}
		}
		
		return result;
	};

/*
	The category computing is no longer used.
			
		if (inRM200({"INTENSITY-INTEGRATION SOFTWARE", "DATA SCALING SOFTWARE", "SOFTWARE USED"}) or
			not m_refinement_software.empty())
			get_category("computing")->emplace({
				{ "entry_id", m_structure_id },
				{ "pdbx_data_reduction_ii", m_remark200["INTENSITY-INTEGRATION SOFTWARE"] },
				{ "pdbx_data_reduction_ds", m_remark200["DATA SCALING SOFTWARE"] },
				{ "structure_solution", m_remark200["SOFTWARE USED"] },
				{ "structure_refinement", m_refinement_software }
			});
*/

	struct { const char* a; const char* b; }
	kSWMap[] = {
		{ "data reduction", "INTENSITY-INTEGRATION SOFTWARE" },
		{ "data scaling", "DATA SCALING SOFTWARE" },
		{ "phasing", "SOFTWARE USED" },
	};
	
	for (auto& sw: kSWMap)
	{
		if (m_remark200[sw.b].empty())
			continue;
		
		get_category("software")->emplace({
			{ "name", m_remark200[sw.b] },
			{ "classification", sw.a },
			{ "version", "." },
			{ "pdbx_ordinal", m_next_software_ord++ }
		});
	}

	string scattering_type;
	if (m_remark200["EXPERIMENT TYPE"] == "X-RAY DIFFRACTION")
		scattering_type = "x-ray";
	else if (m_remark200["EXPERIMENT TYPE"] == "NEUTRON DIFFRACTION")
		scattering_type = "neutron";
	
	set<string> diffrnWaveLengths;
	
	for (int diffrnNr = 1; ; ++diffrnNr)
	{
		string ambientTemp = rm200("TEMPERATURE (KELVIN)", diffrnNr);
		if (ambientTemp.empty())
			break;
		
		get_category("diffrn")->emplace({
			{ "id", diffrnNr },
			{ "ambient_temp", ambientTemp },
	//		{ "ambient_temp_details", seq_id },
			{ "crystal_id", 1 }
		});
		
		get_category("diffrn_detector")->emplace({
			{ "diffrn_id", diffrnNr },
			{ "detector", rm200("DETECTOR TYPE", diffrnNr) },
			{ "type", rm200("DETECTOR MANUFACTURER", diffrnNr) },
			{ "pdbx_collection_date", pdb2cif_date(rm200("DATE OF DATA COLLECTION", diffrnNr)) },
			{ "details", rm200("OPTICS", diffrnNr) }
		});
		
		if (inRM200({"MONOCHROMATIC OR LAUE (M/L)", "MONOCHROMATOR", "DIFFRACTION PROTOCOL"}) or not scattering_type.empty())
			get_category("diffrn_radiation")->emplace({
				{ "diffrn_id", diffrnNr },
				{ "wavelength_id", 1 },
				{ "pdbx_monochromatic_or_laue_m_l", rm200("MONOCHROMATIC OR LAUE (M/L)", diffrnNr) },
				{ "monochromator", rm200("MONOCHROMATOR", diffrnNr) },
				{ "pdbx_diffrn_protocol", rm200("DIFFRACTION PROTOCOL", diffrnNr) },
				{ "pdbx_scattering_type", scattering_type }
			});

		vector<string> wavelengths;
		string wl = rm200("WAVELENGTH OR RANGE (A)", diffrnNr);
		ba::split(wavelengths, wl, ba::is_any_of(", "), ba::token_compress_on);
		
		diffrnWaveLengths.insert(wavelengths.begin(), wavelengths.end());

		string source;
		if (rm200("SYNCHROTRON (Y/N)", diffrnNr) == "Y")
		{
			get_category("diffrn_source")->emplace({
				{ "diffrn_id", diffrnNr },
				{ "source", "SYNCHROTRON" },
				{ "type", rm200("RADIATION SOURCE", diffrnNr) + " BEAMLINE " + rm200("BEAMLINE", diffrnNr) },
				{ "pdbx_synchrotron_site", rm200("RADIATION SOURCE", diffrnNr) },
				{ "pdbx_synchrotron_beamline", rm200("BEAMLINE", diffrnNr) },
				
				{ "pdbx_wavelength", wavelengths.size() == 1 ? wavelengths[0] : "" },
				{ "pdbx_wavelength_list", wavelengths.size() == 1 ? "" : ba::join(wavelengths, ", ") },
			});
		}
		else if (inRM200({"X-RAY GENERATOR MODEL", "RADIATION SOURCE", "BEAMLINE", "WAVELENGTH OR RANGE (A)" }))
		{
			get_category("diffrn_source")->emplace({
				{ "diffrn_id", diffrnNr },
				{ "source", rm200("RADIATION SOURCE", diffrnNr) },
				{ "type", rm200("X-RAY GENERATOR MODEL", diffrnNr) },

				{ "pdbx_wavelength", wavelengths.size() == 1 ? wavelengths[0] : "" },
				{ "pdbx_wavelength_list", wavelengths.size() == 1 ? "" : ba::join(wavelengths, ", ") },
			});
		}
	}

	int wavelengthNr = 1;
	for (auto& wl: diffrnWaveLengths)
	{
		get_category("diffrn_radiation_wavelength")->emplace({
			{ "id", wavelengthNr++ },
			{ "wavelength", wl },
			{ "wt", "1.0" }
		});
	}

	if (inRM200({"METHOD USED TO DETERMINE THE STRUCTURE", "STARTING MODEL"}))
	{
		auto cat = get_category("refine");
		assert(cat->empty());
		
		string resolution = m_remark200["RESOLUTION RANGE HIGH (A)"];
		if (resolution.empty())
			resolution = ".";
		
		cat->emplace({
			{ "pdbx_method_to_determine_struct", m_remark200["METHOD USED TO DETERMINE THE STRUCTURE"] },
			{ "pdbx_starting_model", m_remark200["STARTING MODEL"] },
			{ "ls_d_res_high", resolution },
			{ "pdbx_diffrn_id", 1 },
			{ "pdbx_refine_id", m_exp_method },
			{ "entry_id", m_structure_id }
		});
	}
	
	if (inRM200({"REJECTION CRITERIA (SIGMA(I))", "RESOLUTION RANGE HIGH (A)", "RESOLUTION RANGE LOW (A)", "NUMBER OF UNIQUE REFLECTIONS", "COMPLETENESS FOR RANGE (%)", "<I/SIGMA(I)> FOR THE DATA SET", "R MERGE (I)", "R SYM (I)", "DATA REDUNDANCY"}))
	{
		auto cat = get_category("reflns");
		row r;
		if (cat->empty())
			cat->emplace({});
		r = cat->back();
		r["entry_id"] = m_structure_id;
		r["observed_criterion_sigma_I"] = m_remark200["REJECTION CRITERIA (SIGMA(I))"];
		r["d_resolution_high"] = m_remark200["RESOLUTION RANGE HIGH (A)"];
		r["d_resolution_low"] = m_remark200["RESOLUTION RANGE LOW (A)"];
		r["number_obs"] = m_remark200["NUMBER OF UNIQUE REFLECTIONS"];
		r["percent_possible_obs"] = m_remark200["COMPLETENESS FOR RANGE (%)"];
		r["pdbx_netI_over_sigmaI"] = m_remark200["<I/SIGMA(I)> FOR THE DATA SET"];
		r["pdbx_Rmerge_I_obs"] = m_remark200["R MERGE (I)"];
		r["pdbx_Rsym_value"] = m_remark200["R SYM (I)"];
		r["pdbx_redundancy"] = m_remark200["DATA REDUNDANCY"];
		r["pdbx_ordinal"] = 1;
		r["pdbx_diffrn_id"] = 1;
	}
	
	if (inRM200({ "HIGHEST RESOLUTION SHELL, RANGE HIGH (A)", "HIGHEST RESOLUTION SHELL, RANGE LOW (A)", "COMPLETENESS FOR SHELL (%)",
		"R MERGE FOR SHELL (I)", "R SYM FOR SHELL (I)", "<I/SIGMA(I)> FOR SHELL", "DATA REDUNDANCY IN SHELL" }))
	{
		get_category("reflns_shell")->emplace({
			{ "d_res_high", m_remark200["HIGHEST RESOLUTION SHELL, RANGE HIGH (A)"] },
			{ "d_res_low", m_remark200["HIGHEST RESOLUTION SHELL, RANGE LOW (A)"] },
			{ "percent_possible_all", m_remark200["COMPLETENESS FOR SHELL (%)"] },
			{ "Rmerge_I_obs", m_remark200["R MERGE FOR SHELL (I)"] },
			{ "pdbx_Rsym_value", m_remark200["R SYM FOR SHELL (I)"] },
			{ "meanI_over_sigI_obs", m_remark200["<I/SIGMA(I)> FOR SHELL"] },
			{ "pdbx_redundancy", m_remark200["DATA REDUNDANCY IN SHELL"] },
			{ "pdbx_ordinal", 1},
			{ "pdbx_diffrn_id" , 1}
		});
	}
	
}

void PDBFileParser::ParseRemark350()
{
	auto saved = m_rec;
	
	enum State { eStart, eInfo, eAnd, eApply, eBioMT } state = eStart;
	
	const regex
		kRX1(R"(BIOMOLECULE: (\d+))"),
		kRX2(R"(([^:]+): (.+?)(?: (ANGSTROM\*\*2|KCAL/MOL))?)"),
		kRX8(R"(APPLY THE FOLLOWING TO CHAINS: (.+))"),
		kRX9(R"(AND CHAINS: (.+))"),
		kRX10(R"(BIOMT([123])\s+(\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?))");
	
	int biomolecule = 0, oper_id = 0;
	vector<string> oper_expression;
	map<string,string> values;
	vector<string> asym_id_list;
	smatch m;
	row gen_r;
	
	vector<double> mat, vec;
	
	for (m_rec = FindRecord("REMARK 350"); m_rec != nullptr and m_rec->is("REMARK 350"); GetNextRecord())
	{
		string line = v_s(11);
		
		switch (state)
		{
			case eStart:
				if (regex_match(line, m, kRX1))
				{
					biomolecule = stoi(m[1].str());
					state = eInfo;
				}
				break;
			
			case eInfo:
				if (regex_match(line, m, kRX8))
				{
					state = eApply;

					string value = m[1].str();

					for (auto i = make_split_iterator(value,
						ba::token_finder(ba::is_any_of(", "), ba::token_compress_on)); not i.eof(); ++i)
					{
						string chain = boost::copy_range<string>(*i);
						
						if (chain.empty())	// happens when we have a AND CHAIN line
						{
							state = eAnd;
							break;
						}
						
						if (chain.length() != 1)
							throw runtime_error("Invalid REMARK 350");
						
						MapChainID2AsymIDS(chain[0], asym_id_list);
					}
				}
				else if (regex_match(line, m, kRX2))
					values[m[1].str()] = m[2].str();
				break;
			
			case eAnd:
				if (regex_match(line, m, kRX9))
				{
					state = eApply;
					
					string value = m[1].str();
					for (auto i = make_split_iterator(value,
						ba::token_finder(ba::is_any_of(", "), ba::token_compress_on)); not i.eof(); ++i)
					{
						string chain = boost::copy_range<string>(*i);
						
						if (chain.empty())	// happens when we have another AND CHAIN line
						{
							state = eAnd;
							break;
						}
						
						MapChainID2AsymIDS(chain[0], asym_id_list);
					}
					
					continue;
				}
				// fall through
		
			case eApply:
				if (regex_match(line, m, kRX10))
				{
					int mt = stoi(m[1].str());
					if (mt != 1)
						throw runtime_error("Invalid REMARK 350");

					oper_id = stoi(m[2].str());
					oper_expression.push_back(to_string(oper_id));

					mat.push_back(stod(m[3].str())); 
					mat.push_back(stod(m[4].str())); 
					mat.push_back(stod(m[5].str()));
					vec.push_back(stod(m[6].str()));
					state = eBioMT;
				}
				break;
			
			case eBioMT:
				if (regex_match(line, m, kRX10))
				{
					int mt = stoi(m[1].str());
					
					if (mt == 1)
					{
						oper_id = stoi(m[2].str());
						oper_expression.push_back(to_string(oper_id));
					}
					else if (oper_id != stoi(m[2].str()))
						throw runtime_error("Invalid REMARK 350");
					
					mat.push_back(stod(m[3].str())); 
					mat.push_back(stod(m[4].str())); 
					mat.push_back(stod(m[5].str()));
					vec.push_back(stod(m[6].str()));
					
					if (mt == 3)
					{
						if (vec.size() != 3 or mat.size() != 9)
							throw runtime_error("Invalid REMARK 350");
						
						if (oper_id == 1)
						{
							string oligomer = values["AUTHOR DETERMINED BIOLOGICAL UNIT"];
							if (oligomer.empty())
								oligomer = values["SOFTWARE DETERMINED QUATERNARY STRUCTURE"];
							ba::to_lower(oligomer);
							
							int count = 0;
							smatch m;
							
							if (regex_match(oligomer, m, regex(R"((\d+)-meric)")))
							{
								count = stoi(m[1].str());
							}
							else if (ba::ends_with(oligomer, "meric"))
							{
								string cs = oligomer.substr(0, oligomer.length() - 5);
								if (cs == "mono")			count = 1;
								else if (cs == "di")		count = 2;
								else if (cs == "tri")		count = 3;
								else if (cs == "tetra")		count = 4;
								else if (cs == "hexa")		count = 6;
								else if (cs == "octa")		count = 8;
								else if (cs == "dodeca")	count = 12;
							}
							
							string details;
							if (values["AUTHOR DETERMINED BIOLOGICAL UNIT"].empty())
							{
								if (not values["SOFTWARE DETERMINED QUATERNARY STRUCTURE"].empty())
									details = "software_defined_assembly";
							}
							else if (values["SOFTWARE DETERMINED QUATERNARY STRUCTURE"].empty())
								details = "author_defined_assembly";
							else
								details = "author_and_software_defined_assembly";
							
							get_category("pdbx_struct_assembly")->emplace({
								{ "id",		biomolecule },
								{ "details", details },
								{ "method_details", values["SOFTWARE USED"] },
								{ "oligomeric_details", oligomer },
								{ "oligomeric_count", count > 0 ? to_string(count) : "" }
							});
							
							auto cat = get_category("pdbx_struct_assembly_prop");
							
							if (not values["TOTAL BURIED SURFACE AREA"].empty())
								cat->emplace({
									{ "biol_id",	biomolecule },
									{ "type",		"ABSA (A^2)" },
									{ "value",		values["TOTAL BURIED SURFACE AREA"] }
								});
							
							if (not values["CHANGE IN SOLVENT FREE ENERGY"].empty())
								cat->emplace({
									{ "biol_id",	biomolecule },
									{ "type",		"MORE" },
									{ "value",		values["CHANGE IN SOLVENT FREE ENERGY"] }
								});
							
							if (not values["SURFACE AREA OF THE COMPLEX"].empty())
								cat->emplace({
									{ "biol_id",	biomolecule },
									{ "type",		"SSA (A^2)" },
									{ "value",		values["SURFACE AREA OF THE COMPLEX"] }
								});
							
							values.clear();
						}

						boost::format fmt("%12.10f");
						
						get_category("pdbx_struct_oper_list")->emplace({
							{ "id", oper_id },
							{ "type",
								mat == vector<double>{ 1, 0, 0, 0, 1, 0, 0, 0, 1 } and vec == vector<double>{ 0, 0, 0 }
									? "identity operation" : "crystal symmetry operation" },
//										{ "name", "" }, 
//										{ "symmetry_operation", "" },
							{ "matrix[1][1]", (fmt % mat[0]).str() },
							{ "matrix[1][2]", (fmt % mat[1]).str() },
							{ "matrix[1][3]", (fmt % mat[2]).str() }, 
							{ "vector[1]", (fmt % vec[0]).str() },
							{ "matrix[2][1]", (fmt % mat[3]).str() },
							{ "matrix[2][2]", (fmt % mat[4]).str() },
							{ "matrix[2][3]", (fmt % mat[5]).str() },
							{ "vector[2]", (fmt % vec[1]).str() },
							{ "matrix[3][1]", (fmt % mat[6]).str() },
							{ "matrix[3][2]", (fmt % mat[7]).str() },
							{ "matrix[3][3]", (fmt % mat[8]).str() },
							{ "vector[3]", (fmt % vec[2]).str() }
						});

						mat.clear();
						vec.clear();
					}
				}
				else if (regex_match(line, m, kRX1))
				{
					if (not (vec.empty() and mat.empty()))
						throw runtime_error("Invalid REMARK 350");

					get_category("pdbx_struct_assembly_gen")->emplace({
						{ "assembly_id", biomolecule },
						{ "oper_expression", ba::join(oper_expression, ",") },
						{ "asym_id_list", ba::join(asym_id_list, ",") }
					});

					biomolecule = stoi(m[1].str());
					asym_id_list.clear();
					oper_expression.clear();
					
					state = eInfo;
				}
				break;
				
		}
	}
	
	if (not oper_expression.empty())
	{
		get_category("pdbx_struct_assembly_gen")->emplace({
			{ "assembly_id", biomolecule },
			{ "oper_expression", ba::join(oper_expression, ",") },
			{ "asym_id_list", ba::join(asym_id_list, ",") }
		});
	}
	
	m_rec = saved;
}

void PDBFileParser::ParsePrimaryStructure()
{
	// First locate the DBREF record. Might be missing
	DBREF cur = { m_structure_id };

	while (ba::starts_with(m_rec->m_name, "DBREF"))
	{
		if (m_rec->is("DBREF "))						//	 1 -  6       Record name   "DBREF "                                                    
		{
			cur.PDBIDCode				= v_s(8, 11);	//	 8 - 11       IDcode        idCode             ID code of this datablock.                   
			cur.chainID					= v_c(13);      //	13            Character     chainID            Chain  identifier.                       
			cur.seqBegin				= v_i(15, 18);  //	15 - 18       Integer       seqBegin           Initial sequence number of the           
                                                        //	                                               PDB sequence segment.                    
			cur.insertBegin				= v_c(19);      //	19            AChar         insertBegin        Initial  insertion code of the           
                                                        //	                                               PDB  sequence segment.                   
			cur.seqEnd					= v_i(21, 24);  //	21 - 24       Integer       seqEnd             Ending sequence number of the            
                                                        //	                                               PDB  sequence segment.                   
			cur.insertEnd				= v_c(25);      //	25            AChar         insertEnd          Ending insertion code of the             
                                                        //	                                               PDB  sequence segment.                   
			cur.database				= v_s(27, 32);  //	27 - 32       LString       database           Sequence database name.                  
			cur.dbAccession				= v_s(34, 41);  //	34 - 41       LString       dbAccession        Sequence database accession code.        
			cur.dbIdCode				= v_s(43, 54);  //	43 - 54       LString       dbIdCode           Sequence  database identification code.  
			cur.dbSeqBegin				= v_i(56, 60);  //	56 - 60       Integer       dbseqBegin         Initial sequence number of the           
                                                        //	                                               database seqment.                        
			cur.dbinsBeg				= v_c(61);      //	61            AChar         idbnsBeg           Insertion code of initial residue of the 
                                                        //	                                               segment, if PDB is the reference.        
			cur.dbSeqEnd				= v_i(63, 67);  //	63 - 67       Integer       dbseqEnd           Ending sequence number of the            
                                                        //	                                               database segment.                        
			cur.dbinsEnd				= v_c(68);      //	68            AChar         dbinsEnd           Insertion code of the ending residue of  
			                                            //	                                               the segment, if PDB is the reference.    
			auto& chain = GetChainForID(cur.chainID);
			chain.m_dbref = cur;
		}
		else if (m_rec->is("DBREF1"))					//	 1 -  6        Record name   "DBREF1"                                             
		{
			cur.PDBIDCode				= v_s(8, 11);	//	 8 - 11       IDcode        idCode             ID code of this datablock.                   
			cur.chainID					= v_c(13);      //	13             Character     chainID       Chain identifier.                      
			cur.seqBegin				= v_i(15, 18);  //	15 - 18        Integer       seqBegin      Initial sequence number of the         
                                                        //	                                           PDB sequence segment, right justified. 
			cur.insertBegin				= v_c(19);      //	19             AChar         insertBegin   Initial insertion code of the          
                                                        //	                                           PDB sequence segment.                  
			cur.seqEnd					= v_i(21, 24);  //	21 - 24        Integer       seqEnd        Ending sequence number of the          
                                                        //	                                           PDB sequence segment, right justified. 
			cur.insertEnd				= v_c(25);      //	25             AChar         insertEnd     Ending insertion code of the           
                                                        //	                                           PDB sequence  segment.                 
			cur.database				= v_s(27, 32);  //	27 - 32        LString       database      Sequence database name.                
			cur.dbIdCode				= v_s(48, 67);  //	48 - 67        LString       dbIdCode      Sequence database identification code, 
		}
		else if (m_rec->is("DBREF2"))					//	 1 -  6       Record name   "DBREF2"                                        
		{                                               //	 8 - 11       IDcode        idCode        ID code of this datablock.            
			if (v_c(13) != cur.chainID)			        //	13            Character     chainID       Chain identifier.                 
				Error("Chain ID's for DBREF1/DBREF2 records do not match");
			cur.dbAccession				= v_s(19, 40);  //	19 - 40       LString       dbAccession   Sequence database accession code, 
                                                        //	                                          left justified.                   
			cur.dbSeqBegin				= v_i(46, 55);  //	46 - 55       Integer       seqBegin      Initial sequence number of the    
                                                        //	                                          Database segment, right justified.
			cur.dbSeqEnd				= v_i(58, 67);  //	58 - 67       Integer       seqEnd        Ending sequence number of the     
			                                            //	                                          Database segment, right justified.
			auto& chain = GetChainForID(cur.chainID);
			chain.m_dbref = cur;
		}
		
		GetNextRecord();
	}

	// update chains
	for (auto& chain: m_chains)
	{
		chain.m_next_seq_num = chain.m_dbref.seqBegin;
		chain.m_next_db_seq_num = chain.m_dbref.dbSeqBegin;
	}

	while (m_rec->is("SEQADV"))
	{							//	 1 -  6        Record name   "SEQADV"                                           
		m_seqadvs.push_back({	//	 8 - 11        IDcode        idCode        ID  code of this datablock.              	
			v_s(13, 15),        //	13 - 15        Residue name  resName       Name of the PDB residue in conflict. 
			v_c(17),            //	17             Character     chainID       PDB  chain identifier.               
			v_i(19, 22),        //	19 - 22        Integer       seqNum        PDB  sequence number.                
			v_c(23),            //	23             AChar         iCode         PDB insertion code.                  
			v_s(25, 28),        //	25 - 28        LString       database                                           
			v_s(30, 38),        //	30 - 38        LString       dbAccession   Sequence  database accession number. 
			v_s(40, 42),        //	40 - 42        Residue name  dbRes         Sequence database residue name.      
			v_i(44, 48),        //	44 - 48        Integer       dbSeq         Sequence database sequence number.   
			v_s(50, 70)         //	50 - 70        LString       conflict      Conflict comment.                                            
		});

		GetNextRecord();
	}
	
	while (m_rec->is("SEQRES"))			//	 1 -  6        Record name    "SEQRES"
	{									//	 8 - 10        Integer        serNum       Serial number of the SEQRES record for  the
										//	                                           current  chain. Starts at 1 and increments
										//	                                           by one  each line. Reset to 1 for each chain.
		char chain_id = v_c(12);		//	12             Character      chainID      Chain identifier. This may be any single
										//	                                           legal  character, including a blank which is
										//	                                           is  used if there is only one chain.
		int numRes = v_i(14, 17);		//	14 - 17        Integer        numRes       Number of residues in the chain.
										//	                                           This  value is repeated on every record.
		string monomers = v_s(20, 70);	//	20 - 22        Residue name   resName      Residue name.
										//	 ...

		auto& chain = GetChainForID(chain_id, numRes);
		
		for (auto si = ba::make_split_iterator(monomers, ba::token_finder(ba::is_any_of(" "), ba::token_compress_on)); not si.eof(); ++si)
		{
			string mon_id(si->begin(), si->end());
			if (mon_id.empty())
				continue;
			
			chain.m_seqres.push_back({mon_id, chain.m_next_seq_num++, ' ', chain.m_next_db_seq_num++});
			
			InsertChemComp(mon_id);
		}
		
		GetNextRecord();
	}

	// First pass over MODRES, only store relevant information required in ConstructEntities
	while (m_rec->is("MODRES"))				//	 1 -  6        Record name   "MODRES"                                            												
	{							 			//	 8 - 11        IDcode        idCode      ID code of this datablock.                  
		string resName		= v_s(13, 15);	//	13 - 15        Residue name  resName     Residue name used in this datablock.        
//		char chainID		= v_c(17);		//	17             Character     chainID     Chain identifier.                       
//		int seqNum			= v_i(19, 22);	//	19 - 22        Integer       seqNum      Sequence number.                        
//		char iCode			= v_c(23);		//	23             AChar         iCode       Insertion code.                         
		string stdRes		= v_s(25, 27);	//	25 - 27        Residue name  stdRes      Standard residue name.                  
//		string comment		= v_s(30, 70);	//	30 - 70        String        comment     Description of the residue modification.

		m_mod2parent[resName] = stdRes;

		GetNextRecord();
	}
}

void PDBFileParser::ParseHeterogen()
{
	while (m_rec->is("HET   "))
	{									//	 1 -  6       Record name   "HET   "                                                         
		string hetID = v_s(8, 10);      //	 8 - 10       LString(3)    hetID          Het identifier, right-justified.                  
		char chainID = v_c(13);			//	13            Character     ChainID        Chain  identifier.                                
		int seqNum = v_i(14, 17);		//	14 - 17       Integer       seqNum         Sequence  number.                                 
		char iCode = v_c(18);			//	18            AChar         iCode          Insertion  code.                                  
		int numHetAtoms = v_i(21, 25);	//	21 - 25       Integer       numHetAtoms    Number of HETATM records for the group            
										//	                                           present in the datablock.                             
		string text = v_s(31, 70);		//	31 - 70       String        text           Text describing Het group.                        

		m_hets.push_back({ hetID, chainID, seqNum, iCode, numHetAtoms, text });
		
		GetNextRecord();
	}
	
	while (m_rec->is("HETNAM"))		//	 1 -  6       Record name   "HETNAM"                                                 
	{								//	 9 - 10       Continuation  continuation    Allows concatenation of multiple records.
		string hetID = v_s(12, 14);	//	12 - 14       LString(3)    hetID           Het identifier, right-justified.         
        string text = v_s(16);		//	16 - 70       String        text            Chemical name.                           

		m_hetnams[hetID] = text;
		InsertChemComp(hetID);
		
		GetNextRecord();
	}

	while (m_rec->is("HETSYN"))		 //	 1 -  6       Record name   "HETSYN"                                                           
	{                                //	 9 - 10       Continuation  continuation   Allows concatenation of multiple records.           
         string hetID = v_s(12, 14); //	12 - 14       LString(3)    hetID          Het identifier, right-justified.                    
         string syn = v_s(16);		 //	16 - 70       SList         hetSynonyms    List of synonyms.                                   

		m_hetsyns[hetID] = syn;

		GetNextRecord();
	}
	
	while (m_rec->is("FORMUL"))			//	 1 -  6        Record name   "FORMUL"                          
	{                                   //	 9 - 10        Integer       compNum       Component  number.  
        string hetID = v_s(13, 15);     //	13 - 15        LString(3)    hetID         Het identifier.     
                                        //	17 - 18        Integer       continuation  Continuation number.
		char waterMark = v_c(19);		//	19             Character     asterisk      "*" for water.      
		string formula = v_s(20);		//	20 - 70        String        text          Chemical formula.   
		
		m_formuls[hetID] = formula;
		
		if (waterMark == '*')
			m_water_het_id = hetID;
		
		GetNextRecord();
	}
}

void PDBFileParser::ConstructEntities()
{
	// We parsed the Primary Structure and Heterogen sections, if available.
	// But if we didn't parse anything, we need to fake the data based on residues in ATOM records

	// First iterate all ATOM records and store the residues as found in these records
	int model_nr = 1;
	
	for (auto r = m_data; r != nullptr; r = r->m_next)
	{
		if (r->is("MODEL "))
		{
			model_nr = v_i(11, 14);
			continue;
		}

		if (r->is("ATOM  ") or r->is("HETATM"))		//	 1 -  6        Record name   "ATOM  "
		{											//	 ...
			string name			= r->v_s(13, 16);	//	13 - 16        Atom          name         Atom name.
			string resName		= r->v_s(18, 20);	//	18 - 20        Residue name  resName      Residue name.
			char chainID		= r->v_c(22);		//	22             Character     chainID      Chain identifier.
			int resSeq			= r->v_i(23, 26);	//	23 - 26        Integer       resSeq       Residue sequence number.
			char iCode			= r->v_c(27);		//	27             AChar         iCode        Code for insertion of residues.

			auto& chain = GetChainForID(chainID);
			
			PDBChain::AtomRes ar{ resName, resSeq, iCode };

			if (chain.m_residues_seen.empty() or chain.m_residues_seen.back() != ar)
				chain.m_residues_seen.push_back(ar);

			// now that we're iterating atoms anyway, clean up the m_unobs array
			m_unobs.erase(remove_if(m_unobs.begin(), m_unobs.end(), [=](UNOBS& a)
			{
				bool result = false;
				
				if (model_nr == a.model_nr and
					resName == a.res and
					chainID == a.chain and
					resSeq == a.seq and
					iCode == a.iCode)
				{
					auto i = find(a.atoms.begin(), a.atoms.end(), name);
					if (i != a.atoms.end())
					{
						a.atoms.erase(i);
						result = a.atoms.empty();
					}
				}
				
				return result;
			}), m_unobs.end());

			continue;
		}

		if (r->is("TER   "))						//	 1 -  6 	   Record name	 "TER	"								  
		{											//	 7 - 11 	   Integer		 serial 		 Serial number. 		  
													//	18 - 20 	   Residue name  resName		 Residue name.			  
			char chainID = r->v_c(22);				//	22			   Character	 chainID		 Chain identifier.		  
													//	23 - 26 	   Integer		 resSeq 		 Residue sequence number. 
													//	27			   AChar		 iCode			 Insertion code.		  
			auto& chain = GetChainForID(chainID);
			chain.m_ter_index = chain.m_residues_seen.size();
			continue;
		}
	}
	
	for (auto& chain: m_chains)
	{
		if (not (chain.m_seqres.empty() or chain.m_residues_seen.empty()))
		{
			// seems safe to assume TER record is at the right location...
			// However, some files don't have them at all.
			// When m_ter_index == 0 this is most likely the case. Right?
			
			if (chain.m_ter_index > 0)
				chain.m_residues_seen.erase(chain.m_residues_seen.begin() + chain.m_ter_index, chain.m_residues_seen.end());
			
			chain.AlignResToSeqRes();
		}
		else
		{
			// So, we did not have a SEQRES for this chain. Try to reconstruct it.
			// Problem here is that TER records may be located incorrectly. So
			// first lets shift the ter index until it is past the last known
			// aminoacid or base.
			
			for (size_t ix = chain.m_ter_index; ix < chain.m_residues_seen.size(); ++ix)
			{
				string resName = chain.m_residues_seen[ix].m_mon_id;
				
				if (kAAMap.count(resName) or
					kBaseMap.count(resName) or
					PeptideDB::Instance().IsKnownPeptide(resName) or
					PeptideDB::Instance().IsKnownBase(resName))
				{
					chain.m_ter_index = ix + 1;
				}

				InsertChemComp(resName);
			}
			
			// And now construct our 'SEQRES'...
			for (size_t ix = 0; ix < chain.m_ter_index; ++ix)
			{
				auto& ar = chain.m_residues_seen[ix];
				chain.m_seqres.push_back({ar.m_mon_id, ar.m_seq_num, ar.m_icode, ar.m_seq_num, true});
			}
		}
	}

	for (auto r = m_data; r != nullptr; r = r->m_next)
	{
		if (r->is("ATOM  ") or r->is("HETATM"))
		{											//	 1 -  6        Record name   "ATOM  "
			int serial = r->v_i(7, 11);				//	 7 - 11        Integer       serial       Atom  serial number.
													//	 ...
			char altLoc = v_c(17);					//	17             Character     altLoc       Alternate location indicator.
			string resName		= r->v_s(18, 20);	//	18 - 20        Residue name  resName      Residue name.
			char chainID		= r->v_c(22);		//	22             Character     chainID      Chain identifier.
			int resSeq			= r->v_i(23, 26);	//	23 - 26        Integer       resSeq       Residue sequence number.
			char iCode			= r->v_c(27);		//	27             AChar         iCode        Code for insertion of residues.

			auto& chain = GetChainForID(chainID);

			auto i = find(chain.m_seqres.begin(), chain.m_seqres.end(), PDBSeqRes{resName, resSeq, iCode});

			// might be a hetero
			if (altLoc != ' ' and i == chain.m_seqres.end())
			{
				i = find_if(chain.m_seqres.begin(), chain.m_seqres.end(),
					[resSeq, iCode](const PDBSeqRes& r) -> bool
					{
						return r.m_seq_num == resSeq and r.m_icode == iCode;
					});
			}

			if (i != chain.m_seqres.end())
			{
				i->m_seen = true;
				if (i->m_mon_id != resName)
					i->m_alts.insert(resName);
			}
			else
			{
				auto& residues = chain.m_het;
	
				if (residues.empty() or residues.back().m_seq_num != resSeq)
				{
					i = lower_bound(residues.begin(), residues.end(),
						PDBSeqRes{resName, resSeq, iCode},
						[=](const PDBSeqRes& r1, const PDBSeqRes& r2) -> bool {
							return r1.m_seq_num < r2.m_seq_num;
						});
	
					residues.insert(i, { resName, resSeq, iCode, resSeq, true });

					InsertChemComp(resName);
				}
			}
			
			if (r->is("HETATM"))
			{
				if (is_water(resName))
					m_water_het_id = resName;
				
				auto h = find_if(m_hets.begin(), m_hets.end(), [=](const HET& het) -> bool
					{
						return het.hetID == resName and het.chainID == chainID and
							het.seqNum == resSeq and het.iCode == iCode;
					});
				
				if (h == m_hets.end())
				{
					m_hets.push_back({ resName, chainID, resSeq, iCode, 0 });	// double perhaps, but that does not care
					h = prev(m_hets.end());
				}

				h->atoms.insert(serial);
			}

			continue;
		}
	}

	// Create missing compounds
	for (auto& chain: m_chains)
	{
		if (chain.m_mol_id != 0 or chain.m_seqres.empty())
			continue;

		// now this chain may contain the same residues as another one
		for (auto& other: m_chains)
		{
			if (&other == &chain or other.m_mol_id == 0)
				continue;
			
			if (chain.SameSequence(other))
			{
				chain.m_mol_id = other.m_mol_id;
				break;
			}
		}			
			
		if (chain.m_mol_id != 0)
			continue;
		
		auto& comp = GetOrCreateCompound(m_next_mol_id++);
		comp.m_chains.insert(chain.m_dbref.chainID);

		chain.m_mol_id = comp.m_mol_id;
	}

	set<string> struct_title, struct_description;
	
	// Create poly_scheme and write pdbx_poly_seq_scheme and create mapping table

	auto cat = get_category("pdbx_poly_seq_scheme");
	int asym_nr = 0;
	for (auto& chain: m_chains)
	{
		string asym_id = cif_id_for_int(asym_nr++);
		string entity_id = m_molID2EntityID[chain.m_mol_id];
		
		m_asymID2EntityID[asym_id] = entity_id;
		
		get_category("struct_asym")->emplace({
			{ "id", asym_id },
			{ "pdbx_blank_PDB_chainid_flag", chain.m_dbref.chainID == ' ' ? "Y" : "N" },
//			pdbx_modified 
			{ "entity_id", entity_id },
//			details
		});
		
		int seq_nr = 1;
		for (auto& res: chain.m_seqres)
		{
			string auth_mon_id, auth_seq_num;
			if (res.m_seen)
			{
				auth_mon_id = res.m_mon_id;
				auth_seq_num = to_string(res.m_seq_num);
			}

			m_chainSeq2AsymSeq[make_tuple(chain.m_dbref.chainID, res.m_seq_num, res.m_icode)] = make_tuple(asym_id, seq_nr, true);
			
			string seq_id = to_string(seq_nr);
			++seq_nr;
			
			cat->emplace({
				{ "asym_id", asym_id },
				{ "entity_id", m_molID2EntityID[chain.m_mol_id] },
				{ "seq_id", seq_id },
				{ "mon_id", res.m_mon_id },
				{ "ndb_seq_num", seq_id },
				{ "pdb_seq_num", res.m_seq_num },
				{ "auth_seq_num", auth_seq_num },
				{ "pdb_mon_id", auth_mon_id },
				{ "auth_mon_id", auth_mon_id },
				{ "pdb_strand_id", string{chain.m_dbref.chainID} },
				{ "pdb_ins_code", (res.m_icode == ' ' or res.m_icode == 0) ? "." : string{res.m_icode} },
				{ "hetero", res.m_alts.empty() ? "n" : "y" }
			});
		}
	}
	
	// We have now created all compounds, write them out
	uint32 struct_ref_id = 0, struct_ref_seq_align_id = 0;
	
	for (auto& cmp: m_compounds)
	{
		++struct_ref_id;

		string src_method;
		
		if (not cmp.m_source["SYNTHETIC"].empty())
		{
			src_method = "syn";
			
			get_category("pdbx_entity_src_syn")->emplace({
				{ "entity_id", m_molID2EntityID[cmp.m_mol_id] },
				{ "pdbx_src_id", struct_ref_id },
				{ "organism_scientific", cmp.m_source["ORGANISM_SCIENTIFIC"] },
				{ "ncbi_taxonomy_id", cmp.m_source["ORGANISM_TAXID"] },
			});
		}
		else if (cmp.m_info["ENGINEERED"] == "YES" or
			not cmp.m_source["EXPRESSION_SYSTEM"].empty())
		{
			src_method = "man";
			
			get_category("entity_src_gen")->emplace({
				{ "entity_id", m_molID2EntityID[cmp.m_mol_id] },
				{ "pdbx_src_id", struct_ref_id },
				{ "gene_src_common_name", cmp.m_source["ORGANISM_COMMON"] },
				{ "pdbx_gene_src_gene", cmp.m_source["GENE"] },
				{ "gene_src_strain", cmp.m_source["STRAIN"] },
				{ "pdbx_gene_src_cell_line", cmp.m_source["CELL_LINE"] },
				{ "pdbx_gene_src_organelle", cmp.m_source["ORGANELLE"] },
				{ "pdbx_gene_src_cellular_location", cmp.m_source["CELLULAR_LOCATION"] },
				{ "pdbx_gene_src_scientific_name", cmp.m_source["ORGANISM_SCIENTIFIC"] },
				{ "pdbx_gene_src_ncbi_taxonomy_id", cmp.m_source["ORGANISM_TAXID"] },
				{ "pdbx_host_org_scientific_name", cmp.m_source["EXPRESSION_SYSTEM"] },
				{ "pdbx_host_org_ncbi_taxonomy_id", cmp.m_source["EXPRESSION_SYSTEM_TAXID"] },
				{ "pdbx_host_org_strain", cmp.m_source["EXPRESSION_SYSTEM_STRAIN"] },
				{ "pdbx_host_org_variant", cmp.m_source["EXPRESSION_SYSTEM_VARIANT"] },
				{ "pdbx_host_org_cellular_location", cmp.m_source["EXPRESSION_SYSTEM_CELLULAR_LOCATION"] },
				{ "pdbx_host_org_vector_type", cmp.m_source["EXPRESSION_SYSTEM_VECTOR_TYPE"] },
				{ "pdbx_host_org_vector", cmp.m_source["EXPRESSION_SYSTEM_VECTOR"] },
				{ "pdbx_host_org_gene", cmp.m_source["EXPRESSION_SYSTEM_GENE"] },
				{ "plasmid_name", cmp.m_source["EXPRESSION_SYSTEM_PLASMID"] },
				{ "pdbx_description", cmp.m_source["OTHER_DETAILS"] }
			});
		}
		else if (not cmp.m_source["ORGANISM_SCIENTIFIC"].empty())
		{
			src_method = "nat";
			
			get_category("entity_src_nat")->emplace({
				{ "entity_id", m_molID2EntityID[cmp.m_mol_id] },
				{ "pdbx_src_id", struct_ref_id },
				{ "common_name", cmp.m_source["ORGANISM_COMMON"] },
				{ "strain", cmp.m_source["STRAIN"] },
				{ "pdbx_secretion", cmp.m_source["SECRETION"] },
				{ "pdbx_organism_scientific", cmp.m_source["ORGANISM_SCIENTIFIC"] },
				{ "pdbx_ncbi_taxonomy_id", cmp.m_source["ORGANISM_TAXID"] },
				{ "pdbx_cellular_location", cmp.m_source["CELLULAR_LOCATION"] },
				{ "pdbx_plasmid_name", cmp.m_source["PLASMID"] },
				{ "pdbx_organ", cmp.m_source["ORGAN"] },
			});
		}
		
		get_category("entity")->emplace({
			{ "id", m_molID2EntityID[cmp.m_mol_id] },
			{ "type", "polymer" },
			{ "src_method", src_method },
			{ "pdbx_description", cmp.m_info["MOLECULE"] },
//			{ "pdbx_formula_weight", 		},
			{ "pdbx_number_of_molecules", cmp.m_chains.size() },
			{ "details", cmp.m_info["OTHER_DETAILS"] },
			{ "pdbx_mutation", cmp.m_info["MUTATION"] },
			{ "pdbx_fragment", cmp.m_info["FRAGMENT"] },
			{ "pdbx_ec", cmp.m_info["EC"] }
		});
		
		if (not cmp.m_info["SYNONYM"].empty())
		{
			get_category("entity_name_com")->emplace({
				{ "entity_id", m_molID2EntityID[cmp.m_mol_id] },
				{ "name", cmp.m_info["SYNONYM"] }
			});
		}

		string desc = cmp.m_info["MOLECULE"];
		if (not cmp.m_info["EC"].empty())
			desc += " (E.C." + cmp.m_info["EC"] + ")";

		if (not cmp.m_title.empty())
			struct_title.insert(cmp.m_title);
		
		if (not desc.empty())
			struct_description.insert(desc);
		
		auto ci = find_if(m_chains.begin(), m_chains.end(),
			[cmp](PDBChain& c) -> bool { return cmp.m_chains.count(c.m_dbref.chainID); } );

		if (ci != m_chains.end() and not ci->m_dbref.dbIdCode.empty())
		{
			get_category("struct_ref")->emplace({
				{ "id", struct_ref_id },
				{ "entity_id", m_molID2EntityID[cmp.m_mol_id] },
				{ "db_name", ci->m_dbref.database },
				{ "db_code", ci->m_dbref.dbIdCode },
				{ "pdbx_db_accession", ci->m_dbref.dbAccession },
//				{ "pdbx_align_begin", ci->m_dbref.dbSeqBegin }
			});
		}
		
		bool nstd_monomer = false, nonstandard_linkage = false;
		bool mightBePolyPeptide = true, mightBeDNA = true;
		
		vector<string> chains;
		string seq, seq_can;
		
		// write out the chains for this compound
		for (auto& chain: m_chains)
		{
			if (chain.m_mol_id != cmp.m_mol_id)
				continue;
			
//			chain.m_entity_id = cmp.m_entity_id;

			++struct_ref_seq_align_id;
			DBREF& dbref = chain.m_dbref;
			
			if (not dbref.database.empty())
			{
				auto ins_to_str = [](char i) -> string { return i == ' ' ? "" : string{ i }; };
				
				auto& pdbx_poly_seq_scheme = *get_category("pdbx_poly_seq_scheme");
				
				int seq_align_beg = 0, seq_align_end = 0;
				
				try
				{
					seq_align_beg = pdbx_poly_seq_scheme[
							key("pdb_strand_id") == dbref.chainID and
							key("pdb_seq_num") == dbref.seqBegin and
							key("pdb_ins_code") == ins_to_str(dbref.insertBegin)]
						["seq_id"].as<int>();
	
					seq_align_end = pdbx_poly_seq_scheme[
							key("pdb_strand_id") == dbref.chainID and
							key("pdb_seq_num") == dbref.seqEnd and
							key("pdb_ins_code") == ins_to_str(dbref.insertEnd)]
						["seq_id"].as<int>();
				}
				catch (...) {}
			
				get_category("struct_ref_seq")->emplace({
					{ "align_id", struct_ref_seq_align_id },
					{ "ref_id", struct_ref_id },
					{ "pdbx_PDB_id_code", dbref.PDBIDCode },
					{ "pdbx_strand_id", string{ chain.m_dbref.chainID } },
					{ "seq_align_beg", seq_align_beg },
					{ "pdbx_seq_align_beg_ins_code", ins_to_str(dbref.insertBegin) },
					{ "seq_align_end", seq_align_end },
					{ "pdbx_seq_align_end_ins_code", ins_to_str(dbref.insertEnd) },
					{ "pdbx_db_accession", dbref.dbAccession },
					{ "db_align_beg", dbref.dbSeqBegin },
					{ "pdbx_db_align_beg_ins_code", ins_to_str(dbref.dbinsBeg) },
					{ "db_align_end", dbref.dbSeqEnd },
					{ "pdbx_db_align_end_ins_code", ins_to_str(dbref.dbinsEnd) },
					{ "pdbx_auth_seq_align_beg", dbref.seqBegin },
					{ "pdbx_auth_seq_align_end", dbref.seqEnd }
				});

				// write the struct_ref_seq_dif
				for (auto& seqadv: m_seqadvs)
				{
					if (seqadv.chainID != chain.m_dbref.chainID or seqadv.resName.empty())
						continue;
					
					string asym, seq_num;
					int label_seq = -1;
					boost::system::error_code ec;

					tie(asym, label_seq, ignore) = MapResidue(seqadv.chainID, seqadv.seqNum, seqadv.iCode, ec);
					if (ec)
					{
						if (VERBOSE)
							cerr << "dropping unmatched SEQADV record" << endl;
						continue;
					}
					
					seq_num = to_string(label_seq);
					
					get_category("struct_ref_seq_dif")->emplace({
						{ "align_id", struct_ref_seq_align_id },
						{ "pdbx_PDB_id_code", dbref.PDBIDCode },
						{ "mon_id", seqadv.resName },
						{ "pdbx_pdb_strand_id", seqadv.chainID },
						{ "seq_num", seq_num },
						{ "pdbx_pdb_ins_code", seqadv.iCode == ' ' ? string{} : string{seqadv.iCode} },
						{ "pdbx_seq_db_name", seqadv.database },
						{ "pdbx_seq_db_accession_code", seqadv.dbAccession },
						{ "db_mon_id", seqadv.dbRes },
						{ "pdbx_seq_db_seq_num", seqadv.dbSeq },
						{ "details", seqadv.conflict },
						{ "pdbx_auth_seq_num", seqadv.seqNum },
						{ "pdbx_ordinal", ++m_pdbx_dif_ordinal }						
					});
				}
			}
			
			if (not chains.empty())	// not the first one for this mol_id
			{
				chains.push_back( string{ chain.m_dbref.chainID } );
				continue;
			}
			
			chains.push_back(string{chain.m_dbref.chainID});

			size_t seq_len = 0, seq_can_len = 0;
			
			for (auto& res: chain.m_seqres)
			{
				string letter, stdRes;

				if (m_mod2parent.count(res.m_mon_id))
					stdRes = m_mod2parent.at(res.m_mon_id);

				if (kAAMap.count(res.m_mon_id))
				{
					letter = kAAMap.at(res.m_mon_id);
					mightBeDNA = false;
				}
				else if (kBaseMap.count(res.m_mon_id))
				{
					letter = kBaseMap.at(res.m_mon_id);
					mightBePolyPeptide = false;
				}
				else
				{
					nstd_monomer = true;
					letter = '(' + res.m_mon_id + ')';
					
					// sja...
					auto compound = libcif::compound::create(stdRes.empty() ? res.m_mon_id : stdRes);
					if (compound != nullptr and
						not iequals(compound->type(), "L-peptide linking") and
						not iequals(compound->type(), "RNA linking"))
					{
						nonstandard_linkage = true;
					}
				}
				
				if (seq_len + letter.length() > 80)
				{
					seq += '\n';
					seq_len = 0;
				}

				seq += letter;
				seq_len += letter.length();
	
				if (letter.length() > 1)
				{
					if (not stdRes.empty() and kAAMap.count(stdRes))
						letter = kAAMap.at(stdRes);
					else if (kBaseMap.count(res.m_mon_id))
						letter = kBaseMap.at(res.m_mon_id);
					else
						letter = 'X';
				}
	
				if (seq_can_len + letter.length() > 80)
				{
					seq_can += '\n';
					seq_can_len = 0;
				}
				seq_can += letter;
				seq_can_len += letter.length();
			}

			auto cat = get_category("entity_poly_seq");
			for (size_t i = 0; i < chain.m_seqres.size(); ++i)
			{
				auto& rs = chain.m_seqres[i];
				
				cat->emplace({
					{ "entity_id", m_molID2EntityID[cmp.m_mol_id] },
					{ "num", i + 1 }, 
					{ "mon_id", rs.m_mon_id },
					{ "hetero", rs.m_alts.empty() ? "n" : "y"}
				});
				
				for (auto& a: rs.m_alts)
				{
					cat->emplace({
						{ "entity_id", m_molID2EntityID[cmp.m_mol_id] },
						{ "num", i + 1 }, 
						{ "mon_id", a },
						{ "hetero", "y"}
					});
				}
			}
		}
		
		string type;
		if (mightBePolyPeptide and not mightBeDNA)
			type = "polypeptide(L)";
		else if (mightBeDNA and not mightBePolyPeptide)
			type = "polyribonucleotide";

		get_category("entity_poly")->emplace({
			{ "entity_id", m_molID2EntityID[cmp.m_mol_id] },
			{ "pdbx_seq_one_letter_code", seq },
			{ "pdbx_seq_one_letter_code_can", seq_can },
			{ "nstd_monomer", (nstd_monomer ? "yes" : "no") },
			{ "pdbx_strand_id", ba::join(chains, ",") },
			{ "nstd_linkage", nonstandard_linkage ? "yes" : "no" },
			{ "type", type }
		});
	}

	if (not (struct_title.empty() and struct_description.empty()))
	{
		get_category("struct")->emplace({
			{ "entry_id", m_structure_id },
			{ "title", ba::join(struct_title, ", ") },
			{ "pdbx_descriptor", ba::join(struct_description, ", ") }
		});
	}
	
	map<char,string> water_chains;
	map<tuple<string,string>,int> ndb_seq_num;	// for nonpoly scheme
	
	for (size_t i = 0; i < m_hets.size(); ++i)
	{
		auto& heti = m_hets[i];

		if (not heti.asym_id.empty())
			continue;
	
		if (heti.hetID == m_water_het_id or is_water(heti.hetID))
			continue;

		// See if this residue is part of SEQRES
		auto& chain = GetChainForID(heti.chainID);
		auto ih = find(chain.m_seqres.begin(), chain.m_seqres.end(), PDBSeqRes{heti.hetID, heti.seqNum, heti.iCode}); 
		
		// If so, skip it, it is not an entity then
		if (ih != chain.m_seqres.end())
			continue;
		
		heti.asym_id = cif_id_for_int(asym_nr++);
	}

	set<string> written_asyms;

	map<string,int> het_count;		// for pdbx_number_of_molecules
	for (auto& het: m_hets)
		het_count[het.hetID] += 1;
	
	for (auto& het: m_hets)
	{
		string hetID = het.hetID;
		
		auto& chain = GetChainForID(het.chainID);
		
		// See if this residue is part of SEQRES
		auto i = find(chain.m_seqres.begin(), chain.m_seqres.end(), PDBSeqRes{hetID, het.seqNum, het.iCode}); 
		
		// If so, skip it, it is not an entity then
		if (i != chain.m_seqres.end())
			continue;

		// See if we've already added it to the entities
		if (m_het2EntityID.count(hetID) == 0)
		{
			string entity_id = to_string(m_next_entity_nr++);
			m_het2EntityID[hetID] = entity_id;
			
			if (hetID == m_water_het_id)
			{
				get_category("entity")->emplace({
					{ "id", entity_id },
					{ "type", "water" },
					{ "src_method", "nat" },
					{ "pdbx_description", "water" },
					{ "pdbx_number_of_molecules", het_count[hetID] }
				});
			}
			else
			{
				if (m_hetnams[hetID].empty())
					m_hetnams[hetID] = PeptideDB::Instance().GetNameForResidue(hetID);
				
				get_category("entity")->emplace({
					{ "id", entity_id },
					{ "type", "non-polymer" },
					{ "src_method", "syn" },
					{ "pdbx_description", m_hetnams[hetID] },
					{ "details", m_hetsyns[hetID] },
					{ "pdbx_number_of_molecules", het_count[hetID] }
				});
			}

			// write a pdbx_entity_nonpoly record
			string name = m_hetnams[hetID];
			if (name.empty() and hetID == m_water_het_id)
				name = "water";
			get_category("pdbx_entity_nonpoly")->emplace({
				{ "entity_id", entity_id },
				{ "name", name },
				{ "comp_id", hetID }
			});
		}
		
		// create an asym for this het/chain combo, if needed

		string asym_id = het.asym_id;

		auto k = make_tuple(het.chainID, het.seqNum, het.iCode);
		if (m_chainSeq2AsymSeq.count(k) == 0)
		{
			if (hetID == m_water_het_id or is_water(hetID))
			{
				if (water_chains.count(het.chainID) == 0)
				{
					asym_id = cif_id_for_int(asym_nr++);
					water_chains[het.chainID] = asym_id;
				}
				else
					asym_id = water_chains[het.chainID];
			}
			else
				asym_id = het.asym_id;
			
			assert(asym_id.empty() == false);

			m_asymID2EntityID[asym_id] = m_het2EntityID[hetID];
			
			// NOTE, a nonpoly residue has no label_seq_id
			// but in pdbx_nonpoly_scheme there is such a number.
			// Since this number is not used anywhere else we
			// just use it here and do not store it in the table 
			m_chainSeq2AsymSeq[k] = make_tuple(asym_id, 0, false);

			if (written_asyms.count(asym_id) == 0)
			{
				written_asyms.insert(asym_id);
				get_category("struct_asym")->emplace({
					{ "id", asym_id },
					{ "pdbx_blank_PDB_chainid_flag", het.chainID == ' ' ? "Y" : "N" },
	//					pdbx_modified 
					{ "entity_id", m_het2EntityID[hetID] },
	//					details
				});
			}
		}

		int seq_nr = ++ndb_seq_num[make_tuple(hetID, asym_id)];
		
		string iCode{het.iCode};
		ba::trim(iCode);
		if (iCode.empty())
			iCode = { '.' };
		
		get_category("pdbx_nonpoly_scheme")->emplace({
			{ "asym_id", asym_id },
			{ "entity_id", m_het2EntityID[hetID] },
			{ "mon_id", hetID },
			{ "ndb_seq_num", seq_nr },
			{ "pdb_seq_num", het.seqNum },
			{ "auth_seq_num", het.seqNum },	// ????
			{ "pdb_mon_id", hetID },
			{ "auth_mon_id", hetID },
			{ "pdb_strand_id", string{het.chainID} },
			{ "pdb_ins_code", iCode }
		});

		// mapping needed?
		m_chainSeq2AsymSeq[make_tuple(het.chainID, het.seqNum, het.iCode)] = make_tuple(asym_id, seq_nr, false);
	}

	int mod_res_id = 1;
	set<string> mod_res_set;
	for (auto rec = FindRecord("MODRES"); rec != nullptr and rec->is("MODRES");
			rec = rec->m_next)					//	 1 -  6        Record name   "MODRES"                                            												
	{								 			//	 8 - 11        IDcode        idCode      ID code of this datablock.                  
		string resName		= rec->v_s(13, 15);	//	13 - 15        Residue name  resName     Residue name used in this datablock.        
		char chainID		= rec->v_c(17);		//	17             Character     chainID     Chain identifier.                       
		int seqNum			= rec->v_i(19, 22);	//	19 - 22        Integer       seqNum      Sequence number.                        
		char iCode			= rec->v_c(23);		//	23             AChar         iCode       Insertion code.                         
		string stdRes		= rec->v_s(25, 27);	//	25 - 27        Residue name  stdRes      Standard residue name.                  
		string comment		= rec->v_s(30, 70);	//	30 - 70        String        comment     Description of the residue modification.

		string asym_id;
		int seq;
		boost::system::error_code ec;
		
		tie(asym_id, seq, ignore) = MapResidue(chainID, seqNum, iCode, ec);
		if (ec)	// no need to write a modres if it could not be found
		{
			if (VERBOSE)
				cerr << "dropping unmapped MODRES record" << endl;
			continue;
		}
		
		get_category("pdbx_struct_mod_residue")->emplace({
			{ "id", mod_res_id++ },
			{ "label_asym_id", asym_id },
			{ "label_seq_id", seq },
			{ "label_comp_id", resName },
			{ "auth_asym_id", string(1, chainID) },
			{ "auth_seq_id", seqNum },
			{ "auth_comp_id", resName },
			{ "PDB_ins_code", iCode == ' ' ? "" : string{ iCode } },
			{ "parent_comp_id", stdRes },
			{ "details", comment }
		});
		
		mod_res_set.insert(resName);
	}

//	// chem compounds

	for (auto cc: m_chem_comp)
	{
		auto compound = libcif::compound::create(
			m_mod2parent.count(cc) ? m_mod2parent[cc] : cc
		);
		
		string formula = m_formuls[cc];
		if (formula.empty() and compound != nullptr)
			formula = compound->formula();
		else
		{
			const regex rx(R"(\d+\((.+)\))");
			smatch m;
			if (regex_match(formula, m, rx))
				formula = m[1].str();
		}
		
		string name = m_hetnams[cc];
		if (name.empty() and compound != nullptr)
			name = compound->name();

		string type = "other";
		string nstd = ".";
		
		if (compound != nullptr)
		{
			type = compound->type();
			
			if (type.empty())
				type = "NON-POLYMER";

			if (iequals(type, "L-peptide linking") or iequals(type, "peptide linking"))
				nstd = "y";
		}

		if (mod_res_set.count(cc))
			nstd = "n";

		get_category("chem_comp")->emplace({
			{ "id", cc },
			{ "name", name },
			{ "formula", formula },
			{ "mon_nstd_flag", nstd },
			{ "type", type }
		});
	}
	
	get_category("chem_comp")->reorderByIndex();
	
	// unobserved can now be written as well
	
	int id_res = 0, id_atom = 0;
	sort(m_unobs.begin(), m_unobs.end(), [](const UNOBS& a, const UNOBS& b) -> bool
	{
		int d = a.model_nr - b.model_nr;
		if (d == 0)
			d = a.seq - b.seq;
		return d < 0;
	}); 
	
	for (auto& unobs: m_unobs)
	{
		bool is_polymer = false;
		string asym_id, comp_id = unobs.res;
		int seq_nr = 0;
		boost::system::error_code ec;
		
		tie(asym_id, seq_nr, is_polymer) = MapResidue(unobs.chain, unobs.seq, unobs.iCode, ec);
		if (ec)
		{
			if (VERBOSE)
				cerr << "error mapping unobserved residue" << endl;
			continue;
		}
		
		if (unobs.atoms.empty())
		{
			get_category("pdbx_unobs_or_zero_occ_residues")->emplace({
				{ "id",				to_string(++id_res) },
				{ "polymer_flag",	is_polymer ? "Y" : "N" },
				{ "occupancy_flag",	1 },
				{ "PDB_model_num",	unobs.model_nr ? unobs.model_nr : 1 },
				{ "auth_asym_id",	unobs.chain },
				{ "auth_comp_id",	unobs.res },
				{ "auth_seq_id",	unobs.seq },
				{ "PDB_ins_code",	unobs.iCode == ' ' ? "" : string{ unobs.iCode } },
				{ "label_asym_id",	asym_id },
				{ "label_comp_id",	comp_id },		// TODO: change to correct comp_id
				{ "label_seq_id",	seq_nr > 0 ? to_string(seq_nr) : "" }
			});
		}
		else
		{
			for (auto& atom: unobs.atoms)
			{
				get_category("pdbx_unobs_or_zero_occ_atoms")->emplace({
					{ "id",				to_string(++id_atom) },
					{ "polymer_flag",	is_polymer ? "Y" : "N" },
					{ "occupancy_flag",	1 },
					{ "PDB_model_num",	unobs.model_nr ? unobs.model_nr : 1 },
					{ "auth_asym_id",	unobs.chain },
					{ "auth_comp_id",	unobs.res },
					{ "auth_seq_id",	unobs.seq },
					{ "PDB_ins_code",	unobs.iCode == ' ' ? "" : string{ unobs.iCode } },
					{ "auth_atom_id",	atom },
					{ "label_asym_id",	asym_id },
					{ "label_comp_id",	comp_id },		// TODO: change to correct comp_id
					{ "label_seq_id",	seq_nr > 0 ? to_string(seq_nr) : "" },
					{ "label_atom_id",	atom }
				});
				
			}
		}
	}
}

void PDBFileParser::ParseSecondaryStructure()
{
	bool firstHelix = true;
	
	while (m_rec->is("HELIX "))
	{
		//	 1 -  6        Record name    "HELIX "
		//	 8 - 10        Integer        serNum        Serial number of the helix. This starts
		//	                                            at 1  and increases incrementally.
		//	12 - 14        LString(3)     helixID       Helix  identifier. In addition to a serial
		//	                                            number, each helix is given an 
		//	                                            alphanumeric character helix identifier.
		//	16 - 18        Residue name   initResName   Name of the initial residue.
		//	20             Character      initChainID   Chain identifier for the chain containing
		//	                                            this  helix.
		//	22 - 25        Integer        initSeqNum    Sequence number of the initial residue.
		//	26             AChar          initICode     Insertion code of the initial residue.
		//	28 - 30        Residue  name  endResName    Name of the terminal residue of the helix.
		//	32             Character      endChainID    Chain identifier for the chain containing
		//	                                            this  helix.
		//	34 - 37        Integer        endSeqNum     Sequence number of the terminal residue.
		//	38             AChar          endICode      Insertion code of the terminal residue.
		//	39 - 40        Integer        helixClass    Helix class (see below).
		//	41 - 70        String         comment       Comment about this helix.
		//	72 - 76        Integer        length        Length of this helix.

		string beg_asym_id, end_asym_id;
		int beg_seq, end_seq;
		boost::system::error_code ec;
		
		tie(beg_asym_id, beg_seq, ignore) = MapResidue(v_c(20), v_i(22, 25), v_c(26), ec);
		if (not ec)
			tie(end_asym_id, end_seq, ignore) = MapResidue(v_c(32), v_i(34, 37), v_c(38), ec);
		
		if (ec)
		{
			if (VERBOSE)
				cerr << "Could not map residue for HELIX " << v_i(8, 10) << endl;
		}
		else
		{
			auto cat = get_category("struct_conf");
			cat->emplace({
				{ "conf_type_id", "HELX_P" },
				{ "id", "HELX_P" + to_string(v_i(8, 10)) },
				{ "pdbx_PDB_helix_id", v_s(12, 14) },
				{ "beg_label_comp_id", v_s(16, 18) },
				{ "beg_label_asym_id", beg_asym_id },
				{ "beg_label_seq_id", beg_seq },
				{ "pdbx_beg_PDB_ins_code", v_s(26, 26) },
				{ "end_label_comp_id", v_s(28, 30) },
				{ "end_label_asym_id", end_asym_id },
				{ "end_label_seq_id", end_seq },
				{ "pdbx_end_PDB_ins_code", v_s(38, 38) },
	
				{ "beg_auth_comp_id", v_s(16, 18) },
				{ "beg_auth_asym_id", v_s(20, 20) },
				{ "beg_auth_seq_id", v_i(22, 25) },
				{ "end_auth_comp_id", v_s(28, 30) },
				{ "end_auth_asym_id", v_s(32, 32) },
				{ "end_auth_seq_id", v_i(34, 37) },
	
				{ "pdbx_PDB_helix_class", v_s(39, 40) },
				{ "details", v_s(41, 70) },
				{ "pdbx_PDB_helix_length", v_i(72, 76) }
			});
	
			if (firstHelix)
			{
				cat = get_category("struct_conf_type");
				cat->emplace({
					{ "id", "HELX_P" }
				});
				firstHelix = false;
			}
		}

		GetNextRecord();
	}
	
	set<string> sheetsSeen;
	int rangeID = 1;
	
	while (m_rec->is("SHEET "))
	{
		//	 1 -  6        Record name   "SHEET "
		//	 8 - 10        Integer       strand         Strand  number which starts at 1 for each
		//	                                            strand within a sheet and increases by one.
		//	12 - 14        LString(3)    sheetID        Sheet  identifier.
		//	15 - 16        Integer       numStrands     Number  of strands in sheet.
		//	18 - 20        Residue name  initResName    Residue  name of initial residue.
		//	22             Character     initChainID    Chain identifier of initial residue 
		//	                                            in strand. 
		//	23 - 26        Integer       initSeqNum     Sequence number of initial residue
		//	                                            in strand.
		//	27             AChar         initICode      Insertion code of initial residue
		//	                                            in  strand.
		//	29 - 31        Residue name  endResName     Residue name of terminal residue.
		//	33             Character     endChainID     Chain identifier of terminal residue.
		//	34 - 37        Integer       endSeqNum      Sequence number of terminal residue.
		//	38             AChar         endICode       Insertion code of terminal residue.
		//	39 - 40        Integer       sense          Sense of strand with respect to previous
		//	                                            strand in the sheet. 0 if first strand,
		//	                                            1 if  parallel,and -1 if anti-parallel.
		//	42 - 45        Atom          curAtom        Registration.  Atom name in current strand.
		//	46 - 48        Residue name  curResName     Registration.  Residue name in current strand
		//	50             Character     curChainId     Registration. Chain identifier in
		//	                                            current strand.
		//	51 - 54        Integer       curResSeq      Registration.  Residue sequence number
		//	                                            in current strand.
		//	55             AChar         curICode       Registration. Insertion code in
		//	                                            current strand.
		//	57 - 60        Atom          prevAtom       Registration.  Atom name in previous strand.
		//	61 - 63        Residue name  prevResName    Registration.  Residue name in
		//	                                            previous strand.
		//	65             Character     prevChainId    Registration.  Chain identifier in
		//	                                            previous  strand.
		//	66 - 69        Integer       prevResSeq     Registration. Residue sequence number
		//	                                            in previous strand.
		//	70             AChar         prevICode      Registration.  Insertion code in
		//	                                            previous strand.
		
		string sheetID = ba::trim_copy(v_s(12, 14));
		if (sheetsSeen.count(sheetID) == 0)
		{
			sheetsSeen.insert(sheetID);

			rangeID = 1;

			get_category("struct_sheet")->emplace({
				{ "id", sheetID },
				{ "number_strands", v_i(15, 16) },
			});
		}
		
		int sense = v_i(39, 40);
		
		if (sense != 0)
		{
			get_category("struct_sheet_order")->emplace({
				{ "sheet_id", sheetID },
				{ "range_id_1", rangeID },
				{ "range_id_2", rangeID + 1 },
				{ "sense", sense == -1 ? "anti-parallel" : "parallel" }
			});
		}

		string beg_asym_id, end_asym_id;
		int beg_seq, end_seq;
		boost::system::error_code ec;

		tie(beg_asym_id, beg_seq, ignore) = MapResidue(v_c(22), v_i(23, 26), v_c(27), ec);
		if (not ec)
			tie(end_asym_id, end_seq, ignore) = MapResidue(v_c(33), v_i(34, 37), v_c(38), ec);
		
		if (ec)
		{
			if (VERBOSE)
				cerr << "Dropping SHEET record " << v_i(8, 10) << endl;
		}
		else
		{
			get_category("struct_sheet_range")->emplace({
				{ "sheet_id", sheetID },
				{ "id", v_i(8, 10) },
				{ "beg_label_comp_id", v_s(18, 20) },
				{ "beg_label_asym_id", beg_asym_id },
				{ "beg_label_seq_id", beg_seq },
				{ "pdbx_beg_PDB_ins_code", v_s(27, 27) },
				{ "end_label_comp_id", v_s(29, 31) },
				{ "end_label_asym_id", end_asym_id },
				{ "end_label_seq_id", end_seq },
				{ "pdbx_end_PDB_ins_code", v_s(38, 38) },
				
				{ "beg_auth_comp_id", v_s(18, 20) },
				{ "beg_auth_asym_id", v_s(22, 22) },
				{ "beg_auth_seq_id", v_i(23, 26) },
				{ "end_auth_comp_id", v_s(29, 31) },
				{ "end_auth_asym_id", v_s(33, 33) },
				{ "end_auth_seq_id", v_i(34, 37) },
			});
			
			if (sense != 0 and m_rec->m_vlen > 34)
			{
				string r1_asym_id, r2_asym_id;
				int r1_seq, r2_seq;
				boost::system::error_code ec;
				
				tie(r1_asym_id, r1_seq, ignore) = MapResidue(v_c(65), v_i(66, 69), v_c(70), ec);
				if (not ec)
					tie(r2_asym_id, r2_seq, ignore) = MapResidue(v_c(50), v_i(51, 54), v_c(55), ec);
				
				if (ec)
				{
					if (VERBOSE)
						cerr << "skipping unmatched pdbx_struct_sheet_hbond record" << endl;
				}
				else
					get_category("pdbx_struct_sheet_hbond")->emplace({
						{ "sheet_id", sheetID },
						{ "range_id_1", rangeID },
						{ "range_id_2", rangeID + 1 },
						{ "range_1_label_atom_id", v_s(57, 60) },
						{ "range_1_label_comp_id", v_s(61, 63) },
						{ "range_1_label_asym_id", r1_asym_id },
						{ "range_1_label_seq_id", r1_seq },
						{ "range_1_PDB_ins_code", v_s(70, 70) },
						{ "range_1_auth_atom_id", v_s(57, 60) },
						{ "range_1_auth_comp_id", v_s(61, 63) },
						{ "range_1_auth_asym_id", v_s(65, 65) },
						{ "range_1_auth_seq_id", v_i(66, 69) },
		
						{ "range_2_label_atom_id", v_s(42, 45) },
						{ "range_2_label_comp_id", v_s(46, 48) },
						{ "range_2_label_asym_id", r2_asym_id },
						{ "range_2_label_seq_id", r2_seq },
						{ "range_2_PDB_ins_code", v_s(55, 55) },
						{ "range_2_auth_atom_id", v_s(42, 45) },
						{ "range_2_auth_comp_id", v_s(46, 48) },
						{ "range_2_auth_asym_id", v_s(50, 50) },
						{ "range_2_auth_seq_id", v_i(51, 54) }
					});
			}
			
			if (sense != 0)
				++rangeID;
		}
		
		GetNextRecord();
	}
}

static bool IsMetal(const string& resName, const string& atomID)
{
	bool result = false;

	try
	{
		auto compound = libcif::compound::create(resName);
		if (compound != nullptr)
		{
			auto at = libcif::atom_type_traits(compound->get_atom_by_id(atomID).type_symbol);
			result = at.is_metal();
		}
	}
	catch (...) {}

	return result;
}

void PDBFileParser::ParseConnectivtyAnnotation()
{
	int ssBondNr = 0;
	int linkNr = 0;
	bool firstCovale = true, firstMetalc = true;
	
	// Aaargh... Coot writes the records in the wrong order...
	for (;; GetNextRecord())
	{
		if (m_rec->is("SSBOND"))
		{
			if (ssBondNr == 0)
			{
				get_category("struct_conn_type")->emplace({
					{ "id", "disulf" },
				});
			}
	
			//	 1 -  6        Record name    "SSBOND"
			//	 8 - 10        Integer        serNum           Serial number.
			//	12 - 14        LString(3)     "CYS"            Residue name.
			//	16             Character      chainID1         Chain identifier.
			//	18 - 21        Integer        seqNum1          Residue sequence number.
			//	22             AChar          icode1           Insertion code.
			//	26 - 28        LString(3)     "CYS"            Residue name.
			//	30             Character      chainID2         Chain identifier.
			//	32 - 35        Integer        seqNum2          Residue sequence number.
			//	36             AChar          icode2           Insertion code.
			//	60 - 65        SymOP          sym1             Symmetry operator for residue 1.
			//	67 - 72        SymOP          sym2             Symmetry operator for residue 2.
			//	74 – 78        Real(5.2)      Length           Disulfide bond distance
	
			string p1_asym, p2_asym;
			int p1_seq, p2_seq;
			boost::system::error_code ec;
			
			tie(p1_asym, p1_seq, ignore) = MapResidue(v_c(16), v_i(18, 21), v_c(22), ec);
			if (not ec)
				tie(p2_asym, p2_seq, ignore) = MapResidue(v_c(30), v_i(32, 35), v_c(36), ec);
			
			if (ec)
			{
				if (VERBOSE)
					cerr << "Dropping SSBOND " << v_i(8, 10) << endl;
				continue;
			}

			vector<char> alt1 = alt_locs_for_atom(v_c(16), v_i(18, 21), v_c(22), "SG");
			vector<char> alt2 = alt_locs_for_atom(v_c(30), v_i(32, 35), v_c(36), "SG");
			
			if (alt1.empty())
				alt1.push_back(0);
			if (alt2.empty())
				alt2.push_back(0);
			
			for (auto a1: alt1)
			{
				for (auto a2: alt2)
				{
					get_category("struct_conn")->emplace({
						{ "id", "disulf" + to_string(++ssBondNr) },
						{ "conn_type_id", "disulf" },
			
						{ "ptnr1_label_asym_id", p1_asym },
						{ "pdbx_ptnr1_label_alt_id", a1 ? string{ a1 } : string() },
						{ "ptnr1_label_comp_id", v_s(12, 14) },
						{ "ptnr1_label_seq_id", p1_seq ? to_string(p1_seq) : "." },
						{ "ptnr1_label_atom_id", "SG" },
						{ "ptnr1_symmetry", pdb2cif_symmetry(v_s(60, 65)) },
			
						{ "ptnr2_label_asym_id", p2_asym },
						{ "pdbx_ptnr2_label_alt_id", a2 ? string{ a2 } : string() },
						{ "ptnr2_label_comp_id", v_s(26, 28) },
						{ "ptnr2_label_seq_id", p2_seq ? to_string(p2_seq) : "." },
						{ "ptnr2_label_atom_id", "SG" },
			
						{ "ptnr1_auth_asym_id", v_s(16, 16) },
						{ "ptnr1_auth_comp_id", v_s(12, 14) },
						{ "ptnr1_auth_seq_id", v_i(18, 21) },
						{ "ptnr2_auth_asym_id", v_s(30, 30) },
						{ "ptnr2_auth_comp_id", v_s(26, 28) },
						{ "ptnr2_auth_seq_id", v_i(32, 35) },
			
						{ "ptnr2_symmetry", pdb2cif_symmetry(v_s(67, 72)) },
						
						{ "pdbx_dist_value", v_s(74, 78) },
					});
				}
			}

			continue;
		}
		
		if (m_rec->is("LINK  "))
		{
											//	 1 -  6         Record name    "LINK  "
			string name1 = v_s(13, 16);		//	13 - 16         Atom           name1           Atom name.
											//	17              Character      altLoc1         Alternate location indicator.
			string resName1 = v_s(18,20);	//	18 - 20         Residue name   resName1        Residue  name.
											//	22              Character      chainID1        Chain identifier.
											//	23 - 26         Integer        resSeq1         Residue sequence number.
											//	27              AChar          iCode1          Insertion code.
			string name2 = v_s(43, 46);		//	43 - 46         Atom           name2           Atom name.
											//	47              Character      altLoc2         Alternate location indicator.
			string resName2 = v_s(48, 50);	//	48 - 50         Residue name   resName2        Residue name.
											//	52              Character      chainID2        Chain identifier.
											//	53 - 56         Integer        resSeq2         Residue sequence number.
											//	57              AChar          iCode2          Insertion code.
											//	60 - 65         SymOP          sym1            Symmetry operator atom 1.
											//	67 - 72         SymOP          sym2            Symmetry operator atom 2.
											//	74 – 78         Real(5.2)      Length          Link distance
	
			string type = "covale";
			if (IsMetal(resName1, name1) or IsMetal(resName2, name2))
				type = "metalc";
			
			if (type == "covale" and firstCovale)
			{
				get_category("struct_conn_type")->emplace({
					{ "id", type },
				});
				firstCovale = false;
			}
	
			if (type == "metalc" and firstMetalc)
			{
				get_category("struct_conn_type")->emplace({
					{ "id", type },
				});
				firstMetalc = false;
			}
			
			++linkNr;
	
			string p1_asym, p2_asym;
			int p1_seq, p2_seq;
			bool is_resseq1, is_resseq2;
			boost::system::error_code ec;
			
			tie(p1_asym, p1_seq, is_resseq1) = MapResidue(v_c(22), v_i(23, 26), v_c(27), ec);
			if (not ec)
				tie(p2_asym, p2_seq, is_resseq2) = MapResidue(v_c(52), v_i(53, 56), v_c(57), ec);
			
			if (ec)
			{
				if (VERBOSE)
					cerr << "Dropping LINK record at line " << m_rec->m_line_nr << endl;
				continue;
			}
			
			get_category("struct_conn")->emplace({
				{ "id", type + to_string(linkNr) },
				{ "conn_type_id", type },
	
				{ "ptnr1_label_asym_id", p1_asym },
				{ "ptnr1_label_comp_id", v_s(18, 20) },
				{ "ptnr1_label_seq_id", (is_resseq1 and p1_seq) ? to_string(p1_seq) : "." },
				{ "ptnr1_label_atom_id", v_s(13, 16) },
				{ "pdbx_ptnr1_label_alt_id", v_s(17, 17) },
				{ "pdbx_ptnr1_PDB_ins_code", v_s(27, 27) },
				{ "pdbx_ptnr1_standard_comp_id", "" },
				{ "ptnr1_symmetry", pdb2cif_symmetry(v_s(60, 65)) },
	
				{ "ptnr2_label_asym_id", p2_asym },
				{ "ptnr2_label_comp_id", v_s(48, 50) },
				{ "ptnr2_label_seq_id", (is_resseq2 and p2_seq) ? to_string(p2_seq) : "." },
				{ "ptnr2_label_atom_id", v_s(43, 46) },
				{ "pdbx_ptnr2_label_alt_id", v_s(47, 47) },
				{ "pdbx_ptnr2_PDB_ins_code", v_s(57, 57) },
	
				{ "ptnr1_auth_asym_id", v_s(22, 22) },
				{ "ptnr1_auth_comp_id", v_s(18, 20) },
				{ "ptnr1_auth_seq_id", v_i(23, 26) },
				{ "ptnr2_auth_asym_id", v_s(52, 52) },
				{ "ptnr2_auth_comp_id", v_s(48, 50) },
				{ "ptnr2_auth_seq_id", v_i(53, 56) },
	
				{ "ptnr2_symmetry", pdb2cif_symmetry(v_s(67, 72)) },
				
				{ "pdbx_dist_value", v_s(74, 78) }
			});
			
			continue;
		}
		
		if (m_rec->is("CISPEP"))
		{
											//	 1 -  6       Record name   "CISPEP"
			int serNum = v_i(8, 10);		//	 8 - 10       Integer       serNum        Record serial number.
			string pep1 = v_s(12, 14);		//	12 - 14       LString(3)    pep1          Residue name.
			char chainID1 = v_c(16); 		//	16            Character     chainID1      Chain identifier.
			int seqNum1 = v_i(18, 21);		//	18 - 21       Integer       seqNum1       Residue sequence number.
			char iCode1 = v_c(22);			//	22            AChar         icode1        Insertion code.
			string pep2 = v_s(26, 28);		//	26 - 28       LString(3)    pep2          Residue name.
			char chainID2 = v_c(30); 		//	30            Character     chainID2      Chain identifier.
			int seqNum2 = v_i(32, 35);		//	32 - 35       Integer       seqNum2       Residue sequence number.
			char iCode2 = v_c(36);			//	36            AChar         icode2        Insertion code.
			int modNum = v_i(44, 46);		//	44 - 46       Integer       modNum        Identifies the specific model.
			string measure = v_f(54, 59);	//	54 - 59       Real(6.2)     measure       Angle measurement in degrees.
			
			if (modNum == 0)
				modNum = 1;
	
			string l_asym1, l_asym2;
			int l_resSeq1, l_resSeq2;
			boost::system::error_code ec;
	
			tie(l_asym1, l_resSeq1, ignore) = MapResidue(chainID1, seqNum1, iCode1, ec);
			if (not ec)
				tie(l_asym2, l_resSeq2, ignore) = MapResidue(chainID2, seqNum2, iCode2, ec);
			
			if (ec)
			{
				if (VERBOSE)
					cerr << "Dropping CISPEP record at line " << m_rec->m_line_nr << endl;
				continue;
			}
			
			string iCode1str = iCode1 == ' ' ? string() : string{ iCode1 };
			string iCode2str = iCode2 == ' ' ? string() : string{ iCode2 };
			
			get_category("struct_mon_prot_cis")->emplace({
				{ "pdbx_id", serNum },
				{ "label_comp_id", pep1 },
				{ "label_seq_id", l_resSeq1 },
				{ "label_asym_id", l_asym1 },
				{ "label_alt_id", "." },
				{ "pdbx_PDB_ins_code", iCode1str },
				{ "auth_comp_id", pep1 },
				{ "auth_seq_id", seqNum1 },
				{ "auth_asym_id", string{chainID1} },
				{ "pdbx_label_comp_id_2", pep2 },
				{ "pdbx_label_seq_id_2", l_resSeq2 },
				{ "pdbx_label_asym_id_2", l_asym2 },
				{ "pdbx_PDB_ins_code_2", iCode2str },
				{ "pdbx_auth_comp_id_2", pep2 },
				{ "pdbx_auth_seq_id_2", seqNum2 },
				{ "pdbx_auth_asym_id_2", string{chainID2} },
				{ "pdbx_PDB_model_num", modNum },
				{ "pdbx_omega_angle", measure }
			});
			
			continue;
		}
		
		break;
	}
}

void PDBFileParser::ParseMiscellaneousFeatures()
{
	int struct_site_gen_id = 1;
	
	while (m_rec->is("SITE  "))
	{									//	 1 -  6        Record name   "SITE  "
										//	 8 - 10        Integer       seqNum        Sequence number.
		string siteID = v_s(12, 14);	//	12 - 14        LString(3)    siteID        Site name.
		int numRes = v_i(16, 17);		//	16 - 17        Integer       numRes        Number of residues that compose the site.

		int o = 19;
		
		auto cat = get_category("struct_site_gen");

		for (int i = 0; i < numRes; ++i)
		{
			string resName = v_s(o, o + 2);	//	19 - 21        Residue name  resName1      Residue name for first residue that 
											//	                                           creates the site.
			char chainID = v_c(o + 4);		//	23             Character     chainID1      Chain identifier for first residue of site.
			int seq = v_i(o + 5, o + 8);	//	24 - 27        Integer       seq1          Residue sequence number for first residue
											//	                                           of the  site.
			char iCode = v_c(o + 9);		//	28             AChar         iCode1        Insertion code for first residue of the site.

			int label_seq;
			string asym;
			bool is_resseq;
			boost::system::error_code ec;
			
			tie(asym, label_seq, is_resseq) = MapResidue(chainID, seq, iCode, ec);
			
			if (ec)
			{
				if (VERBOSE)
					cerr << "skipping struct_site_gen record" << endl;
			}
			else
				cat->emplace({
					{ "id", struct_site_gen_id++ },
					{ "site_id", siteID },
					{ "pdbx_num_res", numRes },
					{ "label_comp_id", resName },
					{ "label_asym_id", asym },
					{ "label_seq_id", (label_seq > 0 and is_resseq) ? to_string(label_seq) : string(".") },
					{ "pdbx_auth_ins_code", iCode == ' ' ? "" : string { iCode } },
					{ "auth_comp_id", resName },
					{ "auth_asym_id", string { chainID } },
					{ "auth_seq_id", seq },
					{ "label_atom_id", "." },
					{ "label_alt_id", "." },
				});

			o += 11;
		}

		GetNextRecord();
	}
}

void PDBFileParser::ParseCrystallographic()
{
	Match("CRYST1");

	get_category("cell")->emplace({		
		{ "entry_id", m_structure_id },			//	 1 -  6       Record name   "CRYST1"                          
		{ "length_a", v_f(7, 15) },             //	 7 - 15       Real(9.3)     a              a (Angstroms).     
		{ "length_b", v_f(16, 24) },            //	16 - 24       Real(9.3)     b              b (Angstroms).     
		{ "length_c", v_f(25, 33) },            //	25 - 33       Real(9.3)     c              c (Angstroms).     
		{ "angle_alpha", v_f(34, 40) },         //	34 - 40       Real(7.2)     alpha          alpha (degrees).   
		{ "angle_beta", v_f(41, 47) },          //	41 - 47       Real(7.2)     beta           beta (degrees).    
		{ "angle_gamma", v_f(48, 54) },         //	48 - 54       Real(7.2)     gamma          gamma (degrees).   
		/* goes into symmetry */				//	56 - 66       LString       sGroup         Space  group.      
		{ "Z_PDB", v_f(67, 70) }                //	67 - 70       Integer       z              Z value.           
	});
	
	get_category("symmetry")->emplace({
		{ "entry_id", m_structure_id },
		{ "space_group_name_H-M", v_s(56, 66) }
	});

	GetNextRecord();
}

void PDBFileParser::ParseCoordinateTransformation()
{
	string m[3][3], v[3];
	
	if (ba::starts_with(m_rec->m_name, "ORIGX"))
	{
		for (string n: { "1", "2", "3" })
		{
			int x = stoi(n) - 1;

			Match("ORIGX" + n);			//	 1 -  6         Record name   "ORIGXn"      n=1, 2, or 3  	
			m[x][0] = v_f(11, 20);  	//	11 - 20         Real(10.6)    o[n][1]       On1           
			m[x][1] = v_f(21, 30);  	//	21 - 30         Real(10.6)    o[n][2]       On2           
			m[x][2] = v_f(31, 40);  	//	31 - 40         Real(10.6)    o[n][3]       On3           
			v[x] = v_f(46, 55);     	//	46 - 55         Real(10.5)    t[n]          Tn            
			
			GetNextRecord();
		}

		get_category("database_PDB_matrix")->emplace({
			{ "entry_id", m_structure_id },
			{ "origx[1][1]", m[0][0] },
			{ "origx[1][2]", m[0][1] },
			{ "origx[1][3]", m[0][2] },
			{ "origx[2][1]", m[1][0] },
			{ "origx[2][2]", m[1][1] },
			{ "origx[2][3]", m[1][2] },
			{ "origx[3][1]", m[2][0] },
			{ "origx[3][2]", m[2][1] },
			{ "origx[3][3]", m[2][2] },
			{ "origx_vector[1]", v[0] },
			{ "origx_vector[2]", v[1] },
			{ "origx_vector[3]", v[2] },
		});
	}

	if (ba::starts_with(m_rec->m_name, "SCALE"))
	{
		for (string n: { "1", "2", "3" })
		{
			int x = stoi(n) - 1;

			Match("SCALE" + n);			//	 1 -  6         Record name   "SCALEn" n=1,  2, or 3      	
			m[x][0] = v_f(11, 20);  	//	11 - 20         Real(10.6)    s[n][1]            Sn1      
			m[x][1] = v_f(21, 30);  	//	21 - 30         Real(10.6)    s[n][2]            Sn2      
			m[x][2] = v_f(31, 40);  	//	31 - 40         Real(10.6)    s[n][3]            Sn3      
			v[x] = v_f(46, 55);     	//	46 - 55         Real(10.5)    u[n]               Un       

			GetNextRecord();
		}

		get_category("atom_sites")->emplace({
			{ "entry_id", m_structure_id },
			{ "fract_transf_matrix[1][1]", m[0][0] },
			{ "fract_transf_matrix[1][2]", m[0][1] },
			{ "fract_transf_matrix[1][3]", m[0][2] },
			{ "fract_transf_matrix[2][1]", m[1][0] },
			{ "fract_transf_matrix[2][2]", m[1][1] },
			{ "fract_transf_matrix[2][3]", m[1][2] },
			{ "fract_transf_matrix[3][1]", m[2][0] },
			{ "fract_transf_matrix[3][2]", m[2][1] },
			{ "fract_transf_matrix[3][3]", m[2][2] },
			{ "fract_transf_vector[1]", v[0] },
			{ "fract_transf_vector[2]", v[1] },
			{ "fract_transf_vector[3]", v[2] },
		});
	}

	while (ba::starts_with(m_rec->m_name, "MTRIX1"))
	{
		int serial, igiven;
		
		for (string n: { "1", "2", "3" })
		{
			int x = stoi(n) - 1;

			Match("MTRIX" + n);				//	 1 -  6        Record name   "MTRIXn"      n=1, 2, or 3
			serial = v_i(8, 10);			//	 8 - 10        Integer       serial        Serial number.
			m[x][0] = v_f(11, 20);  		//	11 - 20        Real(10.6)    m[n][1]       Mn1
			m[x][1] = v_f(21, 30);  		//	21 - 30        Real(10.6)    m[n][2]       Mn2
			m[x][2] = v_f(31, 40);  		//	31 - 40        Real(10.6)    m[n][3]       Mn3
			v[x] = v_f(46, 55);     		//	46 - 55        Real(10.5)    v[n]          Vn
			igiven = v_c(60) == '1';		//	60             Integer       iGiven        1 if coordinates for the  representations
											//	                                           which  are approximately related by the 
			GetNextRecord();				//	                                           transformations  of the molecule are
		}									//	                                           contained in the datablock. Otherwise, blank.

		get_category("struct_ncs_oper")->emplace({
			{ "id", serial },
			{ "matrix[1][1]", m[0][0] },
			{ "matrix[1][2]", m[0][1] },
			{ "matrix[1][3]", m[0][2] },
			{ "matrix[2][1]", m[1][0] },
			{ "matrix[2][2]", m[1][1] },
			{ "matrix[2][3]", m[1][2] },
			{ "matrix[3][1]", m[2][0] },
			{ "matrix[3][2]", m[2][1] },
			{ "matrix[3][3]", m[2][2] },
			{ "vector[1]", v[0] },
			{ "vector[2]", v[1] },
			{ "vector[3]", v[2] },
			{ "code", igiven ? "given" : "" }
		});
	}
}

void PDBFileParser::ParseCoordinate(int model_nr)
{
	// oh oh, we have to sort our atom_site records by ascending asym_id
	// This routine used to be so trivial...
	
	typedef tuple<string,int,bool,PDBRecord*,PDBRecord*> atom_rec;
	
	vector<atom_rec> atoms;
	while (m_rec->is("ATOM  ") or m_rec->is("HETATM"))		//	 1 -  6        Record name   "ATOM  "
	{
		char chainID = v_c(22);				//	22             Character     chainID      Chain identifier.
		int resSeq = v_i(23, 26);			//	23 - 26        Integer       resSeq       Residue sequence number.
		char iCode = v_c(27);
		
		string asym_id;
		int seq_id;
		bool is_resseq;
		
		tie(asym_id, seq_id, is_resseq) = MapResidue(chainID, resSeq, iCode);
		
		PDBRecord* atom = m_rec;
		PDBRecord* anisou = nullptr;

		GetNextRecord();
		if (m_rec->is("ANISOU"))
		{
			anisou = m_rec;
			GetNextRecord();
		}

		atoms.emplace_back(asym_id, seq_id, is_resseq, atom, anisou);

		/*if?... */ while (m_rec->is("TER   "))
		{
			Match("TER   ");
			GetNextRecord();
		}
	}
	
	auto last = m_rec;
	
	// use stable sort here
	auto r_less = [](const atom_rec& a, const atom_rec& b) -> bool
	{
		int d;
		
		string chainA = get<0>(a);
		string chainB = get<0>(b);
		
		if (chainA.length() != chainB.length())
			d = chainA.length() - chainB.length();
		else
			d = get<0>(a).compare(get<0>(b));

		if (d == 0)
			d = get<1>(a) - get<1>(b);
		return d < 0;
	};
	
	stable_sort(atoms.begin(), atoms.end(), r_less);
	
	// now reiterate the atoms to reorder alternates
	for (size_t i = 0; i + 1 < atoms.size(); ++i)
	{
		char altLoc = get<3>(atoms[i])->v_c(17);
		
		if (altLoc == ' ' or altLoc == 0)
			continue;
		
		auto b = atoms.begin() + i;
		auto e = b;
		
		map<string,int> atom_index;	// index number of first occurrence of a atom name

		while (e != atoms.end() and r_less(*b, *e) == false)
		{
			string name = get<3>(*e)->v_s(13, 16);
			
			if (atom_index.count(name) == 0)
				atom_index[name] = atom_index.size() + 1;
			
			++e;
		}
		
		auto a_less = [&](atom_rec& a, atom_rec& b) -> bool
		{
			string na = get<3>(a)->v_s(13, 16);
			string nb = get<3>(b)->v_s(13, 16);
			
			int d = atom_index[na] - atom_index[nb];
			if (d == 0)
				d = get<3>(a)->v_c(17) - get<3>(b)->v_c(17);
			assert(d != 0);
			return d < 0;
		};
		
		sort(b, e, a_less);
		
		i += distance(b, e) - 1;
	}

//	while (m_rec->is("ATOM  ") or m_rec->is("HETATM"))		//	 1 -  6        Record name   "ATOM  "
	for (auto& a: atoms)
	{
		string asym_id;
		int seq_id;
		bool is_resseq;
		PDBRecord* atom;
		PDBRecord* anisou;
		tie(asym_id, seq_id, is_resseq, atom, anisou) = a;
		
		m_rec = atom;
		
		++m_atom_id;

		string group_PDB = m_rec->is("ATOM  ") ? "ATOM" : "HETATM";
//		int serial = v_i(7, 11);			//	 7 - 11        Integer       serial       Atom  serial number.
		string name = v_s(13, 16);			//	13 - 16        Atom          name         Atom name.
		char altLoc = v_c(17);				//	17             Character     altLoc       Alternate location indicator.
		string resName = v_s(18, 20);		//	18 - 20        Residue name  resName      Residue name.
		char chainID = v_c(22);				//	22             Character     chainID      Chain identifier.
		int resSeq = v_i(23, 26);			//	23 - 26        Integer       resSeq       Residue sequence number.
		char iCode = v_c(27);				//	27             AChar         iCode        Code for insertion of residues.
		string x = v_f(31, 38);				//	31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		string y = v_f(39, 46);				//	39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
		string z = v_f(47, 54);				//	47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		string occupancy = v_f(55, 60);		//	55 - 60        Real(6.2)     occupancy    Occupancy.
		string tempFactor = v_f(61, 66);	//	61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		string element = v_s(77, 78);		//	77 - 78        LString(2)    element      Element symbol, right-justified.
		string charge = v_s(79, 80);		//	79 - 80        LString(2)    charge       Charge  on the atom.

		string entity_id = m_asymID2EntityID[asym_id];
		
		charge = pdb2cif_charge(charge);

		get_category("atom_site")->emplace({
			{ "group_PDB" , group_PDB },
			{ "id", m_atom_id },
			{ "type_symbol", element },
			{ "label_atom_id", name },
			{ "label_alt_id", altLoc != ' ' ? string { altLoc } : "." },
			{ "label_comp_id", resName },
			{ "label_asym_id", asym_id },
			{ "label_entity_id", entity_id },
			{ "label_seq_id", (is_resseq and seq_id > 0) ? to_string(seq_id) : "." },
			{ "pdbx_PDB_ins_code", iCode == ' ' ? "" : string { iCode } },
			{ "Cartn_x", x },
			{ "Cartn_y", y },
			{ "Cartn_z", z },
			{ "occupancy", occupancy },
			{ "B_iso_or_equiv", tempFactor },
			{ "pdbx_formal_charge", charge },
			{ "auth_seq_id", resSeq },
			{ "auth_comp_id", resName },
			{ "auth_asym_id", string { chainID } },
			{ "auth_atom_id", name },
			{ "pdbx_PDB_model_num", model_nr }
		});
		
		InsertAtomType(element);
		
		string check = v_s(7, 11) + v_s(77, 80);
		
		if (anisou != nullptr)
		{
			m_rec = anisou;						//	 1 - 6        Record name   "ANISOU"
			int u11 = v_i(29, 35);				//	29 - 35       Integer       u[0][0]        U(1,1)
			int u22 = v_i(36, 42);				//	36 - 42       Integer       u[1][1]        U(2,2)
			int u33 = v_i(43, 49);				//	43 - 49       Integer       u[2][2]        U(3,3)
			int u12 = v_i(50, 56);				//	50 - 56       Integer       u[0][1]        U(1,2)
			int u13 = v_i(57, 63);				//	57 - 63       Integer       u[0][2]        U(1,3)
			int u23 = v_i(64, 70);				//	64 - 70       Integer       u[1][2]        U(2,3)
			
			if (v_s(7, 11) + v_s(77, 80) != check)
				Error("ANISOU record should follow corresponding ATOM record");
			
			auto f = [](float f) -> string { return (boost::format("%6.4f") % f).str(); };
			
			get_category("atom_site_anisotrop")->emplace({
				{ "id", m_atom_id },
				{ "type_symbol", element }, 
				{ "pdbx_label_atom_id", name },
				{ "pdbx_label_alt_id", altLoc != ' ' ? string { altLoc } : "." },
				{ "pdbx_label_comp_id", resName },
				{ "pdbx_label_asym_id", asym_id },
				{ "pdbx_label_seq_id", seq_id },
				{ "U[1][1]", f(u11 / 10000.f) },
				{ "U[2][2]", f(u22 / 10000.f) },
				{ "U[3][3]", f(u33 / 10000.f) },
				{ "U[1][2]", f(u12 / 10000.f) },
				{ "U[1][3]", f(u13 / 10000.f) },
				{ "U[2][3]", f(u23 / 10000.f) },
				{ "pdbx_auth_seq_id", resSeq },
				{ "pdbx_auth_comp_id", resName },
				{ "pdbx_auth_asym_id", string { chainID } },
				{ "pdbx_auth_atom_id", name }
			});
		}
	}
	
	m_rec = last;
}

void PDBFileParser::ParseConnectivty()
{
	while (m_rec->is("CONECT"))
		GetNextRecord();
}

void PDBFileParser::ParseBookkeeping()
{
	if (m_rec->is("MASTER"))
	{
		Match("MASTER");
		GetNextRecord();
	}
	Match("END   ");
}

void PDBFileParser::Parse(istream& is, cif::file& result)
{
	try
	{
		PreParseInput(is);

		m_rec = m_data;
	
		ParseTitle();

		result.append(m_datablock);

		ParseRemarks();
		ParsePrimaryStructure();
		ParseHeterogen();
		
		ConstructEntities();
		
		ParseRemark350();
		
		ParseSecondaryStructure();
		ParseConnectivtyAnnotation();
		ParseMiscellaneousFeatures();
		ParseCrystallographic();
		ParseCoordinateTransformation();
	
		uint32 model_nr = 1;
	
		while (m_rec->is("MODEL ") or m_rec->is("ATOM  ") or m_rec->is("HETATM"))
		{
			bool model = false;
			if (m_rec->is("MODEL "))
			{
				model = true;
				
				model_nr = v_i(11, 14);
				
				GetNextRecord();
			}
			
			ParseCoordinate(model_nr);
			
			if (model)
			{
				Match("ENDMDL");
				GetNextRecord();
			}
		}	
	
		for (auto e: m_atom_types)
			get_category("atom_type")->emplace({
				{ "symbol", e }
			});

		// in V5, atom_type is sorted
		get_category("atom_type")->reorderByIndex();
	
		ParseConnectivty();
		ParseBookkeeping();
		
		// almost done, now fix some outstanding issued that could not be done before
		
		try
		{
			auto r = FindRecord("REMARK   3");
			
			if (r != nullptr and Remark3Parser::Parse(m_exp_method, r, *m_datablock))
			{
				// make sure the "exptl" category is created
				auto exptl = get_category("exptl"); 
				if (exptl->empty())
				{
					exptl->emplace({
						{ "entry_id", m_structure_id },
						{ "method", m_exp_method },
						{ "crystals_number", m_remark200["NUMBER OF CRYSTALS USED"] }
					});
				}
			}
		}
		catch (const exception& ex)
		{
			cerr << "Error parsing REMARK 3: " << endl
				 << ex.what() << endl;
		}
//		
//		
//		
//		auto cat = get_category("pdbx_refine_tls_group");
//		for (row r: *cat)
//		{
//			// add the mapped locations
//			
//			try
//			{
//				string asym_id;
//				int resNum;
//				
//				cif::tie(asym_id, resNum) = r.get("beg_auth_asym_id", "beg_auth_seq_id");
//				
//				r["beg_label_asym_id"] = asym_id;
//				r["beg_label_seq_id"] = resNum;
//				
//				cif::tie(asym_id, resNum) = r.get("end_auth_asym_id", "end_auth_seq_id");
//				
//				r["end_label_asym_id"] = asym_id;
//				r["end_label_seq_id"] = resNum;
//			}
//			catch (const exception& ex)
//			{
//				continue;
//			}
//		}
	}
	catch (const exception& ex)
	{
		Error(ex.what());
	}
}

// ----------------------------------------------------------------

void PDBFileParser::PDBChain::AlignResToSeqRes()
{
	// Use dynamic programming to align the found residues (in ATOM records) against
	// the residues in the SEQRES records in order to obtain the residue numbering.
	// sigh...
	
	using namespace boost::numeric::ublas;
	
	auto& rx = m_seqres;
	auto& ry = m_residues_seen;
	
	int dimX = m_seqres.size();
	if (dimX == 0)
		throw runtime_error(string("SEQRES for chain ") + m_dbref.chainID + " is empty");

	int dimY = m_residues_seen.size();
	if (dimY == 0)
		throw runtime_error(string("Number of residues in ATOM records for chain ") + m_dbref.chainID + " is zero");

	matrix<float> B(dimX, dimY), Ix(dimX, dimY), Iy(dimX, dimY);
	matrix<int8> tb(dimX, dimY);
	
	int x, y;
	
	const float kGapOpen = 10, gapExtend = 0.1;

	float high = 0;
	size_t highX = 0, highY = 0;
	
	for (x = 0; x < dimX; ++x)
	{
		for (y = 0; y < dimY; ++y)
		{
			auto& a = rx[x];
			auto& b = ry[y];

			float Ix1 = x > 0 ? Ix(x - 1, y) : 0;
			float Iy1 = y > 0 ? Iy(x, y - 1) : 0;
			
			// score for alignment
			float M;
			if (a.m_mon_id == b.m_mon_id)
				M = 1;
			else
				M = -10000;

			// gap open cost is zero if the PDB ATOM records indicate that a gap
			// should be here.
			float gapOpen = kGapOpen;
			if (y + 1 < dimY and ry[y + 1].m_seq_num > ry[y].m_seq_num + 1)
				gapOpen = 0;

			if (x > 0 and y > 0)
				M += B(x - 1, y - 1);
			
			float s;
			if (M >= Ix1 and M >= Iy1)
			{
				tb(x, y) = 0;
				B(x, y) = s = M;
				
				Ix(x, y) = M - (x < dimX - 1 ? gapOpen : 0);
				Iy(x, y) = M - (y < dimY - 1 ? gapOpen : 0);
			}
			else if (Ix1 >= Iy1)
			{
				tb(x, y) = 1;
				B(x, y) = s = Ix1;

				Ix(x, y) = Ix1 - gapExtend;
				Iy(x, y) = M - (y < dimY - 1 ? gapOpen : 0);
				if (Iy(x, y) < Iy1 - gapExtend)
					Iy(x, y) = Iy1 - gapExtend;
			}
			else
			{
				tb(x, y) = -1;
				B(x, y) = s = Iy1;

				Ix(x, y) = M - (x < dimX - 1 ? gapOpen : 0);
				if (Ix(x, y) < Ix1 - gapExtend)
					Ix(x, y) = Ix1 - gapExtend;
				Iy(x, y) = Iy1 - gapExtend;
			}

			if (/*(x == dimX - 1 or y == dimY - 1) and */high < s)
			{
				high = s;
				highX = x;
				highY = y;
			}
		}
	}
	
	const int kFlagSeqNr = numeric_limits<int>::min();
	
	// reset positions of seqres
	for (auto& sr: rx)
	{
		sr.m_seq_num = kFlagSeqNr;
		sr.m_icode = ' ';
	}
	
	// assign numbers
	x = highX;
	y = highY;
	
	while (x >= 0 and y >= 0)
	{
		switch (tb(x, y))
		{
			case -1:
				throw runtime_error("A residue found in the ATOM records (" + ry[y].m_mon_id + 
					" @ " + string{m_dbref.chainID} + ":" + to_string(ry[y].m_seq_num) +
					") was not found in the SEQRES records");
				--y;
				break;
			
			case 1:
				if (VERBOSE > 3)
					cerr << "Missing residue in ATOM records: " << rx[x].m_mon_id << " at " << rx[x].m_seq_num << endl;

				--x;
				break;
			
			case 0:
				if (VERBOSE > 3 and rx[x].m_mon_id != ry[y].m_mon_id)
					cerr << "Warning, unaligned residues at " << x << "/" << y << "(" << rx[x].m_mon_id << '/' << ry[y].m_mon_id << ')' << endl;
				else if (VERBOSE > 4)
					cerr << rx[x].m_mon_id << " -> " << ry[y].m_mon_id << " (" << ry[y].m_seq_num << ')' << endl;

				rx[x].m_seq_num = ry[y].m_seq_num;
				rx[x].m_icode = ry[y].m_icode;

				--x;
				--y;
		}
	}
	
	// assign numbers to the residues that don't have them yet
	stack<int> unnumbered;
	for (x = 0; x < dimX; ++x)
	{
		if (rx[x].m_seq_num == kFlagSeqNr)
		{
			if (x > 0 and rx[x - 1].m_seq_num != kFlagSeqNr)
				rx[x].m_seq_num = rx[x - 1].m_seq_num + 1;
			else
				unnumbered.push(x);
		}
	}
	
	while (unnumbered.empty() == false)
	{
		x = unnumbered.top();
		if (x >= dimX - 1)
			throw runtime_error("Could not assign sequence numbers");
		rx[x].m_seq_num = rx[x + 1].m_seq_num - 1;
		unnumbered.pop();
	}
}

bool PDBFileParser::PDBChain::SameSequence(const PDBChain& rhs) const
{
	bool result = m_seqres.size() == rhs.m_seqres.size();
	
	for (size_t i = 0; result and i < m_seqres.size(); ++i)
		result = m_seqres[i].m_mon_id == rhs.m_seqres[i].m_mon_id;
	
	return result;
}

// --------------------------------------------------------------------

void ReadPDBFile(istream& pdbFile, cif::file& cifFile)
{
	PDBFileParser p;

	cifFile.load_dictionary("mmcif_pdbx");

	p.Parse(pdbFile, cifFile);
	
	cifFile.validate();
}
