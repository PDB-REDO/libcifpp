#include "cif++/Config.h"

#include <map>
#include <set>
#include <boost/regex.hpp>

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <clipper/core/spacegroup.h>

#include "cif++/PDB2Cif.h"
#include "cif++/AtomType.h"
#include "cif++/Compound.h"
#include "cif++/PDB2CifRemark3.h"
#include "cif++/CifUtils.h"

using namespace std;
namespace ba = boost::algorithm;

using cif::Datablock;
using cif::Category;
using cif::Row;
using cif::Key;
using cif::iequals;
using mmcif::CompoundFactory;

// --------------------------------------------------------------------
// attempt to come up with better error handling

namespace error
{
	enum pdbErrors
	{
		residueNotFound	= 1000,
		invalidDate
	};
	
	namespace detail
	{
		class pdbCategory : public boost::system::error_category
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
					case residueNotFound:
						return "Residue not found";
					
					case invalidDate:
						return "Invalid date";
					
					default:
						return "Error in PDB format";
				}
			}
		};
	}

	boost::system::error_category& pdbCategory()
	{
		static detail::pdbCategory impl;
		return impl;
	}
	
	inline boost::system::error_code make_error_code(pdbErrors e)
	{
		return boost::system::error_code(static_cast<int>(e), pdbCategory());
	}
}

namespace boost {
namespace system {

template<> struct is_error_code_enum<error::pdbErrors>
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
	"MASTER", "END   ",
	
	// bah...
	"LINKR "
};

bool isWater(const string& resname)
{
	return resname == "HOH" or resname == "H2O" or resname == "OH2" or resname == "WAT" or resname == "DOD";
}

// --------------------------------------------------------------------
//	Unfortunately, parsing a PDB file requires several passes over the
//	data. Therefore we first obtain all records where a record has the
//	value flattened out for continuation.

PDBRecord::PDBRecord(uint32 lineNr, const string& name, const string& value)
	: mNext(nullptr), mLineNr(lineNr), mVlen(value.length())
{
	assert(name.length() <= 10);
	
	strcpy(mName, name.c_str());
	strcpy(mValue, value.c_str());
}

PDBRecord::~PDBRecord()
{
}

void* PDBRecord::operator new(size_t size, size_t vLen)
{
	return malloc(size + vLen + 1);
}

void PDBRecord::operator delete(void* p)
{
	free(p);
}

bool PDBRecord::is(const char* name) const
{
	return iequals(mName, name);
}

char PDBRecord::vC(size_t column)
{
	char result = ' ';
	if (column - 7 < mVlen)
		result = mValue[column - 7];
	return result;
}

string PDBRecord::vS(size_t columnFirst, size_t columnLast)
{
	string result;
	
	if (columnLast > mVlen + 6)
		columnLast = mVlen + 6;
	
	if (columnFirst < mVlen + 7)
	{
		result = string{mValue + columnFirst - 7, mValue + columnLast - 7 + 1};
		ba::trim(result);
	}

	return result;
}

int PDBRecord::vI(int columnFirst, int columnLast)
{
	int result = 0;

	const char* e = mValue + mVlen;
	if (e > mValue + columnLast - 7 + 1)
		e = mValue + columnLast - 7 + 1;
	
	enum { start, digit, tail } state = start;
	bool negate = false;

	try
	{
		for (const char* p = mValue + columnFirst - 7; p < e; ++p)
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
	}
	catch (const exception& ex)
	{
		cerr << "Trying to parse '" << string(mValue + columnFirst - 7, mValue + columnLast - 7) << '\'' << endl;
		throw;
	}
	
	if (negate)
		result = -result;
	
	return result;
}

string PDBRecord::vF(size_t columnFirst, size_t columnLast)
{
	// for now... TODO: check format?
	return vS(columnFirst, columnLast);
}

// --------------------------------------------------------------------

class SpecificationListParser
{
  public:
	SpecificationListParser(const string& text)
		: mText(text), mP(mText.begin()) {}
	
	tuple<string,string> GetNextSpecification();
	
  private:
	string				mText;
	string::iterator	mP;
};

tuple<string,string> SpecificationListParser::GetNextSpecification()
{
	string id, value;
	
	string::iterator start = mP, backup;

	enum { eStart, eID, eColon, eValue, eNL, eNL_ID, eSemiColon, eError, eDone } state = eStart;
	
	while (mP != mText.end() and state != eDone)
	{
		char ch = *mP++;
		
		switch (state)
		{
			case eStart:
				if (isalnum(ch) or ch == '_')
				{
					id = { ch };
					value.clear();
					state = eID;
					start = mP;
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
					backup = mP;
					state = eNL;
				}
				else if (ch == ';')
				{
					backup = mP;
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
					value.insert(value.end(), backup, mP);
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
					mP = backup;
					state = eDone;
				}
				else if (ch == ';')
					state = eSemiColon;
				else if (not (isalnum(ch) or ch == '_'))
				{
					value.insert(value.end(), backup, mP);
					state = eValue;
				}
				break;
			
			case eError:
				if (ch == ';')
				{
					if (VERBOSE)
						cerr << "Skipping invalid header line: '" << string(start, mP) << endl;
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
		: mData(nullptr), mRec(nullptr)
	{
	}
	
	~PDBFileParser()
	{
		PDBRecord* r = mData;
		while (r != nullptr)
		{
			PDBRecord* d = r;
			r = d->mNext;
			delete d;
		}	
	}
	
	void Parse(istream& is, cif::File& result);

  private:

	// ----------------------------------------------------------------
	
	struct DBREF
	{
		string	PDBIDCode;
		char	chainID;
		int		seqBegin;
		char	insertBegin = ' ';
		int		seqEnd;
		char	insertEnd = ' ';
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
		string		asymId;
		set<int>	atoms;
	};
	
	struct UNOBS
	{
		int modelNr;
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
		int						mMolId;
		string					mTitle;
		set<char>				mChains;
		map<string,string>		mInfo;
		map<string,string>		mSource;
		int						mCount = 0;
	};
	
	struct PDBSeqRes
	{
		string					mMonId;
		int						mSeqNum;
		char					mIcode;

		int						mDbSeqNum;
		bool					mSeen = false;
		set<string>				mAlts;
		
		bool operator==(const PDBSeqRes& rhs) const
		{
			return mSeqNum == rhs.mSeqNum and mMonId == rhs.mMonId and mIcode == rhs.mIcode;
		}
	};
	
	struct PDBChain
	{
		PDBChain(const string& structureId, char chainId, int molId)
			: mDbref{structureId, chainId}, mWaters(0), mTerIndex(0), mMolId(molId), mNextSeqNum(1), mNextDbSeqNum(1)
		{
		}
		
		DBREF					mDbref;
		vector<PDBSeqRes>		mSeqres, mHet;
		int						mWaters;
		int						mTerIndex;
		
		int						mMolId;

		// scratch values for reading SEQRES records
		int						mNextSeqNum;
		int						mNextDbSeqNum;
		
		// scratch value for aligning
		struct AtomRes
		{
			string				mMonId;
			int					mSeqNum;
			char				mIcode;
			
			bool operator==(const AtomRes& rhs) const		{ return mSeqNum == rhs.mSeqNum and mIcode == rhs.mIcode; }
			bool operator!=(const AtomRes& rhs) const		{ return mSeqNum != rhs.mSeqNum or mIcode != rhs.mIcode; }
		};
		vector<AtomRes>			mResiduesSeen;
		
		int AlignResToSeqRes();
		bool SameSequence(const PDBChain& rhs) const;
	};

	// ----------------------------------------------------------------

	PDBCompound& GetOrCreateCompound(int molId)
	{
		auto i = find_if(mCompounds.begin(), mCompounds.end(), [molId](PDBCompound& comp) -> bool { return comp.mMolId == molId; });
		if (i == mCompounds.end())
		{
			mCompounds.push_back(PDBCompound{ molId });
			
			mMolID2EntityID[molId] = to_string(mNextEntityNr++);
			
			i = prev(mCompounds.end());
		}
		
		return *i;
	}

	// locate the PDBChain record for a chain ID, or create it with dummy data if missing
	PDBChain& GetChainForID(char chainID, int numRes = 0)
	{
		auto i = find_if(mChains.begin(), mChains.end(), [chainID](PDBChain& ch) -> bool { return ch.mDbref.chainID == chainID; });
		
		if (i == mChains.end())
		{
			// locate the compound for this chain, if any (does that happen?)
			int molId = 0;
			for (auto& cmp: mCompounds)
			{
				if (cmp.mChains.count(chainID) > 0)
				{
					molId = cmp.mMolId;
					break;
				}
			}
			
			mChains.emplace_back(mStructureId, chainID, molId);

			i = prev(mChains.end());
		}
		
		return *i;
	};
	
	void InsertChemComp(const string& chemComp)
	{
		if (find(mChemComp.begin(), mChemComp.end(), chemComp) == mChemComp.end())
			mChemComp.push_back(chemComp);
	}

	void InsertAtomType(const string& atomType)
	{
		if (find(mAtomTypes.begin(), mAtomTypes.end(), atomType) == mAtomTypes.end())
			mAtomTypes.push_back(atomType);
	}

	// ----------------------------------------------------------------

	template<typename Predicate>
	PDBRecord* FindRecord(Predicate&& pred)
	{
		PDBRecord* result;
		
		for (result = mData; result != nullptr; result = result->mNext)
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

	char vC(size_t column) const
	{
		return mRec->vC(column);
	}
	
	string vS(size_t columnFirst, size_t columnLast = numeric_limits<size_t>::max()) const
	{
		return mRec->vS(columnFirst, columnLast);
	}

	string vF(size_t columnFirst, size_t columnLast) const
	{
		return mRec->vF(columnFirst, columnLast);
	}

	int vI(int columnFirst, int columnLast) const
	{
		return mRec->vI(columnFirst, columnLast);
	}

	// ----------------------------------------------------------------
	
	// Map a PDB residue location to a seqnum in a struct_asym 
	tuple<string,int,bool> MapResidue(char chainID, int resSeq, char iCode) const
	{
		auto key = make_tuple(chainID, resSeq, iCode);
		
		try
		{
			return mChainSeq2AsymSeq.at(key);
		}
		catch (const exception& ex)
		{
			throw_with_nested(runtime_error(string("Residue ") + chainID + to_string(resSeq) + iCode + " could not be mapped"));
		}
	}

	tuple<string,int,bool> MapResidue(char chainID, int resSeq, char iCode, boost::system::error_code& ec) const
	{
		auto key = make_tuple(chainID, resSeq, iCode);
		
		tuple<string,int,bool> result;
		
		if (not mChainSeq2AsymSeq.count(key))
		{
			ec = error::make_error_code(error::pdbErrors::residueNotFound);
			if (VERBOSE)
				cerr << "Residue " << chainID << resSeq << iCode << " could not be mapped" << endl;
		}
		else
			result = mChainSeq2AsymSeq.at(key);
		
		return result;
	}
	
	tuple<string,int,bool> MapResidue(char chainID, int resSeq, char iCode, const string& resName);

	// ----------------------------------------------------------------
	
	void PreParseInput(istream& is);

	void GetNextRecord();
	void Match(const string& expected, bool throwIfMissing);

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
	void ParseCoordinate(int modelNr);
	void ParseConnectivty();
	void ParseBookkeeping();

	// ----------------------------------------------------------------

	Category* getCategory(string name)
	{
		cif::Datablock::iterator i;
		std::tie(i, ignore) = mDatablock->emplace(name);
		return &*i;
	}
	
	vector<string> SplitCSV(const string& value);

	string pdb2cifDate(string s, boost::system::error_code& ec)
	{
		smatch m;
		const regex
			rx1(R"((\d{2})-(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)-(\d{2}))"),
			rx2(R"((JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)-(\d{2}))");

		try
		{
			if (regex_match(s, m, rx1))
			{
				using namespace boost::gregorian;
				
				int day = stoi(m[1].str());
				auto mi = kMonths.find(m[2].str());
				if (mi == kMonths.end())
					throw runtime_error("Invalid month: '" + m[2].str() + '\'');
				int month = mi->second;
				int year = 1900 + stoi(m[3].str());
				if (year < 1950)
					year += 100;
			
				date dateOriginal(year, month, day);
				
				s = to_iso_extended_string(dateOriginal);
			}
			else if (regex_match(s, m, rx2))
			{
				auto mi = kMonths.find(m[1].str());
				if (mi == kMonths.end())
					throw runtime_error("Invalid month: '" + m[1].str() + '\'');
				int month = mi->second;
				int year = 1900 + stoi(m[2].str());
				if (year < 1950)
					year += 100;
			
				s = (boost::format("%04d-%02d") % year % month).str();
			}
			else
				ec = error::make_error_code(error::pdbErrors::invalidDate);
		}
		catch (const exception& ex)
		{
			if (VERBOSE)
				cerr << ex.what() << endl;
			ec = error::make_error_code(error::pdbErrors::invalidDate);
		}
		
		return s;
	}
	
	string pdb2cifDate(string s)
	{
		boost::system::error_code ec;		
		return pdb2cifDate(s, ec);
	}
	
	string pdb2cifAuth(string author)
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

	string pdb2cifSymmetry(string s)
	{
		static const regex sgRx(R"((\d{1,3})(\d{3}))");
		
		if (not s.empty())
		{
			smatch m;
			if (not regex_match(s, m, sgRx))
				throw runtime_error("invalid symmetry value '" + s + '\'');
	
			s = m[1].str() + "_" + m[2].str();
		}
		
		return s;
	}
	
	string pdb2cifCharge(string c)
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
	
	string cifIdForInt(int nr) const
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
	
	vector<char> altLocsForAtom(char chainID, int seqNum, char iCode, string atomName);
	void MapChainID2AsymIDS(char chainID, vector<string>& asymIds);

	// ----------------------------------------------------------------
	
	PDBRecord*	mData;
	PDBRecord*	mRec;
	cif::Datablock*	mDatablock = nullptr;

	string		mStructureId;
	string		mModelTypeDetails;
	string		mOriginalDate;
	string		mExpMethod = "X-RAY DIFFRACTION";
	int			mCitationAuthorNr = 1, mCitationEditorNr = 1;
	int			mNextMolId = 1, mNextEntityNr = 1;
	int			mNextSoftwareOrd = 1;

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

	vector<SEQADV>						mSeqadvs;
	
	list<PDBCompound>					mCompounds;
	list<PDBChain>						mChains;
	vector<HET> 						mHets;
	map<string,string>					mHetnams;
	map<string,string>					mHetsyns;
	map<string,string>					mFormuls;
	string								mWaterHetId;
	vector<string>						mChemComp, mAtomTypes;
	
	map<string,string>					mRemark200;
	string								mRefinementSoftware;
	int									mAtomId = 0;
	int									mPdbxDifOrdinal = 0;

	vector<UNOBS>						mUnobs;

	// various maps between numbering schemes
	map<tuple<char,int,char>,tuple<string,int,bool>>	mChainSeq2AsymSeq;

	map<int,string>							mMolID2EntityID;
	map<string,string>						mHet2EntityID;
	map<string,string>						mAsymID2EntityID;
	map<string,string>						mMod2parent;
};

// --------------------------------------------------------------------

vector<char> PDBFileParser::altLocsForAtom(char inChainID, int inResSeq, char inICode, string inAtomName)
{
	// well, maybe this could be optimized...
	set<char> result;
	
	for (auto r = mData; r != nullptr; r = r->mNext)
	{
		if (r->is("ATOM  ") or r->is("HETATM"))		//	 1 -  6        Record name   "ATOM  "
		{											//	 ...
			string name			= r->vS(13, 16);	//	13 - 16        Atom          name         Atom name.
			char altLoc			= r->vC(17);		//	17             Character     altLoc       Alternate location indicator.
			char chainID		= r->vC(22);		//	22             Character     chainID      Chain identifier.
			int resSeq			= r->vI(23, 26);	//	23 - 26        Integer       resSeq       Residue sequence number.
			char iCode			= r->vC(27);		//	27             AChar         iCode        Code for insertion of residues.

			if (chainID == inChainID and resSeq == inResSeq and iCode == inICode and name == inAtomName and altLoc != ' ')
				result.insert(altLoc);
		}
	}
	
	return { result.begin(), result.end() };
}

void PDBFileParser::MapChainID2AsymIDS(char chainID, vector<string>& asymIds)
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
	
	set<string,l> asym(asymIds.begin(), asymIds.end());

	for (auto& m: mChainSeq2AsymSeq)
	{
		if (get<0>(m.first) == chainID)
			asym.insert(get<0>(m.second));
	}
	
	asymIds.assign(asym.begin(), asym.end());
}

// --------------------------------------------------------------------

void PDBFileParser::PreParseInput(istream& is)
{
	string lookahead;
	uint32 lineNr = 1;
	getline(is, lookahead);

//	if (ba::starts_with(lookahead, "HEADER") == false)
//		throw runtime_error("This does not look like a PDB file, should start with a HEADER line");

	auto contNr = [&lookahead](int offset, int len) -> int
	{
		string cs = lookahead.substr(offset, len);
		ba::trim(cs);
		int result;
		
		try
		{
			result = cs.empty() ? 0 : stoi(cs);
		}
		catch (...)
		{
			throw runtime_error("Continuation string '" + cs + "' is not valid");
		}
		return result;
	};

	PDBRecord* last = nullptr;
	set<string> dropped;

	for (;;)
	{
		if (lookahead.empty())
		{
			if (is.eof())
				break;

			cerr << "Line number " << lineNr << " is empty!" << endl;

			getline(is, lookahead);
			++lineNr;

			continue;		
		}

		string type = lookahead.substr(0, 6);
		string value;
		if (lookahead.length() > 6)
			value = ba::trim_right_copy(lookahead.substr(6));

		uint32 curLineNr = lineNr;
		getline(is, lookahead);
		++lineNr;
		
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
			while (lookahead.substr(0, 6) == type and contNr(7, 3) == n)
			{
				value += ba::trim_right_copy(lookahead.substr(10));
				getline(is, lookahead);
				++lineNr;
				++n;
			}
		}
		else if (type == "COMPND")
		{
			int n = 2;
			value += '\n';
			while (lookahead.substr(0, 6) == type and contNr(7, 3) == n)
			{
				value += ba::trim_right_copy(lookahead.substr(10));
				value += '\n';
				getline(is, lookahead);
				++lineNr;
				++n;
			}
		}
		else if (type == "REVDAT")
		{
			int revNr = stoi(value.substr(1, 3));
			int n = 2;
			while (lookahead.substr(0, 6) == type and
				stoi(lookahead.substr(7, 3)) == revNr and
				contNr(10, 2) == n)
			{
				value += lookahead.substr(38);
				getline(is, lookahead);
				++lineNr;
				++n;
			}
		}
		else if (type == "CAVEAT")
		{
			int n = 2;
			while (lookahead.substr(0, 6) == type and contNr(7, 3) == n)
			{
				value += ba::trim_right_copy(lookahead.substr(13));
				getline(is, lookahead);
				++lineNr;
				++n;
			}
		}
		else if (type == "OBSLTE")
		{
			while (lookahead.substr(0, 6) == type)
			{
				value += lookahead.substr(31);
				getline(is, lookahead);
				++lineNr;
			}
		}
		else if (type == "SOURCE")
		{
			value += '\n';
			int n = 2;
			while (lookahead.substr(0, 6) == type and contNr(7, 3) == n)
			{
				value += ba::trim_copy(lookahead.substr(10));
				value += '\n';
				getline(is, lookahead);
				++lineNr;
				++n;
			}
		}
		else if (type == "FORMUL")
		{
			try
			{
				int compNr;
				try
				{
					compNr = stoi(value.substr(1, 3));
				}
				catch (const exception& ex)
				{
					throw_with_nested(runtime_error("Invalid component number '" + value.substr(1, 3) + '\''));
				}
	
				int n = 2;
				try
				{
					while (lookahead.substr(0, 6) == type and
						stoi(lookahead.substr(7, 3)) == compNr and
						contNr(16, 2) == n)
					{
						value += lookahead.substr(19);
						getline(is, lookahead);
						++lineNr;
						++n;
					}
				}
				catch (const invalid_argument& ex)
				{
					throw_with_nested(runtime_error("Invalid component number '" + lookahead.substr(7, 3) + '\''));
				}
			}
			catch (const exception& ex)
			{
				cerr << "Error parsing FORMUL at line " << lineNr << endl;
				throw;
			}
		}
		else if (type == "HETNAM" or
				 type == "HETSYN")
		{
			int n = 2;
			while (lookahead.substr(0, 6) == type and contNr(8, 2) == n)
			{
				value += lookahead.substr(16);
				getline(is, lookahead);
				++lineNr;
				++n;
			}
		}
		else if (type == "SITE  ")
		{
			string siteName = value.substr(5, 3);
			ba::trim_right(value);
			size_t n = value.length() - 12;
			value += string(11 - (n % 11), ' ');
			
			while (lookahead.substr(0, 6) == type and lookahead.substr(11, 3) == siteName)
			{
				string s = lookahead.substr(18);
				ba::trim_right(s);
				s += string(11 - (s.length() % 11), ' ');
				value += s;
				
// TODO: improve this... either use numRes or don't lump together all text
//				value += " " + ba::trim_right_copy();
				getline(is, lookahead);
				++lineNr;
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
					
					if (iequals(v, "NONE") or iequals(v, "N/A") or iequals(v, "NAN"))
						mRemark200[k] = ".";
					else if (not iequals(v, "NULL"))
						mRemark200[k] = v;
				}
			}
		}
		
		PDBRecord* cur = new(value.length()) PDBRecord(curLineNr, type, value);
		
		if (last == nullptr)
			last = mData = cur;
		else
			last->mNext = cur;

		last = cur;

		if (type == "END   ")
			break;
	}
	
	if (not dropped.empty())
	{
		cerr << "Dropped unsupported records: " << ba::join(dropped, ", ") << endl;
	}
	
	mRec = mData;
}

void PDBFileParser::GetNextRecord()
{
	if (mRec != nullptr)
		mRec = mRec->mNext;

	if (mRec == nullptr)
	{
		static PDBRecord* end = new(0)PDBRecord({ 0, "END   ", ""});
		mRec = end;
	}
}

void PDBFileParser::Match(const string& expected, bool throwIfMissing)
{
	if (mRec->mName != expected)
	{
		if (throwIfMissing)
			throw runtime_error("Expected record " + expected + " but found " + mRec->mName);
		if (VERBOSE)
			cerr << "Expected record " + expected + " but found " + mRec->mName << endl;
	}
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

	Match("HEADER", false);
	
	string keywords;
	
	if (mRec->is("HEADER"))
	{
		mStructureId		= vS(63, 66);
		keywords			= vS(11, 50);
		mOriginalDate		= pdb2cifDate(vS(51, 59));
	
		ba::trim(keywords);
		
		GetNextRecord();
	}

	ba::trim(mStructureId);
	if (mStructureId.empty())
		mStructureId = "nohd";

	mDatablock = new cif::Datablock(mStructureId);
		
	auto cat = getCategory("entry");
//	cat->addColumn("id");
	cat->emplace({ {"id", mStructureId} });
	
	// OBSLTE
	if (mRec->is("OBSLTE"))
	{
		//	 1 -  6       Record name   "OBSLTE"
		//	 9 - 10       Continuation  continuation  Allows concatenation of multiple records
		//	12 - 20       Date          repDate       Date that this datablock was replaced.
		//	22 - 25       IDcode        idCode        ID code of this datablock.
		//	32 - 35       IDcode        rIdCode       ID code of datablock that replaced this one.
		//	37 - 40       ...


		string old		= vS(22, 25);
		string date		= pdb2cifDate(vS(12, 20));
		cat = getCategory("pdbx_database_PDB_obs");
		
		string value = mRec->vS(32);
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
	Match("TITLE ", false);
	string title;
	if (mRec->is("TITLE "))	//	 1 -  6       Record name    "TITLE "
	{							//	 9 - 10       Continuation   continuation  Allows concatenation of multiple records.
		title = vS(11);		//	11 - 80       String         title         Title of the  experiment.
		GetNextRecord();
	}
	
	// SPLIT
	if (mRec->is("SPLIT "))
	{
		//	 1 -  6        Record  name  "SPLIT "
		//	 9 - 10        Continuation  continuation  Allows concatenation of multiple records.
		//	12 - 15        IDcode        idCode        ID code of related datablock.
		
		throw runtime_error("SPLIT PDB files are not supported");
	}
	
	// CAVEAT
	int caveatID = 1;
	while (mRec->is("CAVEAT"))							//	  1 - 6       Record name   "CAVEAT"
	{
		getCategory("database_PDB_caveat")->emplace({
			{ "id", caveatID++ },
			{ "text", string{mRec->vS(20) } }    		//	20 - 79       String        comment        Free text giving the reason for the  CAVEAT.
		});
		
		GetNextRecord();
	}
	
	// COMPND
	Match("COMPND", false);
	//	 1 -  6       Record name     "COMPND"   
	//	 8 - 10       Continuation    continuation  Allows concatenation of multiple records.
	//	11 - 80       Specification   compound      Description of the molecular components.
	//	              list   

	string value{mRec->vS(11)};
	if (value.find(':') == string::npos)
	{
			// special case for dumb, stripped files		
		auto& comp = GetOrCreateCompound(1);
		comp.mInfo["MOLECULE"] = value;
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
		
			if (not iequals(key, "MOL_ID") and mCompounds.empty())
			{
				cerr << "Ignoring invalid COMPND record" << endl;
				break;
			}
			
			if (key == "MOL_ID")
			{
				auto& comp = GetOrCreateCompound(stoi(value));
				comp.mTitle = title;
			}
			else if (key == "CHAIN")
			{
				vector<string> chains;
				
				ba::split(chains, value, ba::is_any_of(","));
				for (auto& c: chains)
				{
					ba::trim(c);
					mCompounds.back().mChains.insert(c[0]);
				}
			}
			else
				mCompounds.back().mInfo[key] = value;
		}	
	}

	if (mRec->is("COMPND"))
		GetNextRecord();

	// SOURCE
	Match("SOURCE", false);

	if (mRec->is("SOURCE"))
	{
		//	 1 -  6      Record name    "SOURCE"       
		//	 8 - 10      Continuation   continuation   Allows concatenation of multiple records.
		//	11 - 79      Specification  srcName        Identifies the source of the
		//	             List                          macromolecule in a  token: value format.
		
		map<string,string>* source = nullptr;
	
//		value = { mRec->vS(11) };
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
		SpecificationListParser p(vS(11));
		
		for (;;)
		{
			string key, value;
			std::tie(key, value) = p.GetNextSpecification();

			if (key.empty())
				break;
			
			if (key == "MOL_ID")
			{
				for (auto& c: mCompounds)
				{
					if (c.mMolId == stoi(value))
					{
						source = &c.mSource;
						break;
					}
				}
				
				continue;
			}
			
			if (source == nullptr)
				throw runtime_error("At line " + to_string(mRec->mLineNr) + ": missing MOL_ID in SOURCE");
			
			(*source)[key] = value;
		}
		
		GetNextRecord();
	}

	// KEYWDS
	Match("KEYWDS", false);
	string pdbxKeywords;

	if (mRec->is("KEYWDS"))		//	 1 -  6       Record name    "KEYWDS"  
	{								//	 9 - 10       Continuation   continuation  Allows concatenation of records if necessary.
		pdbxKeywords = vS(11);	//	11 - 79       List           keywds        Comma-separated list of keywords relevant
									//	                                           to the datablock.       
		GetNextRecord();
	}
	
	if (not (keywords.empty() and pdbxKeywords.empty()))
	{
		getCategory("struct_keywords")->emplace({
			{ "entry_id",  mStructureId },
			{ "pdbx_keywords", keywords },
			{ "text", pdbxKeywords }
		});
	}

	// EXPDTA
	Match("EXPDTA", false);
	if (mRec->is("EXPDTA"))
	{
		mExpMethod = vS(11);

		cat = getCategory("exptl");
		
		vector<string> crystals;
		ba::split(crystals, mRemark200["NUMBER OF CRYSTALS USED"], ba::is_any_of("; "));
		if (crystals.empty())
			crystals.push_back("");
		auto ci = crystals.begin();
		
		for (auto si = ba::make_split_iterator(mExpMethod, ba::token_finder(ba::is_any_of(";"), ba::token_compress_on)); not si.eof(); ++si, ++ci)
		{
			string expMethod(si->begin(), si->end());
			ba::trim(expMethod);
			
			if (expMethod.empty())
				continue;
			
			cat->emplace({
				{ "entry_id", mStructureId },
				{ "method", expMethod },
				{ "crystals_number", ci != crystals.end() ? *ci : "" }
			});
		}
		
		GetNextRecord();
	}

	// NUMMDL
	if (mRec->is("NUMMDL"))
	{
		if (VERBOSE)	
			cerr << "skipping unimplemented NUMMDL record" << endl;
		GetNextRecord();
	}
	
	// MDLTYP
	if (mRec->is("MDLTYP"))
	{
		mModelTypeDetails = vS(11);
		GetNextRecord();
	}

	// AUTHOR
	Match("AUTHOR", false);
	if (mRec->is("AUTHOR"))
	{
		int n = 1;
		cat = getCategory("audit_author");
		
		value = { mRec->vS(11) };
		for (auto si = ba::make_split_iterator(value, ba::token_finder(ba::is_any_of(","), ba::token_compress_on)); not si.eof(); ++si)
		{
			string author(si->begin(), si->end());
			
			cat->emplace({
				{ "name", pdb2cifAuth(author) },
				{ "pdbx_ordinal", n }
			});
			++n;
		}
	
		GetNextRecord();
	}

	// REVDAT
	bool firstRevDat = true;
	struct RevDat {
		int revNum;
		string date, dateOriginal, replaces;
		int modType;
		vector<string> types;

		bool operator<(const RevDat& rhs) const { return revNum < rhs.revNum; }
	};
	vector<RevDat> revdats;
	
	while (mRec->is("REVDAT"))
	{
													//	 1 -  6       Record name    "REVDAT"                                             
		int revNum = vI(8, 10);						//	 8 - 10       Integer        modNum        Modification number.                   
													//	11 - 12       Continuation   continuation  Allows concatenation of multiple records.
		string date = pdb2cifDate(vS(14, 22));		//	14 - 22       Date           modDate       Date of modification (or release  for   
													//	                                           new entries)  in DD-MMM-YY format. This is
													//	                                           not repeated on continued lines.
		string modId = vS(24, 27);					//	24 - 27       IDCode         modId         ID code of this datablock. This is not repeated on 
													//	                                           continuation lines.    
		int modType = vI(32, 32);					//	32            Integer        modType       An integer identifying the type of    
													//	                                           modification. For all  revisions, the
													//	                                           modification type is listed as 1 
		string value = vS(40);						//	40 - 45       LString(6)     record        Modification detail. 
													//	47 - 52       LString(6)     record        Modification detail. 
													//	54 - 59       LString(6)     record        Modification detail. 
													//	61 - 66       LString(6)     record        Modification detail.
		
		revdats.push_back({revNum, date, modType == 0 ? mOriginalDate : "", modId, modType});
		
		ba::split(revdats.back().types, value, ba::is_any_of(" "));

		if (firstRevDat)
		{
			cat = getCategory("database_2");
			cat->emplace({
				{ "database_id", "PDB" },
				{ "database_code", modId }
			});
		}

		GetNextRecord();
		firstRevDat = false;
	}
	
/*
	This is internal stuff for PDB, don't write it ???
*/
	sort(revdats.begin(), revdats.end());
	for (auto& revdat: revdats)
	{
		getCategory("database_PDB_rev")->emplace({
			{ "num",			revdat.revNum },
			{ "date",			revdat.date },
			{ "date_original",	revdat.dateOriginal },
			{ "replaces", 		revdat.replaces },
			{ "mod_type",		revdat.modType }
		});
		
		for (auto& type: revdat.types)
		{
			if (type.empty())
				continue;
			getCategory("database_PDB_rev_record")->emplace({
				{ "rev_num",	revdat.revNum  },
				{ "type",		type }
			});
		}
	}
//*/

	// SPRSDE
	if (mRec->is("SPRSDE"))
	{
		if (VERBOSE)	
			cerr << "skipping unimplemented SPRSDE record" << endl;
		GetNextRecord();
	}
	
	// JRNL
	if (mRec->is("JRNL  "))
		ParseCitation("primary");
}

void PDBFileParser::ParseCitation(const string& id)
{
	const char* rec = mRec->mName;
	
	string auth, titl, edit, publ, refn, pmid, doi;
	string pubname, volume, astm, country, issn, csd;
	string pageFirst;
	int year = 0;

	auto extend = [](string& s, const string& p)
	{
		if (not s.empty())	
			s += ' ';
		s += ba::trim_copy(p);
	};

	while (mRec->is(rec) and (id == "primary" or vC(12) == ' '))
	{
		string k = vS(13, 16);
		if (k == "AUTH")				extend(auth, vS(20, 79));
		else if (k == "TITL")			extend(titl, vS(20, 79));
		else if (k == "EDIT")			extend(edit, vS(20, 79));
		else if (k == "REF")
		{
			if (pubname.empty())
			{
				extend(pubname, vS(20, 47));
				if (vS(50, 51) == "V.")
					volume = ba::trim_copy(vS(52, 55));
				pageFirst = vS(57, 61);
				year = vI(63, 66);
			}
			else
				extend(pubname, vS(20, 47));
		}
		else if (k == "PUBL")			extend(publ, vS(20, 70));
		else if (k == "REFN")
		{
			if (vS(20, 23) == "ASTN")
				astm = vS(25, 30);
			country = vS(33, 34);
			if (vS(36, 39) == "ISSN")
				issn = vS(41, 65);
		}
		else if (k == "PMID")			pmid = vS(20, 79);
		else if (k == "DOI")			doi = vS(20, 79);
		
		GetNextRecord();
	}
	
	auto cat = getCategory("citation");
	cat->emplace({
		{ "id", id },
		{ "title", titl },
		{ "journal_abbrev", pubname },
		{ "journal_volume", volume },
		{ "page_first", pageFirst },
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
		cat = getCategory("citation_author");
		for (auto si = ba::make_split_iterator(auth, ba::token_finder(ba::is_any_of(","), ba::token_compress_on)); not si.eof(); ++si)
		{
			string author(si->begin(), si->end());
			
			cat->emplace({
				{ "citation_id", id },
				{ "name", pdb2cifAuth(author) },
				{ "ordinal", mCitationAuthorNr }
			});

			++mCitationAuthorNr;
		}
	}

	if (not edit.empty())
	{
		cat = getCategory("citation_editor");
		for (auto si = ba::make_split_iterator(edit, ba::token_finder(ba::is_any_of(","), ba::token_compress_on)); not si.eof(); ++si)
		{
			string editor(si->begin(), si->end());
			
			cat->emplace({
				{ "citation_id", id },
				{ "name", pdb2cifAuth(editor) },
				{ "ordinal", mCitationEditorNr }
			});

			++mCitationEditorNr;
		}
	}
}

void PDBFileParser::ParseRemarks()
{
	string sequenceDetails, compoundDetails, sourceDetails;
	
	while (ba::starts_with(mRec->mName, "REMARK"))
	{
		int remarkNr = vI(8, 10);
		
		try
		{
			switch (remarkNr)
			{
				case 1:
					while (mRec->is("REMARK   1"))
					{
						if (mRec->mVlen > 15 and vS(12, 20) == "REFERENCE")
						{
							string id = vS(22, 70);
							GetNextRecord();
							
							ParseCitation(id);
						}
						else
							GetNextRecord();
					}
					break;
				
				case 3:
					// we skip REMARK 3 until we know the mapping
					while (mRec->is("REMARK   3"))
						GetNextRecord();
					break;
	
				case 4:
					// who cares...
					while (mRec->is("REMARK   4"))
						GetNextRecord();
					break;
	
				case 100:
				{
					const regex rx(R"(THE (\S+) ID CODE IS (\S+?)\.?\s*)");
					smatch m;
					string r = vS(12);
					
					if (regex_match(r, m, rx))
					{
						auto cat = getCategory("database_2");
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
						string r = mRec->vS(12);
	
						if (ba::starts_with(r, "REMARK: "))
						{
							mRemark200["REMARK"] = r.substr(8);
							remark = true;
						}
						else if (remark)
						{
							if (r.empty())
								remark = false;
							else
								mRemark200["REMARK"] += r;
						}
						
						GetNextRecord();
					}
					while (mRec->is("REMARK 200"));
					break;
				}
	
				case 280:
				{
					string density_Matthews, densityPercentSol, conditions;
					
					const regex rx1(R"(SOLVENT CONTENT, VS +\(%\): *(.+))"),
						rx2(R"(MATTHEWS COEFFICIENT, VM \(ANGSTROMS\*\*3/DA\): *(.+))");
					
					smatch m;
	
					do
					{
						string r = vS(12);
	
						if (conditions.empty())
						{
							if (regex_match(r, m, rx1))
								densityPercentSol = m[1].str();
							else if (regex_match(r, m, rx2))
								density_Matthews = m[1].str();
							else if (ba::starts_with(r, "CRYSTALLIZATION CONDITIONS: "))
								conditions = r.substr(28);
						}
						else
							conditions = conditions + ' ' + r;
						
						GetNextRecord();
					}
					while (mRec->is("REMARK 280"));
					
					string desc = mRemark200["REMARK"];
					if (desc == "NULL")
						desc.clear();
					
					getCategory("exptl_crystal")->emplace({
						{ "id", 1 },
						{ "density_Matthews", iequals(density_Matthews, "NULL") ? "" : density_Matthews  },
						{ "density_percent_sol", iequals(densityPercentSol, "NULL") ? "" : densityPercentSol },
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
						getCategory("exptl_crystal_grow")->emplace({
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
					for (; mRec->is("REMARK 350"); GetNextRecord())
						;
					break;
				
				case 400:
				{
					stringstream s;
					GetNextRecord();
					if (vS(12) == "COMPOUND")
						GetNextRecord();
					
					while (mRec->is("REMARK 400"))
					{
						s << vS(12) << endl;
						GetNextRecord();
					}
					
					compoundDetails = s.str();
					break;
				}
				
				case 450:
				{
					stringstream s;
					GetNextRecord();
					if (vS(12) == "SOURCE")
						GetNextRecord();
					
					while (mRec->is("REMARK 450"))
					{
						s << vS(12) << endl;
						GetNextRecord();
					}
					
					sourceDetails = s.str();
					break;
				}
				
				case 465:
				{
					bool headerSeen = false;
					regex rx(R"( *MODELS *(\d+)-(\d+))");
					int models[2] = { -1, -1 };
					
					for (; mRec->is("REMARK 465"); GetNextRecord())
					{
						if (not headerSeen)
						{
							string line = vS(12);
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
							models[0] = models[1] = vI(12, 14);
						
						string res = vS(16, 18);
						char chain = vC(20);
						int seq = vI(22, 26);
						char iCode = vC(27);
						
						for (int modelNr = models[0]; modelNr <= models[1]; ++modelNr)
							mUnobs.push_back({modelNr, res, chain, seq, iCode});
					}
					
					break;
				}
				
				case 470:
				{
					bool headerSeen = false;
					regex rx(R"( *MODELS *(\d+)-(\d+))");
					int models[2] = { -1, -1 };
					
					for (; mRec->is("REMARK 470"); GetNextRecord())
					{
						if (not headerSeen)
						{
							string line = vS(12);
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
							models[0] = models[1] = vI(12, 14);
						
						string res = vS(16, 18);
						char chain = vC(20);
						int seq = vI(21, 24);
						char iCode = vC(25);
						
						vector<string> atoms;
						string atomStr = mRec->vS(29);
						for (auto i = make_split_iterator(atomStr, ba::token_finder(ba::is_any_of(" "), ba::token_compress_on)); not i.eof(); ++i)
							atoms.push_back({ i-> begin(), i->end() });
						
						for (int modelNr = models[0]; modelNr <= models[1]; ++modelNr)
							mUnobs.push_back({modelNr, res, chain, seq, iCode, atoms});
					}
					
					break;
				}
				
				case 500:
				{
					GetNextRecord();
					
					enum State { eStart, eCCinSAU, eCC, eCBL, eCBA, eTA, eCTg, ePG, eMCP, eChC } state = eStart;
					bool headerSeen = false;
					int id = 0;
					
					for (; mRec->is("REMARK 500"); GetNextRecord())
					{
						string line = vS(12);
						
						if (line == "GEOMETRY AND STEREOCHEMISTRY")
							continue;
						
						switch (state)
						{
							case eStart:
							{
								if (line.empty() or not ba::starts_with(line, "SUBTOPIC: "))
									continue;
	
								string subtopic = line.substr(10);
	
								if (subtopic == "CLOSE CONTACTS IN SAME ASYMMETRIC UNIT")	state = eCCinSAU;
								else if (subtopic == "CLOSE CONTACTS")						state = eCC;
								else if (subtopic == "COVALENT BOND LENGTHS")				state = eCBL;
								else if (subtopic == "COVALENT BOND ANGLES")				state = eCBA;
								else if (subtopic == "TORSION ANGLES")						state = eTA;
								else if (subtopic == "NON-CIS, NON-TRANS")					state = eCTg;
								else if (subtopic == "PLANAR GROUPS")						state = ePG;
								else if (subtopic == "MAIN CHAIN PLANARITY")				state = eMCP;
								else if (subtopic == "CHIRAL CENTERS")						state = eChC;
								else if (VERBOSE)
									throw runtime_error("Unknown subtopic in REMARK 500: " + subtopic);
								
								headerSeen = false;
								id = 0;
								break;
							}
							
							case eCCinSAU:
							{
								if (not headerSeen)
									headerSeen =
										line == "ATM1  RES C  SSEQI   ATM2  RES C  SSEQI           DISTANCE";
								else if (line.empty())
									state = eStart;
								else
								{
									string atom1 = vS(13, 16);
									string res1 = vS(19, 21);
									string alt1 = vS(17, 17);
									char chain1 = vC(23);
									int seq1 = vI(25, 29);
									string iCode1 = vS(30, 30);
									
									string atom2 = vS(34, 37);
									string alt2 = vS(38, 38);
									string res2 = vS(40, 42);
									char chain2 = vC(44);
									int seq2 = vI(46, 50);
									string iCode2 = vS(51, 51);
									
									string distance = vF(63, 71);
									
									getCategory("pdbx_validate_close_contact")->emplace({
										{ "id",				to_string(++id) },
										{ "PDB_model_num",	1 },
										{ "auth_atom_id_1",	atom1 },
										{ "auth_asym_id_1", string{ chain1 } },
										{ "auth_comp_id_1", res1 },
										{ "auth_seq_id_1",	seq1 },
										{ "PDB_ins_code_1", iCode1 },
										{ "label_alt_id_1", alt1 },
										{ "auth_atom_id_2",	atom2 },
										{ "auth_asym_id_2", string { chain2 } },
										{ "auth_comp_id_2", res2 },
										{ "auth_seq_id_2",	seq2 },
										{ "PDB_ins_code_2", iCode2 },
										{ "label_alt_id_2", alt2 },
										{ "dist", distance }
									});
								}
								break;
							}
							
							case eCC:
							{
								if (not headerSeen)
									headerSeen = line == "ATM1  RES C  SSEQI   ATM2  RES C  SSEQI  SSYMOP   DISTANCE";
								else if (line.empty())
									state = eStart;
								else
								{
									string atom1 = vS(13, 16);
									string res1 = vS(19, 21);
									char chain1 = vC(23);
									int seq1 = vI(25, 29);
									
									string atom2 = vS(34, 37);
									string res2 = vS(40, 42);
									char chain2 = vC(44);
									int seq2 = vI(46, 50);
									
									string symop = pdb2cifSymmetry(vS(54, 59));
									
									string distance = vF(63, 71);
									
									getCategory("pdbx_validate_symm_contact")->emplace({
										{ "id",				to_string(++id) },
										{ "PDB_model_num",	1 },
										{ "auth_atom_id_1",	atom1 },
										{ "auth_asym_id_1", string{ chain1 } },
										{ "auth_comp_id_1", res1 },
										{ "auth_seq_id_1",	seq1 },
	//									{ "PDB_ins_code_1", "" },
	//									{ "label_alt_id_1", "" },
										{ "site_symmetry_1", "1_555" },
										{ "auth_atom_id_2",	atom2 },
										{ "auth_asym_id_2", string { chain2 } },
										{ "auth_comp_id_2", res2 },
										{ "auth_seq_id_2",	seq2 },
	//									{ "PDB_ins_code_2", "" },
	//									{ "label_alt_id_2", "" },
										{ "site_symmetry_2", symop },
										{ "dist", distance }
									});
								}
								break;
							}
							
							case eCBL:
							{
								if (not headerSeen)
								{
									if (ba::starts_with(line, "FORMAT: ") and line != "FORMAT: (10X,I3,1X,2(A3,1X,A1,I4,A1,1X,A4,3X),1X,F6.3)")
										throw runtime_error("Unexpected format in REMARK 500");
									
									headerSeen = line == "M RES CSSEQI ATM1   RES CSSEQI ATM2   DEVIATION";
								}
								else if (line.empty())
									state = eStart;
								else
								{
									int model = vI(11, 13);
									string resNam1 = vS(15, 17);
									string chainID1 { vC(19) };
									int seqNum1 = vI(20, 23);
									string iCode1 { vC(24) };
									string alt1 = vS(30, 30);
									string atm1 = vS(26, 29);
	
									string resNam2 = vS(33, 35);
									string chainID2 { vC(37) };
									int seqNum2 = vI(38, 41);
									string iCode2 { vC(42) };
									string alt2 = vS(48, 48);
									string atm2 = vS(44, 47);
									
									string deviation = vF(51, 57);
									
									if (iCode1 == " ") iCode1.clear();
									if (iCode2 == " ") iCode2.clear();
									
									getCategory("pdbx_validate_rmsd_bond")->emplace({
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
							
							case eCBA:
								if (not headerSeen)
								{
									if (ba::starts_with(line, "FORMAT: ") and line != "FORMAT: (10X,I3,1X,A3,1X,A1,I4,A1,3(1X,A4,2X),12X,F5.1)")
										throw runtime_error("Unexpected format in REMARK 500");
									
									headerSeen = line == "M RES CSSEQI ATM1   ATM2   ATM3";
								}
								else if (line.empty())
									state = eStart;
								else if (vS(64) == "DEGREES")
								{
									int model = vI(11, 13);
									string resNam = vS(15, 17);
									string chainID { vC(19) };
									int seqNum = vI(20, 23);
									string iCode { vC(24) };
									
									if (iCode == " ")
										iCode.clear();
									
									string atoms[3] = { vS(27, 30), vS(34, 37), vS(41, 44) };
									string deviation = vF(57, 62);
									if (deviation == "*****")
										deviation.clear();
									
									getCategory("pdbx_validate_rmsd_angle")->emplace({
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
							
							case eTA:
								if (not headerSeen)
								{
									if (ba::starts_with(line, "FORMAT: ") and line != "FORMAT:(10X,I3,1X,A3,1X,A1,I4,A1,4X,F7.2,3X,F7.2)")
										throw runtime_error("Unexpected format in REMARK 500");
									
									headerSeen = line == "M RES CSSEQI        PSI       PHI";
								}
								else if (line.empty())
									state = eStart;
								else
								{
									int model = vI(11, 13);
									string resNam = vS(15, 17);
									string chainID { vC(19) };
									int seqNum = vI(20, 23);
									string iCode { vC(24) };
									
									if (iCode == " ")
										iCode.clear();
									
									string psi = vF(27, 35);
									string phi = vF(37, 45);
									
									getCategory("pdbx_validate_torsion")->emplace({
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
							
							case eCTg:
								if (not headerSeen)
									headerSeen = line == "MODEL     OMEGA";
								else if (line.empty())
									state = eStart;
								else
								{
									int model = vI(45, 48);
	
									string resNam1 = vS(12, 14);
									string chainID1 { vC(16) };
									int seqNum1 = vI(17, 21);
									string iCode1 { vC(22) };
									
									if (iCode1 == " ")
										iCode1.clear();
									
									string resNam2 = vS(27, 29);
									string chainID2 { vC(31) };
									int seqNum2 = vI(32, 36);
									string iCode2 { vC(37) };
									
									if (iCode2 == " ")
										iCode2.clear();
									
									string omega = vF(54, 60);
									
									getCategory("pdbx_validate_peptide_omega")->emplace({
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
							
							case ePG:
								if (not headerSeen)
									headerSeen = line == "M RES CSSEQI        RMS     TYPE";
								else if (line.empty())
									state = eStart;
								else
								{
									int model = vI(11, 13);
									string resNam = vS(15, 17);
									string chainID { vC(19) };
									int seqNum = vI(20, 23);
									string iCode { vC(24) };
									
									if (iCode == " ")
										iCode.clear();
									
									string rmsd = vF(32, 36);
									string type = vS(41);
									
									getCategory("pdbx_validate_planes")->emplace({
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
								state = eStart;
								break;
						}
					}
					
					break;
				}
				
				case 610:
				{
					bool headerSeen = false;
					
					for (; mRec->is("REMARK 610"); GetNextRecord())
					{
						if (not headerSeen)
						{
							string line = vS(12);
							headerSeen = ba::contains(line, "RES C SSEQI");
							continue;
						}
						
						int modelNr = vI(12, 14);
						if (modelNr == 0)
							modelNr = 1;
						string res = vS(16, 18);
						char chain = vC(20);
						int seq = vI(22, 25);
						char iCode = vC(26);
						
						auto compound = mmcif::Compound::create(res);
						if (compound == nullptr)
							continue;
						
						vector<string> atoms;
						for (auto atom: compound->atoms())
						{
							if (atom.typeSymbol != mmcif::H)
								atoms.push_back(atom.id);
						}
						
						mUnobs.push_back({modelNr, res, chain, seq, iCode, { atoms }});
					}
					
					break;
				}
				
				case 800:
				{
					const regex rx1(R"(SITE_IDENTIFIER: (.+))"),
						rx2(R"(EVIDENCE_CODE: (.+))"),
						rx3(R"(SITE_DESCRIPTION: (binding site for residue ([[:alnum:]]{1,3}) ([[:alnum:]]) (\d+)|.+))", regex_constants::icase);
					
					string id, evidence, desc;
					string pdbxAuthAsymId, pdbxAuthCompId, pdbxAuthSeqId, pdbxAuthInsCode;
					smatch m;
					
					enum State { sStart, sId, sEvidence, sDesc, sDesc2 } state = sStart;
					
					auto store = [&]()
					{
						// Find the matching SITE record
						auto site = FindRecord([id](PDBRecord& r) -> bool
						{
							return r.is("SITE  ") and r.vS(12, 14) == id;
						});
						
						if (site == nullptr)
							throw runtime_error("Invalid REMARK 800, no SITE record for id " + id);
						
						// next record, store what we have	
						getCategory("struct_site")->emplace({
							{ "id", id },
							{ "details", desc },
							{ "pdbx_auth_asym_id", pdbxAuthAsymId },
							{ "pdbx_auth_comp_id", pdbxAuthCompId },
							{ "pdbx_auth_seq_id", pdbxAuthSeqId },
							{ "pdbx_num_residues", site->vI(16, 17) },
							{ "pdbx_evidence_code", evidence }
						});
					};
	
					for ( ; mRec->is("REMARK 800"); GetNextRecord())
					{
						string s = mRec->vS(12);
						if (s.empty())
							continue;
	
						switch (state)
						{
							case sStart:
	 							if (s == "SITE")
	 								state = sId;
								else if (VERBOSE)
									throw runtime_error("Invalid REMARK 800 record, expected SITE");
	 							break;
							
							case sId:
								if (regex_match(s, m, rx1))
								{
									id = m[1].str();
									state = sEvidence;
								}
								else if (VERBOSE)
									throw runtime_error("Invalid REMARK 800 record, expected SITE_IDENTIFIER");
								break;
							
							case sEvidence:
								if (regex_match(s, m, rx2))
								{
									evidence = m[1].str();
									state = sDesc;
								}
								else if (VERBOSE)
									throw runtime_error("Invalid REMARK 800 record, expected SITE_IDENTIFIER");
								break;
							
							case sDesc:
								if (regex_match(s, m, rx3))
								{
									desc = m[1].str();
									pdbxAuthCompId = m[2].str();
									pdbxAuthAsymId = m[3].str();
									pdbxAuthSeqId = m[4].str();
	
									state = sDesc2;
								}
								break;
							
							case sDesc2:
								if (regex_match(s, m, rx1))
								{
									store();
									
									id = m[1].str();
									state = sEvidence;
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
					if (vS(12) == "SEQUENCE")
						GetNextRecord();
					
					while (mRec->is("REMARK 999"))
					{
						s << vS(12) << endl;
						GetNextRecord();
					}
					
					sequenceDetails = s.str();
					break;
				}
	
				// these are skipped 
	
				case 2:
				case 290:
				case 300:
				case 620:
					GetNextRecord();
					break;
	
				default:
				{
					string skipped = mRec->mName;
					
					stringstream s;
					
					if (not mRec->vS(11).empty())
						s << mRec->vS(11) << endl;
					GetNextRecord();
						
					while (mRec->is(skipped.c_str()))
					{
						s << mRec->vS(11) << endl;
						GetNextRecord();
					}
					
					getCategory("pdbx_database_remark")->emplace({
						{ "id", remarkNr },
						{ "text", s.str() }
					});
					
					break;
				}
			}
		}
		catch (const exception& ex)
		{
			cerr << "Error parsing REMARK " << remarkNr << endl;
			throw;
		}
	}

	if (not (compoundDetails.empty() and sequenceDetails.empty() and sourceDetails.empty()))
	{
		getCategory("pdbx_entry_details")->emplace({
			{ "entry_id",			mStructureId },
			{ "compound_details",	compoundDetails },
			{ "sequence_details",	sequenceDetails },
			{ "source_details",		sourceDetails }
		});
	}

	// store remark 200 info (special case)
	if (not mRemark200.empty())
		ParseRemark200();
}

void PDBFileParser::ParseRemark200()
{
	auto rm200 = [&](const char* name, int diffrnNr) -> string
	{
		int nr = 0;
		string result;
		
		for (auto i = make_split_iterator(mRemark200[name],
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
			if (not this->mRemark200[n].empty())
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
			not mRefinementSoftware.empty())
			getCategory("computing")->emplace({
				{ "entry_id", mStructureId },
				{ "pdbx_data_reduction_ii", mRemark200["INTENSITY-INTEGRATION SOFTWARE"] },
				{ "pdbx_data_reduction_ds", mRemark200["DATA SCALING SOFTWARE"] },
				{ "structure_solution", mRemark200["SOFTWARE USED"] },
				{ "structure_refinement", mRefinementSoftware }
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
		if (mRemark200[sw.b].empty())
			continue;
		
		getCategory("software")->emplace({
			{ "name", mRemark200[sw.b] },
			{ "classification", sw.a },
			{ "version", "." },
			{ "pdbx_ordinal", mNextSoftwareOrd++ }
		});
	}

	string scatteringType;
	if (mRemark200["EXPERIMENT TYPE"] == "X-RAY DIFFRACTION")
		scatteringType = "x-ray";
	else if (mRemark200["EXPERIMENT TYPE"] == "NEUTRON DIFFRACTION")
		scatteringType = "neutron";
	
	set<string> diffrnWaveLengths;
	
	for (int diffrnNr = 1; ; ++diffrnNr)
	{
		string ambientTemp = rm200("TEMPERATURE (KELVIN)", diffrnNr);
		if (ambientTemp.empty())
			break;
		
		if (ba::ends_with(ambientTemp, "K"))
			ambientTemp.erase(ambientTemp.length() - 1, 1);
		
		getCategory("diffrn")->emplace({
			{ "id", diffrnNr },
			{ "ambient_temp", ambientTemp },
	//		{ "ambient_temp_details", seqId },
			{ "crystal_id", 1 }
		});
		
		string collectionDate;
		boost::system::error_code ec;
		collectionDate = pdb2cifDate(rm200("DATE OF DATA COLLECTION", diffrnNr), ec);
		if (ec)
		{
			if (VERBOSE)
				cerr << ec.message() << " for pdbx_collection_date" << endl;
			
			// The date field can become truncated when multiple values are available		
			if (diffrnNr != 1)
				collectionDate.clear();
		}
				
		getCategory("diffrn_detector")->emplace({
			{ "diffrn_id", diffrnNr },
			{ "detector", rm200("DETECTOR TYPE", diffrnNr) },
			{ "type", rm200("DETECTOR MANUFACTURER", diffrnNr) },
			{ "pdbx_collection_date", collectionDate },
			{ "details", rm200("OPTICS", diffrnNr) }
		});
		
		if (inRM200({"MONOCHROMATIC OR LAUE (M/L)", "MONOCHROMATOR", "DIFFRACTION PROTOCOL"}) or not scatteringType.empty())
			getCategory("diffrn_radiation")->emplace({
				{ "diffrn_id", diffrnNr },
				{ "wavelength_id", 1 },
				{ "pdbx_monochromatic_or_laue_m_l", rm200("MONOCHROMATIC OR LAUE (M/L)", diffrnNr) },
				{ "monochromator", rm200("MONOCHROMATOR", diffrnNr) },
				{ "pdbx_diffrn_protocol", rm200("DIFFRACTION PROTOCOL", diffrnNr) },
				{ "pdbx_scattering_type", scatteringType }
			});

		vector<string> wavelengths;
		string wl = rm200("WAVELENGTH OR RANGE (A)", diffrnNr);
		ba::split(wavelengths, wl, ba::is_any_of(", -"), ba::token_compress_on);
		
		diffrnWaveLengths.insert(wavelengths.begin(), wavelengths.end());

		string source;
		if (rm200("SYNCHROTRON (Y/N)", diffrnNr) == "Y")
		{
			getCategory("diffrn_source")->emplace({
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
			getCategory("diffrn_source")->emplace({
				{ "diffrn_id", diffrnNr },
				{ "source", rm200("RADIATION SOURCE", diffrnNr) },
				{ "type", rm200("X-RAY GENERATOR MODEL", diffrnNr) },

				{ "pdbx_wavelength", wavelengths.size() == 1 ? wavelengths[0] : "" },
				{ "pdbx_wavelength_list", wavelengths.size() == 1 ? "" : ba::join(wavelengths, ", ") },
			});
		}
	}

	int wavelengthNr = 1;
	for (auto wl: diffrnWaveLengths)
	{
		if (ba::ends_with(wl, "A"))
			wl.erase(wl.length() - 1, 1);
		
		getCategory("diffrn_radiation_wavelength")->emplace({
			{ "id", wavelengthNr++ },
			{ "wavelength", wl.empty() ? "." : wl },
			{ "wt", "1.0" }
		});
	}

	if (inRM200({"METHOD USED TO DETERMINE THE STRUCTURE", "STARTING MODEL"}))
	{
		auto cat = getCategory("refine");
		assert(cat->empty());
		
		string resolution = mRemark200["RESOLUTION RANGE HIGH (A)"];
		if (resolution.empty())
			resolution = ".";
		
		cat->emplace({
			{ "pdbx_method_to_determine_struct", mRemark200["METHOD USED TO DETERMINE THE STRUCTURE"] },
			{ "pdbx_starting_model", mRemark200["STARTING MODEL"] },
			{ "ls_d_res_high", resolution },
			{ "pdbx_diffrn_id", 1 },
			{ "pdbx_refine_id", mExpMethod },
			{ "entry_id", mStructureId }
		});
	}
	
	if (inRM200({"REJECTION CRITERIA (SIGMA(I))", "RESOLUTION RANGE HIGH (A)", "RESOLUTION RANGE LOW (A)", "NUMBER OF UNIQUE REFLECTIONS", "COMPLETENESS FOR RANGE (%)", "<I/SIGMA(I)> FOR THE DATA SET", "R MERGE (I)", "R SYM (I)", "DATA REDUNDANCY"}))
	{
		auto cat = getCategory("reflns");
		Row r;
		if (cat->empty())
			cat->emplace({});
		r = cat->back();
		r["entry_id"] = mStructureId;
		r["observed_criterion_sigma_I"] = mRemark200["REJECTION CRITERIA (SIGMA(I))"];
		r["d_resolution_high"] = mRemark200["RESOLUTION RANGE HIGH (A)"];
		r["d_resolution_low"] = mRemark200["RESOLUTION RANGE LOW (A)"];
		r["number_obs"] = mRemark200["NUMBER OF UNIQUE REFLECTIONS"];
		r["percent_possible_obs"] = mRemark200["COMPLETENESS FOR RANGE (%)"];
		r["pdbx_netI_over_sigmaI"] = mRemark200["<I/SIGMA(I)> FOR THE DATA SET"];
		r["pdbx_Rmerge_I_obs"] = mRemark200["R MERGE (I)"];
		r["pdbx_Rsym_value"] = mRemark200["R SYM (I)"];
		r["pdbx_redundancy"] = mRemark200["DATA REDUNDANCY"];
		r["pdbx_ordinal"] = 1;
		r["pdbx_diffrn_id"] = 1;
	}
	
	if (inRM200({ "HIGHEST RESOLUTION SHELL, RANGE HIGH (A)"})) // that one field is mandatory...
	{
		getCategory("reflns_shell")->emplace({
			{ "d_res_high", mRemark200["HIGHEST RESOLUTION SHELL, RANGE HIGH (A)"] },
			{ "d_res_low", mRemark200["HIGHEST RESOLUTION SHELL, RANGE LOW (A)"] },
			{ "percent_possible_all", mRemark200["COMPLETENESS FOR SHELL (%)"] },
			{ "Rmerge_I_obs", mRemark200["R MERGE FOR SHELL (I)"] },
			{ "pdbx_Rsym_value", mRemark200["R SYM FOR SHELL (I)"] },
			{ "meanI_over_sigI_obs", mRemark200["<I/SIGMA(I)> FOR SHELL"] },
			{ "pdbx_redundancy", mRemark200["DATA REDUNDANCY IN SHELL"] },
			{ "pdbx_ordinal", 1},
			{ "pdbx_diffrn_id" , 1}
		});
	}
	else if (inRM200({ "HIGHEST RESOLUTION SHELL, RANGE LOW (A)", "COMPLETENESS FOR SHELL (%)",
		"R MERGE FOR SHELL (I)", "R SYM FOR SHELL (I)", "<I/SIGMA(I)> FOR SHELL", "DATA REDUNDANCY IN SHELL" }))
	{
		if (VERBOSE)
			cerr << "Not writing reflns_shell record since d_res_high is missing" << endl;
	}
}

void PDBFileParser::ParseRemark350()
{
	auto saved = mRec;
	
	enum State { eStart, eInfo, eAnd, eApply, eBioMT } state = eStart;
	
	const regex
		kRX1(R"(BIOMOLECULE: (\d+))"),
		kRX2(R"(([^:]+): (.+?)(?: (ANGSTROM\*\*2|KCAL/MOL))?)"),
		kRX8(R"(APPLY THE FOLLOWING TO CHAINS: (.+))"),
		kRX9(R"(AND CHAINS: (.+))"),
		kRX10(R"(BIOMT([123])\s+(\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?))");
	
	int biomolecule = 0, operId = 0;
	vector<string> operExpression;
	map<string,string> values;
	vector<string> asymIdList;
	smatch m;
	Row genR;
	
	vector<double> mat, vec;
	
	for (mRec = FindRecord("REMARK 350"); mRec != nullptr and mRec->is("REMARK 350"); GetNextRecord())
	{
		string line = vS(11);
		
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
						
						MapChainID2AsymIDS(chain[0], asymIdList);
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
						
						MapChainID2AsymIDS(chain[0], asymIdList);
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

					operId = stoi(m[2].str());
					operExpression.push_back(to_string(operId));

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
						operId = stoi(m[2].str());
						operExpression.push_back(to_string(operId));
					}
					else if (operId != stoi(m[2].str()))
						throw runtime_error("Invalid REMARK 350");
					
					mat.push_back(stod(m[3].str())); 
					mat.push_back(stod(m[4].str())); 
					mat.push_back(stod(m[5].str()));
					vec.push_back(stod(m[6].str()));
					
					if (mt == 3)
					{
						if (vec.size() != 3 or mat.size() != 9)
							throw runtime_error("Invalid REMARK 350");
						
						if (operId == 1)
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
							
							getCategory("pdbx_struct_assembly")->emplace({
								{ "id",		biomolecule },
								{ "details", details },
								{ "method_details", values["SOFTWARE USED"] },
								{ "oligomeric_details", oligomer },
								{ "oligomeric_count", count > 0 ? to_string(count) : "" }
							});
							
							auto cat = getCategory("pdbx_struct_assembly_prop");
							
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
						
						getCategory("pdbx_struct_oper_list")->emplace({
							{ "id", operId },
							{ "type",
								mat == vector<double>{ 1, 0, 0, 0, 1, 0, 0, 0, 1 } and vec == vector<double>{ 0, 0, 0 }
									? "identity operation" : "crystal symmetry operation" },
//										{ "name", "" }, 
//										{ "symmetryOperation", "" },
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

					getCategory("pdbx_struct_assembly_gen")->emplace({
						{ "assembly_id", biomolecule },
						{ "oper_expression", ba::join(operExpression, ",") },
						{ "asym_id_list", ba::join(asymIdList, ",") }
					});

					biomolecule = stoi(m[1].str());
					asymIdList.clear();
					operExpression.clear();
					
					state = eInfo;
				}
				break;
				
		}
	}
	
	if (not operExpression.empty())
	{
		getCategory("pdbx_struct_assembly_gen")->emplace({
			{ "assembly_id", biomolecule },
			{ "oper_expression", ba::join(operExpression, ",") },
			{ "asym_id_list", ba::join(asymIdList, ",") }
		});
	}
	
	mRec = saved;
}

void PDBFileParser::ParsePrimaryStructure()
{
	// First locate the DBREF record. Might be missing
	DBREF cur = { mStructureId };

	while (ba::starts_with(mRec->mName, "DBREF"))
	{
		if (mRec->is("DBREF "))						//	 1 -  6       Record name   "DBREF "                                                    
		{
			cur.PDBIDCode				= vS(8, 11);	//	 8 - 11       IDcode        idCode             ID code of this datablock.                   
			cur.chainID					= vC(13);      //	13            Character     chainID            Chain  identifier.                       
			cur.seqBegin				= vI(15, 18);  //	15 - 18       Integer       seqBegin           Initial sequence number of the           
                                                        //	                                               PDB sequence segment.                    
			cur.insertBegin				= vC(19);      //	19            AChar         insertBegin        Initial  insertion code of the           
                                                        //	                                               PDB  sequence segment.                   
			cur.seqEnd					= vI(21, 24);  //	21 - 24       Integer       seqEnd             Ending sequence number of the            
                                                        //	                                               PDB  sequence segment.                   
			cur.insertEnd				= vC(25);      //	25            AChar         insertEnd          Ending insertion code of the             
                                                        //	                                               PDB  sequence segment.                   
			cur.database				= vS(27, 32);  //	27 - 32       LString       database           Sequence database name.                  
			cur.dbAccession				= vS(34, 41);  //	34 - 41       LString       dbAccession        Sequence database accession code.        
			cur.dbIdCode				= vS(43, 54);  //	43 - 54       LString       dbIdCode           Sequence  database identification code.  
			cur.dbSeqBegin				= vI(56, 60);  //	56 - 60       Integer       dbseqBegin         Initial sequence number of the           
                                                        //	                                               database seqment.                        
			cur.dbinsBeg				= vC(61);      //	61            AChar         idbnsBeg           Insertion code of initial residue of the 
                                                        //	                                               segment, if PDB is the reference.        
			cur.dbSeqEnd				= vI(63, 67);  //	63 - 67       Integer       dbseqEnd           Ending sequence number of the            
                                                        //	                                               database segment.                        
			cur.dbinsEnd				= vC(68);      //	68            AChar         dbinsEnd           Insertion code of the ending residue of  
			                                            //	                                               the segment, if PDB is the reference.    
			auto& chain = GetChainForID(cur.chainID);
			chain.mDbref = cur;
		}
		else if (mRec->is("DBREF1"))					//	 1 -  6        Record name   "DBREF1"                                             
		{
			cur.PDBIDCode				= vS(8, 11);	//	 8 - 11       IDcode        idCode             ID code of this datablock.                   
			cur.chainID					= vC(13);      //	13             Character     chainID       Chain identifier.                      
			cur.seqBegin				= vI(15, 18);  //	15 - 18        Integer       seqBegin      Initial sequence number of the         
                                                        //	                                           PDB sequence segment, right justified. 
			cur.insertBegin				= vC(19);      //	19             AChar         insertBegin   Initial insertion code of the          
                                                        //	                                           PDB sequence segment.                  
			cur.seqEnd					= vI(21, 24);  //	21 - 24        Integer       seqEnd        Ending sequence number of the          
                                                        //	                                           PDB sequence segment, right justified. 
			cur.insertEnd				= vC(25);      //	25             AChar         insertEnd     Ending insertion code of the           
                                                        //	                                           PDB sequence  segment.                 
			cur.database				= vS(27, 32);  //	27 - 32        LString       database      Sequence database name.                
			cur.dbIdCode				= vS(48, 67);  //	48 - 67        LString       dbIdCode      Sequence database identification code, 
		}
		else if (mRec->is("DBREF2"))					//	 1 -  6       Record name   "DBREF2"                                        
		{                                               //	 8 - 11       IDcode        idCode        ID code of this datablock.            
			if (vC(13) != cur.chainID)			        //	13            Character     chainID       Chain identifier.                 
				throw runtime_error("Chain ID's for DBREF1/DBREF2 records do not match");
			cur.dbAccession				= vS(19, 40);  //	19 - 40       LString       dbAccession   Sequence database accession code, 
                                                        //	                                          left justified.                   
			cur.dbSeqBegin				= vI(46, 55);  //	46 - 55       Integer       seqBegin      Initial sequence number of the    
                                                        //	                                          Database segment, right justified.
			cur.dbSeqEnd				= vI(58, 67);  //	58 - 67       Integer       seqEnd        Ending sequence number of the     
			                                            //	                                          Database segment, right justified.
			auto& chain = GetChainForID(cur.chainID);
			chain.mDbref = cur;
		}
		
		GetNextRecord();
	}

	// update chains
	for (auto& chain: mChains)
	{
		chain.mNextSeqNum = chain.mDbref.seqBegin;
		chain.mNextDbSeqNum = chain.mDbref.dbSeqBegin;
	}

	while (mRec->is("SEQADV"))
	{							//	 1 -  6        Record name   "SEQADV"                                           
		mSeqadvs.push_back({	//	 8 - 11        IDcode        idCode        ID  code of this datablock.              	
			vS(13, 15),        //	13 - 15        Residue name  resName       Name of the PDB residue in conflict. 
			vC(17),            //	17             Character     chainID       PDB  chain identifier.               
			vI(19, 22),        //	19 - 22        Integer       seqNum        PDB  sequence number.                
			vC(23),            //	23             AChar         iCode         PDB insertion code.                  
			vS(25, 28),        //	25 - 28        LString       database                                           
			vS(30, 38),        //	30 - 38        LString       dbAccession   Sequence  database accession number. 
			vS(40, 42),        //	40 - 42        Residue name  dbRes         Sequence database residue name.      
			vI(44, 48),        //	44 - 48        Integer       dbSeq         Sequence database sequence number.   
			vS(50, 70)         //	50 - 70        LString       conflict      Conflict comment.                                            
		});

		GetNextRecord();
	}
	
	while (mRec->is("SEQRES"))			//	 1 -  6        Record name    "SEQRES"
	{									//	 8 - 10        Integer        serNum       Serial number of the SEQRES record for  the
										//	                                           current  chain. Starts at 1 and increments
										//	                                           by one  each line. Reset to 1 for each chain.
		char chainId = vC(12);		//	12             Character      chainID      Chain identifier. This may be any single
										//	                                           legal  character, including a blank which is
										//	                                           is  used if there is only one chain.
		int numRes = vI(14, 17);		//	14 - 17        Integer        numRes       Number of residues in the chain.
										//	                                           This  value is repeated on every record.
		string monomers = vS(20, 70);	//	20 - 22        Residue name   resName      Residue name.
										//	 ...

		auto& chain = GetChainForID(chainId, numRes);
		
		for (auto si = ba::make_split_iterator(monomers, ba::token_finder(ba::is_any_of(" "), ba::token_compress_on)); not si.eof(); ++si)
		{
			string monId(si->begin(), si->end());
			if (monId.empty())
				continue;
			
			chain.mSeqres.push_back({monId, chain.mNextSeqNum++, ' ', chain.mNextDbSeqNum++});
			
			InsertChemComp(monId);
		}
		
		GetNextRecord();
	}

	// First pass over MODRES, only store relevant information required in ConstructEntities
	while (mRec->is("MODRES"))				//	 1 -  6        Record name   "MODRES"                                            												
	{							 			//	 8 - 11        IDcode        idCode      ID code of this datablock.                  
		string resName		= vS(13, 15);	//	13 - 15        Residue name  resName     Residue name used in this datablock.        
//		char chainID		= vC(17);		//	17             Character     chainID     Chain identifier.                       
//		int seqNum			= vI(19, 22);	//	19 - 22        Integer       seqNum      Sequence number.                        
//		char iCode			= vC(23);		//	23             AChar         iCode       Insertion code.                         
		string stdRes		= vS(25, 27);	//	25 - 27        Residue name  stdRes      Standard residue name.                  
//		string comment		= vS(30, 70);	//	30 - 70        String        comment     Description of the residue modification.

		mMod2parent[resName] = stdRes;

		GetNextRecord();
	}
}

void PDBFileParser::ParseHeterogen()
{
	while (mRec->is("HET   "))
	{									//	 1 -  6       Record name   "HET   "                                                         
		string hetID = vS(8, 10);      //	 8 - 10       LString(3)    hetID          Het identifier, right-justified.                  
		char chainID = vC(13);			//	13            Character     ChainID        Chain  identifier.                                
		int seqNum = vI(14, 17);		//	14 - 17       Integer       seqNum         Sequence  number.                                 
		char iCode = vC(18);			//	18            AChar         iCode          Insertion  code.                                  
		int numHetAtoms = vI(21, 25);	//	21 - 25       Integer       numHetAtoms    Number of HETATM records for the group            
										//	                                           present in the datablock.                             
		string text = vS(31, 70);		//	31 - 70       String        text           Text describing Het group.                        

		mHets.push_back({ hetID, chainID, seqNum, iCode, numHetAtoms, text });
		
		GetNextRecord();
	}
	
	for (;;)
	{
		if (mRec->is("HETNAM"))		//	 1 -  6       Record name   "HETNAM"                                                 
		{								//	 9 - 10       Continuation  continuation    Allows concatenation of multiple records.
			string hetID = vS(12, 14);	//	12 - 14       LString(3)    hetID           Het identifier, right-justified.         
	        string text = vS(16);		//	16 - 70       String        text            Chemical name.                           
	
			mHetnams[hetID] = text;
			InsertChemComp(hetID);
			
			GetNextRecord();
			continue;
		}
	
		if (mRec->is("HETSYN"))		 //	 1 -  6       Record name   "HETSYN"                                                           
		{                                //	 9 - 10       Continuation  continuation   Allows concatenation of multiple records.           
	         string hetID = vS(12, 14); //	12 - 14       LString(3)    hetID          Het identifier, right-justified.                    
	         string syn = vS(16);		 //	16 - 70       SList         hetSynonyms    List of synonyms.                                   
	
			mHetsyns[hetID] = syn;
	
			GetNextRecord();
			continue;
		}
		
		break;
	}
	
	while (mRec->is("FORMUL"))			//	 1 -  6        Record name   "FORMUL"                          
	{                                   //	 9 - 10        Integer       compNum       Component  number.  
        string hetID = vS(13, 15);     //	13 - 15        LString(3)    hetID         Het identifier.     
                                        //	17 - 18        Integer       continuation  Continuation number.
		char waterMark = vC(19);		//	19             Character     asterisk      "*" for water.      
		string formula = vS(20);		//	20 - 70        String        text          Chemical formula.   
		
		mFormuls[hetID] = formula;
		
		if (waterMark == '*')
			mWaterHetId = hetID;
		
		GetNextRecord();
	}
}

void PDBFileParser::ConstructEntities()
{
	// We parsed the Primary Structure and Heterogen sections, if available.
	// But if we didn't parse anything, we need to fake the data based on residues in ATOM records

	// First iterate all ATOM records and store the residues as found in these records
	int modelNr = 1;
	
	typedef map<tuple<char,int,char,char>,string> CompTypeMap;
	CompTypeMap residuesSeen;	// used to validate PDB files...
	
	for (auto r = mData; r != nullptr; r = r->mNext)
	{
		if (r->is("MODEL "))
		{
			modelNr = r->vI(11, 14);
			if (modelNr != 1)
				break;
			continue;
		}

		if (r->is("ATOM  ") or r->is("HETATM"))		//	 1 -  6        Record name   "ATOM  "
		{											//	 ...
			string name			= r->vS(13, 16);	//	13 - 16        Atom          name         Atom name.
			char altLoc			= r->vC(17);		//	17             Character     altLoc       Alternate location indicator.
			string resName		= r->vS(18, 20);	//	18 - 20        Residue name  resName      Residue name.
			char chainID		= r->vC(22);		//	22             Character     chainID      Chain identifier.
			int resSeq			= r->vI(23, 26);	//	23 - 26        Integer       resSeq       Residue sequence number.
			char iCode			= r->vC(27);		//	27             AChar         iCode        Code for insertion of residues.

			// first validate, too sad this is required...
			CompTypeMap::key_type k = make_tuple(chainID, resSeq, iCode, altLoc);
			if (residuesSeen.count(k) == 0)
				residuesSeen[k] = resName;
			else if (residuesSeen[k] != resName)
				throw runtime_error("inconsistent residue type for " + string{chainID} + to_string(resSeq) + iCode + altLoc + "\n" +
					"  (" + residuesSeen[k] + " != " + resName + ")");

			auto& chain = GetChainForID(chainID);
			
			PDBChain::AtomRes ar{ resName, resSeq, iCode };

			if ((chain.mResiduesSeen.empty() or chain.mResiduesSeen.back() != ar) and
				(CompoundFactory::instance().isKnownPeptide(resName) or CompoundFactory::instance().isKnownBase(resName)))
			{
				chain.mResiduesSeen.push_back(ar);
			}

			// now that we're iterating atoms anyway, clean up the mUnobs array
			mUnobs.erase(remove_if(mUnobs.begin(), mUnobs.end(), [=](UNOBS& a)
			{
				bool result = false;
				
				if (modelNr == a.modelNr and
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
			}), mUnobs.end());

			continue;
		}

		if (r->is("TER   "))						//	 1 -  6 	   Record name	 "TER	"								  
		{											//	 7 - 11 	   Integer		 serial 		 Serial number. 		  
													//	18 - 20 	   Residue name  resName		 Residue name.			  
			char chainID = r->vC(22);				//	22			   Character	 chainID		 Chain identifier.		  
													//	23 - 26 	   Integer		 resSeq 		 Residue sequence number. 
													//	27			   AChar		 iCode			 Insertion code.		  
			auto& chain = GetChainForID(chainID);
			if (chain.mTerIndex == 0)				// Is this the first TER record? (Refmac writes out multiple TER records...)
				chain.mTerIndex = chain.mResiduesSeen.size();
			continue;
		}
	}
	
	for (auto& chain: mChains)
	{
		if (not (chain.mSeqres.empty() or chain.mResiduesSeen.empty()))
		{
			// seems safe to assume TER record is at the right location...
			// However, some files don't have them at all.
			// When mTerIndex == 0 this is most likely the case. Right?
			
			if (chain.mTerIndex > 0)
				chain.mResiduesSeen.erase(chain.mResiduesSeen.begin() + chain.mTerIndex, chain.mResiduesSeen.end());
			
			int lastResidueIndex = chain.AlignResToSeqRes();
			
			if (lastResidueIndex > 0 and lastResidueIndex + 1 < static_cast<int>(chain.mResiduesSeen.size()))
			{
				auto& r = chain.mResiduesSeen[lastResidueIndex + 1];

				if (VERBOSE)
				{
					cerr << "Detected residues that cannot be aligned to SEQRES" << endl
						 << "First residue is " << chain.mDbref.chainID << ':' << r.mSeqNum << r.mIcode << endl;
				}
				
				chain.mTerIndex = lastResidueIndex + 1;
			}
		}
		else
		{
			// So, we did not have a SEQRES for this chain. Try to reconstruct it.
			// Problem here is that TER records may be located incorrectly. So
			// first lets shift the ter index until it is past the last known
			// aminoacid or base.
			
			for (int ix = chain.mTerIndex; ix < static_cast<int>(chain.mResiduesSeen.size()); ++ix)
			{
				string resName = chain.mResiduesSeen[ix].mMonId;
				
				if (mmcif::kAAMap.count(resName) or
					mmcif::kBaseMap.count(resName) or
					CompoundFactory::instance().isKnownPeptide(resName) or
					CompoundFactory::instance().isKnownBase(resName))
				{
					chain.mTerIndex = ix + 1;
				}

				InsertChemComp(resName);
			}
			
			// And now construct our 'SEQRES'...
			for (int ix = 0; ix < chain.mTerIndex; ++ix)
			{
				auto& ar = chain.mResiduesSeen[ix];
				chain.mSeqres.push_back({ar.mMonId, ar.mSeqNum, ar.mIcode, ar.mSeqNum, true});
			}
		}
	}

	set<char> terminatedChains;
	map<char,int> residuePerChainCounter;

	for (auto r = mData; r != nullptr; r = r->mNext)
	{
		if (r->is("MODEL "))
		{
			modelNr = r->vI(11, 14);
			if (modelNr != 1)
				break;
			continue;
		}

		if (r->is("ATOM  ") or r->is("HETATM"))
		{											//	 1 -  6        Record name   "ATOM  "
			int serial = r->vI(7, 11);				//	 7 - 11        Integer       serial       Atom  serial number.
													//	 ...
			char altLoc = vC(17);					//	17             Character     altLoc       Alternate location indicator.
			string resName		= r->vS(18, 20);	//	18 - 20        Residue name  resName      Residue name.
			char chainID		= r->vC(22);		//	22             Character     chainID      Chain identifier.
			int resSeq			= r->vI(23, 26);	//	23 - 26        Integer       resSeq       Residue sequence number.
			char iCode			= r->vC(27);		//	27             AChar         iCode        Code for insertion of residues.

			auto& chain = GetChainForID(chainID);

			auto i = find(chain.mSeqres.begin(), chain.mSeqres.end(), PDBSeqRes{resName, resSeq, iCode});

			// might be a hetero
			if (altLoc != ' ' and i == chain.mSeqres.end())
			{
				i = find_if(chain.mSeqres.begin(), chain.mSeqres.end(),
					[resSeq, iCode](const PDBSeqRes& r) -> bool
					{
						return r.mSeqNum == resSeq and r.mIcode == iCode;
					});
			}

			if (i != chain.mSeqres.end())
			{
				i->mSeen = true;
				if (i->mMonId != resName)
					i->mAlts.insert(resName);
			}
			else
			{
				auto& residues = chain.mHet;
	
				if (residues.empty() or residues.back().mSeqNum != resSeq)
				{
					i = lower_bound(residues.begin(), residues.end(),
						PDBSeqRes{resName, resSeq, iCode},
						[=](const PDBSeqRes& r1, const PDBSeqRes& r2) -> bool {
							return r1.mSeqNum < r2.mSeqNum;
						});
	
					residues.insert(i, { resName, resSeq, iCode, resSeq, true });

					InsertChemComp(resName);
				}
			}
			
			int residueCount = (residuePerChainCounter[chainID] += 1);
			
			// There appears to be a program that writes out HETATM records as ATOM records....
			if (r->is("HETATM") or
				terminatedChains.count(chainID) or
				(chain.mTerIndex > 0 and residueCount >= chain.mTerIndex))
			{
				if (isWater(resName))
					mWaterHetId = resName;
				
				auto h = find_if(mHets.begin(), mHets.end(), [=](const HET& het) -> bool
					{
						return het.hetID == resName and het.chainID == chainID and
							het.seqNum == resSeq and het.iCode == iCode;
					});
				
				if (h == mHets.end())
				{
					mHets.push_back({ resName, chainID, resSeq, iCode, 0 });	// double perhaps, but that does not care
					h = prev(mHets.end());
				}

				h->atoms.insert(serial);
			}

			continue;
		}

		if (r->is("TER   "))
		{
			char chainID		= r->vC(22);		//	22             Character     chainID      Chain identifier.
			terminatedChains.insert(chainID);
		}
	}

	// Create missing compounds
	for (auto& chain: mChains)
	{
		if (chain.mMolId != 0 or chain.mSeqres.empty())
			continue;

		// now this chain may contain the same residues as another one
		for (auto& other: mChains)
		{
			if (&other == &chain or other.mMolId == 0)
				continue;
			
			if (chain.SameSequence(other))
			{
				chain.mMolId = other.mMolId;
				break;
			}
		}			
			
		if (chain.mMolId != 0)
			continue;
		
		auto& comp = GetOrCreateCompound(mNextMolId++);
		comp.mChains.insert(chain.mDbref.chainID);

		chain.mMolId = comp.mMolId;
	}

	set<string> structTitle, structDescription;
	
	// Create poly_scheme and write pdbx_poly_seq_scheme and create mapping table

	auto cat = getCategory("pdbx_poly_seq_scheme");
	int asymNr = 0;
	for (auto& chain: mChains)
	{
		string asymId = cifIdForInt(asymNr++);
		string entityId = mMolID2EntityID[chain.mMolId];
		
		mAsymID2EntityID[asymId] = entityId;
		
		getCategory("struct_asym")->emplace({
			{ "id", asymId },
			{ "pdbx_blank_PDB_chainid_flag", chain.mDbref.chainID == ' ' ? "Y" : "N" },
//			pdbx_modified 
			{ "entity_id", entityId },
//			details
		});
		
		int seqNr = 1;
		for (auto& res: chain.mSeqres)
		{
			mChainSeq2AsymSeq[make_tuple(chain.mDbref.chainID, res.mSeqNum, res.mIcode)] = make_tuple(asymId, seqNr, true);
			
			string seqId = to_string(seqNr);
			++seqNr;
			
			set<string> monIds = { res.mMonId };
			monIds.insert(res.mAlts.begin(), res.mAlts.end());
			
			for (string monId: monIds)
			{
				string authMonId, authSeqNum;
				if (res.mSeen)
				{
					authMonId = monId;
					authSeqNum = to_string(res.mSeqNum);
				}

				cat->emplace({
					{ "asym_id", asymId },
					{ "entity_id", mMolID2EntityID[chain.mMolId] },
					{ "seq_id", seqId },
					{ "mon_id", monId },
					{ "ndb_seq_num", seqId },
					{ "pdb_seq_num", res.mSeqNum },
					{ "auth_seq_num", authSeqNum },
					{ "pdb_mon_id", authMonId },
					{ "auth_mon_id", authMonId },
					{ "pdb_strand_id", string{chain.mDbref.chainID} },
					{ "pdb_ins_code", (res.mIcode == ' ' or res.mIcode == 0) ? "." : string{res.mIcode} },
					{ "hetero", res.mAlts.empty() ? "n" : "y" }
				});
			}
		}
	}
	
	// We have now created all compounds, write them out
	uint32 structRefId = 0, structRefSeqAlignId = 0;
	
	for (auto& cmp: mCompounds)
	{
		++structRefId;

		string srcMethod;
		
		if (not cmp.mSource["SYNTHETIC"].empty())
		{
			srcMethod = "syn";
			
			getCategory("pdbx_entity_src_syn")->emplace({
				{ "entity_id", mMolID2EntityID[cmp.mMolId] },
				{ "pdbx_src_id", structRefId },
				{ "organism_scientific", cmp.mSource["ORGANISM_SCIENTIFIC"] },
				{ "ncbi_taxonomy_id", cmp.mSource["ORGANISM_TAXID"] },
			});
		}
		else if (cmp.mInfo["ENGINEERED"] == "YES" or
			not cmp.mSource["EXPRESSION_SYSTEM"].empty())
		{
			srcMethod = "man";
			
			getCategory("entity_src_gen")->emplace({
				{ "entity_id", mMolID2EntityID[cmp.mMolId] },
				{ "pdbx_src_id", structRefId },
				{ "gene_src_common_name", cmp.mSource["ORGANISM_COMMON"] },
				{ "pdbx_gene_src_gene", cmp.mSource["GENE"] },
				{ "gene_src_strain", cmp.mSource["STRAIN"] },
				{ "gene_src_tissue", cmp.mSource["TISSUE"] },
				{ "gene_src_tissue_fraction", cmp.mSource["TISSUE_FRACTION"] },
				{ "pdbx_gene_src_cell_line", cmp.mSource["CELL_LINE"] },
				{ "pdbx_gene_src_organelle", cmp.mSource["ORGANELLE"] },
				{ "pdbx_gene_src_cell", cmp.mSource["CELL"] },
				{ "pdbx_gene_src_cellular_location", cmp.mSource["CELLULAR_LOCATION"] },
				{ "host_org_common_name", cmp.mSource["EXPRESSION_SYSTEM_COMMON"] },
				{ "pdbx_gene_src_scientific_name", cmp.mSource["ORGANISM_SCIENTIFIC"] },
				{ "pdbx_gene_src_ncbi_taxonomy_id", cmp.mSource["ORGANISM_TAXID"] },
				{ "pdbx_host_org_scientific_name", cmp.mSource["EXPRESSION_SYSTEM"] },
				{ "pdbx_host_org_ncbi_taxonomy_id", cmp.mSource["EXPRESSION_SYSTEM_TAXID"] },
				{ "pdbx_host_org_strain", cmp.mSource["EXPRESSION_SYSTEM_STRAIN"] },
				{ "pdbx_host_org_variant", cmp.mSource["EXPRESSION_SYSTEM_VARIANT"] },
				{ "pdbx_host_org_cell_line", cmp.mSource["EXPRESSION_SYSTEM_CELL_LINE"] },
				{ "pdbx_host_org_cellular_location", cmp.mSource["EXPRESSION_SYSTEM_CELLULAR_LOCATION"] },
				{ "pdbx_host_org_vector_type", cmp.mSource["EXPRESSION_SYSTEM_VECTOR_TYPE"] },
				{ "pdbx_host_org_vector", cmp.mSource["EXPRESSION_SYSTEM_VECTOR"] },
				{ "pdbx_host_org_gene", cmp.mSource["EXPRESSION_SYSTEM_GENE"] },
				{ "plasmid_name", cmp.mSource["EXPRESSION_SYSTEM_PLASMID"] },
				{ "pdbx_description", cmp.mSource["OTHER_DETAILS"] }
			});
		}
		else if (not cmp.mSource["ORGANISM_SCIENTIFIC"].empty())
		{
			srcMethod = "nat";
			
			getCategory("entity_src_nat")->emplace({
				{ "entity_id", mMolID2EntityID[cmp.mMolId] },
				{ "pdbx_src_id", structRefId },
				{ "common_name", cmp.mSource["ORGANISM_COMMON"] },
				{ "strain", cmp.mSource["STRAIN"] },
				{ "pdbx_secretion", cmp.mSource["SECRETION"] },
				{ "pdbx_organism_scientific", cmp.mSource["ORGANISM_SCIENTIFIC"] },
				{ "pdbx_ncbi_taxonomy_id", cmp.mSource["ORGANISM_TAXID"] },
				{ "pdbx_cellular_location", cmp.mSource["CELLULAR_LOCATION"] },
				{ "pdbx_plasmid_name", cmp.mSource["PLASMID"] },
				{ "pdbx_organ", cmp.mSource["ORGAN"] },
			});
		}
		
		getCategory("entity")->emplace({
			{ "id", mMolID2EntityID[cmp.mMolId] },
			{ "type", "polymer" },
			{ "src_method", srcMethod },
			{ "pdbx_description", cmp.mInfo["MOLECULE"] },
//			{ "pdbx_formula_weight", 		},
			{ "pdbx_number_of_molecules", cmp.mChains.size() },
			{ "details", cmp.mInfo["OTHER_DETAILS"] },
			{ "pdbx_mutation", cmp.mInfo["MUTATION"] },
			{ "pdbx_fragment", cmp.mInfo["FRAGMENT"] },
			{ "pdbx_ec", cmp.mInfo["EC"] }
		});
		
		if (not cmp.mInfo["SYNONYM"].empty())
		{
			getCategory("entity_name_com")->emplace({
				{ "entity_id", mMolID2EntityID[cmp.mMolId] },
				{ "name", cmp.mInfo["SYNONYM"] }
			});
		}

		string desc = cmp.mInfo["MOLECULE"];
		if (not cmp.mInfo["EC"].empty())
			desc += " (E.C." + cmp.mInfo["EC"] + ")";

		if (not cmp.mTitle.empty())
			structTitle.insert(cmp.mTitle);
		
		if (not desc.empty())
			structDescription.insert(desc);
		
		auto ci = find_if(mChains.begin(), mChains.end(),
			[cmp](PDBChain& c) -> bool { return cmp.mChains.count(c.mDbref.chainID); } );

		if (ci != mChains.end() and not ci->mDbref.dbIdCode.empty())
		{
			getCategory("struct_ref")->emplace({
				{ "id", structRefId },
				{ "entity_id", mMolID2EntityID[cmp.mMolId] },
				{ "db_name", ci->mDbref.database },
				{ "db_code", ci->mDbref.dbIdCode },
				{ "pdbx_db_accession", ci->mDbref.dbAccession },
//				{ "pdbx_align_begin", ci->mDbref.dbSeqBegin }
			});
		}
		
		bool nstdMonomer = false, nonstandardLinkage = false;
		bool mightBePolyPeptide = true, mightBeDNA = true;
		
		vector<string> chains;
		string seq, seqCan;
		
		// write out the chains for this compound
		for (auto& chain: mChains)
		{
			if (chain.mMolId != cmp.mMolId)
				continue;
			
//			chain.mEntityId = cmp.mEntityId;

			++structRefSeqAlignId;
			DBREF& dbref = chain.mDbref;
			
			if (not dbref.database.empty())
			{
				auto insToStr = [](char i) -> string { return i == ' ' or not isprint(i) ? "" : string{ i }; };
				
				auto& pdbxPolySeqScheme = *getCategory("pdbx_poly_seq_scheme");
				
				int seqAlignBeg = 0, seqAlignEnd = 0;
				
				try
				{
					seqAlignBeg = pdbxPolySeqScheme[
							Key("pdb_strand_id") == dbref.chainID and
							Key("pdb_seq_num") == dbref.seqBegin and
							Key("pdb_ins_code") == insToStr(dbref.insertBegin)]
						["seq_id"].as<int>();
	
					seqAlignEnd = pdbxPolySeqScheme[
							Key("pdb_strand_id") == dbref.chainID and
							Key("pdb_seq_num") == dbref.seqEnd and
							Key("pdb_ins_code") == insToStr(dbref.insertEnd)]
						["seq_id"].as<int>();
				}
				catch (...) {}
			
				getCategory("struct_ref_seq")->emplace({
					{ "align_id", structRefSeqAlignId },
					{ "ref_id", structRefId },
					{ "pdbx_PDB_id_code", dbref.PDBIDCode },
					{ "pdbx_strand_id", string{ chain.mDbref.chainID } },
					{ "seq_align_beg", seqAlignBeg },
					{ "pdbx_seq_align_beg_ins_code", insToStr(dbref.insertBegin) },
					{ "seq_align_end", seqAlignEnd },
					{ "pdbx_seq_align_end_ins_code", insToStr(dbref.insertEnd) },
					{ "pdbx_db_accession", dbref.dbAccession },
					{ "db_align_beg", dbref.dbSeqBegin },
					{ "pdbx_db_align_beg_ins_code", insToStr(dbref.dbinsBeg) },
					{ "db_align_end", dbref.dbSeqEnd },
					{ "pdbx_db_align_end_ins_code", insToStr(dbref.dbinsEnd) },
					{ "pdbx_auth_seq_align_beg", dbref.seqBegin },
					{ "pdbx_auth_seq_align_end", dbref.seqEnd }
				});

				// write the struct_ref_seq_dif
				for (auto& seqadv: mSeqadvs)
				{
					if (seqadv.chainID != chain.mDbref.chainID or seqadv.resName.empty())
						continue;
					
					string asym, seqNum;
					int labelSeq = -1;
					boost::system::error_code ec;

					tie(asym, labelSeq, ignore) = MapResidue(seqadv.chainID, seqadv.seqNum, seqadv.iCode, ec);
					if (ec)
					{
						if (VERBOSE)
							cerr << "dropping unmatched SEQADV record" << endl;
						continue;
					}
					
					seqNum = to_string(labelSeq);
					
					getCategory("struct_ref_seq_dif")->emplace({
						{ "align_id", structRefSeqAlignId },
						{ "pdbx_PDB_id_code", dbref.PDBIDCode },
						{ "mon_id", seqadv.resName },
						{ "pdbx_pdb_strand_id", seqadv.chainID },
						{ "seq_num", seqNum },
						{ "pdbx_pdb_ins_code", seqadv.iCode == ' ' ? string{} : string{seqadv.iCode} },
						{ "pdbx_seq_db_name", seqadv.database },
						{ "pdbx_seq_db_accession_code", seqadv.dbAccession },
						{ "db_mon_id", seqadv.dbRes },
						{ "pdbx_seq_db_seq_num", seqadv.dbSeq },
						{ "details", seqadv.conflict },
						{ "pdbx_auth_seq_num", seqadv.seqNum },
						{ "pdbx_ordinal", ++mPdbxDifOrdinal }						
					});
				}
			}
			
			if (not chains.empty())	// not the first one for this molId
			{
				chains.push_back( string{ chain.mDbref.chainID } );
				continue;
			}
			
			chains.push_back(string{chain.mDbref.chainID});

			size_t seqLen = 0, seqCanLen = 0;
			
			for (auto& res: chain.mSeqres)
			{
				string letter, stdRes;

				if (mMod2parent.count(res.mMonId))
					stdRes = mMod2parent.at(res.mMonId);

				if (mmcif::kAAMap.count(res.mMonId))
				{
					letter = mmcif::kAAMap.at(res.mMonId);
					mightBeDNA = false;
				}
				else if (mmcif::kBaseMap.count(res.mMonId))
				{
					letter = mmcif::kBaseMap.at(res.mMonId);
					mightBePolyPeptide = false;
				}
				else
				{
					nstdMonomer = true;
					letter = '(' + res.mMonId + ')';
					
					// sja...
					auto compound = mmcif::Compound::create(stdRes.empty() ? res.mMonId : stdRes);
					if (compound != nullptr and
						not iequals(compound->type(), "L-peptide linking") and
						not iequals(compound->type(), "RNA linking"))
					{
						nonstandardLinkage = true;
					}
				}
				
				if (seqLen + letter.length() > 80)
				{
					seq += '\n';
					seqLen = 0;
				}

				seq += letter;
				seqLen += letter.length();
	
				if (letter.length() > 1)
				{
					if (not stdRes.empty() and mmcif::kAAMap.count(stdRes))
						letter = mmcif::kAAMap.at(stdRes);
					else if (mmcif::kBaseMap.count(res.mMonId))
						letter = mmcif::kBaseMap.at(res.mMonId);
					else
						letter = 'X';
				}
	
				if (seqCanLen + letter.length() > 80)
				{
					seqCan += '\n';
					seqCanLen = 0;
				}
				seqCan += letter;
				seqCanLen += letter.length();
			}

			auto cat = getCategory("entity_poly_seq");
			for (size_t i = 0; i < chain.mSeqres.size(); ++i)
			{
				auto& rs = chain.mSeqres[i];
				
				cat->emplace({
					{ "entity_id", mMolID2EntityID[cmp.mMolId] },
					{ "num", i + 1 }, 
					{ "mon_id", rs.mMonId },
					{ "hetero", rs.mAlts.empty() ? "n" : "y"}
				});
				
				for (auto& a: rs.mAlts)
				{
					cat->emplace({
						{ "entity_id", mMolID2EntityID[cmp.mMolId] },
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

		getCategory("entity_poly")->emplace({
			{ "entity_id", mMolID2EntityID[cmp.mMolId] },
			{ "pdbx_seq_one_letter_code", seq },
			{ "pdbx_seq_one_letter_code_can", seqCan },
			{ "nstd_monomer", (nstdMonomer ? "yes" : "no") },
			{ "pdbx_strand_id", ba::join(chains, ",") },
			{ "nstd_linkage", nonstandardLinkage ? "yes" : "no" },
			{ "type", type }
		});
	}

	if (not (structTitle.empty() and structDescription.empty()))
	{
		getCategory("struct")->emplace({
			{ "entry_id", mStructureId },
			{ "title", ba::join(structTitle, ", ") },
			{ "pdbx_descriptor", ba::join(structDescription, ", ") },
			{ "pdbx_model_type_details", mModelTypeDetails }
		});
	}
	
	map<char,string> waterChains;
	map<tuple<string,string>,int> ndbSeqNum;	// for nonpoly scheme
	
	for (size_t i = 0; i < mHets.size(); ++i)
	{
		auto& heti = mHets[i];

		if (not heti.asymId.empty())
			continue;
	
		if (heti.hetID == mWaterHetId or isWater(heti.hetID))
			continue;

		// See if this residue is part of SEQRES
		auto& chain = GetChainForID(heti.chainID);
		auto ih = find(chain.mSeqres.begin(), chain.mSeqres.end(), PDBSeqRes{heti.hetID, heti.seqNum, heti.iCode}); 
		
		// If so, skip it, it is not an entity then
		if (ih != chain.mSeqres.end())
			continue;
		
		heti.asymId = cifIdForInt(asymNr++);
	}

	set<string> writtenAsyms;

	map<string,int> hetCount;		// for pdbx_number_of_molecules
	for (auto& het: mHets)
		hetCount[het.hetID] += 1;
	
	for (auto& het: mHets)
	{
		string hetID = het.hetID;
		
		auto& chain = GetChainForID(het.chainID);
		
		// See if this residue is part of SEQRES
		auto i = find(chain.mSeqres.begin(), chain.mSeqres.end(), PDBSeqRes{hetID, het.seqNum, het.iCode}); 
		
		// If so, skip it, it is not an entity then
		if (i != chain.mSeqres.end())
			continue;

		// See if we've already added it to the entities
		if (mHet2EntityID.count(hetID) == 0)
		{
			string entityId = to_string(mNextEntityNr++);
			mHet2EntityID[hetID] = entityId;
			
			if (hetID == mWaterHetId)
			{
				getCategory("entity")->emplace({
					{ "id", entityId },
					{ "type", "water" },
					{ "src_method", "nat" },
					{ "pdbx_description", "water" },
					{ "pdbx_number_of_molecules", hetCount[hetID] }
				});
			}
			else
			{
				if (mHetnams[hetID].empty())
				{
					auto compound = mmcif::Compound::create(hetID);
					if (compound != nullptr)
						mHetnams[hetID] = compound->name();
				}
				
				getCategory("entity")->emplace({
					{ "id", entityId },
					{ "type", "non-polymer" },
					{ "src_method", "syn" },
					{ "pdbx_description", mHetnams[hetID] },
					{ "details", mHetsyns[hetID] },
					{ "pdbx_number_of_molecules", hetCount[hetID] }
				});
			}

			// write a pdbx_entity_nonpoly record
			string name = mHetnams[hetID];
			if (name.empty() and hetID == mWaterHetId)
				name = "water";
			getCategory("pdbx_entity_nonpoly")->emplace({
				{ "entity_id", entityId },
				{ "name", name },
				{ "comp_id", hetID }
			});
		}
		
		// create an asym for this het/chain combo, if needed

		string asymId = het.asymId;

		auto k = make_tuple(het.chainID, het.seqNum, het.iCode);
		if (mChainSeq2AsymSeq.count(k) == 0)
		{
			if (hetID == mWaterHetId or isWater(hetID))
			{
				if (waterChains.count(het.chainID) == 0)
				{
					asymId = cifIdForInt(asymNr++);
					waterChains[het.chainID] = asymId;
				}
				else
					asymId = waterChains[het.chainID];
			}
			else
				asymId = het.asymId;
			
			assert(asymId.empty() == false);

			mAsymID2EntityID[asymId] = mHet2EntityID[hetID];
			
			// NOTE, a nonpoly residue has no label_seq_id
			// but in pdbx_nonpoly_scheme there is such a number.
			// Since this number is not used anywhere else we
			// just use it here and do not store it in the table 
			mChainSeq2AsymSeq[k] = make_tuple(asymId, 0, false);

			if (writtenAsyms.count(asymId) == 0)
			{
				writtenAsyms.insert(asymId);
				getCategory("struct_asym")->emplace({
					{ "id", asymId },
					{ "pdbx_blank_PDB_chainid_flag", het.chainID == ' ' ? "Y" : "N" },
	//					pdbx_modified 
					{ "entity_id", mHet2EntityID[hetID] },
	//					details
				});
			}
		}

		int seqNr = ++ndbSeqNum[make_tuple(hetID, asymId)];
		
		string iCode{het.iCode};
		ba::trim(iCode);
		if (iCode.empty())
			iCode = { '.' };
		
		getCategory("pdbx_nonpoly_scheme")->emplace({
			{ "asym_id", asymId },
			{ "entity_id", mHet2EntityID[hetID] },
			{ "mon_id", hetID },
			{ "ndb_seq_num", seqNr },
			{ "pdb_seq_num", het.seqNum },
//			{ "auth_seq_num", het.seqNum },	// ????
			{ "pdb_mon_id", hetID },
//			{ "auth_mon_id", hetID },
			{ "pdb_strand_id", string{het.chainID} },
			{ "pdb_ins_code", iCode }
		});

		// mapping needed?
		mChainSeq2AsymSeq[make_tuple(het.chainID, het.seqNum, het.iCode)] = make_tuple(asymId, seqNr, false);
	}

	int modResId = 1;
	set<string> modResSet;
	for (auto rec = FindRecord("MODRES"); rec != nullptr and rec->is("MODRES");
			rec = rec->mNext)					//	 1 -  6        Record name   "MODRES"                                            												
	{								 			//	 8 - 11        IDcode        idCode      ID code of this datablock.                  
		string resName		= rec->vS(13, 15);	//	13 - 15        Residue name  resName     Residue name used in this datablock.        
		char chainID		= rec->vC(17);		//	17             Character     chainID     Chain identifier.                       
		int seqNum			= rec->vI(19, 22);	//	19 - 22        Integer       seqNum      Sequence number.                        
		char iCode			= rec->vC(23);		//	23             AChar         iCode       Insertion code.                         
		string stdRes		= rec->vS(25, 27);	//	25 - 27        Residue name  stdRes      Standard residue name.                  
		string comment		= rec->vS(30, 70);	//	30 - 70        String        comment     Description of the residue modification.

		string asymId;
		int seq;
		boost::system::error_code ec;
		
		tie(asymId, seq, ignore) = MapResidue(chainID, seqNum, iCode, ec);
		if (ec)	// no need to write a modres if it could not be found
		{
			if (VERBOSE)
				cerr << "dropping unmapped MODRES record" << endl;
			continue;
		}
		
		getCategory("pdbx_struct_mod_residue")->emplace({
			{ "id", modResId++ },
			{ "label_asym_id", asymId },
			{ "label_seq_id", seq },
			{ "label_comp_id", resName },
			{ "auth_asym_id", string(1, chainID) },
			{ "auth_seq_id", seqNum },
			{ "auth_comp_id", resName },
			{ "PDB_ins_code", iCode == ' ' ? "" : string{ iCode } },
			{ "parent_comp_id", stdRes },
			{ "details", comment }
		});
		
		modResSet.insert(resName);
	}

//	// chem compounds

	for (auto cc: mChemComp)
	{
		auto compound = mmcif::Compound::create(
			mMod2parent.count(cc) ? mMod2parent[cc] : cc
		);
		
		string formula = mFormuls[cc];
		if (formula.empty() and compound != nullptr)
			formula = compound->formula();
		else
		{
			const regex rx(R"(\d+\((.+)\))");
			smatch m;
			if (regex_match(formula, m, rx))
				formula = m[1].str();
		}
		
		string name = mHetnams[cc];
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

		if (modResSet.count(cc))
			nstd = "n";

		getCategory("chem_comp")->emplace({
			{ "id", cc },
			{ "name", name },
			{ "formula", formula },
			{ "mon_nstd_flag", nstd },
			{ "type", type }
		});
	}
	
	getCategory("chem_comp")->reorderByIndex();
	
	// unobserved can now be written as well
	
	int idRes = 0, idAtom = 0;
	sort(mUnobs.begin(), mUnobs.end(), [](const UNOBS& a, const UNOBS& b) -> bool
	{
		int d = a.modelNr - b.modelNr;
		if (d == 0)
			d = a.seq - b.seq;
		return d < 0;
	}); 
	
	for (auto& unobs: mUnobs)
	{
		bool isPolymer = false;
		string asymId, compId = unobs.res;
		int seqNr = 0;
		boost::system::error_code ec;
		
		tie(asymId, seqNr, isPolymer) = MapResidue(unobs.chain, unobs.seq, unobs.iCode, ec);
		if (ec)
		{
			if (VERBOSE)
				cerr << "error mapping unobserved residue" << endl;
			continue;
		}
		
		if (unobs.atoms.empty())
		{
			getCategory("pdbx_unobs_or_zero_occ_residues")->emplace({
				{ "id",				to_string(++idRes) },
				{ "polymer_flag",	isPolymer ? "Y" : "N" },
				{ "occupancy_flag",	1 },
				{ "PDB_model_num",	unobs.modelNr ? unobs.modelNr : 1 },
				{ "auth_asym_id",	unobs.chain },
				{ "auth_comp_id",	unobs.res },
				{ "auth_seq_id",	unobs.seq },
				{ "PDB_ins_code",	unobs.iCode == ' ' ? "" : string{ unobs.iCode } },
				{ "label_asym_id",	asymId },
				{ "label_comp_id",	compId },		// TODO: change to correct comp_id
				{ "label_seq_id",	seqNr > 0 ? to_string(seqNr) : "" }
			});
		}
		else
		{
			for (auto& atom: unobs.atoms)
			{
				getCategory("pdbx_unobs_or_zero_occ_atoms")->emplace({
					{ "id",				to_string(++idAtom) },
					{ "polymer_flag",	isPolymer ? "Y" : "N" },
					{ "occupancy_flag",	1 },
					{ "PDB_model_num",	unobs.modelNr ? unobs.modelNr : 1 },
					{ "auth_asym_id",	unobs.chain },
					{ "auth_comp_id",	unobs.res },
					{ "auth_seq_id",	unobs.seq },
					{ "PDB_ins_code",	unobs.iCode == ' ' ? "" : string{ unobs.iCode } },
					{ "auth_atom_id",	atom },
					{ "label_asym_id",	asymId },
					{ "label_comp_id",	compId },		// TODO: change to correct comp_id
					{ "label_seq_id",	seqNr > 0 ? to_string(seqNr) : "" },
					{ "label_atom_id",	atom }
				});
				
			}
		}
	}
}

void PDBFileParser::ParseSecondaryStructure()
{
	bool firstHelix = true;
	
	while (mRec->is("HELIX "))
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

		string begAsymId, endAsymId;
		int begSeq, endSeq;
		boost::system::error_code ec;
		
		tie(begAsymId, begSeq, ignore) = MapResidue(vC(20), vI(22, 25), vC(26), ec);
		if (not ec)
			tie(endAsymId, endSeq, ignore) = MapResidue(vC(32), vI(34, 37), vC(38), ec);
		
		if (ec)
		{
			if (VERBOSE)
				cerr << "Could not map residue for HELIX " << vI(8, 10) << endl;
		}
		else
		{
			auto cat = getCategory("struct_conf");
			cat->emplace({
				{ "conf_type_id", "HELX_P" },
				{ "id", "HELX_P" + to_string(vI(8, 10)) },
				{ "pdbx_PDB_helix_id", vS(12, 14) },
				{ "beg_label_comp_id", vS(16, 18) },
				{ "beg_label_asym_id", begAsymId },
				{ "beg_label_seq_id", begSeq },
				{ "pdbx_beg_PDB_ins_code", vS(26, 26) },
				{ "end_label_comp_id", vS(28, 30) },
				{ "end_label_asym_id", endAsymId },
				{ "end_label_seq_id", endSeq },
				{ "pdbx_end_PDB_ins_code", vS(38, 38) },
	
				{ "beg_auth_comp_id", vS(16, 18) },
				{ "beg_auth_asym_id", vS(20, 20) },
				{ "beg_auth_seq_id", vI(22, 25) },
				{ "end_auth_comp_id", vS(28, 30) },
				{ "end_auth_asym_id", vS(32, 32) },
				{ "end_auth_seq_id", vI(34, 37) },
	
				{ "pdbx_PDB_helix_class", vS(39, 40) },
				{ "details", vS(41, 70) },
				{ "pdbx_PDB_helix_length", vI(72, 76) }
			});
	
			if (firstHelix)
			{
				cat = getCategory("struct_conf_type");
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
	
	while (mRec->is("SHEET "))
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
		
		string sheetID = ba::trim_copy(vS(12, 14));
		if (sheetsSeen.count(sheetID) == 0)
		{
			sheetsSeen.insert(sheetID);

			rangeID = 1;

			getCategory("struct_sheet")->emplace({
				{ "id", sheetID },
				{ "number_strands", vI(15, 16) },
			});
		}
		
		int sense = vI(39, 40);
		
		if (sense != 0)
		{
			getCategory("struct_sheet_order")->emplace({
				{ "sheet_id", sheetID },
				{ "range_id_1", rangeID },
				{ "range_id_2", rangeID + 1 },
				{ "sense", sense == -1 ? "anti-parallel" : "parallel" }
			});
		}

		string begAsymId, endAsymId;
		int begSeq, endSeq;
		boost::system::error_code ec;

		tie(begAsymId, begSeq, ignore) = MapResidue(vC(22), vI(23, 26), vC(27), ec);
		if (not ec)
			tie(endAsymId, endSeq, ignore) = MapResidue(vC(33), vI(34, 37), vC(38), ec);
		
		if (ec)
		{
			if (VERBOSE)
				cerr << "Dropping SHEET record " << vI(8, 10) << endl;
		}
		else
		{
			getCategory("struct_sheet_range")->emplace({
				{ "sheet_id", sheetID },
				{ "id", vI(8, 10) },
				{ "beg_label_comp_id", vS(18, 20) },
				{ "beg_label_asym_id", begAsymId },
				{ "beg_label_seq_id", begSeq },
				{ "pdbx_beg_PDB_ins_code", vS(27, 27) },
				{ "end_label_comp_id", vS(29, 31) },
				{ "end_label_asym_id", endAsymId },
				{ "end_label_seq_id", endSeq },
				{ "pdbx_end_PDB_ins_code", vS(38, 38) },
				
				{ "beg_auth_comp_id", vS(18, 20) },
				{ "beg_auth_asym_id", vS(22, 22) },
				{ "beg_auth_seq_id", vI(23, 26) },
				{ "end_auth_comp_id", vS(29, 31) },
				{ "end_auth_asym_id", vS(33, 33) },
				{ "end_auth_seq_id", vI(34, 37) },
			});
			
			if (sense != 0 and mRec->mVlen > 34)
			{
				string r1AsymId, r2AsymId;
				int r1Seq, r2Seq;
				boost::system::error_code ec;
				
				tie(r1AsymId, r1Seq, ignore) = MapResidue(vC(65), vI(66, 69), vC(70), ec);
				if (not ec)
					tie(r2AsymId, r2Seq, ignore) = MapResidue(vC(50), vI(51, 54), vC(55), ec);
				
				if (ec)
				{
					if (VERBOSE)
						cerr << "skipping unmatched pdbx_struct_sheet_hbond record" << endl;
				}
				else
					getCategory("pdbx_struct_sheet_hbond")->emplace({
						{ "sheet_id", sheetID },
						{ "range_id_1", rangeID },
						{ "range_id_2", rangeID + 1 },
						{ "range_1_label_atom_id", vS(57, 60) },
						{ "range_1_label_comp_id", vS(61, 63) },
						{ "range_1_label_asym_id", r1AsymId },
						{ "range_1_label_seq_id", r1Seq },
						{ "range_1_PDB_ins_code", vS(70, 70) },
						{ "range_1_auth_atom_id", vS(57, 60) },
						{ "range_1_auth_comp_id", vS(61, 63) },
						{ "range_1_auth_asym_id", vS(65, 65) },
						{ "range_1_auth_seq_id", vI(66, 69) },
		
						{ "range_2_label_atom_id", vS(42, 45) },
						{ "range_2_label_comp_id", vS(46, 48) },
						{ "range_2_label_asym_id", r2AsymId },
						{ "range_2_label_seq_id", r2Seq },
						{ "range_2_PDB_ins_code", vS(55, 55) },
						{ "range_2_auth_atom_id", vS(42, 45) },
						{ "range_2_auth_comp_id", vS(46, 48) },
						{ "range_2_auth_asym_id", vS(50, 50) },
						{ "range_2_auth_seq_id", vI(51, 54) }
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
		auto compound = mmcif::Compound::create(resName);
		if (compound != nullptr)
		{
			auto at = mmcif::AtomTypeTraits(compound->getAtomById(atomID).typeSymbol);
			result = at.isMetal();
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
		if (mRec->is("SSBOND"))
		{
			if (ssBondNr == 0)
			{
				getCategory("struct_conn_type")->emplace({
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
	
			string p1Asym, p2Asym;
			int p1Seq, p2Seq;
			boost::system::error_code ec;
			
			tie(p1Asym, p1Seq, ignore) = MapResidue(vC(16), vI(18, 21), vC(22), ec);
			if (not ec)
				tie(p2Asym, p2Seq, ignore) = MapResidue(vC(30), vI(32, 35), vC(36), ec);
			
			if (ec)
			{
				if (VERBOSE)
					cerr << "Dropping SSBOND " << vI(8, 10) << endl;
				continue;
			}

			vector<char> alt1 = altLocsForAtom(vC(16), vI(18, 21), vC(22), "SG");
			vector<char> alt2 = altLocsForAtom(vC(30), vI(32, 35), vC(36), "SG");
			
			if (alt1.empty())
				alt1.push_back(0);
			if (alt2.empty())
				alt2.push_back(0);
			
			for (auto a1: alt1)
			{
				for (auto a2: alt2)
				{
					getCategory("struct_conn")->emplace({
						{ "id", "disulf" + to_string(++ssBondNr) },
						{ "conn_type_id", "disulf" },
			
						{ "ptnr1_label_asym_id", p1Asym },
						{ "pdbx_ptnr1_label_alt_id", a1 ? string{ a1 } : string() },
						{ "ptnr1_label_comp_id", vS(12, 14) },
						{ "ptnr1_label_seq_id", p1Seq ? to_string(p1Seq) : "." },
						{ "ptnr1_label_atom_id", "SG" },
						{ "ptnr1_symmetry", pdb2cifSymmetry(vS(60, 65)) },
			
						{ "ptnr2_label_asym_id", p2Asym },
						{ "pdbx_ptnr2_label_alt_id", a2 ? string{ a2 } : string() },
						{ "ptnr2_label_comp_id", vS(26, 28) },
						{ "ptnr2_label_seq_id", p2Seq ? to_string(p2Seq) : "." },
						{ "ptnr2_label_atom_id", "SG" },
			
						{ "ptnr1_auth_asym_id", vS(16, 16) },
						{ "ptnr1_auth_comp_id", vS(12, 14) },
						{ "ptnr1_auth_seq_id", vI(18, 21) },
						{ "ptnr2_auth_asym_id", vS(30, 30) },
						{ "ptnr2_auth_comp_id", vS(26, 28) },
						{ "ptnr2_auth_seq_id", vI(32, 35) },
			
						{ "ptnr2_symmetry", pdb2cifSymmetry(vS(67, 72)) },
						
						{ "pdbx_dist_value", vS(74, 78) },
					});
				}
			}

			continue;
		}
		
		if (mRec->is("LINK  ") or mRec->is("LINKR "))
		{
			if (VERBOSE and mRec->is("LINKR "))
				cerr << "Accepting non-standard LINKR record, but ignoring extra information" << endl;
			
											//	 1 -  6         Record name    "LINK  "
			string name1 = vS(13, 16);		//	13 - 16         Atom           name1           Atom name.
											//	17              Character      altLoc1         Alternate location indicator.
			string resName1 = vS(18,20);	//	18 - 20         Residue name   resName1        Residue  name.
											//	22              Character      chainID1        Chain identifier.
											//	23 - 26         Integer        resSeq1         Residue sequence number.
											//	27              AChar          iCode1          Insertion code.
			string name2 = vS(43, 46);		//	43 - 46         Atom           name2           Atom name.
											//	47              Character      altLoc2         Alternate location indicator.
			string resName2 = vS(48, 50);	//	48 - 50         Residue name   resName2        Residue name.
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
				getCategory("struct_conn_type")->emplace({
					{ "id", type },
				});
				firstCovale = false;
			}
	
			if (type == "metalc" and firstMetalc)
			{
				getCategory("struct_conn_type")->emplace({
					{ "id", type },
				});
				firstMetalc = false;
			}
			
			++linkNr;
	
			string p1Asym, p2Asym;
			int p1Seq, p2Seq;
			bool isResseq1, isResseq2;
			boost::system::error_code ec;
			
			tie(p1Asym, p1Seq, isResseq1) = MapResidue(vC(22), vI(23, 26), vC(27), ec);
			if (not ec)
				tie(p2Asym, p2Seq, isResseq2) = MapResidue(vC(52), vI(53, 56), vC(57), ec);
			
			if (ec)
			{
				if (VERBOSE)
					cerr << "Dropping LINK record at line " << mRec->mLineNr << endl;
				continue;
			}
			
			string distance, details;
			
			if (mRec->is("LINK  "))
			{
				distance = vS(74, 78);
				try
				{
					stod(distance);
				}
				catch (const invalid_argument&)
				{
					if (VERBOSE)
						cerr << "Distance value '" << distance << "' is not a valid float in LINK record" << endl;
					distance.clear();
				}
			}
			else	// LINKR
			{
				details = vS(74, 78);	// the link ID 
			}
			
			getCategory("struct_conn")->emplace({
				{ "id", type + to_string(linkNr) },
				{ "conn_type_id", type },
				
				{ "details", details },
	
				{ "ptnr1_label_asym_id", p1Asym },
				{ "ptnr1_label_comp_id", vS(18, 20) },
				{ "ptnr1_label_seq_id", (isResseq1 and p1Seq) ? to_string(p1Seq) : "." },
				{ "ptnr1_label_atom_id", vS(13, 16) },
				{ "pdbx_ptnr1_label_alt_id", vS(17, 17) },
				{ "pdbx_ptnr1_PDB_ins_code", vS(27, 27) },
				{ "pdbx_ptnr1_standard_comp_id", "" },
				{ "ptnr1_symmetry", pdb2cifSymmetry(vS(60, 65)) },
	
				{ "ptnr2_label_asym_id", p2Asym },
				{ "ptnr2_label_comp_id", vS(48, 50) },
				{ "ptnr2_label_seq_id", (isResseq2 and p2Seq) ? to_string(p2Seq) : "." },
				{ "ptnr2_label_atom_id", vS(43, 46) },
				{ "pdbx_ptnr2_label_alt_id", vS(47, 47) },
				{ "pdbx_ptnr2_PDB_ins_code", vS(57, 57) },
	
				{ "ptnr1_auth_asym_id", vS(22, 22) },
				{ "ptnr1_auth_comp_id", vS(18, 20) },
				{ "ptnr1_auth_seq_id", vI(23, 26) },
				{ "ptnr2_auth_asym_id", vS(52, 52) },
				{ "ptnr2_auth_comp_id", vS(48, 50) },
				{ "ptnr2_auth_seq_id", vI(53, 56) },
	
				{ "ptnr2_symmetry", pdb2cifSymmetry(vS(67, 72)) },
				
				{ "pdbx_dist_value", distance }
			});
			
			continue;
		}
		
		if (mRec->is("CISPEP"))
		{
											//	 1 -  6       Record name   "CISPEP"
			int serNum = vI(8, 10);			//	 8 - 10       Integer       serNum        Record serial number.
			string pep1 = vS(12, 14);		//	12 - 14       LString(3)    pep1          Residue name.
			char chainID1 = vC(16); 		//	16            Character     chainID1      Chain identifier.
			int seqNum1 = vI(18, 21);		//	18 - 21       Integer       seqNum1       Residue sequence number.
			char iCode1 = vC(22);			//	22            AChar         icode1        Insertion code.
			string pep2 = vS(26, 28);		//	26 - 28       LString(3)    pep2          Residue name.
			char chainID2 = vC(30); 		//	30            Character     chainID2      Chain identifier.
			int seqNum2 = vI(32, 35);		//	32 - 35       Integer       seqNum2       Residue sequence number.
			char iCode2 = vC(36);			//	36            AChar         icode2        Insertion code.
			int modNum = vI(44, 46);		//	44 - 46       Integer       modNum        Identifies the specific model.
			string measure = vF(54, 59);	//	54 - 59       Real(6.2)     measure       Angle measurement in degrees.
			
			if (modNum == 0)
				modNum = 1;
	
			string lAsym1, lAsym2;
			int lResSeq1, lResSeq2;
			boost::system::error_code ec;
	
			tie(lAsym1, lResSeq1, ignore) = MapResidue(chainID1, seqNum1, iCode1, ec);
			if (not ec)
				tie(lAsym2, lResSeq2, ignore) = MapResidue(chainID2, seqNum2, iCode2, ec);
			
			if (ec)
			{
				if (VERBOSE)
					cerr << "Dropping CISPEP record at line " << mRec->mLineNr << endl;
				continue;
			}
			
			string iCode1str = iCode1 == ' ' ? string() : string{ iCode1 };
			string iCode2str = iCode2 == ' ' ? string() : string{ iCode2 };
			
			getCategory("struct_mon_prot_cis")->emplace({
				{ "pdbx_id", serNum },
				{ "label_comp_id", pep1 },
				{ "label_seq_id", lResSeq1 },
				{ "label_asym_id", lAsym1 },
				{ "label_alt_id", "." },
				{ "pdbx_PDB_ins_code", iCode1str },
				{ "auth_comp_id", pep1 },
				{ "auth_seq_id", seqNum1 },
				{ "auth_asym_id", string{chainID1} },
				{ "pdbx_label_comp_id_2", pep2 },
				{ "pdbx_label_seq_id_2", lResSeq2 },
				{ "pdbx_label_asym_id_2", lAsym2 },
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
	int structSiteGenId = 1;
	
	while (mRec->is("SITE  "))
	{									//	 1 -  6        Record name   "SITE  "
										//	 8 - 10        Integer       seqNum        Sequence number.
		string siteID = vS(12, 14);	//	12 - 14        LString(3)    siteID        Site name.
		int numRes = vI(16, 17);		//	16 - 17        Integer       numRes        Number of residues that compose the site.

		int o = 19;
		
		auto cat = getCategory("struct_site_gen");

		for (int i = 0; i < numRes; ++i)
		{
			string resName = vS(o, o + 2);	//	19 - 21        Residue name  resName1      Residue name for first residue that 
											//	                                           creates the site.
			char chainID = vC(o + 4);		//	23             Character     chainID1      Chain identifier for first residue of site.
			int seq = vI(o + 5, o + 8);	//	24 - 27        Integer       seq1          Residue sequence number for first residue
											//	                                           of the  site.
			char iCode = vC(o + 9);		//	28             AChar         iCode1        Insertion code for first residue of the site.

			int labelSeq;
			string asym;
			bool isResseq;
			boost::system::error_code ec;
			
			tie(asym, labelSeq, isResseq) = MapResidue(chainID, seq, iCode, ec);
			
			if (ec)
			{
				if (VERBOSE)
					cerr << "skipping struct_site_gen record" << endl;
			}
			else
				cat->emplace({
					{ "id", structSiteGenId++ },
					{ "site_id", siteID },
					{ "pdbx_num_res", numRes },
					{ "label_comp_id", resName },
					{ "label_asym_id", asym },
					{ "label_seq_id", (labelSeq > 0 and isResseq) ? to_string(labelSeq) : string(".") },
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
	Match("CRYST1", true);

	getCategory("cell")->emplace({		
		{ "entry_id", mStructureId },			//	 1 -  6       Record name   "CRYST1"                          
		{ "length_a", vF(7, 15) },             //	 7 - 15       Real(9.3)     a              a (Angstroms).     
		{ "length_b", vF(16, 24) },            //	16 - 24       Real(9.3)     b              b (Angstroms).     
		{ "length_c", vF(25, 33) },            //	25 - 33       Real(9.3)     c              c (Angstroms).     
		{ "angle_alpha", vF(34, 40) },         //	34 - 40       Real(7.2)     alpha          alpha (degrees).   
		{ "angle_beta", vF(41, 47) },          //	41 - 47       Real(7.2)     beta           beta (degrees).    
		{ "angle_gamma", vF(48, 54) },         //	48 - 54       Real(7.2)     gamma          gamma (degrees).   
		/* goes into symmetry */				//	56 - 66       LString       sGroup         Space  group.      
		{ "Z_PDB", vF(67, 70) }                //	67 - 70       Integer       z              Z value.           
	});
	
	string spageGroup, intTablesNr;
	try
	{
		spageGroup = vS(56, 66);
		clipper::Spacegroup sg(clipper::Spgr_descr{spageGroup});
		intTablesNr = to_string(sg.spacegroup_number());
	}
	catch (...)
	{
	}

	getCategory("symmetry")->emplace({
		{ "entry_id", mStructureId },
		{ "space_group_name_H-M", spageGroup },
		{ "Int_Tables_number", intTablesNr }
	});

	GetNextRecord();
}

void PDBFileParser::ParseCoordinateTransformation()
{
	string m[3][3], v[3];
	
	if (ba::starts_with(mRec->mName, "ORIGX"))
	{
		for (string n: { "1", "2", "3" })
		{
			int x = stoi(n) - 1;

			Match("ORIGX" + n, true);	//	 1 -  6         Record name   "ORIGXn"      n=1, 2, or 3  	
			m[x][0] = vF(11, 20);   	//	11 - 20         Real(10.6)    o[n][1]       On1           
			m[x][1] = vF(21, 30);   	//	21 - 30         Real(10.6)    o[n][2]       On2           
			m[x][2] = vF(31, 40);   	//	31 - 40         Real(10.6)    o[n][3]       On3           
			v[x] = vF(46, 55);      	//	46 - 55         Real(10.5)    t[n]          Tn            
			
			GetNextRecord();
		}

		getCategory("database_PDB_matrix")->emplace({
			{ "entry_id", mStructureId },
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

	if (ba::starts_with(mRec->mName, "SCALE"))
	{
		for (string n: { "1", "2", "3" })
		{
			int x = stoi(n) - 1;

			Match("SCALE" + n, true);	//	 1 -  6         Record name   "SCALEn" n=1,  2, or 3      	
			m[x][0] = vF(11, 20);  		//	11 - 20         Real(10.6)    s[n][1]            Sn1      
			m[x][1] = vF(21, 30);	  	//	21 - 30         Real(10.6)    s[n][2]            Sn2      
			m[x][2] = vF(31, 40);  		//	31 - 40         Real(10.6)    s[n][3]            Sn3      
			v[x] = vF(46, 55);     		//	46 - 55         Real(10.5)    u[n]               Un       

			GetNextRecord();
		}

		getCategory("atom_sites")->emplace({
			{ "entry_id", mStructureId },
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

	while (ba::starts_with(mRec->mName, "MTRIX1"))
	{
		int serial, igiven;
		
		for (string n: { "1", "2", "3" })
		{
			int x = stoi(n) - 1;

			Match("MTRIX" + n, true);	//	 1 -  6        Record name   "MTRIXn"      n=1, 2, or 3
			serial = vI(8, 10);			//	 8 - 10        Integer       serial        Serial number.
			m[x][0] = vF(11, 20);  		//	11 - 20        Real(10.6)    m[n][1]       Mn1
			m[x][1] = vF(21, 30);  		//	21 - 30        Real(10.6)    m[n][2]       Mn2
			m[x][2] = vF(31, 40);  		//	31 - 40        Real(10.6)    m[n][3]       Mn3
			v[x] = vF(46, 55);     		//	46 - 55        Real(10.5)    v[n]          Vn
			igiven = vC(60) == '1';		//	60             Integer       iGiven        1 if coordinates for the  representations
											//	                                           which  are approximately related by the 
			GetNextRecord();				//	                                           transformations  of the molecule are
		}									//	                                           contained in the datablock. Otherwise, blank.

		getCategory("struct_ncs_oper")->emplace({
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

void PDBFileParser::ParseCoordinate(int modelNr)
{
	// oh oh, we have to sort our atom_site records by ascending asym_id
	// This routine used to be so trivial...
	
	typedef tuple<string,int,bool,PDBRecord*,PDBRecord*> atomRec;
	
	vector<atomRec> atoms;
	while (mRec->is("ATOM  ") or mRec->is("HETATM"))		//	 1 -  6        Record name   "ATOM  "
	{
		char chainID = vC(22);				//	22             Character     chainID      Chain identifier.
		int resSeq = vI(23, 26);			//	23 - 26        Integer       resSeq       Residue sequence number.
		char iCode = vC(27);
		
		string asymId;
		int seqId;
		bool isResseq;
		
		tie(asymId, seqId, isResseq) = MapResidue(chainID, resSeq, iCode);
		
		PDBRecord* atom = mRec;
		PDBRecord* anisou = nullptr;

		GetNextRecord();
		if (mRec->is("ANISOU"))
		{
			anisou = mRec;
			GetNextRecord();
		}

		atoms.emplace_back(asymId, seqId, isResseq, atom, anisou);

		/*if?... */ while (mRec->is("TER   "))
		{
			Match("TER   ", true);
			GetNextRecord();
		}
	}
	
	auto last = mRec;
	
	// use stable sort here
	auto rLess = [](const atomRec& a, const atomRec& b) -> bool
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
	
	stable_sort(atoms.begin(), atoms.end(), rLess);
	
	// now reiterate the atoms to reorder alternates
	for (size_t i = 0; i + 1 < atoms.size(); ++i)
	{
		char altLoc = get<3>(atoms[i])->vC(17);
		
		if (altLoc == ' ' or altLoc == 0)
			continue;
		
		auto b = atoms.begin() + i;
		auto e = b;
		
		map<string,int> atomIndex;	// index number of first occurrence of a atom name

		while (e != atoms.end() and rLess(*b, *e) == false)
		{
			string name = get<3>(*e)->vS(13, 16);
			
			if (atomIndex.count(name) == 0)
				atomIndex[name] = atomIndex.size() + 1;
			
			++e;
		}
		
		auto aLess = [&](atomRec& a, atomRec& b) -> bool
		{
			string na = get<3>(a)->vS(13, 16);
			string nb = get<3>(b)->vS(13, 16);
			
			int d = atomIndex[na] - atomIndex[nb];
			if (d == 0)
				d = get<3>(a)->vC(17) - get<3>(b)->vC(17);
			assert(d != 0);
			return d < 0;
		};
		
		sort(b, e, aLess);
		
		i += distance(b, e) - 1;
	}

//	while (mRec->is("ATOM  ") or mRec->is("HETATM"))		//	 1 -  6        Record name   "ATOM  "
	for (auto& a: atoms)
	{
		string asymId;
		int seqId;
		bool isResseq;
		PDBRecord* atom;
		PDBRecord* anisou;
		tie(asymId, seqId, isResseq, atom, anisou) = a;
		
		mRec = atom;
		
		++mAtomId;

		string groupPDB = mRec->is("ATOM  ") ? "ATOM" : "HETATM";
//		int serial = vI(7, 11);			//	 7 - 11        Integer       serial       Atom  serial number.
		string name = vS(13, 16);			//	13 - 16        Atom          name         Atom name.
		char altLoc = vC(17);				//	17             Character     altLoc       Alternate location indicator.
		string resName = vS(18, 20);		//	18 - 20        Residue name  resName      Residue name.
		char chainID = vC(22);				//	22             Character     chainID      Chain identifier.
		int resSeq = vI(23, 26);			//	23 - 26        Integer       resSeq       Residue sequence number.
		char iCode = vC(27);				//	27             AChar         iCode        Code for insertion of residues.
		string x = vF(31, 38);				//	31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		string y = vF(39, 46);				//	39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
		string z = vF(47, 54);				//	47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		string occupancy = vF(55, 60);		//	55 - 60        Real(6.2)     occupancy    Occupancy.
		string tempFactor = vF(61, 66);	//	61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		string element = vS(77, 78);		//	77 - 78        LString(2)    element      Element symbol, right-justified.
		string charge = vS(79, 80);		//	79 - 80        LString(2)    charge       Charge  on the atom.

		string entityId = mAsymID2EntityID[asymId];
		
		charge = pdb2cifCharge(charge);

		getCategory("atom_site")->emplace({
			{ "group_PDB" , groupPDB },
			{ "id", mAtomId },
			{ "type_symbol", element },
			{ "label_atom_id", name },
			{ "label_alt_id", altLoc != ' ' ? string { altLoc } : "." },
			{ "label_comp_id", resName },
			{ "label_asym_id", asymId },
			{ "label_entity_id", entityId },
			{ "label_seq_id", (isResseq and seqId > 0) ? to_string(seqId) : "." },
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
			{ "pdbx_PDB_model_num", modelNr }
		});
		
		InsertAtomType(element);
		
		string check = vS(7, 11) + vS(77, 80);
		
		if (anisou != nullptr)
		{
			mRec = anisou;						//	 1 - 6        Record name   "ANISOU"
			int u11 = vI(29, 35);				//	29 - 35       Integer       u[0][0]        U(1,1)
			int u22 = vI(36, 42);				//	36 - 42       Integer       u[1][1]        U(2,2)
			int u33 = vI(43, 49);				//	43 - 49       Integer       u[2][2]        U(3,3)
			int u12 = vI(50, 56);				//	50 - 56       Integer       u[0][1]        U(1,2)
			int u13 = vI(57, 63);				//	57 - 63       Integer       u[0][2]        U(1,3)
			int u23 = vI(64, 70);				//	64 - 70       Integer       u[1][2]        U(2,3)
			
			if (vS(7, 11) + vS(77, 80) != check)
				throw runtime_error("ANISOU record should follow corresponding ATOM record");
			
			auto f = [](float f) -> string { return (boost::format("%6.4f") % f).str(); };
			
			getCategory("atom_site_anisotrop")->emplace({
				{ "id", mAtomId },
				{ "type_symbol", element }, 
				{ "pdbx_label_atom_id", name },
				{ "pdbx_label_alt_id", altLoc != ' ' ? string { altLoc } : "." },
				{ "pdbx_label_comp_id", resName },
				{ "pdbx_label_asym_id", asymId },
				{ "pdbx_label_seq_id", (isResseq and seqId > 0) ? to_string(seqId) : "." },
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
	
	mRec = last;
}

void PDBFileParser::ParseConnectivty()
{
	while (mRec->is("CONECT"))
		GetNextRecord();
}

void PDBFileParser::ParseBookkeeping()
{
	if (mRec->is("MASTER"))
	{
		Match("MASTER", false);
		GetNextRecord();
	}
	Match("END   ", false);
}

void PDBFileParser::Parse(istream& is, cif::File& result)
{
	try
	{
		PreParseInput(is);

		mRec = mData;
	
		ParseTitle();

		result.append(mDatablock);

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
	
		uint32 modelNr = 1;
		bool hasAtoms = false;
	
		while (mRec->is("MODEL ") or mRec->is("ATOM  ") or mRec->is("HETATM"))
		{
			bool model = false;
			if (mRec->is("MODEL "))
			{
				model = true;
				
				modelNr = vI(11, 14);
				
				GetNextRecord();
			}
			
			hasAtoms = hasAtoms or mRec->is("ATOM  ") or mRec->is("HETATM");
			
			ParseCoordinate(modelNr);
			
			if (model)
			{
				Match("ENDMDL", true);
				GetNextRecord();
			}
		}
		
		if (not hasAtoms)
			throw runtime_error("Either the PDB file has no atom records, or the field " + string(mRec->mName) + " is not at the correct location");
	
		for (auto e: mAtomTypes)
			getCategory("atom_type")->emplace({
				{ "symbol", e }
			});

		// in V5, atom_type is sorted
		getCategory("atom_type")->reorderByIndex();
	
		ParseConnectivty();
		ParseBookkeeping();
		
		// almost done, now fix some outstanding issued that could not be done before
		
		try
		{
			auto r = FindRecord("REMARK   3");
			
			if (r != nullptr and Remark3Parser::parse(mExpMethod, r, *mDatablock))
			{
				// make sure the "exptl" category is created
				auto exptl = getCategory("exptl"); 
				if (exptl->empty())
				{
					exptl->emplace({
						{ "entry_id", mStructureId },
						{ "method", mExpMethod },
						{ "crystals_number", mRemark200["NUMBER OF CRYSTALS USED"] }
					});
				}
			}
		}
		catch (const exception& ex)
		{
			cerr << "Error parsing REMARK 3" << endl;
			throw;
		}
//		
//		auto cat = getCategory("pdbx_refine_tls_group");
//		for (Row r: *cat)
//		{
//			// add the mapped locations
//			
//			try
//			{
//				string asymId;
//				int resNum;
//				
//				cif::tie(asymId, resNum) = r.get("beg_auth_asym_id", "beg_auth_seq_id");
//				
//				r["beg_label_asym_id"] = asymId;
//				r["beg_label_seq_id"] = resNum;
//				
//				cif::tie(asymId, resNum) = r.get("end_auth_asym_id", "end_auth_seq_id");
//				
//				r["end_label_asym_id"] = asymId;
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
		cerr << "Error parsing PDB at line " << mRec->mLineNr << endl;
		throw;
	}
}

// ----------------------------------------------------------------
// A blast like alignment. Returns index of last aligned residue.

int PDBFileParser::PDBChain::AlignResToSeqRes()
{
	// Use dynamic programming to align the found residues (in ATOM records) against
	// the residues in the SEQRES records in order to obtain the residue numbering.
	// sigh...
	
	using namespace boost::numeric::ublas;
	
	auto& rx = mSeqres;
	auto& ry = mResiduesSeen;
	
	int dimX = mSeqres.size();
	if (dimX == 0)
		throw runtime_error(string("SEQRES for chain ") + mDbref.chainID + " is empty");

	int dimY = mResiduesSeen.size();
	if (dimY == 0)
		throw runtime_error(string("Number of residues in ATOM records for chain ") + mDbref.chainID + " is zero");

	matrix<float> B(dimX, dimY), Ix(dimX, dimY), Iy(dimX, dimY);
	matrix<int8> tb(dimX, dimY);
	
	int x, y;
	
	const float
		kMatchReward = 5,
		kMismatchCost = -10,
		kGapOpen = 10, gapExtend = 0.1;

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
			if (a.mMonId == b.mMonId)
				M = kMatchReward;
			else
				M = kMismatchCost;

			// gap open cost is zero if the PDB ATOM records indicate that a gap
			// should be here.
			float gapOpen = kGapOpen;
			if (y == 0 or (y + 1 < dimY and ry[y + 1].mSeqNum > ry[y].mSeqNum + 1))
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
		sr.mSeqNum = kFlagSeqNr;
		sr.mIcode = ' ';
	}
	
	// assign numbers
	x = highX;
	y = highY;

	// C++ is getting closer to Pascal :-)
	auto printAlignment = [=]()
	{
		cerr << string(cif::get_terminal_width(), '-') << endl
			 << "Alignment for chain " << mDbref.chainID << endl
			 << endl;
		std::vector<pair<string,string>> alignment;

		int x = highX;
		int y = highY;
		
		for (x = highX, y = highY; x >= 0 and y >= 0; )
		{
			switch (tb(x, y))
			{
				case -1:
					alignment.push_back(make_pair("...", ry[y].mMonId));
					--y;
					break;
				
				case 1:
					alignment.push_back(make_pair(rx[x].mMonId, "..."));
					--x;
					break;
				
				case 0:
					alignment.push_back(make_pair(rx[x].mMonId, ry[y].mMonId));
					--x;
					--y;
					break;
			}
		}
		
		while (x >= 0)
		{
			alignment.push_back(make_pair(rx[x].mMonId, "..."));
			--x;				
		}

		while (y >= 0)
		{
			alignment.push_back(make_pair("...", ry[y].mMonId));
			--y;
		}
			
		reverse(alignment.begin(), alignment.end());
		for (auto a: alignment)
			cerr << "  " << a.first << " -- " << a.second << endl;
		
		cerr << endl;
	};
	
	if (VERBOSE > 1)
		printAlignment();
	
	try
	{
		while (x >= 0 and y >= 0)
		{
			switch (tb(x, y))
			{
				case -1:
//					if (VERBOSE)
//						cerr << "A residue found in the ATOM records "
//							 << "(" << ry[y].mMonId << " @ " << mDbref.chainID << ":" << ry[y].mSeqNum
//							 	<<  ((ry[y].mIcode == ' ' or ry[y].mIcode == 0) ? "" : string{ ry[y].mIcode }) << ")"
//							 << " was not found in the SEQRES records" << endl;

					throw runtime_error("A residue found in the ATOM records (" + ry[y].mMonId + 
						" @ " + string{mDbref.chainID} + ":" + to_string(ry[y].mSeqNum) +
						((ry[y].mIcode == ' ' or ry[y].mIcode == 0) ? "" : string{ ry[y].mIcode })+
						") was not found in the SEQRES records");
					--y;
					break;
				
				case 1:
					if (VERBOSE > 3)
						cerr << "Missing residue in ATOM records: " << rx[x].mMonId << " at " << rx[x].mSeqNum << endl;
	
					--x;
					break;
				
				case 0:
					if (VERBOSE > 3 and rx[x].mMonId != ry[y].mMonId)
						cerr << "Warning, unaligned residues at " << x << "/" << y << "(" << rx[x].mMonId << '/' << ry[y].mMonId << ')' << endl;
					else if (VERBOSE > 4)
						cerr << rx[x].mMonId << " -> " << ry[y].mMonId << " (" << ry[y].mSeqNum << ')' << endl;
	
					rx[x].mSeqNum = ry[y].mSeqNum;
					rx[x].mIcode = ry[y].mIcode;
	
					--x;
					--y;
			}
		}
	}
	catch (const exception& ex)
	{
		if (VERBOSE == 1)
			printAlignment();
		
		throw;
	}

	// assign numbers to the residues that don't have them yet
	stack<int> unnumbered;
	for (x = 0; x < dimX; ++x)
	{
		if (rx[x].mSeqNum == kFlagSeqNr)
		{
			if (x > 0 and rx[x - 1].mSeqNum != kFlagSeqNr)
				rx[x].mSeqNum = rx[x - 1].mSeqNum + 1;
			else
				unnumbered.push(x);
		}
	}
	
	while (unnumbered.empty() == false)
	{
		x = unnumbered.top();
		if (x >= dimX - 1)
			throw runtime_error("Could not assign sequence numbers");
		rx[x].mSeqNum = rx[x + 1].mSeqNum - 1;
		unnumbered.pop();
	}
	
	return highY;
}

bool PDBFileParser::PDBChain::SameSequence(const PDBChain& rhs) const
{
	bool result = mSeqres.size() == rhs.mSeqres.size();
	
	for (size_t i = 0; result and i < mSeqres.size(); ++i)
		result = mSeqres[i].mMonId == rhs.mSeqres[i].mMonId;
	
	return result;
}

// --------------------------------------------------------------------

void ReadPDBFile(istream& pdbFile, cif::File& cifFile)
{
	PDBFileParser p;

	cifFile.loadDictionary("mmcif_pdbx");

	p.Parse(pdbFile, cifFile);
	
	if (not cifFile.isValid())
//		throw runtime_error("Resulting mmCIF file is invalid");
		cerr << "Resulting mmCIF file is not valid!" << endl;
}
