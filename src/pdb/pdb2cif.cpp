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

#include "pdb2cif_remark_3.hpp"

#include "cif++.hpp"

#include <iomanip>
#include <map>
#include <set>
#include <stack>

using cif::category;
using cif::datablock;
using cif::iequals;
using cif::key;
using cif::to_lower;
using cif::to_lower_copy;

// --------------------------------------------------------------------
// attempt to come up with better error handling

namespace error
{
enum pdbErrors
{
	residueNotFound = 1000,
	invalidDate
};

namespace detail
{
	class pdbCategory : public std::error_category
	{
	  public:
		const char *name() const noexcept
		{
			return "pdb";
		}

		std::string message(int value) const
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
} // namespace detail

std::error_category &pdbCategory()
{
	static detail::pdbCategory impl;
	return impl;
}

inline std::error_code make_error_code(pdbErrors e)
{
	return std::error_code(static_cast<int>(e), pdbCategory());
}
} // namespace error

namespace std
{

template <>
struct is_error_code_enum<error::pdbErrors>
{
	static const bool value = true;
};

} // namespace std

namespace cif::pdb
{

// --------------------------------------------------------------------

const std::map<std::string, int> kMonths{
	{ "JAN", 1 },
	{ "FEB", 2 },
	{ "MAR", 3 },
	{ "APR", 4 },
	{ "MAY", 5 },
	{ "JUN", 6 },
	{ "JUL", 7 },
	{ "AUG", 8 },
	{ "SEP", 9 },
	{ "OCT", 10 },
	{ "NOV", 11 },
	{ "DEC", 12 },
};

const std::set<std::string> kSupportedRecords{
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

bool isWater(const std::string &resname)
{
	return resname == "HOH" or resname == "H2O" or resname == "OH2" or resname == "WAT" or resname == "DOD" or resname == "WAT";
}

// --------------------------------------------------------------------
//	Unfortunately, parsing a PDB file requires several passes over the
//	data. Therefore we first obtain all records where a record has the
//	value flattened out for continuation.

PDBRecord::PDBRecord(uint32_t lineNr, const std::string &name, const std::string &value)
	: mNext(nullptr)
	, mLineNr(lineNr)
	, mVlen(value.length())
{
	assert(name.length() <= 10);

	strcpy(mName, name.c_str());
	strcpy(mValue, value.c_str());
}

PDBRecord::~PDBRecord()
{
}

void *PDBRecord::operator new(std::size_t size, std::size_t vLen)
{
	return malloc(size + vLen + 1);
}

void PDBRecord::operator delete(void *p)
{
	free(p);
}

void PDBRecord::operator delete(void *p, std::size_t vLen)
{
	free(p);
}

bool PDBRecord::is(const char *name) const
{
	return iequals(mName, name);
}

char PDBRecord::vC(std::size_t column)
{
	char result = ' ';
	if (column - 7 < mVlen)
		result = mValue[column - 7];
	return result;
}

std::string PDBRecord::vS(std::size_t columnFirst, std::size_t columnLast)
{
	std::string result;

	if (columnLast > mVlen + 6)
		columnLast = mVlen + 6;

	if (columnFirst < mVlen + 7)
	{
		result = std::string{ mValue + columnFirst - 7, mValue + columnLast - 7 + 1 };
		cif::trim(result);
	}

	return result;
}

int PDBRecord::vI(int columnFirst, int columnLast)
{
	int result = 0;

	const char *e = mValue + mVlen;
	if (e > mValue + columnLast - 7 + 1)
		e = mValue + columnLast - 7 + 1;

	enum
	{
		start,
		digit,
		tail
	} state = start;
	bool negate = false;

	try
	{
		for (const char *p = mValue + columnFirst - 7; p < e; ++p)
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
						throw std::runtime_error("Not a valid integer in PDB record");
					break;

				case digit:
					if (isspace(*p))
						state = tail;
					else if (not isdigit(*p))
						throw std::runtime_error("Not a valid integer in PDB record");
					else
						result = result * 10 + *p - '0';
					break;

				case tail:
					if (not isspace(*p))
						throw std::runtime_error("Not a valid integer in PDB record");
					break;
			}
		}
	}
	catch (const std::exception &ex)
	{
		if (cif::VERBOSE >= 0)
			std::cerr << "Trying to parse '" << std::string(mValue + columnFirst - 7, mValue + columnLast - 7) << '\'' << '\n';
		throw;
	}

	if (negate)
		result = -result;

	return result;
}

std::string PDBRecord::vF(std::size_t columnFirst, std::size_t columnLast)
{
	// for now... TODO: check format?
	return vS(columnFirst, columnLast);
}

// --------------------------------------------------------------------

class SpecificationListParser
{
  public:
	SpecificationListParser(const std::string &text)
		: mText(text)
		, mP(mText.begin())
	{
	}

	std::tuple<std::string, std::string> GetNextSpecification();

  private:
	std::string mText;
	std::string::iterator mP;
};

std::tuple<std::string, std::string> SpecificationListParser::GetNextSpecification()
{
	std::string id, value;

	std::string::iterator start = mP, backup;

	enum
	{
		eStart,
		eID,
		eColon,
		eValue,
		eNL,
		eNL_ID,
		eSemiColon,
		eError,
		eDone
	} state = eStart;

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
					if (cif::VERBOSE > 0)
						std::cerr << "skipping invalid character in SOURCE ID: " << ch << '\n';
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
					if (cif::VERBOSE > 0)
						std::cerr << "Empty value for SOURCE: " << id << '\n';
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
				else if (not(isalnum(ch) or ch == '_'))
				{
					value.insert(value.end(), backup, mP);
					state = eValue;
				}
				break;

			case eError:
				if (ch == ';')
				{
					if (cif::VERBOSE > 0)
						std::cerr << "Skipping invalid header line: '" << std::string(start, mP) << '\n';
					state = eStart;
				}
				break;

			case eDone: break; // keep compiler happy
		}
	}

	cif::trim(value);

	return std::make_tuple(id, value);
}

// --------------------------------------------------------------------

class PDBFileParser
{
  public:
	PDBFileParser()
		: mData(nullptr)
		, mRec(nullptr)
	{
	}

	~PDBFileParser()
	{
		PDBRecord *r = mData;
		while (r != nullptr)
		{
			PDBRecord *d = r;
			r = d->mNext;
			delete d;
		}
	}

	void Parse(std::istream &is, cif::file &result);

  private:
	// ----------------------------------------------------------------

	struct DBREF
	{
		std::string PDBIDCode;
		char chainID;
		int seqBegin;
		char insertBegin = ' ';
		int seqEnd;
		char insertEnd = ' ';
		std::string database;
		std::string dbAccession;
		std::string dbIdCode;
		int dbSeqBegin;
		char dbinsBeg;
		int dbSeqEnd;
		char dbinsEnd;
	};

	struct HET
	{
		std::string hetID;
		char chainID;
		int seqNum;
		char iCode;
		int numHetAtoms = 0;
		std::string text;
		std::string asymID;
		std::vector<PDBRecord *> atoms;
		bool processed = false;
		bool branch = false;
		PDBRecord *asn = nullptr;

		HET(const std::string &hetID, char chainID, int seqNum, char iCode, int numHetAtoms = 0, const std::string &text = {})
			: hetID(hetID)
			, chainID(chainID)
			, seqNum(seqNum)
			, iCode(iCode)
			, numHetAtoms(numHetAtoms)
			, text(text)
		{
		}
	};

	struct UNOBS
	{
		int modelNr;
		std::string res;
		char chain;
		int seq;
		char iCode;
		std::vector<std::string> atoms;
	};

	struct ATOM_REF
	{
		std::string name;
		std::string resName;
		int resSeq;
		char chainID;
		char iCode;
		char altLoc;

		bool operator==(const ATOM_REF &rhs) const
		{
			return name == rhs.name and
			       resName == rhs.resName and
			       resSeq == rhs.resSeq and
			       (altLoc == rhs.altLoc or altLoc == ' ' or rhs.altLoc == ' ') and
			       chainID == rhs.chainID and
			       iCode == rhs.iCode;
		}

		bool operator!=(const ATOM_REF &rhs) const
		{
			return not operator==(rhs);
		}

		bool operator<(const ATOM_REF &rhs) const
		{
			int d = chainID - rhs.chainID;
			if (d == 0)
				d = resSeq - rhs.resSeq;
			if (d == 0)
				d = iCode - rhs.iCode;
			// if (d == 0) d = resName.compare(rhs.resName);
			if (d == 0)
				d = name.compare(rhs.name);
			if (d == 0 and altLoc != ' ' and rhs.altLoc != ' ')
				d = altLoc - rhs.altLoc;
			return d < 0;
		}

		friend std::ostream &operator<<(std::ostream &os, const ATOM_REF &a)
		{
			os << a.name << ' ' << a.resName << ' ' << a.chainID << ' ' << a.resSeq << (a.iCode == ' ' ? "" : std::string{ a.iCode }) << (a.altLoc != ' ' ? std::string{ ' ', a.altLoc } : "");
			return os;
		}
	};

	struct LINK
	{
		ATOM_REF a, b;
		std::string symOpA, symOpB;
		float distance;
	};

	struct SUGAR
	{
		ATOM_REF c1;
		int leaving_o;
		ATOM_REF next;
	};

	class SUGAR_TREE : public std::vector<SUGAR>
	{
	  public:
		std::string entityName() const
		{
			return empty() ? "" : entityName(begin());
		}

	  private:
		std::string entityName(const_iterator sugar) const
		{
			std::string result;

			for (auto i = begin(); i != end(); ++i)
			{
				if (i->next != sugar->c1)
					continue;

				auto n = entityName(i) + "-(1-" + std::to_string(i->leaving_o) + ")";

				if (result.empty())
					result = n;
				else
					result += "-[" + n + ']';
			}

			if (not result.empty() and result.back() != ']')
				result += '-';

			auto compound = cif::compound_factory::instance().create(sugar->c1.resName);
			if (compound)
				result += compound->name();
			else if (sugar->c1.resName == "MAN")
				result += "alpha-D-mannopyranose";
			else if (sugar->c1.resName == "BMA")
				result += "beta-D-mannopyranose";
			else if (sugar->c1.resName == "NAG")
				result += "2-acetamido-2-deoxy-beta-D-glucopyranose";
			else if (sugar->c1.resName == "NDG")
				result += "2-acetamido-2-deoxy-alpha-D-glucopyranose";
			else if (sugar->c1.resName == "FUC")
				result += "alpha-L-fucopyranose";
			else if (sugar->c1.resName == "FUL")
				result += "beta-L-fucopyranose";
			else
				result += sugar->c1.resName;

			return result;
		}
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
		int mMolID;
		std::string mTitle;
		std::set<char> mChains;
		std::map<std::string, std::string> mInfo;
		std::map<std::string, std::string> mSource;
		int mCount = 0;
	};

	struct PDBSeqRes
	{
		std::string mMonID;
		int mSeqNum;
		char mIcode;

		int mDbSeqNum = 0;
		bool mSeen = false;
		std::set<std::string> mAlts;

		bool operator==(const PDBSeqRes &rhs) const
		{
			return mSeqNum == rhs.mSeqNum and mMonID == rhs.mMonID and mIcode == rhs.mIcode;
		}
	};

	struct PDBChain
	{
		PDBChain(const std::string &structureID, char chainID, int molID)
			: mDbref{ structureID, chainID }
			, mWaters(0)
			, mTerIndex(0)
			, mMolID(molID)
			, mNextSeqNum(1)
			, mNextDbSeqNum(1)
		{
		}

		DBREF mDbref;
		std::vector<PDBSeqRes> mSeqres, mHet;
		int mWaters;
		int mTerIndex;

		int mMolID;

		// scratch values for reading SEQRES records
		int mNextSeqNum;
		int mNextDbSeqNum;

		// scratch value for aligning
		struct AtomRes
		{
			std::string mMonID;
			int mSeqNum;
			char mIcode;

			bool operator==(const AtomRes &rhs) const { return mSeqNum == rhs.mSeqNum and mIcode == rhs.mIcode; }
			bool operator!=(const AtomRes &rhs) const { return mSeqNum != rhs.mSeqNum or mIcode != rhs.mIcode; }
		};
		std::vector<AtomRes> mResiduesSeen;

		int AlignResToSeqRes();
		bool SameSequence(const PDBChain &rhs) const;
	};

	// ----------------------------------------------------------------

	PDBCompound &GetOrCreateCompound(int molID)
	{
		auto i = std::find_if(mCompounds.begin(), mCompounds.end(), [molID](PDBCompound &comp) -> bool
			{ return comp.mMolID == molID; });
		if (i == mCompounds.end())
		{
			mCompounds.push_back(PDBCompound{ molID });

			mMolID2EntityID[molID] = std::to_string(mNextEntityNr++);

			i = prev(mCompounds.end());
		}

		return *i;
	}

	// locate the PDBChain record for a chain ID, or create it with dummy data if missing
	PDBChain &GetChainForID(char chainID, int numRes = 0)
	{
		auto i = std::find_if(mChains.begin(), mChains.end(), [chainID](PDBChain &ch) -> bool
			{ return ch.mDbref.chainID == chainID; });

		if (i == mChains.end())
		{
			// locate the compound for this chain, if any (does that happen?)
			int molID = 0;
			for (auto &cmp : mCompounds)
			{
				if (cmp.mChains.count(chainID) > 0)
				{
					molID = cmp.mMolID;
					break;
				}
			}

			mChains.emplace_back(mStructureID, chainID, molID);

			i = prev(mChains.end());
		}

		return *i;
	};

	void InsertChemComp(const std::string &chemComp)
	{
		if (find(mChemComp.begin(), mChemComp.end(), chemComp) == mChemComp.end())
			mChemComp.push_back(chemComp);
	}

	void InsertAtomType(const std::string &atomType)
	{
		if (find(mAtomTypes.begin(), mAtomTypes.end(), atomType) == mAtomTypes.end())
			mAtomTypes.push_back(atomType);
	}

	// ----------------------------------------------------------------

	template <typename Predicate>
	PDBRecord *FindRecord(Predicate &&pred)
	{
		PDBRecord *result;

		for (result = mData; result != nullptr; result = result->mNext)
		{
			if (pred(*result))
				break;
		}

		return result;
	}

	PDBRecord *FindRecord(const char *name)
	{
		return FindRecord([name](PDBRecord &rec) -> bool
			{ return rec.is(name); });
	}

	// ----------------------------------------------------------------

	char vC(std::size_t column) const
	{
		return mRec->vC(column);
	}

	std::string vS(std::size_t columnFirst, std::size_t columnLast = std::numeric_limits<std::size_t>::max()) const
	{
		return mRec->vS(columnFirst, columnLast);
	}

	std::string vF(std::size_t columnFirst, std::size_t columnLast) const
	{
		return mRec->vF(columnFirst, columnLast);
	}

	int vI(int columnFirst, int columnLast) const
	{
		return mRec->vI(columnFirst, columnLast);
	}

	// ----------------------------------------------------------------

	// Map a PDB residue location to a seqnum in a struct_asym
	std::tuple<std::string, int, bool> MapResidue(char chainID, int resSeq, char iCode) const
	{
		auto key = std::make_tuple(chainID, resSeq, iCode);

		try
		{
			return mChainSeq2AsymSeq.at(key);
		}
		catch (const std::exception &ex)
		{
			throw_with_nested(std::runtime_error(std::string("Residue ") + chainID + std::to_string(resSeq) + iCode + " could not be mapped"));
		}
	}

	std::tuple<std::string, int, bool> MapResidue(char chainID, int resSeq, char iCode, std::error_code &ec) const
	{
		auto key = std::make_tuple(chainID, resSeq, iCode);

		std::tuple<std::string, int, bool> result;

		if (not mChainSeq2AsymSeq.count(key))
		{
			ec = error::make_error_code(error::pdbErrors::residueNotFound);
			if (cif::VERBOSE > 0)
				std::cerr << "Residue " << chainID << resSeq << iCode << " could not be mapped\n";
		}
		else
			result = mChainSeq2AsymSeq.at(key);

		return result;
	}

	// ----------------------------------------------------------------

	void PreParseInput(std::istream &is);

	void GetNextRecord();
	void Match(const std::string &expected, bool throwIfMissing);

	void ParseTitle();
	void ParseCitation(const std::string &id);
	void ParseRemarks();

	//	void ParseRemark3();
	//	std::size_t ParseRemark3(const std::string& program, const Remark3Template templ[], std::size_t N);
	//	std::string NextRemark3Line();

	void ParseRemark200();
	void ParseRemark350();

	void ParsePrimaryStructure();
	void ParseHeterogen();
	void ConstructEntities();
	void ConstructSugarTrees(int &asymNr);
	void ParseSecondaryStructure();
	void ParseConnectivtyAnnotation();
	void ParseMiscellaneousFeatures();
	void ParseCrystallographic();
	void ParseCoordinateTransformation();
	void ParseCoordinate(int modelNr);
	void ParseConnectivty();
	void ParseBookkeeping();

	// ----------------------------------------------------------------

	category *getCategory(std::string name)
	{
		return &mDatablock[name];
	}

	std::vector<std::string> SplitCSV(const std::string &value);

	std::string pdb2cifDate(std::string s, std::error_code &ec)
	{
		std::smatch m;
		const std::regex
			rx1(R"((\d{2})-(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)-(\d{2}))"),
			rx2(R"((JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)-(\d{2}))");

		try
		{
			if (regex_match(s, m, rx1))
			{
				int day = stoi(m[1].str());
				auto mi = kMonths.find(m[2].str());
				if (mi == kMonths.end())
					throw std::runtime_error("Invalid month: '" + m[2].str() + '\'');
				int month = mi->second;
				int year = 1900 + stoi(m[3].str());
				if (year < 1950)
					year += 100;

				std::stringstream ss;
				ss << std::setw(4) << std::setfill('0') << year << '-'
				   << std::setw(2) << std::setfill('0') << month << '-'
				   << std::setw(2) << std::setfill('0') << day;

				s = ss.str();
			}
			else if (regex_match(s, m, rx2))
			{
				auto mi = kMonths.find(m[1].str());
				if (mi == kMonths.end())
					throw std::runtime_error("Invalid month: '" + m[1].str() + '\'');
				int month = mi->second;
				int year = 1900 + stoi(m[2].str());
				if (year < 1950)
					year += 100;

				s = cif::format("%04d-%02d", year, month).str();
			}
			else
				ec = error::make_error_code(error::pdbErrors::invalidDate);
		}
		catch (const std::exception &ex)
		{
			if (cif::VERBOSE > 0)
				std::cerr << ex.what() << '\n';
			ec = error::make_error_code(error::pdbErrors::invalidDate);
		}

		return s;
	}

	std::string pdb2cifDate(std::string s)
	{
		std::error_code ec;
		auto result = pdb2cifDate(s, ec);
		if (ec and cif::VERBOSE > 0)
			std::cerr << "Invalid date(" << s << "): " << ec.message() << '\n';
		return result;
	}

	std::string pdb2cifAuth(std::string author)
	{
		cif::trim(author);

		const std::regex rx(R"(((?:[A-Z]+\.)+)(.+))");
		std::smatch m;
		if (regex_match(author, m, rx))
			author = m[2].str() + ", " + m[1].str();

		bool upper = true;
		for (auto &c : author)
		{
			if (ispunct(c) or isspace(c))
				upper = true;
			else if (upper)
				upper = false;
			else
				c = cif::tolower(c);
		}

		return author;
	}

	std::string pdb2cifSymmetry(std::string s)
	{
		static const std::regex sgRx(R"((\d{1,3})(\d{3}))");

		if (not s.empty())
		{
			std::smatch m;
			if (not std::regex_match(s, m, sgRx))
				throw std::runtime_error("invalid symmetry value '" + s + '\'');

			s = m[1].str() + "_" + m[2].str();
		}

		return s;
	}

	std::string pdb2cifCharge(std::string c)
	{
		std::regex rx(R"((\d+)(\+|-))");
		std::smatch m;

		if (std::regex_match(c, m, rx))
		{
			if (m[2].str() == "-")
				c = '-' + m[1].str();
			else
				c = m[1].str();
		}

		return c;
	}

	std::vector<char> altLocsForAtom(char chainID, int seqNum, char iCode, std::string atomName);
	void MapChainID2AsymIDS(char chainID, std::vector<std::string> &asymIds);

	std::tuple<ATOM_REF, bool> FindLink(const std::string &name1, const std::string &resName1, int resSeq1, char altLoc1, char chainID1, char iCode1,
		const std::string &name2, const std::string &resName2 = "")
	{
		return FindLink(ATOM_REF{ name1, resName1, resSeq1, altLoc1, chainID1, iCode1 }, name2, resName2);
	}

	std::tuple<ATOM_REF, bool> FindLink(const ATOM_REF &atom, const std::string &name2, const std::string &resName2 = "") const
	{
		auto i = std::find_if(mLinks.begin(), mLinks.end(), [&](const LINK &link)
			{ return (link.a == atom and link.b.name == name2 and (resName2.empty() or link.b.resName == resName2)) or
			         (link.b == atom and link.a.name == name2 and (resName2.empty() or link.a.resName == resName2)); });

		if (i != mLinks.end())
			return { i->a == atom ? i->b : i->a, true };

		return {};
	}

	// ----------------------------------------------------------------

	PDBRecord *mData;
	PDBRecord *mRec;
	cif::datablock mDatablock;

	std::string mStructureID;
	std::string mModelTypeDetails;
	std::string mOriginalDate;
	std::string mExpMethod = "X-RAY DIFFRACTION";
	int mCitationAuthorNr = 1, mCitationEditorNr = 1;
	int mNextMolID = 1, mNextEntityNr = 1;
	int mNextSoftwareOrd = 1;

	struct SEQADV
	{
		std::string resName;
		char chainID;
		int seqNum;
		char iCode;
		std::string database;
		std::string dbAccession;
		std::string dbRes;
		int dbSeq;
		std::string conflict;
	};

	std::vector<SEQADV> mSeqadvs;

	std::list<PDBCompound> mCompounds;
	std::list<PDBChain> mChains;
	std::vector<HET> mHets;
	std::map<std::string, std::string> mHetnams;
	std::map<std::string, std::string> mHetsyns;
	std::map<std::string, std::string> mFormuls;
	std::string mWaterHetID;
	std::vector<std::string> mChemComp, mAtomTypes;

	std::map<std::string, std::string> mRemark200;
	std::string mRefinementSoftware;
	int mAtomID = 0;
	int mPdbxDifOrdinal = 0;

	std::vector<UNOBS> mUnobs;
	std::vector<LINK> mLinks;

	// various maps between numbering schemes
	std::map<std::tuple<char, int, char>, std::tuple<std::string, int, bool>> mChainSeq2AsymSeq;

	std::map<int, std::string> mMolID2EntityID;
	std::map<std::string, std::string> mHet2EntityID;
	std::map<std::string, std::string> mBranch2EntityID;
	std::map<std::string, std::string> mAsymID2EntityID;
	std::map<std::string, std::string> mMod2parent;
	std::set<std::string> mSugarEntities;
};

// --------------------------------------------------------------------

std::vector<char> PDBFileParser::altLocsForAtom(char inChainID, int inResSeq, char inICode, std::string inAtomName)
{
	// well, maybe this could be optimized...
	std::set<char> result;

	for (auto r = mData; r != nullptr; r = r->mNext)
	{
		if (r->is("ATOM  ") or r->is("HETATM")) //	 1 -  6        Record name   "ATOM  "
		{                                       //	 ...
			std::string name = r->vS(13, 16);   //	13 - 16        Atom          name         Atom name.
			char altLoc = r->vC(17);            //	17             Character     altLoc       Alternate location indicator.
			char chainID = r->vC(22);           //	22             Character     chainID      Chain identifier.
			int resSeq = r->vI(23, 26);         //	23 - 26        Integer       resSeq       Residue sequence number.
			char iCode = r->vC(27);             //	27             AChar         iCode        Code for insertion of residues.

			if (chainID == inChainID and resSeq == inResSeq and iCode == inICode and name == inAtomName and altLoc != ' ')
				result.insert(altLoc);
		}
	}

	return { result.begin(), result.end() };
}

void PDBFileParser::MapChainID2AsymIDS(char chainID, std::vector<std::string> &asymIds)
{
	for (const auto &[key, value] : mChainSeq2AsymSeq)
	{
		if (std::get<0>(key) == chainID)
			asymIds.push_back(std::get<0>(value));
	}

	std::sort(asymIds.begin(), asymIds.end(), [](const std::string &a, const std::string &b)
		{
			int d = static_cast<int>(a.length() - b.length());
			if (d == 0)
				d = a.compare(b);
			return d < 0; });

	asymIds.erase(std::unique(asymIds.begin(), asymIds.end()), asymIds.end());
}

// --------------------------------------------------------------------

void PDBFileParser::PreParseInput(std::istream &is)
{
	std::string lookahead;
	uint32_t lineNr = 1;
	getline(is, lookahead);

	if (lookahead.back() == '\r')
		lookahead.pop_back();

	auto contNr = [&lookahead](int offset, int len) -> int
	{
		std::string cs = lookahead.substr(offset, len);
		cif::trim(cs);
		int result = 0;

		if (not cs.empty())
		{
			auto r = std::from_chars(cs.data(), cs.data() + cs.length(), result);
			if ((bool)r.ec)
				throw std::runtime_error("Continuation std::string '" + cs + "' is not valid");
		}

		return result;
	};

	PDBRecord *last = nullptr;
	std::set<std::string> dropped;

	for (;;)
	{
		if (lookahead.empty())
		{
			if (is.eof())
				break;

			if (cif::VERBOSE > 0)
				std::cerr << "Line number " << lineNr << " is empty!\n";

			getline(is, lookahead);
			++lineNr;

			continue;
		}

		std::string type = lookahead.substr(0, 6);
		std::string value;
		if (lookahead.length() > 6)
			value = cif::trim_right_copy(lookahead.substr(6));

		lookahead.clear();

		uint32_t curLineNr = lineNr;
		getline(is, lookahead);
		++lineNr;

		if (kSupportedRecords.count(type) == 0)
		{
			cif::trim(type);

			if (type != "END") // special case
				dropped.insert(type);

			lookahead.clear();

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
				value += cif::trim_right_copy(lookahead.substr(10));
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
				value += cif::trim_right_copy(lookahead.substr(10));
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
				value += cif::trim_right_copy(lookahead.substr(13));
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
				value += cif::trim_copy(lookahead.substr(10));
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
				catch (const std::exception &ex)
				{
					if (cif::VERBOSE >= 0)
						std::cerr << "Dropping FORMUL line (" << (lineNr - 1) << ") with invalid component number '" << value.substr(1, 3) << '\'' << '\n';
					continue;
					// throw_with_nested(std::runtime_error("Invalid component number '" + value.substr(1, 3) + '\''));
				}

				int n = 2;
				try
				{
					while (lookahead.substr(0, 6) == type and
						   stoi(lookahead.substr(7, 3)) == compNr and
						   contNr(16, 2) == n)
					{
						value += cif::trim_right_copy(lookahead.substr(19));
						;
						getline(is, lookahead);
						++lineNr;
						++n;
					}
				}
				catch (const std::invalid_argument &ex)
				{
					continue;
					// throw_with_nested(std::runtime_error("Invalid component number '" + lookahead.substr(7, 3) + '\''));
				}
			}
			catch (const std::exception &ex)
			{
				if (cif::VERBOSE >= 0)
					std::cerr << "Error parsing FORMUL at line " << lineNr << '\n';
				throw;
			}
		}
		else if (type == "HETNAM" or
				 type == "HETSYN")
		{
			int n = 2;
			while (lookahead.substr(0, 6) == type and contNr(8, 2) == n)
			{
				value += cif::trim_right_copy(lookahead.substr(16));
				;
				getline(is, lookahead);
				++lineNr;
				++n;
			}
		}
		else if (type == "SITE  ")
		{
			std::string siteName = value.substr(5, 3);
			cif::trim_right(value);
			std::size_t n = value.length() - 12;
			value += std::string(11 - (n % 11), ' ');

			while (lookahead.substr(0, 6) == type and lookahead.substr(11, 3) == siteName)
			{
				std::string s = lookahead.substr(18);
				cif::trim_right(s);
				s += std::string(11 - (s.length() % 11), ' ');
				value += s;

				// TODO: improve this... either use numRes or don't lump together all text
				//				value += " " + cif::trim_right_copy();
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

				if (i != std::string::npos)
				{
					std::string k = value.substr(4, i - 4);
					std::string v = value.substr(i + 1);

					cif::trim(k);
					while (k.find("  ") != std::string::npos)
						cif::replace_all(k, "  ", " ");
					cif::trim(v);

					if (iequals(v, "NONE") or iequals(v, "N/A") or iequals(v, "NAN"))
						mRemark200[k] = ".";
					else if (not iequals(v, "NULL"))
						mRemark200[k] = v;
				}
			}
		}

		PDBRecord *cur = new (value.length()) PDBRecord(curLineNr, type, value);

		if (last == nullptr)
			last = mData = cur;
		else
			last->mNext = cur;

		last = cur;

		cif::trim(type);

		if (type == "LINK" or type == "LINKR")
		{
			LINK link = {};

			link.a.name = cur->vS(13, 16);    //	13 - 16         Atom           name1           Atom name.
			link.a.altLoc = cur->vC(17);      //	17              Character      altLoc1         Alternate location indicator.
			link.a.resName = cur->vS(18, 20); //	18 - 20         Residue name   resName1        Residue  name.
			link.a.chainID = cur->vC(22);     //	22              Character      chainID1        Chain identifier.
			link.a.resSeq = cur->vI(23, 26);  //	23 - 26         Integer        resSeq1         Residue sequence number.
			link.a.iCode = cur->vC(27);       //	27              AChar          iCode1          Insertion code.
			link.b.name = cur->vS(43, 46);    //	43 - 46         Atom           name2           Atom name.
			link.b.altLoc = cur->vC(47);      //	47              Character      altLoc2         Alternate location indicator.
			link.b.resName = cur->vS(48, 50); //	48 - 50         Residue name   resName2        Residue name.
			link.b.chainID = cur->vC(52);     //	52              Character      chainID2        Chain identifier.
			link.b.resSeq = cur->vI(53, 56);  //	53 - 56         Integer        resSeq2         Residue sequence number.
			link.b.iCode = cur->vC(57);       //	57              AChar          iCode2          Insertion code.
			link.symOpA = cur->vS(60, 65);    //	60 - 65         SymOP          sym1            Symmetry operator atom 1.
			link.symOpB = cur->vS(67, 72);    //	67 - 72         SymOP          sym2            Symmetry operator atom 2.

			if (type == "LINK") //	 1 -  6         Record name    "LINK  "
			{
				auto f = cur->vF(74, 78);
				auto r = cif::from_chars(f.data(), f.data() + f.length(), link.distance);
				if ((bool)r.ec and cif::VERBOSE > 0)
					std::cerr << "Error parsing link distance at line " << cur->mLineNr << '\n';
			}
			//	74 – 78         Real(5.2)      Length          Link distance

			mLinks.push_back(link);
		}

		if (type == "END")
			break;
	}

	if (not dropped.empty())
	{
		if (cif::VERBOSE >= 0)
			std::cerr << "Dropped unsupported records: " << cif::join(dropped, ", ") << '\n';
	}

	if (mData == nullptr)
		throw std::runtime_error("Empty file?");

	mRec = mData;
}

void PDBFileParser::GetNextRecord()
{
	if (mRec != nullptr)
		mRec = mRec->mNext;

	if (mRec == nullptr)
	{
		static PDBRecord *end = new (0) PDBRecord({ 0, "END   ", "" });
		mRec = end;
	}
}

void PDBFileParser::Match(const std::string &expected, bool throwIfMissing)
{
	assert(mRec);
	if (mRec->mName != expected)
	{
		if (throwIfMissing)
			throw std::runtime_error("Expected record " + expected + " but found " + mRec->mName);
		if (cif::VERBOSE > 0)
			std::cerr << "Expected record " + expected + " but found " + mRec->mName << '\n';
	}
}

std::vector<std::string> PDBFileParser::SplitCSV(const std::string &value)
{
	auto vs = cif::split<std::string>(value, ",");
	for (auto &v : vs)
		cif::trim(v);
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

	std::string keywords;

	if (mRec->is("HEADER"))
	{
		mStructureID = vS(63, 66);
		keywords = vS(11, 50);
		mOriginalDate = pdb2cifDate(vS(51, 59));

		cif::trim(keywords);

		GetNextRecord();
	}

	cif::trim(mStructureID);
	if (mStructureID.empty())
		mStructureID = "nohd";

	mDatablock.set_name(mStructureID);

	auto cat = getCategory("entry");
	//	cat->addColumn("id");
	cat->emplace({ { "id", mStructureID } });

	// OBSLTE
	if (mRec->is("OBSLTE"))
	{
		//	 1 -  6       Record name   "OBSLTE"
		//	 9 - 10       Continuation  continuation  Allows concatenation of multiple records
		//	12 - 20       Date          repDate       Date that this datablock was replaced.
		//	22 - 25       IDcode        idCode        ID code of this datablock.
		//	32 - 35       IDcode        rIdCode       ID code of datablock that replaced this one.
		//	37 - 40       ...

		std::string old = vS(22, 25);
		std::string date = pdb2cifDate(vS(12, 20));
		cat = getCategory("pdbx_database_PDB_obs");

		std::string value = mRec->vS(32);
		for (auto i : cif::split<std::string>(value, " ", true))
		{
			cat->emplace({ { "id", "OBSLTE" },
				{ "date", date },
				{ "replace_pdb_id", old },
				{ "pdb_id", i } });
		}

		GetNextRecord();
	}

	// TITLE
	Match("TITLE ", false);
	std::string title;
	if (mRec->is("TITLE ")) //	 1 -  6       Record name    "TITLE "
	{                       //	 9 - 10       Continuation   continuation  Allows concatenation of multiple records.
		title = vS(11);     //	11 - 80       String         title         Title of the  experiment.
		GetNextRecord();
	}

	// SPLIT
	if (mRec->is("SPLIT "))
	{
		//	 1 -  6        Record  name  "SPLIT "
		//	 9 - 10        Continuation  continuation  Allows concatenation of multiple records.
		//	12 - 15        IDcode        idCode        ID code of related datablock.

		throw std::runtime_error("SPLIT PDB files are not supported");
	}

	// CAVEAT
	int caveatID = 1;
	while (mRec->is("CAVEAT")) //	  1 - 6       Record name   "CAVEAT"
	{
		// clang-format off
		getCategory("database_PDB_caveat")->emplace({
			{ "id", caveatID++ },
			{ "text", std::string{ mRec->vS(20) } } //	20 - 79       String        comment        Free text giving the reason for the  CAVEAT.
		});
		// clang-format on

		GetNextRecord();
	}

	// COMPND
	Match("COMPND", false);
	//	 1 -  6       Record name     "COMPND"
	//	 8 - 10       Continuation    continuation  Allows concatenation of multiple records.
	//	11 - 80       Specification   compound      Description of the molecular components.
	//	              list

	if (mRec->is("COMPND"))
	{
		std::string value{ mRec->vS(11) };
		if (value.find(':') == std::string::npos)
		{
			// special case for dumb, stripped files
			auto &comp = GetOrCreateCompound(1);
			comp.mInfo["MOLECULE"] = value;
		}
		else
		{
			SpecificationListParser p(value);

			for (;;)
			{
				std::string key, val;
				std::tie(key, val) = p.GetNextSpecification();

				if (key.empty())
					break;

				if (not iequals(key, "MOL_ID") and mCompounds.empty())
				{
					if (cif::VERBOSE > 0)
						std::cerr << "Ignoring invalid COMPND record\n";
					break;
				}

				if (key == "MOL_ID")
				{
					auto &comp = GetOrCreateCompound(stoi(val));
					comp.mTitle = title;
				}
				else if (key == "CHAIN")
				{
					for (auto c : cif::split<std::string>(val, ","))
					{
						cif::trim(c);
						mCompounds.back().mChains.insert(c[0]);
					}
				}
				else
					mCompounds.back().mInfo[key] = val;
			}
		}

		GetNextRecord();
	}

	// SOURCE
	Match("SOURCE", false);

	if (mRec->is("SOURCE"))
	{
		//	 1 -  6      Record name    "SOURCE"
		//	 8 - 10      Continuation   continuation   Allows concatenation of multiple records.
		//	11 - 79      Specification  srcName        Identifies the source of the
		//	             List                          macromolecule in a  token: value format.

		std::map<std::string, std::string> *source = nullptr;

		//		value = { mRec->vS(11) };
		//		for (auto si = ba::make_split_iterator(value, ba::token_finder(ba::is_any_of(";"), ba::token_compress_on)); not si.eof(); ++si)
		//		{
		//			std::string s(si->begin(), si->end());
		//			if (s.empty())
		//				continue;
		//
		//			auto colon = s.find(": ");
		//			if (colon == std::string::npos)
		//			{
		//				if (cif::VERBOSE > 0)
		//					std::cerr << "invalid source field, missing colon (" << s << ')' << '\n';
		//				continue;
		//			}
		SpecificationListParser p(vS(11));

		for (;;)
		{
			std::string key, val;
			std::tie(key, val) = p.GetNextSpecification();

			if (key.empty())
				break;

			if (key == "MOL_ID")
			{
				for (auto &c : mCompounds)
				{
					if (c.mMolID == stoi(val))
					{
						source = &c.mSource;
						break;
					}
				}

				continue;
			}

			if (source == nullptr)
				throw std::runtime_error("At line " + std::to_string(mRec->mLineNr) + ": missing MOL_ID in SOURCE");

			(*source)[key] = val;
		}

		GetNextRecord();
	}

	// KEYWDS
	Match("KEYWDS", false);
	std::string pdbxKeywords;

	if (mRec->is("KEYWDS"))    //	 1 -  6       Record name    "KEYWDS"
	{                          //	 9 - 10       Continuation   continuation  Allows concatenation of records if necessary.
		pdbxKeywords = vS(11); //	11 - 79       List           keywds        Comma-separated list of keywords relevant
		                       //	                                           to the datablock.
		GetNextRecord();
	}

	if (not(keywords.empty() and pdbxKeywords.empty()))
	{
		// clang-format off
		getCategory("struct_keywords")->emplace({
			{ "entry_id", mStructureID },
			{ "pdbx_keywords", keywords },
			{ "text", pdbxKeywords }
		});
		// clang-format on
	}

	// EXPDTA
	Match("EXPDTA", false);
	if (mRec->is("EXPDTA"))
	{
		mExpMethod = vS(11);

		cat = getCategory("exptl");

		auto crystals = cif::split<std::string>(mRemark200["NUMBER OF CRYSTALS USED"], "; ");
		if (crystals.empty())
			crystals.push_back("");
		auto ci = crystals.begin();

		for (auto expMethod : cif::split<std::string>(mExpMethod, ";"))
		{
			cif::trim(expMethod);

			if (expMethod.empty())
				continue;

			// clang-format off
			cat->emplace({
				{ "entry_id", mStructureID },
				{ "method", expMethod },
				{ "crystals_number", ci != crystals.end() ? *ci : "" }
			});
		// clang-format ob
		}

		GetNextRecord();
	}

	// NUMMDL
	if (mRec->is("NUMMDL"))
	{
		if (cif::VERBOSE > 0)
			std::cerr << "skipping unimplemented NUMMDL record\n";
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

		std::string value = { mRec->vS(11) };
		for (auto author : cif::split<std::string>(value, ",", true))
		{
			// clang-format off
			cat->emplace({
				{ "name", pdb2cifAuth(author) },
				{ "pdbx_ordinal", n }
			});
			// clang-format on
			++n;
		}

		GetNextRecord();
	}

	// REVDAT
	bool firstRevDat = true;
	struct RevDat
	{
		int revNum;
		std::string date, dateOriginal, replaces;
		int modType;
		std::vector<std::string> types;

		bool operator<(const RevDat &rhs) const { return revNum < rhs.revNum; }
	};
	std::vector<RevDat> revdats;

	while (mRec->is("REVDAT"))
	{
		//	 1 -  6       Record name    "REVDAT"
		int revNum = vI(8, 10);                     //	 8 - 10       Integer        modNum        Modification number.
		                                            //	11 - 12       Continuation   continuation  Allows concatenation of multiple records.
		std::string date = pdb2cifDate(vS(14, 22)); //	14 - 22       Date           modDate       Date of modification (or release  for
		                                            //	                                           new entries)  in DD-MMM-YY format. This is
		                                            //	                                           not repeated on continued lines.
		std::string modID = vS(24, 27);             //	24 - 27       IDCode         modID         ID code of this datablock. This is not repeated on
		                                            //	                                           continuation lines.
		int modType = vI(32, 32);                   //	32            Integer        modType       An integer identifying the type of
		                                            //	                                           modification. For all  revisions, the
		                                            //	                                           modification type is listed as 1
		std::string detail = vS(40);                //	40 - 45       LString(6)     record        Modification detail.
		                                            //	47 - 52       LString(6)     record        Modification detail.
		                                            //	54 - 59       LString(6)     record        Modification detail.
		                                            //	61 - 66       LString(6)     record        Modification detail.

		revdats.push_back({ revNum, date, modType == 0 ? mOriginalDate : "", modID, modType });

		revdats.back().types = cif::split<std::string>(detail, " ");

		if (firstRevDat)
		{
			// clang-format off
			getCategory("database_2")->emplace({
				{ "database_id", "PDB" },
				{ "database_code", modID }
			});
			// clang-format on
		}

		GetNextRecord();
		firstRevDat = false;
	}

	/*
	This is internal stuff for PDB, don't write it ???
*/
	sort(revdats.begin(), revdats.end());
	for (auto &revdat : revdats)
	{
		// clang-format off
		getCategory("database_PDB_rev")->emplace({
			{ "num", revdat.revNum },
			{ "date", revdat.date },
			{ "date_original", revdat.dateOriginal },
			{ "replaces", revdat.replaces },
			{ "mod_type", revdat.modType }
		});
		// clang-format on

		for (auto &type : revdat.types)
		{
			if (type.empty())
				continue;

			// clang-format off
			getCategory("database_PDB_rev_record")->emplace({
				{ "rev_num", revdat.revNum },
				{ "type", type }
			});
			// clang-format on
		}
	}
	//*/

	// SPRSDE
	if (mRec->is("SPRSDE"))
	{
		if (cif::VERBOSE > 0)
			std::cerr << "skipping unimplemented SPRSDE record\n";
		GetNextRecord();
	}

	// JRNL
	if (mRec->is("JRNL  "))
		ParseCitation("primary");
}

void PDBFileParser::ParseCitation(const std::string &id)
{
	const char *rec = mRec->mName;

	std::string auth, titl, edit, publ, refn, pmid, doi;
	std::string pubname, volume, astm, country, issn, csd;
	std::string pageFirst;
	int year = 0;

	auto extend = [](std::string &s, const std::string &p)
	{
		if (not s.empty())
			s += ' ';
		s += cif::trim_copy(p);
	};

	while (mRec->is(rec) and (id == "primary" or vC(12) == ' '))
	{
		std::string k = vS(13, 16);
		if (k == "AUTH")
			extend(auth, vS(20, 79));
		else if (k == "TITL")
			extend(titl, vS(20, 79));
		else if (k == "EDIT")
			extend(edit, vS(20, 79));
		else if (k == "REF")
		{
			if (pubname.empty())
			{
				extend(pubname, vS(20, 47));
				if (vS(50, 51) == "V.")
					volume = cif::trim_copy(vS(52, 55));
				pageFirst = vS(57, 61);
				year = vI(63, 66);
			}
			else
				extend(pubname, vS(20, 47));
		}
		else if (k == "PUBL")
			extend(publ, vS(20, 70));
		else if (k == "REFN")
		{
			if (vS(20, 23) == "ASTN")
				astm = vS(25, 30);
			country = vS(33, 34);
			if (vS(36, 39) == "ISSN")
				issn = vS(41, 65);
		}
		else if (k == "PMID")
			pmid = vS(20, 79);
		else if (k == "DOI")
			doi = vS(20, 79);

		GetNextRecord();
	}

	auto cat = getCategory("citation");
	// clang-format off
	cat->emplace({
		{ "id", id },
		{ "title", titl },
		{ "journal_abbrev", pubname },
		{ "journal_volume", volume },
		{ "page_first", pageFirst },
		{ "year", year > 0 ? std::to_string(year) : "" },
		{ "journal_id_ASTM", astm },
		{ "country", country },
		{ "journal_id_ISSN", issn },
		{ "journal_id_CSD", csd },
		{ "book_publisher", publ },
		{ "pdbx_database_id_PubMed", pmid },
		{ "pdbx_database_id_DOI", doi }
	});
	// clang-format on

	if (not auth.empty())
	{
		cat = getCategory("citation_author");
		for (auto author : cif::split<std::string>(auth, ",", true))
		{
			cat->emplace({ { "citation_id", id },
				{ "name", pdb2cifAuth(author) },
				{ "ordinal", mCitationAuthorNr } });

			++mCitationAuthorNr;
		}
	}

	if (not edit.empty())
	{
		cat = getCategory("citation_editor");
		for (auto editor : cif::split<std::string>(edit, ",", true))
		{
			cat->emplace({ { "citation_id", id },
				{ "name", pdb2cifAuth(editor) },
				{ "ordinal", mCitationEditorNr } });

			++mCitationEditorNr;
		}
	}
}

void PDBFileParser::ParseRemarks()
{
	std::string sequenceDetails, compoundDetails, sourceDetails;

	while (cif::starts_with(mRec->mName, "REMARK"))
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
							std::string id = vS(22, 70);
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
					const std::regex rx(R"(THE (\S+) ID CODE IS (\S+?)\.?\s*)");
					std::smatch m;
					std::string r = vS(12);

					if (std::regex_match(r, m, rx))
					{
						auto cat = getCategory("database_2");
						cat->emplace({ { "database_id", m[1].str() },
							{ "database_code", m[2].str() } });
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
						std::string r = mRec->vS(12);

						if (cif::starts_with(r, "REMARK: "))
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
					} while (mRec->is("REMARK 200"));
					break;
				}

				case 280:
				{
					std::string density_Matthews, densityPercentSol, conditions;

					const std::regex rx1(R"(SOLVENT CONTENT, VS +\(%\): *(.+))"),
						rx2(R"(MATTHEWS COEFFICIENT, VM \(ANGSTROMS\*\*3/DA\): *(.+))");

					std::smatch m;

					do
					{
						std::string r = vS(12);

						if (conditions.empty())
						{
							if (std::regex_match(r, m, rx1))
								densityPercentSol = m[1].str();
							else if (std::regex_match(r, m, rx2))
								density_Matthews = m[1].str();
							else if (cif::starts_with(r, "CRYSTALLIZATION CONDITIONS: "))
								conditions = r.substr(28);
						}
						else
							conditions = conditions + ' ' + r;

						GetNextRecord();
					} while (mRec->is("REMARK 280"));

					std::string desc = mRemark200["REMARK"];
					if (desc == "NULL")
						desc.clear();

					// clang-format off
					getCategory("exptl_crystal")->emplace({
						{ "id", 1 },
						{ "density_Matthews", iequals(density_Matthews, "NULL") ? "" : density_Matthews },
						{ "density_percent_sol", iequals(densityPercentSol, "NULL") ? "" : densityPercentSol },
						{ "description", desc }
					});
					// clang-format on

					// now try to parse the conditions
					const std::regex rx3(R"(TEMPERATURE +(\d+)K)"), rx4(R"(PH *(?:: *)?(\d+(?:\.\d+)?))") /*, rx5(R"(\b(\d+)C\b)")*/;

					std::string temp, ph, method;

					for (auto s : cif::split<std::string>(conditions, ",", true))
					{
						cif::trim(s);

						if (std::regex_search(s, m, rx3))
							temp = m[1].str();
						if (std::regex_search(s, m, rx4))
							ph = m[1].str();
						if (s.length() < 60 and
							(cif::icontains(s, "drop") or cif::icontains(s, "vapor") or cif::icontains(s, "batch")))
						{
							if (not method.empty())
								method = method + ", " + s;
							else
								method = s;
						}
					}

					if (not(method.empty() and temp.empty() and ph.empty() and (conditions.empty() or conditions == "NULL")))
					{
						// clang-format off
						getCategory("exptl_crystal_grow")->emplace({
							{ "crystal_id", 1 },
							{ "method", method },
							{ "temp", temp },
							{ "pH", ph },
							{ "pdbx_details", conditions }
						});
						// clang-format on
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
					std::stringstream s;
					GetNextRecord();
					if (vS(12) == "COMPOUND")
						GetNextRecord();

					while (mRec->is("REMARK 400"))
					{
						s << vS(12) << '\n';
						GetNextRecord();
					}

					compoundDetails = s.str();
					break;
				}

				case 450:
				{
					std::stringstream s;
					GetNextRecord();
					if (vS(12) == "SOURCE")
						GetNextRecord();

					while (mRec->is("REMARK 450"))
					{
						s << vS(12) << '\n';
						GetNextRecord();
					}

					sourceDetails = s.str();
					break;
				}

				case 465:
				{
					bool headerSeen = false;
					std::regex rx(R"( *MODELS *(\d+)-(\d+))");
					int models[2] = { -1, -1 };

					for (; mRec->is("REMARK 465"); GetNextRecord())
					{
						if (not headerSeen)
						{
							std::string line = vS(12);
							std::smatch m;

							if (std::regex_match(line, m, rx))
							{
								models[0] = std::stoi(m[1].str());
								models[1] = stoi(m[2].str());
							}
							else
								headerSeen = cif::contains(line, "RES C SSSEQI");
							continue;
						}

						if (models[0] == models[1])
							models[0] = models[1] = vI(12, 14);

						std::string res = vS(16, 18);
						char chain = vC(20);
						int seq = vI(22, 26);
						char iCode = vC(27);

						for (int modelNr = models[0]; modelNr <= models[1]; ++modelNr)
							mUnobs.push_back({ modelNr, res, chain, seq, iCode });
					}

					break;
				}

				case 470:
				{
					bool headerSeen = false;
					std::regex rx(R"( *MODELS *(\d+)-(\d+))");
					int models[2] = { -1, -1 };

					for (; mRec->is("REMARK 470"); GetNextRecord())
					{
						if (not headerSeen)
						{
							std::string line = vS(12);
							std::smatch m;

							if (std::regex_match(line, m, rx))
							{
								models[0] = stoi(m[1].str());
								models[1] = stoi(m[2].str());
							}
							else
								headerSeen = cif::contains(line, "RES CSSEQI  ATOMS");
							continue;
						}

						if (models[0] == models[1])
							models[0] = models[1] = vI(12, 14);

						std::string res = vS(16, 18);
						char chain = vC(20);
						int seq = vI(21, 24);
						char iCode = vC(25);

						std::string atomStr = mRec->vS(29);
						auto atoms = cif::split<std::string>(atomStr, " ", true);

						for (int modelNr = models[0]; modelNr <= models[1]; ++modelNr)
							mUnobs.push_back({ modelNr, res, chain, seq, iCode, atoms });
					}

					break;
				}

				case 500:
				{
					GetNextRecord();

					enum State
					{
						eStart,
						eCCinSAU,
						eCC,
						eCBL,
						eCBA,
						eTA,
						eCTg,
						ePG,
						eMCP,
						eChC
					} state = eStart;
					bool headerSeen = false;
					int id = 0;

					for (; mRec->is("REMARK 500"); GetNextRecord())
					{
						std::string line = vS(12);

						if (line == "GEOMETRY AND STEREOCHEMISTRY")
							continue;

						switch (state)
						{
							case eStart:
							{
								if (line.empty() or not cif::starts_with(line, "SUBTOPIC: "))
									continue;

								std::string subtopic = line.substr(10);

								if (subtopic == "CLOSE CONTACTS IN SAME ASYMMETRIC UNIT")
									state = eCCinSAU;
								else if (subtopic == "CLOSE CONTACTS")
									state = eCC;
								else if (subtopic == "COVALENT BOND LENGTHS")
									state = eCBL;
								else if (subtopic == "COVALENT BOND ANGLES")
									state = eCBA;
								else if (subtopic == "TORSION ANGLES")
									state = eTA;
								else if (subtopic == "NON-CIS, NON-TRANS")
									state = eCTg;
								else if (subtopic == "PLANAR GROUPS")
									state = ePG;
								else if (subtopic == "MAIN CHAIN PLANARITY")
									state = eMCP;
								else if (subtopic == "CHIRAL CENTERS")
									state = eChC;
								else if (cif::VERBOSE > 0)
									throw std::runtime_error("Unknown subtopic in REMARK 500: " + subtopic);

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
									std::string atom1 = vS(13, 16);
									std::string res1 = vS(19, 21);
									std::string alt1 = vS(17, 17);
									char chain1 = vC(23);
									int seq1 = vI(25, 29);
									std::string iCode1 = vS(30, 30);

									std::string atom2 = vS(34, 37);
									std::string alt2 = vS(38, 38);
									std::string res2 = vS(40, 42);
									char chain2 = vC(44);
									int seq2 = vI(46, 50);
									std::string iCode2 = vS(51, 51);

									std::string distance = vF(63, 71);

									// clang-format off
									getCategory("pdbx_validate_close_contact")->emplace({
										{ "id", std::to_string(++id) },
										{ "PDB_model_num", 1 },
										{ "auth_atom_id_1", atom1 },
										{ "auth_asym_id_1", std::string{ chain1 } },
										{ "auth_comp_id_1", res1 },
										{ "auth_seq_id_1", seq1 },
										{ "PDB_ins_code_1", iCode1 },
										{ "label_alt_id_1", alt1 },
										{ "auth_atom_id_2", atom2 },
										{ "auth_asym_id_2", std::string{ chain2 } },
										{ "auth_comp_id_2", res2 },
										{ "auth_seq_id_2", seq2 },
										{ "PDB_ins_code_2", iCode2 },
										{ "label_alt_id_2", alt2 },
										{ "dist", distance }
									});
									// clang-format on
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
									std::string atom1 = vS(13, 16);
									std::string res1 = vS(19, 21);
									char chain1 = vC(23);
									int seq1 = vI(25, 29);

									std::string atom2 = vS(34, 37);
									std::string res2 = vS(40, 42);
									char chain2 = vC(44);
									int seq2 = vI(46, 50);

									std::string symop;
									try
									{
										symop = pdb2cifSymmetry(vS(54, 59));
									}
									catch (const std::exception &ex)
									{
										if (cif::VERBOSE > 0)
											std::cerr << "Dropping REMARK 500 at line " << mRec->mLineNr << " due to invalid symmetry operation\n";
										continue;
									}

									std::string distance = vF(63, 71);

									// clang-format off
									getCategory("pdbx_validate_symm_contact")->emplace({
										{ "id", std::to_string(++id) },
										{ "PDB_model_num", 1 },
										{ "auth_atom_id_1", atom1 },
										{ "auth_asym_id_1", std::string{ chain1 } },
										{ "auth_comp_id_1", res1 },
										{ "auth_seq_id_1", seq1 },
//										{ "PDB_ins_code_1", "" },
//										{ "label_alt_id_1", "" },
										{ "site_symmetry_1", "1_555" },
										{ "auth_atom_id_2", atom2 },
										{ "auth_asym_id_2", std::string{ chain2 } },
										{ "auth_comp_id_2", res2 },
										{ "auth_seq_id_2", seq2 },
//										{ "PDB_ins_code_2", "" },
//										{ "label_alt_id_2", "" },
										{ "site_symmetry_2", symop },
										{ "dist", distance }
									});
									// clang-format on
								}
								break;
							}

							case eCBL:
							{
								if (not headerSeen)
								{
									if (cif::starts_with(line, "FORMAT: ") and line != "FORMAT: (10X,I3,1X,2(A3,1X,A1,I4,A1,1X,A4,3X),1X,F6.3)")
										throw std::runtime_error("Unexpected format in REMARK 500");

									headerSeen = line == "M RES CSSEQI ATM1   RES CSSEQI ATM2   DEVIATION";
								}
								else if (line.empty())
									state = eStart;
								else
								{
									int model = vI(11, 13);
									std::string resNam1 = vS(15, 17);
									std::string chainID1{ vC(19) };
									int seqNum1 = vI(20, 23);
									std::string iCode1{ vC(24) };
									std::string alt1 = vS(30, 30);
									std::string atm1 = vS(26, 29);

									std::string resNam2 = vS(33, 35);
									std::string chainID2{ vC(37) };
									int seqNum2 = vI(38, 41);
									std::string iCode2{ vC(42) };
									std::string alt2 = vS(48, 48);
									std::string atm2 = vS(44, 47);

									std::string deviation = vF(51, 57);

									if (iCode1 == " ")
										iCode1.clear();
									if (iCode2 == " ")
										iCode2.clear();

									// clang-format off
									getCategory("pdbx_validate_rmsd_bond")->emplace({
										{ "id", std::to_string(++id) },
										{ "PDB_model_num", model ? model : 1 },
										{ "auth_atom_id_1", atm1 },
										{ "auth_asym_id_1", chainID1 },
										{ "auth_comp_id_1", resNam1 },
										{ "auth_seq_id_1", seqNum1 },
										{ "PDB_ins_code_1", iCode1 },
										{ "label_alt_id_1", alt1 },
										{ "auth_atom_id_2", atm2 },
										{ "auth_asym_id_2", chainID2 },
										{ "auth_comp_id_2", resNam2 },
										{ "auth_seq_id_2", seqNum2 },
										{ "PDB_ins_code_2", iCode2 },
										{ "label_alt_id_2", alt2 },
										{ "bond_deviation", deviation }
									});
									// clang-format on
								}

								break;
							}

							case eCBA:
								if (not headerSeen)
								{
									if (cif::starts_with(line, "FORMAT: ") and line != "FORMAT: (10X,I3,1X,A3,1X,A1,I4,A1,3(1X,A4,2X),12X,F5.1)")
										throw std::runtime_error("Unexpected format in REMARK 500");

									headerSeen = line == "M RES CSSEQI ATM1   ATM2   ATM3";
								}
								else if (line.empty())
									state = eStart;
								else if (vS(64) == "DEGREES")
								{
									int model = vI(11, 13);
									std::string resNam = vS(15, 17);
									std::string chainID{ vC(19) };
									int seqNum = vI(20, 23);
									std::string iCode{ vC(24) };

									if (iCode == " ")
										iCode.clear();

									std::string atoms[3] = { vS(27, 30), vS(34, 37), vS(41, 44) };
									std::string deviation = vF(57, 62);
									if (deviation == "*****")
										deviation.clear();

									// clang-format off
									getCategory("pdbx_validate_rmsd_angle")->emplace({
										{ "id", std::to_string(++id) },
										{ "PDB_model_num", model ? model : 1 },
										{ "auth_atom_id_1", atoms[0] },
										{ "auth_asym_id_1", chainID },
										{ "auth_comp_id_1", resNam },
										{ "auth_seq_id_1", seqNum },
										{ "PDB_ins_code_1", iCode },
										{ "auth_atom_id_2", atoms[1] },
										{ "auth_asym_id_2", chainID },
										{ "auth_comp_id_2", resNam },
										{ "auth_seq_id_2", seqNum },
										{ "PDB_ins_code_2", iCode },
										{ "auth_atom_id_3", atoms[2] },
										{ "auth_asym_id_3", chainID },
										{ "auth_comp_id_3", resNam },
										{ "auth_seq_id_3", seqNum },
										{ "PDB_ins_code_3", iCode },
										{ "angle_deviation", deviation }
									});
									// clang-format on
								}

								break;

							case eTA:
								if (not headerSeen)
								{
									if (cif::starts_with(line, "FORMAT: ") and line != "FORMAT:(10X,I3,1X,A3,1X,A1,I4,A1,4X,F7.2,3X,F7.2)")
										throw std::runtime_error("Unexpected format in REMARK 500");

									headerSeen = line == "M RES CSSEQI        PSI       PHI";
								}
								else if (line.empty())
									state = eStart;
								else
								{
									int model = vI(11, 13);
									std::string resNam = vS(15, 17);
									std::string chainID{ vC(19) };
									int seqNum = vI(20, 23);
									std::string iCode{ vC(24) };

									if (iCode == " ")
										iCode.clear();

									std::string psi = vF(27, 35);
									std::string phi = vF(37, 45);

									// clang-format off
									getCategory("pdbx_validate_torsion")->emplace({
										{ "id", std::to_string(++id) },
										{ "PDB_model_num", model ? model : 1 },
										{ "auth_comp_id", resNam },
										{ "auth_asym_id", chainID },
										{ "auth_seq_id", seqNum },
										{ "PDB_ins_code", iCode },
										{ "phi", phi },
										{ "psi", psi }
									});
									// clang-format on
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

									std::string resNam1 = vS(12, 14);
									std::string chainID1{ vC(16) };
									int seqNum1 = vI(17, 21);
									std::string iCode1{ vC(22) };

									if (iCode1 == " ")
										iCode1.clear();

									std::string resNam2 = vS(27, 29);
									std::string chainID2{ vC(31) };
									int seqNum2 = vI(32, 36);
									std::string iCode2{ vC(37) };

									if (iCode2 == " ")
										iCode2.clear();

									std::string omega = vF(54, 60);

									// clang-format off
									getCategory("pdbx_validate_peptide_omega")->emplace({
										{ "id", std::to_string(++id) },
										{ "PDB_model_num", model ? model : 1 },
										{ "auth_comp_id_1", resNam1 },
										{ "auth_asym_id_1", chainID1 },
										{ "auth_seq_id_1", seqNum1 },
										{ "PDB_ins_code_1", iCode1 },
										{ "auth_comp_id_2", resNam2 },
										{ "auth_asym_id_2", chainID2 },
										{ "auth_seq_id_2", seqNum2 },
										{ "PDB_ins_code_2", iCode2 },
										{ "omega", omega }
									});
									// clang-format on
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
									std::string resNam = vS(15, 17);
									std::string chainID{ vC(19) };
									int seqNum = vI(20, 23);
									std::string iCode{ vC(24) };

									if (iCode == " ")
										iCode.clear();

									std::string rmsd = vF(32, 36);
									std::string type = vS(41);

									// clang-format off
									getCategory("pdbx_validate_planes")->emplace({
										{ "id", std::to_string(++id) },
										{ "PDB_model_num", model ? model : 1 },
										{ "auth_comp_id", resNam },
										{ "auth_asym_id", chainID },
										{ "auth_seq_id", seqNum },
										{ "PDB_ins_code", iCode },
										{ "rmsd", rmsd },
										{ "type", type }
									});
									// clang-format on
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
							std::string line = vS(12);
							headerSeen = cif::contains(line, "RES C SSEQI");
							continue;
						}

						int modelNr = vI(12, 14);
						if (modelNr == 0)
							modelNr = 1;
						std::string res = vS(16, 18);
						char chain = vC(20);
						int seq = vI(22, 25);
						char iCode = vC(26);

						auto compound = cif::compound_factory::instance().create(res);
						if (compound == nullptr)
							continue;

						std::vector<std::string> atoms;
						for (auto atom : compound->atoms())
						{
							if (atom.type_symbol != cif::H)
								atoms.push_back(atom.id);
						}

						mUnobs.push_back({ modelNr, res, chain, seq, iCode, { atoms } });
					}

					break;
				}

				case 800:
				{
					const std::regex rx1(R"(SITE_IDENTIFIER: (.+))"),
						rx2(R"(EVIDENCE_CODE: (.+))"),
						rx3(R"(SITE_DESCRIPTION: (binding site for residue ([[:alnum:]]{1,3}) ([[:alnum:]]) (\d+)|.+))", std::regex_constants::icase);

					std::string id, evidence, desc;
					std::string pdbxAuthAsymID, pdbxAuthCompID, pdbxAuthSeqID, pdbxAuthInsCode;
					std::smatch m;

					enum State
					{
						sStart,
						sID,
						sEvidence,
						sDesc,
						sDesc2
					} state = sStart;

					auto store = [&]()
					{
						// Find the matching SITE record
						auto site = FindRecord([id](PDBRecord &r) -> bool
							{ return r.is("SITE  ") and r.vS(12, 14) == id; });

						if (site == nullptr)
							throw std::runtime_error("Invalid REMARK 800, no SITE record for id " + id);

						// next record, store what we have
						// clang-format off
						getCategory("struct_site")->emplace({
							{ "id", id },
							{ "details", desc },
							{ "pdbx_auth_asym_id", pdbxAuthAsymID },
							{ "pdbx_auth_comp_id", pdbxAuthCompID },
							{ "pdbx_auth_seq_id", pdbxAuthSeqID },
							{ "pdbx_num_residues", site->vI(16, 17) },
							{ "pdbx_evidence_code", evidence }
						});
						// clang-format on
					};

					for (; mRec->is("REMARK 800"); GetNextRecord())
					{
						std::string s = mRec->vS(12);
						if (s.empty())
							continue;

						switch (state)
						{
							case sStart:
								if (s == "SITE")
									state = sID;
								else if (cif::VERBOSE > 0)
									throw std::runtime_error("Invalid REMARK 800 record, expected SITE");
								break;

							case sID:
								if (std::regex_match(s, m, rx1))
								{
									id = m[1].str();
									state = sEvidence;
								}
								else if (cif::VERBOSE > 0)
									throw std::runtime_error("Invalid REMARK 800 record, expected SITE_IDENTIFIER");
								break;

							case sEvidence:
								if (regex_match(s, m, rx2))
								{
									evidence = m[1].str();
									state = sDesc;
								}
								else if (cif::VERBOSE > 0)
									throw std::runtime_error("Invalid REMARK 800 record, expected SITE_IDENTIFIER");
								break;

							case sDesc:
								if (regex_match(s, m, rx3))
								{
									desc = m[1].str();
									pdbxAuthCompID = m[2].str();
									pdbxAuthAsymID = m[3].str();
									pdbxAuthSeqID = m[4].str();

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
					std::stringstream s;
					GetNextRecord();
					if (vS(12) == "SEQUENCE")
						GetNextRecord();

					while (mRec->is("REMARK 999"))
					{
						s << vS(12) << '\n';
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
					std::string skipped = mRec->mName;

					std::stringstream s;

					if (not mRec->vS(11).empty())
						s << mRec->vS(11) << '\n';
					GetNextRecord();

					while (mRec->is(skipped.c_str()))
					{
						s << mRec->vS(11) << '\n';
						GetNextRecord();
					}

					// clang-format off
					getCategory("pdbx_database_remark")->emplace({
						{ "id", remarkNr },
						{ "text", s.str() }
					});
					// clang-format on

					break;
				}
			}
		}
		catch (const std::exception &ex)
		{
			std::throw_with_nested(std::runtime_error("Error parsing REMARK " + std::to_string(remarkNr)));
		}
	}

	if (not(compoundDetails.empty() and sequenceDetails.empty() and sourceDetails.empty()))
	{
		// clang-format off
		getCategory("pdbx_entry_details")->emplace({
			{ "entry_id", mStructureID },
			{ "compound_details", compoundDetails },
			{ "sequence_details", sequenceDetails },
			{ "source_details", sourceDetails }
		});
		// clang-format on
	}

	// store remark 200 info (special case)
	if (not mRemark200.empty())
		ParseRemark200();
}

void PDBFileParser::ParseRemark200()
{
	auto rm200 = [&](const char *name, int diffrnNr) -> std::string
	{
		int nr = 0;
		std::string result;

		for (auto s : cif::split<std::string>(mRemark200[name], ";"))
		{
			if (++nr != diffrnNr)
				continue;

			cif::trim(s);

			if (s == "NULL")
				s.clear();

			result = std::move(s);
			break;
		}

		return result;
	};

	auto inRM200 = [this](std::initializer_list<const char *> s) -> bool
	{
		bool result = false;

		for (auto *n : s)
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
	            { "entry_id", mStructureID },
	            { "pdbx_data_reduction_ii", mRemark200["INTENSITY-INTEGRATION SOFTWARE"] },
	            { "pdbx_data_reduction_ds", mRemark200["DATA SCALING SOFTWARE"] },
	            { "structure_solution", mRemark200["SOFTWARE USED"] },
	            { "structure_refinement", mRefinementSoftware }
	        });
*/

	struct
	{
		const char *a;
		const char *b;
	} kSWMap[] = {
		{ "data reduction", "INTENSITY-INTEGRATION SOFTWARE" },
		{ "data scaling", "DATA SCALING SOFTWARE" },
		{ "phasing", "SOFTWARE USED" },
	};

	for (auto &sw : kSWMap)
	{
		if (mRemark200[sw.b].empty())
			continue;

		// clang-format off
		getCategory("software")->emplace({
			{ "name", mRemark200[sw.b] },
			{ "classification", sw.a },
			{ "version", "." },
			{ "pdbx_ordinal", mNextSoftwareOrd++ }
		});
		// clang-format on
	}

	std::string scatteringType;
	if (mRemark200["EXPERIMENT TYPE"] == "X-RAY DIFFRACTION")
		scatteringType = "x-ray";
	else if (mRemark200["EXPERIMENT TYPE"] == "NEUTRON DIFFRACTION")
		scatteringType = "neutron";

	std::set<std::string> diffrnWaveLengths;

	for (int diffrnNr = 1;; ++diffrnNr)
	{
		std::string ambientTemp = rm200("TEMPERATURE (KELVIN)", diffrnNr);
		if (ambientTemp.empty())
			break;

		if (cif::ends_with(ambientTemp, "K"))
			ambientTemp.erase(ambientTemp.length() - 1, 1);

		// clang-format off
		getCategory("diffrn")->emplace({
			{ "id", diffrnNr },
			{ "ambient_temp", ambientTemp },
//			{ "ambient_temp_details", seqID },
			{ "crystal_id", 1 } });
		// clang-format on

		std::string collectionDate;
		std::error_code ec;
		collectionDate = pdb2cifDate(rm200("DATE OF DATA COLLECTION", diffrnNr), ec);
		if (ec)
		{
			if (cif::VERBOSE > 0)
				std::cerr << ec.message() << " for pdbx_collection_date\n";

			// The date field can become truncated when multiple values are available
			if (diffrnNr != 1)
				collectionDate.clear();
		}

		// clang-format off
		getCategory("diffrn_detector")->emplace({
			{ "diffrn_id", diffrnNr },
			{ "detector", rm200("DETECTOR TYPE", diffrnNr) },
			{ "type", rm200("DETECTOR MANUFACTURER", diffrnNr) },
			{ "pdbx_collection_date", collectionDate },
			{ "details", rm200("OPTICS", diffrnNr) }
		});
		// clang-format on

		if (inRM200({ "MONOCHROMATIC OR LAUE (M/L)", "MONOCHROMATOR", "DIFFRACTION PROTOCOL" }) or not scatteringType.empty())
			// clang-format off
			getCategory("diffrn_radiation")->emplace({
				{ "diffrn_id", diffrnNr },
				{ "wavelength_id", 1 },
				{ "pdbx_monochromatic_or_laue_m_l", rm200("MONOCHROMATIC OR LAUE (M/L)", diffrnNr) },
				{ "monochromator", rm200("MONOCHROMATOR", diffrnNr) },
				{ "pdbx_diffrn_protocol", rm200("DIFFRACTION PROTOCOL", diffrnNr) },
				{ "pdbx_scattering_type", scatteringType }
			});
		// clang-format on

		std::string wl = rm200("WAVELENGTH OR RANGE (A)", diffrnNr);
		auto wavelengths = cif::split<std::string>(wl, ", -", true);

		diffrnWaveLengths.insert(wavelengths.begin(), wavelengths.end());

		std::string source;
		if (rm200("SYNCHROTRON (Y/N)", diffrnNr) == "Y")
		{
			// clang-format off
			getCategory("diffrn_source")->emplace({
				{ "diffrn_id", diffrnNr },
				{ "source", "SYNCHROTRON" },
				{ "type", rm200("RADIATION SOURCE", diffrnNr) + " BEAMLINE " + rm200("BEAMLINE", diffrnNr) },
				{ "pdbx_synchrotron_site", rm200("RADIATION SOURCE", diffrnNr) },
				{ "pdbx_synchrotron_beamline", rm200("BEAMLINE", diffrnNr) },

				{ "pdbx_wavelength", wavelengths.size() == 1 ? wavelengths[0] : "" },
				{ "pdbx_wavelength_list", wavelengths.size() == 1 ? "" : cif::join(wavelengths, ", ") },
			});
			// clang-format on
		}
		else if (inRM200({ "X-RAY GENERATOR MODEL", "RADIATION SOURCE", "BEAMLINE", "WAVELENGTH OR RANGE (A)" }))
		{
			// clang-format off
			getCategory("diffrn_source")->emplace({
				{ "diffrn_id", diffrnNr },
				{ "source", rm200("RADIATION SOURCE", diffrnNr) },
				{ "type", rm200("X-RAY GENERATOR MODEL", diffrnNr) },

				{ "pdbx_wavelength", wavelengths.size() == 1 ? wavelengths[0] : "" },
				{ "pdbx_wavelength_list", wavelengths.size() == 1 ? "" : cif::join(wavelengths, ", ") },
			});
			// clang-format on
		}
	}

	int wavelengthNr = 1;
	for (auto wl : diffrnWaveLengths)
	{
		if (cif::ends_with(wl, "A"))
			wl.erase(wl.length() - 1, 1);

		// clang-format off
		getCategory("diffrn_radiation_wavelength")->emplace({
			{ "id", wavelengthNr++ },
			{ "wavelength", wl.empty() ? "." : wl },
			{ "wt", "1.0" }
		});
		// clang-format on
	}

	if (inRM200({ "METHOD USED TO DETERMINE THE STRUCTURE", "STARTING MODEL" }))
	{
		auto cat = getCategory("refine");
		assert(cat->empty());

		std::string resolution = mRemark200["RESOLUTION RANGE HIGH (A)"];
		if (resolution.empty())
			resolution = ".";

		// clang-format off
		cat->emplace({
			{ "pdbx_method_to_determine_struct", mRemark200["METHOD USED TO DETERMINE THE STRUCTURE"] },
			{ "pdbx_starting_model", mRemark200["STARTING MODEL"] },
			{ "ls_d_res_high", resolution },
			{ "pdbx_diffrn_id", 1 },
			{ "pdbx_refine_id", mExpMethod },
			{ "entry_id", mStructureID } });
		// clang-format on
	}

	if (inRM200({ "REJECTION CRITERIA (SIGMA(I))", "RESOLUTION RANGE HIGH (A)", "RESOLUTION RANGE LOW (A)", "NUMBER OF UNIQUE REFLECTIONS", "COMPLETENESS FOR RANGE (%)", "<I/SIGMA(I)> FOR THE DATA SET", "R MERGE (I)", "R SYM (I)", "DATA REDUNDANCY" }))
	{
		auto cat = getCategory("reflns");
		// clang-format off
		cat->emplace({
			{ "entry_id", mStructureID },
			{ "observed_criterion_sigma_I", mRemark200["REJECTION CRITERIA (SIGMA(I))"] },
			{ "d_resolution_high", mRemark200["RESOLUTION RANGE HIGH (A)"] },
			{ "d_resolution_low", mRemark200["RESOLUTION RANGE LOW (A)"] },
			{ "number_obs", mRemark200["NUMBER OF UNIQUE REFLECTIONS"] },
			{ "percent_possible_obs", mRemark200["COMPLETENESS FOR RANGE (%)"] },
			{ "pdbx_netI_over_sigmaI", mRemark200["<I/SIGMA(I)> FOR THE DATA SET"] },
			{ "pdbx_Rmerge_I_obs", mRemark200["R MERGE (I)"] },
			{ "pdbx_Rsym_value", mRemark200["R SYM (I)"] },
			{ "pdbx_redundancy", mRemark200["DATA REDUNDANCY"] },
			{ "pdbx_ordinal", 1 },
			{ "pdbx_diffrn_id", 1 }
		});
		// clang-format on
	}

	if (inRM200({ "HIGHEST RESOLUTION SHELL, RANGE HIGH (A)" })) // that one field is mandatory...
	{
		// clang-format off
		getCategory("reflns_shell")->emplace({
			{ "d_res_high", mRemark200["HIGHEST RESOLUTION SHELL, RANGE HIGH (A)"] },
			{ "d_res_low", mRemark200["HIGHEST RESOLUTION SHELL, RANGE LOW (A)"] },
			{ "percent_possible_all", mRemark200["COMPLETENESS FOR SHELL (%)"] },
			{ "Rmerge_I_obs", mRemark200["R MERGE FOR SHELL (I)"] },
			{ "pdbx_Rsym_value", mRemark200["R SYM FOR SHELL (I)"] },
			{ "meanI_over_sigI_obs", mRemark200["<I/SIGMA(I)> FOR SHELL"] },
			{ "pdbx_redundancy", mRemark200["DATA REDUNDANCY IN SHELL"] },
			{ "pdbx_ordinal", 1 },
			{ "pdbx_diffrn_id", 1 }
		});
		// clang-format on
	}
	else if (inRM200({ "HIGHEST RESOLUTION SHELL, RANGE LOW (A)", "COMPLETENESS FOR SHELL (%)",
				 "R MERGE FOR SHELL (I)", "R SYM FOR SHELL (I)", "<I/SIGMA(I)> FOR SHELL", "DATA REDUNDANCY IN SHELL" }))
	{
		if (cif::VERBOSE > 0)
			std::cerr << "Not writing reflns_shell record since d_res_high is missing\n";
	}
}

void PDBFileParser::ParseRemark350()
{
	auto saved = mRec;

	enum State
	{
		eStart,
		eInfo,
		eAnd,
		eApply,
		eBioMT
	} state = eStart;

	const std::regex
		kRX1(R"(BIOMOLECULE: (\d+))"),
		kRX2(R"(([^:]+): (.+?)(?: (ANGSTROM\*\*2|KCAL/MOL))?)"),
		kRX8(R"(APPLY THE FOLLOWING TO CHAINS: (.+))"),
		kRX9(R"(AND CHAINS: (.+))"),
		kRX10(R"(BIOMT([123])\s+(\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?))");

	int biomolecule = 0, operID = 0;
	std::vector<std::string> operExpression;
	std::map<std::string, std::string> values;
	std::vector<std::string> asymIdList;
	std::smatch m;
	cif::row_handle genR;

	std::vector<double> mat, vec;

	for (mRec = FindRecord("REMARK 350"); mRec != nullptr and mRec->is("REMARK 350"); GetNextRecord())
	{
		std::string line = vS(11);

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

					std::string value = m[1].str();

					for (auto chain : cif::split<std::string>(value, ", ", true))
					{
						if (chain.empty()) // happens when we have a AND CHAIN line
						{
							state = eAnd;
							break;
						}

						if (chain.length() != 1)
							throw std::runtime_error("Invalid REMARK 350");

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

					std::string value = m[1].str();

					for (auto chain : cif::split<std::string>(value, ", ", true))
					{
						if (chain.empty()) // happens when we have another AND CHAIN line
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
						throw std::runtime_error("Invalid REMARK 350");

					operID = stoi(m[2].str());
					operExpression.push_back(std::to_string(operID));

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
						operID = stoi(m[2].str());
						operExpression.push_back(std::to_string(operID));
					}
					else if (operID != stoi(m[2].str()))
						throw std::runtime_error("Invalid REMARK 350");

					mat.push_back(stod(m[3].str()));
					mat.push_back(stod(m[4].str()));
					mat.push_back(stod(m[5].str()));
					vec.push_back(stod(m[6].str()));

					if (mt == 3)
					{
						if (vec.size() != 3 or mat.size() != 9)
							throw std::runtime_error("Invalid REMARK 350");

						if (operID == 1)
						{
							std::string oligomer = values["AUTHOR DETERMINED BIOLOGICAL UNIT"];
							if (oligomer.empty())
								oligomer = values["SOFTWARE DETERMINED QUATERNARY STRUCTURE"];
							to_lower(oligomer);

							int count = 0;
							std::smatch m2;

							if (std::regex_match(oligomer, m2, std::regex(R"((\d+)-meric)")))
							{
								count = stoi(m2[1].str());
							}
							else if (cif::ends_with(oligomer, "meric"))
							{
								std::string cs = oligomer.substr(0, oligomer.length() - 5);
								if (cs == "mono")
									count = 1;
								else if (cs == "di")
									count = 2;
								else if (cs == "tri")
									count = 3;
								else if (cs == "tetra")
									count = 4;
								else if (cs == "hexa")
									count = 6;
								else if (cs == "octa")
									count = 8;
								else if (cs == "dodeca")
									count = 12;
							}

							std::string details;
							if (values["AUTHOR DETERMINED BIOLOGICAL UNIT"].empty())
							{
								if (not values["SOFTWARE DETERMINED QUATERNARY STRUCTURE"].empty())
									details = "software_defined_assembly";
							}
							else if (values["SOFTWARE DETERMINED QUATERNARY STRUCTURE"].empty())
								details = "author_defined_assembly";
							else
								details = "author_and_software_defined_assembly";

							// clang-format off
							getCategory("pdbx_struct_assembly")->emplace({
								{ "id", biomolecule },
								{ "details", details },
								{ "method_details", values["SOFTWARE USED"] },
								{ "oligomeric_details", oligomer },
								{ "oligomeric_count", count > 0 ? std::to_string(count) : "" }
							});

							auto cat = getCategory("pdbx_struct_assembly_prop");

							if (not values["TOTAL BURIED SURFACE AREA"].empty())
								cat->emplace({
									{ "biol_id", biomolecule },
									{ "type", "ABSA (A^2)" },
									{ "value", values["TOTAL BURIED SURFACE AREA"] }
								});

							if (not values["CHANGE IN SOLVENT FREE ENERGY"].empty())
								cat->emplace({
									{ "biol_id", biomolecule },
									{ "type", "MORE" },
									{ "value", values["CHANGE IN SOLVENT FREE ENERGY"] }
								});

							if (not values["SURFACE AREA OF THE COMPLEX"].empty())
								cat->emplace({
									{ "biol_id", biomolecule },
									{ "type", "SSA (A^2)" },
									{ "value", values["SURFACE AREA OF THE COMPLEX"] }
								});
							// clang-format on

							values.clear();
						}

						std::string type = mat == std::vector<double>{ 1, 0, 0, 0, 1, 0, 0, 0, 1 } and vec == std::vector<double>{ 0, 0, 0 } ? "identity operation" : "crystal symmetry operation";

						// if (type == "identity operation")
						// {

						// }
						// else
						try
						{
							// clang-format off
							getCategory("pdbx_struct_oper_list")->emplace({
								{ "id", operID },
								{ "type", type },
								// { "name", "" },
							    // { "symmetryOperation", "" },
								{ "matrix[1][1]", cif::format("%12.10f", mat[0]).str() },
								{ "matrix[1][2]", cif::format("%12.10f", mat[1]).str() },
								{ "matrix[1][3]", cif::format("%12.10f", mat[2]).str() },
								{ "vector[1]", cif::format("%12.10f", vec[0]).str() },
								{ "matrix[2][1]", cif::format("%12.10f", mat[3]).str() },
								{ "matrix[2][2]", cif::format("%12.10f", mat[4]).str() },
								{ "matrix[2][3]", cif::format("%12.10f", mat[5]).str() },
								{ "vector[2]", cif::format("%12.10f", vec[1]).str() },
								{ "matrix[3][1]", cif::format("%12.10f", mat[6]).str() },
								{ "matrix[3][2]", cif::format("%12.10f", mat[7]).str() },
								{ "matrix[3][3]", cif::format("%12.10f", mat[8]).str() },
								{ "vector[3]", cif::format("%12.10f", vec[2]).str() }
							});
							// clang-format on
						}
						catch (duplicate_key_error &ex)
						{
							// so what?
						}

						mat.clear();
						vec.clear();
					}
				}
				else if (regex_match(line, m, kRX1))
				{
					if (not(vec.empty() and mat.empty()))
						throw std::runtime_error("Invalid REMARK 350");

					// clang-format off
					getCategory("pdbx_struct_assembly_gen")->emplace({
						{ "assembly_id", biomolecule },
						{ "oper_expression", cif::join(operExpression, ",") },
						{ "asym_id_list", cif::join(asymIdList, ",") }
					});
					// clang-format on

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
		// clang-format off
		getCategory("pdbx_struct_assembly_gen")->emplace({
			{ "assembly_id", biomolecule },
			{ "oper_expression", cif::join(operExpression, ",") },
			{ "asym_id_list", cif::join(asymIdList, ",") }
		});
		// clang-format on
	}

	mRec = saved;
}

void PDBFileParser::ParsePrimaryStructure()
{
	// First locate the DBREF record. Might be missing
	DBREF cur = { mStructureID };

	while (cif::starts_with(mRec->mName, "DBREF"))
	{
		if (mRec->is("DBREF ")) //	 1 -  6       Record name   "DBREF "
		{
			cur.PDBIDCode = vS(8, 11);    //	 8 - 11       IDcode        idCode             ID code of this datablock.
			cur.chainID = vC(13);         //	13            Character     chainID            Chain  identifier.
			cur.seqBegin = vI(15, 18);    //	15 - 18       Integer       seqBegin           Initial sequence number of the
			                              //	                                               PDB sequence segment.
			cur.insertBegin = vC(19);     //	19            AChar         insertBegin        Initial  insertion code of the
			                              //	                                               PDB  sequence segment.
			cur.seqEnd = vI(21, 24);      //	21 - 24       Integer       seqEnd             Ending sequence number of the
			                              //	                                               PDB  sequence segment.
			cur.insertEnd = vC(25);       //	25            AChar         insertEnd          Ending insertion code of the
			                              //	                                               PDB  sequence segment.
			cur.database = vS(27, 32);    //	27 - 32       LString       database           Sequence database name.
			cur.dbAccession = vS(34, 41); //	34 - 41       LString       dbAccession        Sequence database accession code.
			cur.dbIdCode = vS(43, 54);    //	43 - 54       LString       dbIdCode           Sequence  database identification code.
			cur.dbSeqBegin = vI(56, 60);  //	56 - 60       Integer       dbseqBegin         Initial sequence number of the
			                              //	                                               database seqment.
			cur.dbinsBeg = vC(61);        //	61            AChar         idbnsBeg           Insertion code of initial residue of the
			                              //	                                               segment, if PDB is the reference.
			cur.dbSeqEnd = vI(63, 67);    //	63 - 67       Integer       dbseqEnd           Ending sequence number of the
			                              //	                                               database segment.
			cur.dbinsEnd = vC(68);        //	68            AChar         dbinsEnd           Insertion code of the ending residue of
			                              //	                                               the segment, if PDB is the reference.
			auto &chain = GetChainForID(cur.chainID);
			chain.mDbref = cur;
		}
		else if (mRec->is("DBREF1")) //	 1 -  6        Record name   "DBREF1"
		{
			cur.PDBIDCode = vS(8, 11); //	 8 - 11       IDcode        idCode             ID code of this datablock.
			cur.chainID = vC(13);      //	13             Character     chainID       Chain identifier.
			cur.seqBegin = vI(15, 18); //	15 - 18        Integer       seqBegin      Initial sequence number of the
			                           //	                                           PDB sequence segment, right justified.
			cur.insertBegin = vC(19);  //	19             AChar         insertBegin   Initial insertion code of the
			                           //	                                           PDB sequence segment.
			cur.seqEnd = vI(21, 24);   //	21 - 24        Integer       seqEnd        Ending sequence number of the
			                           //	                                           PDB sequence segment, right justified.
			cur.insertEnd = vC(25);    //	25             AChar         insertEnd     Ending insertion code of the
			                           //	                                           PDB sequence  segment.
			cur.database = vS(27, 32); //	27 - 32        LString       database      Sequence database name.
			cur.dbIdCode = vS(48, 67); //	48 - 67        LString       dbIdCode      Sequence database identification code,
		}
		else if (mRec->is("DBREF2"))   //	 1 -  6       Record name   "DBREF2"
		{                              //	 8 - 11       IDcode        idCode        ID code of this datablock.
			if (vC(13) != cur.chainID) //	13            Character     chainID       Chain identifier.
				throw std::runtime_error("Chain ID's for DBREF1/DBREF2 records do not match");
			cur.dbAccession = vS(19, 40); //	19 - 40       LString       dbAccession   Sequence database accession code,
			                              //	                                          left justified.
			cur.dbSeqBegin = vI(46, 55);  //	46 - 55       Integer       seqBegin      Initial sequence number of the
			                              //	                                          Database segment, right justified.
			cur.dbSeqEnd = vI(58, 67);    //	58 - 67       Integer       seqEnd        Ending sequence number of the
			                              //	                                          Database segment, right justified.
			auto &chain = GetChainForID(cur.chainID);
			chain.mDbref = cur;
		}

		GetNextRecord();
	}

	// update chains
	for (auto &chain : mChains)
	{
		chain.mNextSeqNum = chain.mDbref.seqBegin;
		chain.mNextDbSeqNum = chain.mDbref.dbSeqBegin;
	}

	while (mRec->is("SEQADV"))
	{ //	 1 -  6        Record name   "SEQADV"
		mSeqadvs.push_back({
			//	 8 - 11        IDcode        idCode        ID  code of this datablock.
			vS(13, 15), //	13 - 15        Residue name  resName       Name of the PDB residue in conflict.
			vC(17),     //	17             Character     chainID       PDB  chain identifier.
			vI(19, 22), //	19 - 22        Integer       seqNum        PDB  sequence number.
			vC(23),     //	23             AChar         iCode         PDB insertion code.
			vS(25, 28), //	25 - 28        LString       database
			vS(30, 38), //	30 - 38        LString       dbAccession   Sequence  database accession number.
			vS(40, 42), //	40 - 42        Residue name  dbRes         Sequence database residue name.
			vI(44, 48), //	44 - 48        Integer       dbSeq         Sequence database sequence number.
			vS(50, 70)  //	50 - 70        LString       conflict      Conflict comment.
		});

		GetNextRecord();
	}

	while (mRec->is("SEQRES"))             //	 1 -  6        Record name    "SEQRES"
	{                                      //	 8 - 10        Integer        serNum       Serial number of the SEQRES record for  the
		                                   //	                                           current  chain. Starts at 1 and increments
		                                   //	                                           by one  each line. Reset to 1 for each chain.
		char chainID = vC(12);             //	12             Character      chainID      Chain identifier. This may be any single
		                                   //	                                           legal  character, including a blank which is
		                                   //	                                           is  used if there is only one chain.
		int numRes = vI(14, 17);           //	14 - 17        Integer        numRes       Number of residues in the chain.
		                                   //	                                           This  value is repeated on every record.
		std::string monomers = vS(20, 70); //	20 - 22        Residue name   resName      Residue name.
		                                   //	 ...

		auto &chain = GetChainForID(chainID, numRes);

		for (auto monID : cif::split<std::string>(monomers, " ", true))
		{
			if (monID.empty())
				continue;

			chain.mSeqres.push_back({ monID, chain.mNextSeqNum++, ' ', chain.mNextDbSeqNum++ });

			InsertChemComp(monID);
		}

		GetNextRecord();
	}

	// First pass over MODRES, only store relevant information required in ConstructEntities
	while (mRec->is("MODRES"))            //	 1 -  6        Record name   "MODRES"
	{                                     //	 8 - 11        IDcode        idCode      ID code of this datablock.
		std::string resName = vS(13, 15); //	13 - 15        Residue name  resName     Residue name used in this datablock.
		                                  //		char chainID		= vC(17);			//	17             Character     chainID     Chain identifier.
		                                  //		int seqNum			= vI(19, 22);		//	19 - 22        Integer       seqNum      Sequence number.
		                                  //		char iCode			= vC(23);			//	23             AChar         iCode       Insertion code.
		std::string stdRes = vS(25, 27);  //	25 - 27        Residue name  stdRes      Standard residue name.
		                                  //		std::string comment		= vS(30, 70);	//	30 - 70        String        comment     Description of the residue modification.

		mMod2parent[resName] = stdRes;

		GetNextRecord();
	}
}

void PDBFileParser::ParseHeterogen()
{
	while (mRec->is("HET   "))
	{                                  //	 1 -  6       Record name   "HET   "
		std::string hetID = vS(8, 10); //	 8 - 10       LString(3)    hetID          Het identifier, right-justified.
		char chainID = vC(13);         //	13            Character     ChainID        Chain  identifier.
		int seqNum = vI(14, 17);       //	14 - 17       Integer       seqNum         Sequence  number.
		char iCode = vC(18);           //	18            AChar         iCode          Insertion  code.
		int numHetAtoms = vI(21, 25);  //	21 - 25       Integer       numHetAtoms    Number of HETATM records for the group
		                               //	                                           present in the datablock.
		std::string text = vS(31, 70); //	31 - 70       String        text           Text describing Het group.

		mHets.emplace_back(hetID, chainID, seqNum, iCode, numHetAtoms, text);

		GetNextRecord();
	}

	for (;;)
	{
		if (mRec->is("HETNAM"))             //	 1 -  6       Record name   "HETNAM"
		{                                   //	 9 - 10       Continuation  continuation    Allows concatenation of multiple records.
			std::string hetID = vS(12, 14); //	12 - 14       LString(3)    hetID           Het identifier, right-justified.
			std::string text = vS(16);      //	16 - 70       String        text            Chemical name.

			mHetnams[hetID] = text;
			InsertChemComp(hetID);

			GetNextRecord();
			continue;
		}

		if (mRec->is("HETSYN"))             //	 1 -  6       Record name   "HETSYN"
		{                                   //	 9 - 10       Continuation  continuation   Allows concatenation of multiple records.
			std::string hetID = vS(12, 14); //	12 - 14       LString(3)    hetID          Het identifier, right-justified.
			std::string syn = vS(16);       //	16 - 70       SList         hetSynonyms    List of synonyms.

			mHetsyns[hetID] = syn;

			GetNextRecord();
			continue;
		}

		break;
	}

	while (mRec->is("FORMUL"))          //	 1 -  6        Record name   "FORMUL"
	{                                   //	 9 - 10        Integer       compNum       Component  number.
		std::string hetID = vS(13, 15); //	13 - 15        LString(3)    hetID         Het identifier.
		                                //	17 - 18        Integer       continuation  Continuation number.
		char waterMark = vC(19);        //	19             Character     asterisk      "*" for water.
		std::string formula = vS(20);   //	20 - 70        String        text          Chemical formula.

		mFormuls[hetID] = formula;

		if (waterMark == '*')
			mWaterHetID = hetID;

		GetNextRecord();
	}
}

void PDBFileParser::ConstructEntities()
{
	// We parsed the Primary Structure and Heterogen sections, if available.
	// But if we didn't parse anything, we need to fake the data based on residues in ATOM records

	// First iterate all ATOM records and store the residues as found in these records
	int modelNr = 1;

	typedef std::map<std::tuple<char, int, char, char>, std::string> CompTypeMap;
	CompTypeMap residuesSeen; // used to validate PDB files...

	for (auto r = mData; r != nullptr; r = r->mNext)
	{
		if (r->is("MODEL "))
		{
			modelNr = r->vI(11, 14);
			if (modelNr != 1)
				break;
			continue;
		}

		if (r->is("ATOM  ") or r->is("HETATM"))  //	 1 -  6        Record name   "ATOM  "
		{                                        //	 ...
			std::string name = r->vS(13, 16);    //	13 - 16        Atom          name         Atom name.
			char altLoc = r->vC(17);             //	17             Character     altLoc       Alternate location indicator.
			std::string resName = r->vS(18, 20); //	18 - 20        Residue name  resName      Residue name.
			char chainID = r->vC(22);            //	22             Character     chainID      Chain identifier.
			int resSeq = r->vI(23, 26);          //	23 - 26        Integer       resSeq       Residue sequence number.
			char iCode = r->vC(27);              //	27             AChar         iCode        Code for insertion of residues.

			// first validate, too sad this is required...
			CompTypeMap::key_type k = std::make_tuple(chainID, resSeq, iCode, altLoc);
			if (residuesSeen.count(k) == 0)
				residuesSeen[k] = resName;
			else if (residuesSeen[k] != resName)
				throw std::runtime_error("inconsistent residue type for " + std::string{ chainID } + std::to_string(resSeq) + iCode + altLoc + "\n" +
										 "  (" + residuesSeen[k] + " != " + resName + ")");

			auto &chain = GetChainForID(chainID);

			PDBChain::AtomRes ar{ resName, resSeq, iCode };

			if ((chain.mResiduesSeen.empty() or chain.mResiduesSeen.back() != ar) and
				cif::compound_factory::instance().is_monomer(resName))
			{
				chain.mResiduesSeen.push_back(ar);
			}

			// now that we're iterating atoms anyway, clean up the mUnobs array
			mUnobs.erase(remove_if(mUnobs.begin(), mUnobs.end(), [=](UNOBS &a)
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

				return result; }),
				mUnobs.end());

			continue;
		}

		if (r->is("TER   "))          //	 1 -  6 	   Record name	 "TER	"
		{                             //	 7 - 11 	   Integer		 serial 		 Serial number.
			                          //	18 - 20 	   Residue name  resName		 Residue name.
			char chainID = r->vC(22); //	22			   Character	 chainID		 Chain identifier.
			                          //	23 - 26 	   Integer		 resSeq 		 Residue sequence number.
			                          //	27			   AChar		 iCode			 Insertion code.
			auto &chain = GetChainForID(chainID);
			if (chain.mTerIndex == 0) // Is this the first TER record? (Refmac writes out multiple TER records...)
				chain.mTerIndex = static_cast<int>(chain.mResiduesSeen.size());
			continue;
		}
	}

	// prune completely empty chains?
	mChains.erase(remove_if(mChains.begin(), mChains.end(), [](auto &chain)
					  { return chain.mResiduesSeen.empty() and chain.mSeqres.empty(); }),
		mChains.end());

	for (auto &chain : mChains)
	{
		if (not(chain.mSeqres.empty() or chain.mResiduesSeen.empty()))
		{
			// seems safe to assume TER record is at the right location...
			// However, some files don't have them at all.
			// When mTerIndex == 0 this is most likely the case. Right?

			if (chain.mTerIndex > 0)
				chain.mResiduesSeen.erase(chain.mResiduesSeen.begin() + chain.mTerIndex, chain.mResiduesSeen.end());

			int lastResidueIndex = chain.AlignResToSeqRes();

			if (lastResidueIndex > 0 and lastResidueIndex + 1 < static_cast<int>(chain.mResiduesSeen.size()))
			{
				auto &r = chain.mResiduesSeen[lastResidueIndex + 1];

				if (cif::VERBOSE > 0)
				{
					std::cerr << "Detected residues that cannot be aligned to SEQRES\n"
							  << "First residue is " << chain.mDbref.chainID << ':' << r.mSeqNum << r.mIcode << '\n';
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
				std::string resName = chain.mResiduesSeen[ix].mMonID;

				if (cif::compound_factory::instance().is_monomer(resName))
					chain.mTerIndex = ix + 1;

				InsertChemComp(resName);
			}

			// And now construct our 'SEQRES'...
			for (int ix = 0; ix < chain.mTerIndex; ++ix)
			{
				auto &ar = chain.mResiduesSeen[ix];
				chain.mSeqres.push_back({ ar.mMonID, ar.mSeqNum, ar.mIcode, ar.mSeqNum, true });
			}
		}
	}

	std::set<char> terminatedChains;
	std::map<char, int> residuePerChainCounter;

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
		{ //	 1 -  6        Record name   "ATOM  "
			// int serial = r->vI(7, 11);			//	 7 - 11        Integer       serial       Atom  serial number.
			//	 ...
			char altLoc = vC(17);                //	17             Character     altLoc       Alternate location indicator.
			std::string resName = r->vS(18, 20); //	18 - 20        Residue name  resName      Residue name.
			char chainID = r->vC(22);            //	22             Character     chainID      Chain identifier.
			int resSeq = r->vI(23, 26);          //	23 - 26        Integer       resSeq       Residue sequence number.
			char iCode = r->vC(27);              //	27             AChar         iCode        Code for insertion of residues.

			auto &chain = GetChainForID(chainID);

			auto i = find(chain.mSeqres.begin(), chain.mSeqres.end(), PDBSeqRes{ resName, resSeq, iCode });

			// might be a hetero
			if (altLoc != ' ' and i == chain.mSeqres.end())
			{
				i = find_if(chain.mSeqres.begin(), chain.mSeqres.end(),
					[resSeq, iCode](const PDBSeqRes &r) -> bool
					{
						return r.mSeqNum == resSeq and r.mIcode == iCode;
					});
			}

			if (i != chain.mSeqres.end())
			{
				i->mSeen = true;
				if (i->mMonID != resName)
					i->mAlts.insert(resName);
			}
			else
			{
				auto &residues = chain.mHet;

				if (residues.empty() or residues.back().mSeqNum != resSeq)
				{
					i = lower_bound(residues.begin(), residues.end(),
						PDBSeqRes{ resName, resSeq, iCode },
						[=](const PDBSeqRes &r1, const PDBSeqRes &r2) -> bool
						{
							return r1.mSeqNum < r2.mSeqNum;
						});

					residues.insert(i, { resName, resSeq, iCode, resSeq, true });

					InsertChemComp(resName);
				}
			}

			int residueCount = (residuePerChainCounter[chainID] += 1);

			// There appears to be a program that writes out HETATM records as ATOM records....
			if (not cif::compound_factory::instance().is_monomer(resName) or
				terminatedChains.count(chainID) or
				(chain.mTerIndex > 0 and residueCount >= chain.mTerIndex))
			{
				if (isWater(resName))
					mWaterHetID = resName;

				auto h = find_if(mHets.begin(), mHets.end(), [=](const HET &het) -> bool
					{ return het.hetID == resName and het.chainID == chainID and
					         het.seqNum == resSeq and het.iCode == iCode; });

				if (h == mHets.end())
				{
					mHets.push_back({ resName, chainID, resSeq, iCode, 0 }); // double perhaps, but that does not care
					h = prev(mHets.end());
				}

				h->atoms.push_back(r);
			}

			continue;
		}

		if (r->is("TER   "))
		{
			char chainID = r->vC(22); //	22             Character     chainID      Chain identifier.
			terminatedChains.insert(chainID);
		}
	}

	// Create missing compounds
	for (auto &chain : mChains)
	{
		if (chain.mMolID != 0 or chain.mSeqres.empty())
			continue;

		// now this chain may contain the same residues as another one
		for (auto &other : mChains)
		{
			if (&other == &chain or other.mMolID == 0)
				continue;

			if (chain.SameSequence(other))
			{
				chain.mMolID = other.mMolID;
				break;
			}
		}

		if (chain.mMolID != 0)
			continue;

		auto &comp = GetOrCreateCompound(mNextMolID++);
		comp.mChains.insert(chain.mDbref.chainID);

		chain.mMolID = comp.mMolID;
	}

	std::set<std::string> structTitle, structDescription;

	// Create poly_scheme and write pdbx_poly_seq_scheme and create mapping table

	auto cat = getCategory("pdbx_poly_seq_scheme");
	int asymNr = 0;
	for (auto &chain : mChains)
	{
		std::string asymID = cif::cif_id_for_number(asymNr++);

		if (mMolID2EntityID.count(chain.mMolID) == 0)
			continue;

		std::string entityID = mMolID2EntityID[chain.mMolID];

		mAsymID2EntityID[asymID] = entityID;

		// clang-format off
		getCategory("struct_asym")->emplace({
			{ "id", asymID },
			{ "pdbx_blank_PDB_chainid_flag", chain.mDbref.chainID == ' ' ? "Y" : "N" },
			// pdbx_modified
			{ "entity_id", entityID },
			// details
		});
		// clang-format on

		int seqNr = 1;
		for (auto &res : chain.mSeqres)
		{
			mChainSeq2AsymSeq[std::make_tuple(chain.mDbref.chainID, res.mSeqNum, res.mIcode)] = std::make_tuple(asymID, seqNr, true);

			std::string seqID = std::to_string(seqNr);
			++seqNr;

			std::set<std::string> monIds = { res.mMonID };
			monIds.insert(res.mAlts.begin(), res.mAlts.end());

			for (std::string monID : monIds)
			{
				std::string authMonID, authSeqNum, authInsCode{ '.' };

				if (res.mSeen)
				{
					authMonID = monID;
					authSeqNum = std::to_string(res.mSeqNum);
					if (res.mIcode != ' ' and res.mIcode != 0)
						authInsCode = std::string{ res.mIcode };

					// clang-format off
					cat->emplace({
						{ "asym_id", asymID },
						{ "entity_id", mMolID2EntityID[chain.mMolID] },
						{ "seq_id", seqID },
						{ "mon_id", monID },
						{ "ndb_seq_num", seqID },
						{ "pdb_seq_num", res.mSeqNum },
						{ "auth_seq_num", authSeqNum },
						{ "pdb_mon_id", authMonID },
						{ "auth_mon_id", authMonID },
						{ "pdb_strand_id", std::string{ chain.mDbref.chainID } },
						{ "pdb_ins_code", authInsCode },
						{ "hetero", res.mAlts.empty() ? "n" : "y" }
					});
					// clang-format on
				}
				else
				{
					if (res.mIcode != ' ' and res.mIcode != 0)
						authInsCode = std::string{ res.mIcode } + "A";

					// clang-format off
					cat->emplace({
						{ "asym_id", asymID },
						{ "entity_id", mMolID2EntityID[chain.mMolID] },
						{ "seq_id", seqID },
						{ "mon_id", monID },
						{ "ndb_seq_num", seqID },
						{ "pdb_seq_num", res.mSeqNum },
						{ "auth_seq_num", "." },
						{ "pdb_mon_id", "." },
						{ "auth_mon_id", "." },
						{ "pdb_strand_id", std::string{ chain.mDbref.chainID } },
						{ "pdb_ins_code", authInsCode },
						{ "hetero", res.mAlts.empty() ? "n" : "y" }
					});
					// clang-format on
				}
			}
		}
	}

	// We have now created all compounds, write them out
	uint32_t structRefID = 0, structRefSeqAlignID = 0;

	for (auto &cmp : mCompounds)
	{
		++structRefID;

		std::string srcMethod;

		if (not cmp.mSource["SYNTHETIC"].empty())
		{
			srcMethod = "syn";

			// clang-format off
			getCategory("pdbx_entity_src_syn")->emplace({
				{ "entity_id", mMolID2EntityID[cmp.mMolID] },
				{ "pdbx_src_id", structRefID },
				{ "organism_scientific", cmp.mSource["ORGANISM_SCIENTIFIC"] },
				{ "ncbi_taxonomy_id", cmp.mSource["ORGANISM_TAXID"] },
			});
			// clang-format on
		}
		else if (cmp.mInfo["ENGINEERED"] == "YES" or
				 not cmp.mSource["EXPRESSION_SYSTEM"].empty())
		{
			srcMethod = "man";

			// clang-format off
			getCategory("entity_src_gen")->emplace({
				{ "entity_id", mMolID2EntityID[cmp.mMolID] },
				{ "pdbx_src_id", structRefID },
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
			// clang-format on
		}
		else if (not cmp.mSource["ORGANISM_SCIENTIFIC"].empty())
		{
			srcMethod = "nat";

			// clang-format off
			getCategory("entity_src_nat")->emplace({
				{ "entity_id", mMolID2EntityID[cmp.mMolID] },
				{ "pdbx_src_id", structRefID },
				{ "common_name", cmp.mSource["ORGANISM_COMMON"] },
				{ "strain", cmp.mSource["STRAIN"] },
				{ "pdbx_secretion", cmp.mSource["SECRETION"] },
				{ "pdbx_organism_scientific", cmp.mSource["ORGANISM_SCIENTIFIC"] },
				{ "pdbx_ncbi_taxonomy_id", cmp.mSource["ORGANISM_TAXID"] },
				{ "pdbx_cellular_location", cmp.mSource["CELLULAR_LOCATION"] },
				{ "pdbx_plasmid_name", cmp.mSource["PLASMID"] },
				{ "pdbx_organ", cmp.mSource["ORGAN"] },
			});
			// clang-format on
		}

		// clang-format off
		getCategory("entity")->emplace({
			{ "id", mMolID2EntityID[cmp.mMolID] },
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
		// clang-format on

		if (not cmp.mInfo["SYNONYM"].empty())
		{
			// clang-format off
			getCategory("entity_name_com")->emplace({
				{ "entity_id", mMolID2EntityID[cmp.mMolID] },
				{ "name", cmp.mInfo["SYNONYM"] }
			});
			// clang-format on
		}

		std::string desc = cmp.mInfo["MOLECULE"];
		if (not cmp.mInfo["EC"].empty())
			desc += " (E.C." + cmp.mInfo["EC"] + ")";

		if (not cmp.mTitle.empty())
			structTitle.insert(cmp.mTitle);

		if (not desc.empty())
			structDescription.insert(desc);

		auto ci = find_if(mChains.begin(), mChains.end(),
			[cmp](PDBChain &c) -> bool
			{ return cmp.mChains.count(c.mDbref.chainID); });

		if (ci != mChains.end() and not ci->mDbref.dbIdCode.empty())
		{
			// clang-format off
			getCategory("struct_ref")->emplace({
				{ "id", structRefID },
				{ "entity_id", mMolID2EntityID[cmp.mMolID] },
				{ "db_name", ci->mDbref.database },
				{ "db_code", ci->mDbref.dbIdCode },
				{ "pdbx_db_accession", ci->mDbref.dbAccession },
//				{ "pdbx_align_begin", ci->mDbref.dbSeqBegin }
			});
			// clang-format on
		}

		bool nstdMonomer = false, nonstandardLinkage = false;
		bool mightBePolyPeptide = true, mightBeDNA = true;

		std::vector<std::string> chains;
		std::string seq, seqCan;

		// write out the chains for this compound
		for (auto &chain : mChains)
		{
			if (chain.mMolID != cmp.mMolID)
				continue;

			//			chain.mEntityID = cmp.mEntityID;

			++structRefSeqAlignID;
			DBREF &dbref = chain.mDbref;

			if (not dbref.database.empty())
			{
				auto insToStr = [](char i) -> std::string
				{
					return i == ' ' or not isprint(i) ? "" : std::string{ i };
				};

				auto &pdbxPolySeqScheme = *getCategory("pdbx_poly_seq_scheme");

				int seqAlignBeg = 0, seqAlignEnd = 0;

				try
				{
					seqAlignBeg = pdbxPolySeqScheme.find1<int>(key("pdb_strand_id") == std::string{ dbref.chainID } and
																   key("pdb_seq_num") == dbref.seqBegin and
																   (key("pdb_ins_code") == insToStr(dbref.insertBegin) or key("pdb_ins_code") == cif::null),
						"seq_id");

					seqAlignEnd = pdbxPolySeqScheme.find1<int>(key("pdb_strand_id") == std::string{ dbref.chainID } and
																   key("pdb_seq_num") == dbref.seqEnd and
																   (key("pdb_ins_code") == insToStr(dbref.insertEnd) or key("pdb_ins_code") == cif::null),
						"seq_id");
				}
				catch (...)
				{
				}

				// clang-format off
				getCategory("struct_ref_seq")->emplace({
					{ "align_id", structRefSeqAlignID },
					{ "ref_id", structRefID },
					{ "pdbx_PDB_id_code", dbref.PDBIDCode },
					{ "pdbx_strand_id", std::string{ chain.mDbref.chainID } },
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
				// clang-format on

				// write the struct_ref_seq_dif
				for (auto &seqadv : mSeqadvs)
				{
					if (seqadv.chainID != chain.mDbref.chainID or seqadv.resName.empty())
						continue;

					std::string asym, seqNum;
					int labelSeq = -1;
					std::error_code ec;

					std::tie(asym, labelSeq, std::ignore) = MapResidue(seqadv.chainID, seqadv.seqNum, seqadv.iCode, ec);
					if (ec)
					{
						if (cif::VERBOSE > 0)
							std::cerr << "dropping unmatched SEQADV record\n";
						continue;
					}

					seqNum = std::to_string(labelSeq);

					// clang-format off
					getCategory("struct_ref_seq_dif")->emplace({
						{ "align_id", structRefSeqAlignID },
						{ "pdbx_PDB_id_code", dbref.PDBIDCode },
						{ "mon_id", seqadv.resName },
						{ "pdbx_pdb_strand_id", seqadv.chainID },
						{ "seq_num", seqNum },
						{ "pdbx_pdb_ins_code", seqadv.iCode == ' ' ? std::string{} : std::string{ seqadv.iCode } },
						{ "pdbx_seq_db_name", seqadv.database },
						{ "pdbx_seq_db_accession_code", seqadv.dbAccession },
						{ "db_mon_id", seqadv.dbRes },
						{ "pdbx_seq_db_seq_num", seqadv.dbSeq },
						{ "details", seqadv.conflict },
						{ "pdbx_auth_seq_num", seqadv.seqNum },
						{ "pdbx_ordinal", ++mPdbxDifOrdinal }
					});
					// clang-format on
				}
			}

			if (not chains.empty()) // not the first one for this molID
			{
				chains.push_back(std::string{ chain.mDbref.chainID });
				continue;
			}

			chains.push_back(std::string{ chain.mDbref.chainID });

			std::size_t seqLen = 0, seqCanLen = 0;

			for (auto &res : chain.mSeqres)
			{
				std::string letter, stdRes;

				if (mMod2parent.count(res.mMonID))
					stdRes = mMod2parent.at(res.mMonID);

				if (cif::compound_factory::kAAMap.count(res.mMonID))
				{
					letter = cif::compound_factory::kAAMap.at(res.mMonID);
					mightBeDNA = false;
				}
				else if (cif::compound_factory::kBaseMap.count(res.mMonID))
				{
					letter = cif::compound_factory::kBaseMap.at(res.mMonID);
					mightBePolyPeptide = false;
				}
				else
				{
					nstdMonomer = true;
					letter = '(' + res.mMonID + ')';

					// sja...
					auto compound = cif::compound_factory::instance().create(stdRes.empty() ? res.mMonID : stdRes);
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
					if (not stdRes.empty() and cif::compound_factory::kAAMap.count(stdRes))
						letter = cif::compound_factory::kAAMap.at(stdRes);
					else if (cif::compound_factory::kBaseMap.count(res.mMonID))
						letter = cif::compound_factory::kBaseMap.at(res.mMonID);
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

			auto cat_ps = getCategory("entity_poly_seq");
			for (std::size_t i = 0; i < chain.mSeqres.size(); ++i)
			{
				auto &rs = chain.mSeqres[i];

				if (std::find(mChemComp.begin(), mChemComp.end(), rs.mMonID) == mChemComp.end())
					mChemComp.emplace_back(rs.mMonID);

				// clang-format off
				cat_ps->emplace({
					{ "entity_id", mMolID2EntityID[cmp.mMolID] },
					{ "num", i + 1 },
					{ "mon_id", rs.mMonID },
					{ "hetero", rs.mAlts.empty() ? "n" : "y" }
				});
				// clang-format on

				for (auto &a : rs.mAlts)
				{
					// clang-format off
					cat_ps->emplace({
						{ "entity_id", mMolID2EntityID[cmp.mMolID] },
						{ "num", i + 1 },
						{ "mon_id", a },
						{ "hetero", "y" }
					});
					// clang-format on
				}
			}
		}

		std::string type;
		if (mightBePolyPeptide and not mightBeDNA)
			type = "polypeptide(L)";
		else if (mightBeDNA and not mightBePolyPeptide)
			type = "polyribonucleotide";

		// clang-format off
		getCategory("entity_poly")->emplace({
			{ "entity_id", mMolID2EntityID[cmp.mMolID] },
			{ "pdbx_seq_one_letter_code", seq },
			{ "pdbx_seq_one_letter_code_can", seqCan },
			{ "nstd_monomer", (nstdMonomer ? "yes" : "no") },
			{ "pdbx_strand_id", cif::join(chains, ",") },
			{ "nstd_linkage", nonstandardLinkage ? "yes" : "no" },
			{ "type", type }
		});
		// clang-format on
	}

	if (not(structTitle.empty() and structDescription.empty()))
	{
		// clang-format off
		getCategory("struct")->emplace({
			{ "entry_id", mStructureID },
			{ "title", cif::join(structTitle, ", ") },
			{ "pdbx_descriptor", cif::join(structDescription, ", ") },
			{ "pdbx_model_type_details", mModelTypeDetails }
		});
		// clang-format on
	}

	// build sugar trees first
	ConstructSugarTrees(asymNr);

	// done with the sugar, resume operation as before

	std::map<char, std::string> waterChains;
	std::map<std::tuple<std::string, std::string>, int> ndbSeqNum; // for nonpoly scheme
	std::map<std::string, int> entityAuthSeqNum;                   // for nonpoly scheme too

	for (std::size_t i = 0; i < mHets.size(); ++i)
	{
		auto &heti = mHets[i];

		if (not heti.asymID.empty())
			continue;

		if (heti.hetID == mWaterHetID or isWater(heti.hetID))
			continue;

		// See if this residue is part of SEQRES
		auto &chain = GetChainForID(heti.chainID);
		auto ih = find(chain.mSeqres.begin(), chain.mSeqres.end(), PDBSeqRes{ heti.hetID, heti.seqNum, heti.iCode });

		// If so, skip it, it is not an entity then
		if (ih != chain.mSeqres.end())
			continue;

		heti.asymID = cif::cif_id_for_number(asymNr++);
	}

	std::set<std::string> writtenAsyms;

	std::map<std::string, int> hetCount; // for pdbx_number_of_molecules
	for (auto &het : mHets)
		hetCount[het.hetID] += 1;

	for (auto &het : mHets)
	{
		std::string hetID = het.hetID;

		auto &chain = GetChainForID(het.chainID);

		// See if this residue is part of SEQRES
		auto i = find(chain.mSeqres.begin(), chain.mSeqres.end(), PDBSeqRes{ hetID, het.seqNum, het.iCode });

		// If so, skip it, it is not an entity then
		if (i != chain.mSeqres.end())
			continue;

		// See if we've already added it to the entities
		if (mHet2EntityID.count(hetID) == 0)
		{
			std::string entityID = std::to_string(mNextEntityNr++);
			mHet2EntityID[hetID] = entityID;

			if (hetID == mWaterHetID)
			{
				// clang-format off
				getCategory("entity")->emplace({
					{ "id", entityID },
					{ "type", "water" },
					{ "src_method", "nat" },
					{ "pdbx_description", "water" },
					{ "pdbx_number_of_molecules", hetCount[hetID] }
				});
				// clang-format on
			}
			else
			{
				if (mHetnams[hetID].empty())
				{
					auto compound = cif::compound_factory::instance().create(hetID);
					if (compound != nullptr)
						mHetnams[hetID] = compound->name();
				}

				// clang-format off
				getCategory("entity")->emplace({
					{ "id", entityID },
					{ "type", "non-polymer" },
					{ "src_method", "syn" },
					{ "pdbx_description", mHetnams[hetID] },
					{ "details", mHetsyns[hetID] },
					{ "pdbx_number_of_molecules", hetCount[hetID] }
				});
				// clang-format on
			}

			// write a pdbx_entity_nonpoly record
			std::string name = mHetnams[hetID];
			if (name.empty() and hetID == mWaterHetID)
				name = "water";

			// clang-format off
			getCategory("pdbx_entity_nonpoly")->emplace({
				{ "entity_id", entityID },
				{ "name", name },
				{ "comp_id", hetID }
			});
			// clang-format on
		}

		// create an asym for this het/chain combo, if needed

		std::string asymID = het.asymID;

		auto k = std::make_tuple(het.chainID, het.seqNum, het.iCode);
		if (mChainSeq2AsymSeq.count(k) == 0)
		{
			if (hetID == mWaterHetID or isWater(hetID))
			{
				if (waterChains.count(het.chainID) == 0)
				{
					asymID = cif::cif_id_for_number(asymNr++);
					waterChains[het.chainID] = asymID;
				}
				else
					asymID = waterChains[het.chainID];
			}
			else
				asymID = het.asymID;

			assert(asymID.empty() == false);

			mAsymID2EntityID[asymID] = mHet2EntityID[hetID];

			// NOTE, a nonpoly residue has no label_seq_id
			// but in pdbx_nonpoly_scheme there is such a number.
			// Since this number is not used anywhere else we
			// just use it here and do not store it in the table
			mChainSeq2AsymSeq[k] = std::make_tuple(asymID, 0, false);

			if (writtenAsyms.count(asymID) == 0)
			{
				writtenAsyms.insert(asymID);

				// clang-format off
				getCategory("struct_asym")->emplace({
					{ "id", asymID },
					{ "pdbx_blank_PDB_chainid_flag", het.chainID == ' ' ? "Y" : "N" },
					//					pdbx_modified
					{ "entity_id", mHet2EntityID[hetID] },
					//					details
				});

				// clang-format on
			}
		}

		int seqNr = ++ndbSeqNum[std::make_tuple(hetID, asymID)];
		int authSeqNr = ++entityAuthSeqNum[hetID];

		std::string iCode{ het.iCode };
		cif::trim(iCode);
		if (iCode.empty())
			iCode = { '.' };

		// clang-format off
		getCategory("pdbx_nonpoly_scheme")->emplace({
			{ "asym_id", asymID },
			{ "entity_id", mHet2EntityID[hetID] },
			{ "mon_id", hetID },
			{ "ndb_seq_num", seqNr },
			{ "pdb_seq_num", het.seqNum },
			{ "auth_seq_num", authSeqNr }, // Yes
			{ "pdb_mon_id", hetID },
			{ "auth_mon_id", hetID },
			{ "pdb_strand_id", std::string{ het.chainID } },
			{ "pdb_ins_code", iCode }
		});
		// clang-format on

		// mapping needed?
		mChainSeq2AsymSeq[std::make_tuple(het.chainID, het.seqNum, het.iCode)] = std::make_tuple(asymID, seqNr, false);
	}

	int modResID = 1;
	std::set<std::string> modResSet;
	for (auto rec = FindRecord("MODRES"); rec != nullptr and rec->is("MODRES");
		 rec = rec->mNext)                     //	 1 -  6        Record name   "MODRES"
	{                                          //	 8 - 11        IDcode        idCode      ID code of this datablock.
		std::string resName = rec->vS(13, 15); //	13 - 15        Residue name  resName     Residue name used in this datablock.
		char chainID = rec->vC(17);            //	17             Character     chainID     Chain identifier.
		int seqNum = rec->vI(19, 22);          //	19 - 22        Integer       seqNum      Sequence number.
		char iCode = rec->vC(23);              //	23             AChar         iCode       Insertion code.
		std::string stdRes = rec->vS(25, 27);  //	25 - 27        Residue name  stdRes      Standard residue name.
		std::string comment = rec->vS(30, 70); //	30 - 70        String        comment     Description of the residue modification.

		std::string asymID;
		int seq;
		std::error_code ec;

		std::tie(asymID, seq, std::ignore) = MapResidue(chainID, seqNum, iCode, ec);
		if (ec) // no need to write a modres if it could not be found
		{
			if (cif::VERBOSE > 0)
				std::cerr << "dropping unmapped MODRES record\n";
			continue;
		}

		// clang-format off
		getCategory("pdbx_struct_mod_residue")->emplace({
			{ "id", modResID++ },
			{ "label_asym_id", asymID },
			{ "label_seq_id", seq },
			{ "label_comp_id", resName },
			{ "auth_asym_id", std::string(1, chainID) },
			{ "auth_seq_id", seqNum },
			{ "auth_comp_id", resName },
			{ "PDB_ins_code", iCode == ' ' ? "" : std::string{ iCode } },
			{ "parent_comp_id", stdRes },
			{ "details", comment }
		});
		// clang-format on

		modResSet.insert(resName);
	}

	//	// chem compounds

	for (auto cc : mChemComp)
	{
		auto compound = cif::compound_factory::instance().create(
			mMod2parent.count(cc) ? mMod2parent[cc] : cc);

		std::string name;
		std::string formula;
		std::string type;
		std::string nstd = ".";
		std::optional<float> formulaWeight;

		if (compound != nullptr)
		{
			name = compound->name();
			type = compound->type();

			if (iequals(type, "L-peptide linking") or iequals(type, "peptide linking"))
				nstd = "y";

			formula = compound->formula();
			formulaWeight = compound->formula_weight();
		}

		if (name.empty())
			name = mHetnams[cc];

		if (type.empty())
			type = "NON-POLYMER";

		if (formula.empty())
		{
			formula = mFormuls[cc];

			const std::regex rx(R"(\d+\((.+)\))");
			std::smatch m;
			if (std::regex_match(formula, m, rx))
				formula = m[1].str();
		}

		if (modResSet.count(cc))
			nstd = "n";

		// clang-format off
		getCategory("chem_comp")->emplace({
			{ "id", cc },
			{ "name", name },
			{ "formula", formula },
			{ "formula_weight", formulaWeight, 3 },
			{ "mon_nstd_flag", nstd },
			{ "type", type }
		});
		// clang-format on
	}

	getCategory("chem_comp")->reorder_by_index();

	// unobserved can now be written as well

	int idRes = 0, idAtom = 0;
	sort(mUnobs.begin(), mUnobs.end(), [](const UNOBS &a, const UNOBS &b) -> bool
		{
			 int d = a.modelNr - b.modelNr;
			 if (d == 0)
				 d = a.seq - b.seq;
			 return d < 0; });

	for (auto &unobs : mUnobs)
	{
		bool isPolymer = false;
		std::string asymID, compID = unobs.res;
		int seqNr = 0;
		std::error_code ec;

		std::tie(asymID, seqNr, isPolymer) = MapResidue(unobs.chain, unobs.seq, unobs.iCode, ec);
		if (ec)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "error mapping unobserved residue\n";
			continue;
		}

		if (unobs.atoms.empty())
		{
			// clang-format off
			getCategory("pdbx_unobs_or_zero_occ_residues")->emplace({
				{ "id", std::to_string(++idRes) },
				{ "polymer_flag", isPolymer ? "Y" : "N" },
				{ "occupancy_flag", 1 },
				{ "PDB_model_num", unobs.modelNr ? unobs.modelNr : 1 },
				{ "auth_asym_id", std::string{ unobs.chain } },
				{ "auth_comp_id", unobs.res },
				{ "auth_seq_id", unobs.seq },
				{ "PDB_ins_code", unobs.iCode == ' ' ? "" : std::string{ unobs.iCode } },
				{ "label_asym_id", asymID },
				{ "label_comp_id", compID }, // TODO: change to correct comp_id
				{ "label_seq_id", seqNr > 0 ? std::to_string(seqNr) : "" }
			});
			// clang-format on
		}
		else
		{
			for (auto &atom : unobs.atoms)
			{
				// clang-format off
				getCategory("pdbx_unobs_or_zero_occ_atoms")->emplace({
					{ "id", std::to_string(++idAtom) },
					{ "polymer_flag", isPolymer ? "Y" : "N" },
					{ "occupancy_flag", 1 },
					{ "PDB_model_num", unobs.modelNr ? unobs.modelNr : 1 },
					{ "auth_asym_id", std::string{ unobs.chain } },
					{ "auth_comp_id", unobs.res },
					{ "auth_seq_id", unobs.seq },
					{ "PDB_ins_code", unobs.iCode == ' ' ? "" : std::string{ unobs.iCode } },
					{ "auth_atom_id", atom },
					{ "label_asym_id", asymID },
					{ "label_comp_id", compID }, // TODO: change to correct comp_id
					{ "label_seq_id", seqNr > 0 ? std::to_string(seqNr) : "" },
					{ "label_atom_id", atom }
				});
				// clang-format on
			}
		}
	}
}

void PDBFileParser::ConstructSugarTrees(int &asymNr)
{
	for (;;)
	{
		// find a first NAG/NDG
		auto si = std::find_if(mHets.begin(), mHets.end(), [](const HET &h)
			{ return (h.hetID == "NAG" or h.hetID == "NDG") and not(h.processed or h.branch); });
		if (si != mHets.end())
		{
			si->processed = true;

			// take the location of the C1 atom(s?)
			std::set<char> ci;

			for (auto a : si->atoms)
			{
				std::string name = a->vS(13, 16); //	13 - 16        Atom          name         Atom name.

				if (name != "C1")
					continue;

				ci.insert(a->vC(17)); //	17             Character     altLoc       Alternate location indicator.
			}

			if (ci.empty())
				continue;

			for (auto alt : ci)
			{
				ATOM_REF c1{ "C1", si->hetID, si->seqNum, si->chainID, si->iCode, alt };

				const auto &[asn, linked] = FindLink(c1, "ND2", "ASN");
				if (not linked)
					continue;

				std::stack<ATOM_REF> c1s;
				c1s.push(c1);

				SUGAR_TREE sugarTree;
				sugarTree.push_back({ c1 });

				// naive implementation
				while (not c1s.empty())
				{
					c1 = c1s.top();
					c1s.pop();

					for (auto o : { "O1", "O2", "O3", "O4", "O5", "O6" })
					{
						ATOM_REF leaving = c1;
						leaving.name = o;

						const auto &[nc1, linked_c1] = FindLink(leaving, "C1");
						if (linked_c1)
						{
							sugarTree.push_back({ nc1, o[1] - '0', c1 });
							c1s.push(nc1);
						}
					}
				}

				if (sugarTree.size() < 2) // not really a tree
					continue;

				auto branchName = sugarTree.entityName();
				auto entityID = mBranch2EntityID[branchName];

				// See if we've already added it to the entities
				if (entityID.empty())
				{
					entityID = std::to_string(mNextEntityNr++);
					mBranch2EntityID[branchName] = entityID;

					// clang-format off
					getCategory("entity")->emplace({
						{ "id", entityID },
						{ "type", "branched" },
						{ "src_method", "man" },
						{ "pdbx_description", branchName }
					});

					getCategory("pdbx_entity_branch")->emplace({
						{ "entity_id", entityID },
						{ "type", "oligosaccharide" }
					});
					// clang-format on

					int num = 0;
					std::map<ATOM_REF, int> branch_list;

					for (auto &s : sugarTree)
					{
						// clang-format off
						getCategory("pdbx_entity_branch_list")->emplace({
							{ "entity_id", entityID },
							{ "comp_id", s.c1.resName },
							{ "num", ++num },
							{ "hetero", ci.size() == 1 ? "n" : "y" }
						});
						// clang-format on

						branch_list[s.c1] = num;
					}

					auto &branch_link = *getCategory("pdbx_entity_branch_link");

					for (auto &s : sugarTree)
					{
						if (s.leaving_o == 0)
							continue;

						// clang-format off
						branch_link.emplace({
							{ "link_id", branch_link.size() + 1 },
							{ "entity_id", entityID },
							{ "entity_branch_list_num_1", branch_list[s.c1] },
							{ "comp_id_1", s.c1.resName },
							{ "atom_id_1", s.c1.name },
							{ "leaving_atom_id_1", "O1" },
							{ "entity_branch_list_num_2", branch_list[s.next] },
							{ "comp_id_2", s.next.resName },
							{ "atom_id_2", "O" + std::to_string(s.leaving_o) },
							{ "leaving_atom_id_2", "HO" + std::to_string(s.leaving_o) },
							{ "value_order", "sing" } /// ??
						});
						// clang-format on
					}
				}

				mSugarEntities.insert(entityID);

				// create an asym for this sugar tree

				std::string asymID = cif::cif_id_for_number(asymNr++);

				mAsymID2EntityID[asymID] = entityID;

				// clang-format off
				getCategory("struct_asym")->emplace({
					{ "id", asymID },
					{ "pdbx_blank_PDB_chainid_flag", si->chainID == ' ' ? "Y" : "N" },
					{ "pdbx_modified", "N" },
					{ "entity_id", entityID }
				});
				// clang-format on

				std::string iCode{ si->iCode };
				cif::trim(iCode);
				if (iCode.empty())
					iCode = { '.' };

				int num = 0;
				for (auto s : sugarTree)
				{
					// clang-format off
					getCategory("pdbx_branch_scheme")->emplace({
						{ "asym_id", asymID },
						{ "entity_id", entityID },
						{ "mon_id", s.c1.resName },
						{ "num", ++num },
						{ "pdb_asym_id", asymID },
						{ "pdb_mon_id", s.c1.resName },
						{ "pdb_seq_num", num },
						{ "auth_asym_id", std::string{ s.c1.chainID } },
						{ "auth_mon_id", s.next.resName },
						{ "auth_seq_num", s.c1.resSeq },
						{ "hetero", ci.size() == 1 ? "n" : "y" }
					});
					// clang-format on

					auto k = std::make_tuple(s.c1.chainID, s.c1.resSeq, s.c1.iCode);
					assert(mChainSeq2AsymSeq.count(k) == 0);

					mChainSeq2AsymSeq[k] = std::make_tuple(asymID, num, false);

					// mark all hets as part of tree

					for (auto &h : mHets)
					{
						if (h.hetID == s.c1.resName and h.chainID == s.c1.chainID and h.seqNum == s.c1.resSeq and h.iCode == s.c1.iCode)
						{
							h.branch = true;
							break; // should be only one of course... right?
						}
					}
				}

				break;
			}

			continue;
		}

		break;
	}

	// remove the branched HET's
	mHets.erase(std::remove_if(mHets.begin(), mHets.end(), [](auto &h)
					{ return h.branch; }),
		mHets.end());
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

		std::string begAsymID, endAsymID;
		int begSeq, endSeq;
		std::error_code ec;

		std::tie(begAsymID, begSeq, std::ignore) = MapResidue(vC(20), vI(22, 25), vC(26), ec);
		if (not ec)
			std::tie(endAsymID, endSeq, std::ignore) = MapResidue(vC(32), vI(34, 37), vC(38), ec);

		if (ec)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Could not map residue for HELIX " << vI(8, 10) << '\n';
		}
		else
		{
			auto cat = getCategory("struct_conf");
			// clang-format off
			cat->emplace({
				{ "conf_type_id", "HELX_P" },
				{ "id", "HELX_P" + std::to_string(vI(8, 10)) },
				{ "pdbx_PDB_helix_id", vS(12, 14) },
				{ "beg_label_comp_id", vS(16, 18) },
				{ "beg_label_asym_id", begAsymID },
				{ "beg_label_seq_id", begSeq },
				{ "pdbx_beg_PDB_ins_code", vS(26, 26) },
				{ "end_label_comp_id", vS(28, 30) },
				{ "end_label_asym_id", endAsymID },
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
			// clang-format off

			if (firstHelix)
			{
				cat = getCategory("struct_conf_type");
				cat->emplace({ { "id", "HELX_P" } });
				firstHelix = false;
			}
		}

		GetNextRecord();
	}

	std::set<std::string> sheetsSeen;
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
		//	50             Character     curChainID     Registration. Chain identifier in
		//	                                            current strand.
		//	51 - 54        Integer       curResSeq      Registration.  Residue sequence number
		//	                                            in current strand.
		//	55             AChar         curICode       Registration. Insertion code in
		//	                                            current strand.
		//	57 - 60        Atom          prevAtom       Registration.  Atom name in previous strand.
		//	61 - 63        Residue name  prevResName    Registration.  Residue name in
		//	                                            previous strand.
		//	65             Character     prevChainID    Registration.  Chain identifier in
		//	                                            previous  strand.
		//	66 - 69        Integer       prevResSeq     Registration. Residue sequence number
		//	                                            in previous strand.
		//	70             AChar         prevICode      Registration.  Insertion code in
		//	                                            previous strand.

		std::string sheetID = cif::trim_copy(vS(12, 14));
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
			// clang-format off
			getCategory("struct_sheet_order")->emplace({
				{ "sheet_id", sheetID },
				{ "range_id_1", rangeID },
				{ "range_id_2", rangeID + 1 },
				{ "sense", sense == -1 ? "anti-parallel" : "parallel" }
			});
			// clang-format on
		}

		std::string begAsymID, endAsymID;
		int begSeq, endSeq;
		std::error_code ec;

		std::tie(begAsymID, begSeq, std::ignore) = MapResidue(vC(22), vI(23, 26), vC(27), ec);
		if (not ec)
			std::tie(endAsymID, endSeq, std::ignore) = MapResidue(vC(33), vI(34, 37), vC(38), ec);

		if (ec)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Dropping SHEET record " << vI(8, 10) << '\n';
		}
		else
		{
			// clang-format off
			getCategory("struct_sheet_range")->emplace({
				{ "sheet_id", sheetID },
				{ "id", vI(8, 10) },
				{ "beg_label_comp_id", vS(18, 20) },
				{ "beg_label_asym_id", begAsymID },
				{ "beg_label_seq_id", begSeq },
				{ "pdbx_beg_PDB_ins_code", vS(27, 27) },
				{ "end_label_comp_id", vS(29, 31) },
				{ "end_label_asym_id", endAsymID },
				{ "end_label_seq_id", endSeq },
				{ "pdbx_end_PDB_ins_code", vS(38, 38) },

				{ "beg_auth_comp_id", vS(18, 20) },
				{ "beg_auth_asym_id", vS(22, 22) },
				{ "beg_auth_seq_id", vI(23, 26) },
				{ "end_auth_comp_id", vS(29, 31) },
				{ "end_auth_asym_id", vS(33, 33) },
				{ "end_auth_seq_id", vI(34, 37) },
			});
			// clang-format on

			if (sense != 0 and mRec->mVlen > 34)
			{
				std::string r1AsymID, r2AsymID;
				int r1Seq, r2Seq;

				std::tie(r1AsymID, r1Seq, std::ignore) = MapResidue(vC(65), vI(66, 69), vC(70), ec);
				if (not ec)
					std::tie(r2AsymID, r2Seq, std::ignore) = MapResidue(vC(50), vI(51, 54), vC(55), ec);

				if (ec)
				{
					if (cif::VERBOSE > 0)
						std::cerr << "skipping unmatched pdbx_struct_sheet_hbond record\n";
				}
				else
					// clang-format off
					getCategory("pdbx_struct_sheet_hbond")->emplace({
						{ "sheet_id", sheetID },
						{ "range_id_1", rangeID },
						{ "range_id_2", rangeID + 1 },
						{ "range_1_label_atom_id", vS(57, 60) },
						{ "range_1_label_comp_id", vS(61, 63) },
						{ "range_1_label_asym_id", r1AsymID },
						{ "range_1_label_seq_id", r1Seq },
						{ "range_1_PDB_ins_code", vS(70, 70) },
						{ "range_1_auth_atom_id", vS(57, 60) },
						{ "range_1_auth_comp_id", vS(61, 63) },
						{ "range_1_auth_asym_id", vS(65, 65) },
						{ "range_1_auth_seq_id", vI(66, 69) },

						{ "range_2_label_atom_id", vS(42, 45) },
						{ "range_2_label_comp_id", vS(46, 48) },
						{ "range_2_label_asym_id", r2AsymID },
						{ "range_2_label_seq_id", r2Seq },
						{ "range_2_PDB_ins_code", vS(55, 55) },
						{ "range_2_auth_atom_id", vS(42, 45) },
						{ "range_2_auth_comp_id", vS(46, 48) },
						{ "range_2_auth_asym_id", vS(50, 50) },
						{ "range_2_auth_seq_id", vI(51, 54) }
					});
				// clang-format on
			}

			if (sense != 0)
				++rangeID;
		}

		GetNextRecord();
	}
}

static bool IsMetal(const std::string &resName, const std::string &atomID)
{
	bool result = false;

	try
	{
		auto compound = cif::compound_factory::instance().create(resName);
		if (compound != nullptr)
		{
			auto at = cif::atom_type_traits(compound->get_atom_by_atom_id(atomID).type_symbol);
			result = at.is_metal();
		}
	}
	catch (...)
	{
	}

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

			std::string p1Asym, p2Asym;
			int p1Seq = 0, p2Seq = 0;
			std::error_code ec;

			std::tie(p1Asym, p1Seq, std::ignore) = MapResidue(vC(16), vI(18, 21), vC(22), ec);
			if (not ec)
				std::tie(p2Asym, p2Seq, std::ignore) = MapResidue(vC(30), vI(32, 35), vC(36), ec);

			if (ec)
			{
				if (cif::VERBOSE > 0)
					std::cerr << "Dropping SSBOND " << vI(8, 10) << '\n';
				continue;
			}

			std::vector<char> alt1 = altLocsForAtom(vC(16), vI(18, 21), vC(22), "SG");
			std::vector<char> alt2 = altLocsForAtom(vC(30), vI(32, 35), vC(36), "SG");

			if (alt1.empty())
				alt1.push_back(0);
			if (alt2.empty())
				alt2.push_back(0);

			std::string sym1, sym2;
			try
			{
				sym1 = pdb2cifSymmetry(vS(60, 65));
				sym2 = pdb2cifSymmetry(vS(67, 72));
			}
			catch (const std::exception &ex)
			{
				if (cif::VERBOSE > 0)
					std::cerr << "Dropping SSBOND " << vI(8, 10) << " due to invalid symmetry operation\n";
				continue;
			}

			for (auto a1 : alt1)
			{
				for (auto a2 : alt2)
				{
					// clang-format off
					getCategory("struct_conn")->emplace({
						{ "id", "disulf" + std::to_string(++ssBondNr) },
						{ "conn_type_id", "disulf" },

						{ "ptnr1_label_asym_id", p1Asym },
						{ "pdbx_ptnr1_label_alt_id", a1 ? std::string{ a1 } : std::string() },
						{ "ptnr1_label_comp_id", vS(12, 14) },
						{ "ptnr1_label_seq_id", p1Seq ? std::to_string(p1Seq) : "." },
						{ "ptnr1_label_atom_id", "SG" },
						{ "ptnr1_symmetry", sym1 },

						{ "ptnr2_label_asym_id", p2Asym },
						{ "pdbx_ptnr2_label_alt_id", a2 ? std::string{ a2 } : std::string() },
						{ "ptnr2_label_comp_id", vS(26, 28) },
						{ "ptnr2_label_seq_id", p2Seq ? std::to_string(p2Seq) : "." },
						{ "ptnr2_label_atom_id", "SG" },

						{ "ptnr1_auth_asym_id", vS(16, 16) },
						{ "ptnr1_auth_comp_id", vS(12, 14) },
						{ "ptnr1_auth_seq_id", vI(18, 21) },
						{ "ptnr2_auth_asym_id", vS(30, 30) },
						{ "ptnr2_auth_comp_id", vS(26, 28) },
						{ "ptnr2_auth_seq_id", vI(32, 35) },

						{ "ptnr2_symmetry", sym2 },

						{ "pdbx_dist_value", vS(74, 78) },
					});
					// clang-format on
				}
			}

			continue;
		}

		if (mRec->is("LINK  ") or mRec->is("LINKR "))
		{
			if (cif::VERBOSE > 0 and mRec->is("LINKR "))
				std::cerr << "Accepting non-standard LINKR record, but ignoring extra information\n";

			//	 1 -  6         Record name    "LINK  "
			std::string name1 = vS(13, 16);    //	13 - 16         Atom           name1           Atom name.
			                                   //	17              Character      altLoc1         Alternate location indicator.
			std::string resName1 = vS(18, 20); //	18 - 20         Residue name   resName1        Residue  name.
			                                   //	22              Character      chainID1        Chain identifier.
			                                   //	23 - 26         Integer        resSeq1         Residue sequence number.
			                                   //	27              AChar          iCode1          Insertion code.
			std::string name2 = vS(43, 46);    //	43 - 46         Atom           name2           Atom name.
			                                   //	47              Character      altLoc2         Alternate location indicator.
			std::string resName2 = vS(48, 50); //	48 - 50         Residue name   resName2        Residue name.
			                                   //	52              Character      chainID2        Chain identifier.
			                                   //	53 - 56         Integer        resSeq2         Residue sequence number.
			                                   //	57              AChar          iCode2          Insertion code.
			                                   //	60 - 65         SymOP          sym1            Symmetry operator atom 1.
			                                   //	67 - 72         SymOP          sym2            Symmetry operator atom 2.
			                                   //	74 – 78         Real(5.2)      Length          Link distance

			std::string type = "covale";
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

			std::string p1Asym, p2Asym;
			int p1Seq = 0, p2Seq = 0;
			bool isResseq1 = false, isResseq2 = false;
			std::error_code ec;

			std::tie(p1Asym, p1Seq, isResseq1) = MapResidue(vC(22), vI(23, 26), vC(27), ec);
			if (not ec)
				std::tie(p2Asym, p2Seq, isResseq2) = MapResidue(vC(52), vI(53, 56), vC(57), ec);

			if (ec)
			{
				if (cif::VERBOSE > 0)
					std::cerr << "Dropping LINK record at line " << mRec->mLineNr << '\n';
				continue;
			}

			std::string distance, ccp4LinkID;

			if (mRec->is("LINK  "))
			{
				distance = vS(74, 78);

				double d;
				auto r = cif::from_chars(distance.data(), distance.data() + distance.length(), d);
				if ((bool)r.ec)
				{
					if (cif::VERBOSE > 0)
						std::cerr << "Distance value '" << distance << "' is not a valid float in LINK record\n";
					swap(ccp4LinkID, distance); // assume this is a ccp4_link_id... oh really?
				}
			}
			else                         // LINKR
				ccp4LinkID = vS(74, 78); // the link ID

			std::string sym1, sym2;
			try
			{
				sym1 = pdb2cifSymmetry(vS(60, 65));
				sym2 = pdb2cifSymmetry(vS(67, 72));
			}
			catch (const std::exception &ex)
			{
				if (cif::VERBOSE > 0)
					std::cerr << "Dropping LINK record at line " << mRec->mLineNr << " due to invalid symmetry operation\n";
				continue;
			}

			// clang-format off
			getCategory("struct_conn")->emplace({
				{ "id", type + std::to_string(linkNr) },
				{ "conn_type_id", type },

				// { "ccp4_link_id", ccp4LinkID },

				{ "ptnr1_label_asym_id", p1Asym },
				{ "ptnr1_label_comp_id", vS(18, 20) },
				{ "ptnr1_label_seq_id", (isResseq1 and p1Seq) ? std::to_string(p1Seq) : "." },
				{ "ptnr1_label_atom_id", vS(13, 16) },
				{ "pdbx_ptnr1_label_alt_id", vS(17, 17) },
				{ "pdbx_ptnr1_PDB_ins_code", vS(27, 27) },
				{ "pdbx_ptnr1_standard_comp_id", "" },
				{ "ptnr1_symmetry", sym1 },

				{ "ptnr2_label_asym_id", p2Asym },
				{ "ptnr2_label_comp_id", vS(48, 50) },
				{ "ptnr2_label_seq_id", (isResseq2 and p2Seq) ? std::to_string(p2Seq) : "." },
				{ "ptnr2_label_atom_id", vS(43, 46) },
				{ "pdbx_ptnr2_label_alt_id", vS(47, 47) },
				{ "pdbx_ptnr2_PDB_ins_code", vS(57, 57) },

				{ "ptnr1_auth_asym_id", vS(22, 22) },
				{ "ptnr1_auth_comp_id", vS(18, 20) },
				{ "ptnr1_auth_seq_id", vI(23, 26) },
				{ "ptnr2_auth_asym_id", vS(52, 52) },
				{ "ptnr2_auth_comp_id", vS(48, 50) },
				{ "ptnr2_auth_seq_id", vI(53, 56) },

				// { "ptnr1_auth_atom_id", vS(13, 16) },
			    // { "ptnr2_auth_atom_id", vS(43, 46) },

				{ "ptnr2_symmetry", sym2 },

				{ "pdbx_dist_value", distance }
			});
			// clang-format on

			continue;
		}

		if (mRec->is("CISPEP"))
		{
			//	 1 -  6       Record name   "CISPEP"
			int serNum = vI(8, 10);           //	 8 - 10       Integer       serNum        Record serial number.
			std::string pep1 = vS(12, 14);    //	12 - 14       LString(3)    pep1          Residue name.
			char chainID1 = vC(16);           //	16            Character     chainID1      Chain identifier.
			int seqNum1 = vI(18, 21);         //	18 - 21       Integer       seqNum1       Residue sequence number.
			char iCode1 = vC(22);             //	22            AChar         icode1        Insertion code.
			std::string pep2 = vS(26, 28);    //	26 - 28       LString(3)    pep2          Residue name.
			char chainID2 = vC(30);           //	30            Character     chainID2      Chain identifier.
			int seqNum2 = vI(32, 35);         //	32 - 35       Integer       seqNum2       Residue sequence number.
			char iCode2 = vC(36);             //	36            AChar         icode2        Insertion code.
			int modNum = vI(44, 46);          //	44 - 46       Integer       modNum        Identifies the specific model.
			std::string measure = vF(54, 59); //	54 - 59       Real(6.2)     measure       Angle measurement in degrees.

			if (modNum == 0)
				modNum = 1;

			std::string lAsym1, lAsym2;
			int lResSeq1, lResSeq2;
			std::error_code ec;

			std::tie(lAsym1, lResSeq1, std::ignore) = MapResidue(chainID1, seqNum1, iCode1, ec);
			if (not ec)
				std::tie(lAsym2, lResSeq2, std::ignore) = MapResidue(chainID2, seqNum2, iCode2, ec);

			if (ec)
			{
				if (cif::VERBOSE > 0)
					std::cerr << "Dropping CISPEP record at line " << mRec->mLineNr << '\n';
				continue;
			}

			std::string iCode1str = iCode1 == ' ' ? std::string() : std::string{ iCode1 };
			std::string iCode2str = iCode2 == ' ' ? std::string() : std::string{ iCode2 };

			// clang-format off
			getCategory("struct_mon_prot_cis")->emplace({
				{ "pdbx_id", serNum },
				{ "label_comp_id", pep1 },
				{ "label_seq_id", lResSeq1 },
				{ "label_asym_id", lAsym1 },
				{ "label_alt_id", "." },
				{ "pdbx_PDB_ins_code", iCode1str },
				{ "auth_comp_id", pep1 },
				{ "auth_seq_id", seqNum1 },
				{ "auth_asym_id", std::string{ chainID1 } },
				{ "pdbx_label_comp_id_2", pep2 },
				{ "pdbx_label_seq_id_2", lResSeq2 },
				{ "pdbx_label_asym_id_2", lAsym2 },
				{ "pdbx_PDB_ins_code_2", iCode2str },
				{ "pdbx_auth_comp_id_2", pep2 },
				{ "pdbx_auth_seq_id_2", seqNum2 },
				{ "pdbx_auth_asym_id_2", std::string{ chainID2 } },
				{ "pdbx_PDB_model_num", modNum },
				{ "pdbx_omega_angle", measure }
			});
			// clang-format on

			continue;
		}

		break;
	}
}

void PDBFileParser::ParseMiscellaneousFeatures()
{
	int structSiteGenID = 1;

	while (mRec->is("SITE  "))
	{                                    //	 1 -  6        Record name   "SITE  "
		                                 //	 8 - 10        Integer       seqNum        Sequence number.
		std::string siteID = vS(12, 14); //	12 - 14        LString(3)    siteID        Site name.
		int numRes = vI(16, 17);         //	16 - 17        Integer       numRes        Number of residues that compose the site.

		int o = 19;

		auto cat = getCategory("struct_site_gen");

		for (int i = 0; i < numRes; ++i)
		{
			std::string resName = vS(o, o + 2); //	19 - 21        Residue name  resName1      Residue name for first residue that
			                                    //	                                           creates the site.
			char chainID = vC(o + 4);           //	23             Character     chainID1      Chain identifier for first residue of site.
			int seq = vI(o + 5, o + 8);         //	24 - 27        Integer       seq1          Residue sequence number for first residue
			                                    //	                                           of the  site.
			char iCode = vC(o + 9);             //	28             AChar         iCode1        Insertion code for first residue of the site.

			int labelSeq;
			std::string asym;
			bool isResseq;
			std::error_code ec;

			std::tie(asym, labelSeq, isResseq) = MapResidue(chainID, seq, iCode, ec);

			if (ec)
			{
				if (cif::VERBOSE > 0)
					std::cerr << "skipping struct_site_gen record\n";
			}
			else
				// clang-format off
				cat->emplace({
					{ "id", structSiteGenID++ },
					{ "site_id", siteID },
					{ "pdbx_num_res", numRes },
					{ "label_comp_id", resName },
					{ "label_asym_id", asym },
					{ "label_seq_id", (labelSeq > 0 and isResseq) ? std::to_string(labelSeq) : std::string(".") },
					{ "pdbx_auth_ins_code", iCode == ' ' ? "" : std::string{ iCode } },
					{ "auth_comp_id", resName },
					{ "auth_asym_id", std::string{ chainID } },
					{ "auth_seq_id", seq },
					{ "label_atom_id", "." },
					{ "label_alt_id", "." },
				});
			// clang-format on

			o += 11;
		}

		GetNextRecord();
	}
}

void PDBFileParser::ParseCrystallographic()
{
	if (mRec->is("CRYST1"))
	{
		Match("CRYST1", true);

		// clang-format off
		getCategory("cell")->emplace({
			{ "entry_id", mStructureID },  //	 1 -  6       Record name   "CRYST1"
			{ "length_a", vF(7, 15) },     //	 7 - 15       Real(9.3)     a              a (Angstroms).
			{ "length_b", vF(16, 24) },    //	16 - 24       Real(9.3)     b              b (Angstroms).
			{ "length_c", vF(25, 33) },    //	25 - 33       Real(9.3)     c              c (Angstroms).
			{ "angle_alpha", vF(34, 40) }, //	34 - 40       Real(7.2)     alpha          alpha (degrees).
			{ "angle_beta", vF(41, 47) },  //	41 - 47       Real(7.2)     beta           beta (degrees).
			{ "angle_gamma", vF(48, 54) }, //	48 - 54       Real(7.2)     gamma          gamma (degrees).
			/* goes into symmetry */       //	56 - 66       LString       sGroup         Space  group.
			{ "Z_PDB", vF(67, 70) }        //	67 - 70       Integer       z              Z value.
		});
		// clang-format on

		std::string spaceGroup, intTablesNr;
		try
		{
			spaceGroup = vS(56, 66);
			intTablesNr = std::to_string(get_space_group_number(spaceGroup));
		}
		catch (...)
		{
		}

		// clang-format off
		getCategory("symmetry")->emplace({
			{ "entry_id", mStructureID },
			{ "space_group_name_H-M", spaceGroup },
			{ "Int_Tables_number", intTablesNr }
		});

		GetNextRecord();
	}
}

void PDBFileParser::ParseCoordinateTransformation()
{
	std::string m[3][3], v[3];

	if (cif::starts_with(mRec->mName, "ORIGX"))
	{
		for (std::string n : { "1", "2", "3" })
		{
			int x = stoi(n) - 1;

			Match("ORIGX" + n, true); //	 1 -  6         Record name   "ORIGXn"      n=1, 2, or 3
			m[x][0] = vF(11, 20);     //	11 - 20         Real(10.6)    o[n][1]       On1
			m[x][1] = vF(21, 30);     //	21 - 30         Real(10.6)    o[n][2]       On2
			m[x][2] = vF(31, 40);     //	31 - 40         Real(10.6)    o[n][3]       On3
			v[x] = vF(46, 55);        //	46 - 55         Real(10.5)    t[n]          Tn

			GetNextRecord();
		}

		// clang-format off
		getCategory("database_PDB_matrix")->emplace({
			{ "entry_id", mStructureID },
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
		// clang-format on
	}

	if (cif::starts_with(mRec->mName, "SCALE"))
	{
		for (std::string n : { "1", "2", "3" })
		{
			int x = stoi(n) - 1;

			Match("SCALE" + n, true); //	 1 -  6         Record name   "SCALEn" n=1,  2, or 3
			m[x][0] = vF(11, 20);     //	11 - 20         Real(10.6)    s[n][1]            Sn1
			m[x][1] = vF(21, 30);     //	21 - 30         Real(10.6)    s[n][2]            Sn2
			m[x][2] = vF(31, 40);     //	31 - 40         Real(10.6)    s[n][3]            Sn3
			v[x] = vF(46, 55);        //	46 - 55         Real(10.5)    u[n]               Un

			GetNextRecord();
		}

		// clang-format off
		getCategory("atom_sites")->emplace({
			{ "entry_id", mStructureID },
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
		// clang-format on
	}

	while (cif::starts_with(mRec->mName, "MTRIX1"))
	{
		int serial = 0, igiven = 0;

		for (std::string n : { "1", "2", "3" })
		{
			int x = stoi(n) - 1;

			Match("MTRIX" + n, true); //	 1 -  6        Record name   "MTRIXn"      n=1, 2, or 3
			serial = vI(8, 10);       //	 8 - 10        Integer       serial        Serial number.
			m[x][0] = vF(11, 20);     //	11 - 20        Real(10.6)    m[n][1]       Mn1
			m[x][1] = vF(21, 30);     //	21 - 30        Real(10.6)    m[n][2]       Mn2
			m[x][2] = vF(31, 40);     //	31 - 40        Real(10.6)    m[n][3]       Mn3
			v[x] = vF(46, 55);        //	46 - 55        Real(10.5)    v[n]          Vn
			igiven = vC(60) == '1';   //	60             Integer       iGiven        1 if coordinates for the  representations
			                          //	                                           which  are approximately related by the
			GetNextRecord();          //	                                           transformations  of the molecule are
		}                             //	                                           contained in the datablock. Otherwise, blank.

		// clang-format off
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
		// clang-format on
	}
}

void PDBFileParser::ParseCoordinate(int modelNr)
{
	// oh oh, we have to sort our atom_site records by ascending asym_id
	// This routine used to be so trivial...

	typedef std::tuple<std::string, int, bool, PDBRecord *, PDBRecord *> atomRec;

	std::vector<atomRec> atoms;
	while (mRec->is("ATOM  ") or mRec->is("HETATM")) //	 1 -  6        Record name   "ATOM  "
	{
		char chainID = vC(22);   //	22             Character     chainID      Chain identifier.
		int resSeq = vI(23, 26); //	23 - 26        Integer       resSeq       Residue sequence number.
		char iCode = vC(27);

		std::string asymID;
		int seqID;
		bool isResseq;

		std::tie(asymID, seqID, isResseq) = MapResidue(chainID, resSeq, iCode);

		PDBRecord *atom = mRec;
		PDBRecord *anisou = nullptr;

		GetNextRecord();
		if (mRec->is("ANISOU"))
		{
			anisou = mRec;
			GetNextRecord();
		}

		atoms.emplace_back(asymID, seqID, isResseq, atom, anisou);

		/*if?... */ while (mRec->is("TER   "))
		{
			Match("TER   ", true);
			GetNextRecord();
		}
	}

	auto last = mRec;

	// use stable sort here
	auto rLess = [](const atomRec &a, const atomRec &b) -> bool
	{
		int d;

		std::string chainA = std::get<0>(a);
		std::string chainB = std::get<0>(b);

		if (chainA.length() != chainB.length())
			d = static_cast<int>(chainA.length() - chainB.length());
		else
			d = std::get<0>(a).compare(std::get<0>(b));

		if (d == 0)
			d = std::get<1>(a) - std::get<1>(b);
		return d < 0;
	};

	stable_sort(atoms.begin(), atoms.end(), rLess);

	// now reiterate the atoms to reorder alternates
	for (std::size_t i = 0; i + 1 < atoms.size(); ++i)
	{
		char altLoc = std::get<3>(atoms[i])->vC(17);

		if (altLoc == ' ' or altLoc == 0)
			continue;

		auto b = atoms.begin() + i;
		auto e = b;

		std::map<std::string, int> atomIndex; // index number of first occurrence of a atom name

		while (e != atoms.end() and rLess(*b, *e) == false)
		{
			std::string name = std::get<3>(*e)->vS(13, 16);

			if (atomIndex.count(name) == 0)
				atomIndex[name] = static_cast<int>(atomIndex.size() + 1);

			++e;
		}

		auto aLess = [&](atomRec &a, atomRec &b) -> bool
		{
			std::string na = std::get<3>(a)->vS(13, 16);
			std::string nb = std::get<3>(b)->vS(13, 16);

			int d = atomIndex[na] - atomIndex[nb];
			if (d == 0)
				d = std::get<3>(a)->vC(17) - std::get<3>(b)->vC(17);
			assert(d != 0);
			return d < 0;
		};

		sort(b, e, aLess);

		i += distance(b, e) - 1;
	}

	//	while (mRec->is("ATOM  ") or mRec->is("HETATM"))		//	 1 -  6        Record name   "ATOM  "
	for (auto &a : atoms)
	{
		std::string asymID;
		int seqID;
		bool isResseq;
		PDBRecord *atom;
		PDBRecord *anisou;
		std::tie(asymID, seqID, isResseq, atom, anisou) = a;

		mRec = atom;

		++mAtomID;

		std::string groupPDB = mRec->is("ATOM  ") ? "ATOM" : "HETATM";
		//		int serial = vI(7, 11);				//	 7 - 11        Integer       serial       Atom  serial number.
		std::string name = vS(13, 16);       //	13 - 16        Atom          name         Atom name.
		char altLoc = vC(17);                //	17             Character     altLoc       Alternate location indicator.
		std::string resName = vS(18, 20);    //	18 - 20        Residue name  resName      Residue name.
		char chainID = vC(22);               //	22             Character     chainID      Chain identifier.
		int resSeq = vI(23, 26);             //	23 - 26        Integer       resSeq       Residue sequence number.
		char iCode = vC(27);                 //	27             AChar         iCode        Code for insertion of residues.
		std::string x = vF(31, 38);          //	31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		std::string y = vF(39, 46);          //	39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
		std::string z = vF(47, 54);          //	47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		std::string occupancy = vF(55, 60);  //	55 - 60        Real(6.2)     occupancy    Occupancy.
		std::string tempFactor = vF(61, 66); //	61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		std::string element = vS(77, 78);    //	77 - 78        LString(2)    element      Element symbol, right-justified.
		std::string charge = vS(79, 80);     //	79 - 80        LString(2)    charge       Charge  on the atom.

		std::string entityID = mAsymID2EntityID[asymID];

		charge = pdb2cifCharge(charge);

		// if (cif::compound_factory::instance().is_known_peptide(resName) or cif::compound_factory::instance().is_known_base(resName))
		if (resName == "UNK" or cif::compound_factory::kAAMap.count(resName) or cif::compound_factory::kBaseMap.count(resName))
		{
			if (groupPDB == "HETATM")
			{
				if (cif::VERBOSE > 0)
					std::cerr << "Changing atom from HETATM to ATOM at line " << mRec->mLineNr << '\n';
				groupPDB = "ATOM";
			}
		}
		else
		{
			if (groupPDB == "ATOM")
			{
				if (cif::VERBOSE > 0)
					std::cerr << "Changing atom from ATOM to HETATM at line " << mRec->mLineNr << '\n';
				groupPDB = "HETATM";
			}
		}

		// if the atom is part of a sugar, we need to replace the auth_seq_id/resSeq
		if (mSugarEntities.count(entityID))
		{
			using namespace cif::literals;

			auto &branch_scheme = *getCategory("pdbx_branch_scheme");
			resSeq = branch_scheme.find1<int>("asym_id"_key == asymID and "auth_seq_num"_key == resSeq, "pdb_seq_num");
		}

		// clang-format off
		getCategory("atom_site")->emplace({
			{ "group_PDB", groupPDB },
			{ "id", mAtomID },
			{ "type_symbol", element },
			{ "label_atom_id", name },
			{ "label_alt_id", altLoc != ' ' ? std::string{ altLoc } : "." },
			{ "label_comp_id", resName },
			{ "label_asym_id", asymID },
			{ "label_entity_id", entityID },
			{ "label_seq_id", (isResseq and seqID > 0) ? std::to_string(seqID) : "." },
			{ "pdbx_PDB_ins_code", iCode == ' ' ? "" : std::string{ iCode } },
			{ "Cartn_x", x },
			{ "Cartn_y", y },
			{ "Cartn_z", z },
			{ "occupancy", occupancy },
			{ "B_iso_or_equiv", tempFactor },
			{ "pdbx_formal_charge", charge },
			{ "auth_seq_id", resSeq },
			{ "auth_comp_id", resName },
			{ "auth_asym_id", std::string{ chainID } },
			{ "auth_atom_id", name },
			{ "pdbx_PDB_model_num", modelNr }
		});
		// clang-format on

		InsertAtomType(element);

		std::string check = vS(7, 11) + vS(77, 80);

		if (anisou != nullptr)
		{
			mRec = anisou;        //	 1 - 6        Record name   "ANISOU"
			int u11 = vI(29, 35); //	29 - 35       Integer       u[0][0]        U(1,1)
			int u22 = vI(36, 42); //	36 - 42       Integer       u[1][1]        U(2,2)
			int u33 = vI(43, 49); //	43 - 49       Integer       u[2][2]        U(3,3)
			int u12 = vI(50, 56); //	50 - 56       Integer       u[0][1]        U(1,2)
			int u13 = vI(57, 63); //	57 - 63       Integer       u[0][2]        U(1,3)
			int u23 = vI(64, 70); //	64 - 70       Integer       u[1][2]        U(2,3)

			if (vS(7, 11) + vS(77, 80) != check)
				throw std::runtime_error("ANISOU record should follow corresponding ATOM record");

			auto f = [](float f) -> std::string
			{
				return cif::format("%6.4f", f).str();
			};

			// clang-format off
			getCategory("atom_site_anisotrop")->emplace({
				{ "id", mAtomID },
				{ "type_symbol", element },
				{ "pdbx_label_atom_id", name },
				{ "pdbx_label_alt_id", altLoc != ' ' ? std::string{ altLoc } : "." },
				{ "pdbx_label_comp_id", resName },
				{ "pdbx_label_asym_id", asymID },
				{ "pdbx_label_seq_id", (isResseq and seqID > 0) ? std::to_string(seqID) : "." },
				{ "U[1][1]", f(u11 / 10000.f) },
				{ "U[2][2]", f(u22 / 10000.f) },
				{ "U[3][3]", f(u33 / 10000.f) },
				{ "U[1][2]", f(u12 / 10000.f) },
				{ "U[1][3]", f(u13 / 10000.f) },
				{ "U[2][3]", f(u23 / 10000.f) },
				{ "pdbx_auth_seq_id", resSeq },
				{ "pdbx_auth_comp_id", resName },
				{ "pdbx_auth_asym_id", std::string{ chainID } },
				{ "pdbx_auth_atom_id", name }
			});
			// clang-format on
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

void PDBFileParser::Parse(std::istream &is, cif::file &result)
{
	try
	{
		mDatablock.set_validator(result.get_validator());

		PreParseInput(is);

		mRec = mData;

		ParseTitle();

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

		uint32_t modelNr = 1;
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
			throw std::runtime_error("Either the PDB file has no atom records, or the field " + std::string(mRec->mName) + " is not at the correct location");

		for (auto e : mAtomTypes)
			getCategory("atom_type")->emplace({ { "symbol", e } });

		// in V5, atom_type is sorted
		getCategory("atom_type")->reorder_by_index();

		ParseConnectivty();
		ParseBookkeeping();

		// almost done, now fix some outstanding issued that could not be done before

		try
		{
			auto r = FindRecord("REMARK   3");

			if (r != nullptr and Remark3Parser::parse(mExpMethod, r, mDatablock))
			{
				// make sure the "exptl" category is created
				auto exptl = getCategory("exptl");
				if (exptl->empty())
				{
					exptl->emplace({ { "entry_id", mStructureID },
						{ "method", mExpMethod },
						{ "crystals_number", mRemark200["NUMBER OF CRYSTALS USED"] } });
				}
			}
		}
		catch (const std::exception &ex)
		{
			if (cif::VERBOSE >= 0)
				std::cerr << "Error parsing REMARK 3\n";
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
		//				std::string asymID;
		//				int resNum;
		//
		//				cif::tie(asymID, resNum) = r.get("beg_auth_asym_id", "beg_auth_seq_id");
		//
		//				r["beg_label_asym_id"] = asymID;
		//				r["beg_label_seq_id"] = resNum;
		//
		//				cif::tie(asymID, resNum) = r.get("end_auth_asym_id", "end_auth_seq_id");
		//
		//				r["end_label_asym_id"] = asymID;
		//				r["end_label_seq_id"] = resNum;
		//			}
		//			catch (const std::exception& ex)
		//			{
		//				continue;
		//			}
		//		}

		using namespace cif::literals;

		auto &atom_site = *getCategory("atom_site");

		for (auto r : getCategory("struct_conn")->find("pdbx_dist_value"_key == 0 or "pdbx_dist_value"_key == cif::null))
		{
			const auto &[asym1, seq1, atom1, symm1, asym2, seq2, atom2, symm2] = r.get<std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string>(
				"ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_label_atom_id", "ptnr1_symmetry",
				"ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr2_label_atom_id", "ptnr2_symmetry");

			float distance = 1.0f;

			try
			{
				auto a1 = atom_site.find1("label_asym_id"_key == asym1 and "label_seq_id"_key == seq1 and "label_atom_id"_key == atom1);
				auto a2 = atom_site.find1("label_asym_id"_key == asym2 and "label_seq_id"_key == seq2 and "label_atom_id"_key == atom2);

				if (not a1 or not a2)
					throw std::runtime_error("cannot find atom");

				const auto &[x1, y1, z1] = a1.get<float, float, float>("cartn_x", "cartn_y", "cartn_z");
				const auto &[x2, y2, z2] = a2.get<float, float, float>("cartn_x", "cartn_y", "cartn_z");

				if ((symm1.empty() or symm1 == "1_555") and (symm2.empty() or symm2 == "1_555"))
					distance = std::sqrt(
						(x1 - x2) * (x1 - x2) +
						(y1 - y2) * (y1 - y2) +
						(z1 - z2) * (z1 - z2));
				else if (cif::VERBOSE > 0)
					std::cerr << "Cannot calculate distance for link since one of the atoms is in another dimension\n";
			}
			catch (std::exception &ex)
			{
				if (cif::VERBOSE > 0)
					std::cerr << "Error finding atom for LINK distance calculation: " << ex.what() << '\n';
			}

			r["pdbx_dist_value"] = distance;
		}

		result.emplace_back(std::move(mDatablock));
	}
	catch (const std::exception &ex)
	{
		if (cif::VERBOSE >= 0)
		{
			std::cerr << "Error parsing PDB";
			if (mRec != nullptr)
				std::cerr << " at line " << mRec->mLineNr;
			std::cerr << '\n';
		}
		throw;
	}
}

// ----------------------------------------------------------------
// A blast like alignment. Returns index of last aligned residue.

// matrix is m x n, addressing i,j is 0 <= i < m and 0 <= j < n
// element m i,j is mapped to [i * n + j] and thus storage is row major

template <typename T>
class matrix
{
  public:
	using value_type = T;

	matrix() = delete;
	matrix(const matrix &) = delete;
	matrix &operator=(const matrix &) = delete;

	matrix(uint32_t m, uint32_t n, T v = T())
		: m_m(m)
		, m_n(n)
	{
		m_data = new value_type[m_m * m_n];
		std::fill(m_data, m_data + (m_m * m_n), v);
	}

	~matrix()
	{
		delete[] m_data;
	}

	uint32_t dim_m() const { return m_m; }
	uint32_t dim_n() const { return m_n; }

	value_type operator()(uint32_t i, uint32_t j) const
	{
		assert(i < m_m);
		assert(j < m_n);
		return m_data[i * m_n + j];
	}

	value_type &operator()(uint32_t i, uint32_t j)
	{
		assert(i < m_m);
		assert(j < m_n);
		return m_data[i * m_n + j];
	}

  private:
	value_type *m_data;
	uint32_t m_m, m_n;
};

int PDBFileParser::PDBChain::AlignResToSeqRes()
{
	// Use dynamic programming to align the found residues (in ATOM records) against
	// the residues in the SEQRES records in order to obtain the residue numbering.
	// sigh...

	auto &rx = mSeqres;
	auto &ry = mResiduesSeen;

	int dimX = static_cast<int>(mSeqres.size());
	if (dimX == 0)
		throw std::runtime_error(std::string("SEQRES for chain ") + mDbref.chainID + " is empty");

	int dimY = static_cast<int>(mResiduesSeen.size());
	if (dimY == 0)
		throw std::runtime_error(std::string("Number of residues in ATOM records for chain ") + mDbref.chainID + " is zero");

	matrix<float> B(dimX, dimY), Ix(dimX, dimY), Iy(dimX, dimY);
	matrix<int8_t> tb(dimX, dimY);

	int x, y;

	const float
		kMatchReward = 5,
		kMismatchCost = -10,
		kGapOpen = 10, gapExtend = 0.1f;

	float high = 0;
	int highX = 0, highY = 0;

	for (x = 0; x < dimX; ++x)
	{
		for (y = 0; y < dimY; ++y)
		{
			auto &a = rx[x];
			auto &b = ry[y];

			float Ix1 = x > 0 ? Ix(x - 1, y) : 0;
			float Iy1 = y > 0 ? Iy(x, y - 1) : 0;

			// score for alignment
			float M;
			if (a.mMonID == b.mMonID)
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

			if (/*(x == dimX - 1 or y == dimY - 1) and */ high < s)
			{
				high = s;
				highX = x;
				highY = y;
			}
		}
	}

	const int kFlagSeqNr = std::numeric_limits<int>::min();

	// reset positions of seqres
	for (auto &sr : rx)
	{
		sr.mSeqNum = kFlagSeqNr;
		sr.mIcode = ' ';
	}

	// assign numbers
	x = highX;
	y = highY;

	// C++ is getting closer to Pascal :-)
	auto printAlignment = [&tb, highX, highY, &rx, &ry, this]()
	{
		std::cerr << std::string(22, '-') << '\n'
				  << "Alignment for chain " << mDbref.chainID << '\n'
				  << '\n';
		std::vector<std::pair<std::string, std::string>> alignment;

		int x = highX;
		int y = highY;

		for (x = highX, y = highY; x >= 0 and y >= 0;)
		{
			switch (tb(x, y))
			{
				case -1:
					alignment.push_back(make_pair("...", ry[y].mMonID));
					--y;
					break;

				case 1:
					alignment.push_back(make_pair(rx[x].mMonID, "..."));
					--x;
					break;

				case 0:
					alignment.push_back(make_pair(rx[x].mMonID, ry[y].mMonID));
					--x;
					--y;
					break;
			}
		}

		while (x >= 0)
		{
			alignment.push_back(make_pair(rx[x].mMonID, "..."));
			--x;
		}

		while (y >= 0)
		{
			alignment.push_back(make_pair("...", ry[y].mMonID));
			--y;
		}

		reverse(alignment.begin(), alignment.end());
		for (auto a : alignment)
			std::cerr << "  " << a.first << " -- " << a.second << '\n';

		std::cerr << '\n';
	};

	if (cif::VERBOSE > 1)
		printAlignment();

	try
	{
		while (x >= 0 and y >= 0)
		{
			switch (tb(x, y))
			{
				case -1:
					throw std::runtime_error("A residue found in the ATOM records (" + ry[y].mMonID +
											 " @ " + std::string{ mDbref.chainID } + ":" + std::to_string(ry[y].mSeqNum) +
											 ((ry[y].mIcode == ' ' or ry[y].mIcode == 0) ? "" : std::string{ ry[y].mIcode }) +
											 ") was not found in the SEQRES records");
					break;

				case 1:
					if (cif::VERBOSE > 3)
						std::cerr << "Missing residue in ATOM records: " << rx[x].mMonID << " at " << rx[x].mSeqNum << '\n';

					--x;
					break;

				case 0:
					if (rx[x].mMonID != ry[y].mMonID)
					{
						std::cerr << "Warning, unaligned residues at " << x << "/" << y << "(" << rx[x].mMonID << '/' << ry[y].mMonID << ") SEQRES does not agree with ATOM records\n";
						rx[x].mMonID = ry[y].mMonID;
					}

					rx[x].mSeqNum = ry[y].mSeqNum;
					rx[x].mIcode = ry[y].mIcode;

					--x;
					--y;
			}
		}
	}
	catch (const std::exception &ex)
	{
		if (cif::VERBOSE == 1)
			printAlignment();

		throw;
	}

	// assign numbers to the residues that don't have them yet
	std::stack<int> unnumbered;
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
			throw std::runtime_error("Could not assign sequence numbers");
		rx[x].mSeqNum = rx[x + 1].mSeqNum - 1;
		unnumbered.pop();
	}

	return highY;
}

bool PDBFileParser::PDBChain::SameSequence(const PDBChain &rhs) const
{
	bool result = mSeqres.size() == rhs.mSeqres.size();

	for (std::size_t i = 0; result and i < mSeqres.size(); ++i)
		result = mSeqres[i].mMonID == rhs.mSeqres[i].mMonID;

	return result;
}

// --------------------------------------------------------------------

void read_pdb_file(std::istream &pdbFile, cif::file &cifFile)
{
	PDBFileParser p;

	cifFile.load_dictionary("mmcif_pdbx.dic");

	p.Parse(pdbFile, cifFile);

	if (not cifFile.is_valid() and cif::VERBOSE >= 0)
		std::cerr << "Resulting mmCIF file is not valid!\n";
}

// --------------------------------------------------------------------

file read(std::istream &is)
{
	file result;

	auto *buffer = is.rdbuf();
	if (buffer)
	{
		char ch = std::char_traits<char>::to_char_type(buffer->sgetc());

		// All PDB files should always start with a HEADER line
		// and so the very first character in a valid PDB file
		// is 'H'. It is as simple as that.

		// Well, not quite, Unfortunately... People insisted that
		// having only ATOM records also makes up a valid PDB file...
		// Since mmCIF files cannot validly start with a letter character
		// apart from the letter 'd', the test has changed into the following:

		if (std::isalpha(ch) and std::toupper(ch) != 'D')
			read_pdb_file(is, result);
		else
		{
			try
			{
				result.load(is);
			}
			catch (const std::exception &ex)
			{
				std::throw_with_nested(std::runtime_error("Since the file did not start with a valid PDB HEADER line mmCIF was assumed, but that failed."));
			}
		}

		// Since we're using the cif::pdb way of reading the file, the data may need
		// reconstruction
		reconstruct_pdbx(result);
	}

	// Must be a PDB like file, right?
	if (result.get_validator() == nullptr)
		result.load_dictionary("mmcif_pdbx.dic");

	return result;
}

file read(const std::filesystem::path &file)
{
	try
	{
		gzio::ifstream in(file);
		if (not in.is_open())
			throw std::runtime_error("Could not open file " + file.string() + " for input");

		return read(in);
	}
	catch (const std::exception &ex)
	{
		throw_with_nested(std::runtime_error("Error reading file " + file.string()));
	}
}

} // namespace cif::pdb
