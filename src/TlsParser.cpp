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

#include <boost/algorithm/string.hpp>

#include <cif++/TlsParser.hpp>

namespace ba = boost::algorithm;

namespace cif
{

const int
	kResidueNrWildcard = std::numeric_limits<int>::min(),
	kNoSeqNum = std::numeric_limits<int>::max() - 1;

// --------------------------------------------------------------------
// We parse selection statements and create a selection expression tree
// which is then interpreted by setting the selected flag for the
// residues. After that, the selected ranges are collected and printed.

struct TLSResidue
{
	std::string	chainID;
	int			seqNr;
	char		iCode;
	std::string	name;
	bool		selected;
	
	std::string	asymID;
	int			seqID;
	
	bool operator==(const TLSResidue& rhs) const
	{
		return chainID == rhs.chainID and
			seqNr == rhs.seqNr and
			iCode == rhs.iCode and
			iequals(name, rhs.name) and
			selected == rhs.selected;
	}
};

void DumpSelection(const std::vector<TLSResidue>& selected, std::size_t indentLevel)
{
	std::string indent(indentLevel * 2, ' ');
	
	auto i = selected.begin();
	bool first = true;

	// First print in PDB space	
	while (i != selected.end())
	{
		auto b = std::find_if(i, selected.end(), [](auto s) -> bool { return s.selected; });
		if (b == selected.end())
			break;
		
		if (first)
			std::cout << indent << "PDB:" << std::endl;
		first = false;
		
		auto e = std::find_if(b, selected.end(), [b](auto s) -> bool { return s.chainID != b->chainID or not s.selected; });
		
		std::cout << indent << " >> " << b->chainID << ' ' << b->seqNr << ':' << (e - 1)->seqNr << std::endl;
		i = e;
	}
	
	// Then in mmCIF space
	
	if (not first)
		std::cout << indent << "mmCIF:" << std::endl;
	
	i = selected.begin();
	while (i != selected.end())
	{
		auto b = std::find_if(i, selected.end(), [](auto s) -> bool { return s.selected; });
		if (b == selected.end())
			break;
		
		auto e = std::find_if(b, selected.end(), [b](auto s) -> bool { return s.asymID != b->asymID or not s.selected; });
	
		std::string asymID = b->asymID;
		int from = b->seqID, to = from;

		for (auto j = b + 1; j != e; ++j)
		{
			if (j->seqID == to + 1)
				to = j->seqID;
			else if (j->seqID != to)		// probably an insertion code
			{
				if (from == kNoSeqNum or to == kNoSeqNum)
					std::cout << indent << " >> " << asymID << std::endl;
				else
					std::cout << indent << " >> " << asymID << ' ' << from << ':' << to << std::endl;
				asymID = b->asymID;
				from = to = b->seqID;
			}
		}

		if (from == kNoSeqNum or to == kNoSeqNum)
			std::cout << indent << " >> " << asymID << std::endl;
		else
			std::cout << indent << " >> " << asymID << ' ' << from << ':' << to << std::endl;

		i = e;
	}
	
	if (first)
	{
		if (isatty(STDOUT_FILENO))
			std::cout << indent << cif::coloured("Empty selection") << std::endl;
		else
			std::cout << indent << cif::coloured("Empty selection") << std::endl;
	}
}

std::vector<std::tuple<std::string,int,int>> TLSSelection::GetRanges(Datablock& db, bool pdbNamespace) const
{
	std::vector<TLSResidue> selected;

	// Collect the residues from poly seq scheme...		
	for (auto r: db["pdbx_poly_seq_scheme"])
	{
		std::string chain, seqNr, iCode, name;

		std::string asymID;
		int seqID;
		
		if (pdbNamespace)
			cif::tie(chain, seqNr, iCode, name, asymID, seqID) = r.get("pdb_strand_id", "pdb_seq_num", "pdb_ins_code", "pdb_comp_id", "asym_id", "seq_id");
		else
		{
			cif::tie(chain, seqNr, name) = r.get("asym_id", "seq_id", "mon_id");
			asymID = chain;
			seqID = stoi(seqNr);
		}
		
		if (seqNr.empty())
			continue;

		if (iCode.length() > 1)
			throw std::runtime_error("invalid iCode");
		
		selected.push_back({chain, stoi(seqNr), iCode[0], name, false, asymID, seqID});
	}

	// ... those from the nonpoly scheme
	for (auto r: db["pdbx_nonpoly_scheme"])
	{
		std::string chain, seqNr, iCode, name, asymID;

		if (pdbNamespace)
		{
			cif::tie(chain, seqNr, iCode, name, asymID) = r.get("pdb_strand_id", "pdb_seq_num", "pdb_ins_code", "pdb_mon_id", "asym_id");
			if (seqNr.empty())
				continue;
		}
		else
		{
			cif::tie(chain, name) = r.get("asym_id", "mon_id");
			asymID = chain;
			seqNr = "0";
		}
		
		if (iequals(name, "HOH") or iequals(name, "H2O"))
			continue;
		
		if (iCode.length() > 1)
			throw std::runtime_error("invalid iCode");
		
		selected.push_back({chain, stoi(seqNr), iCode[0], name, false, asymID, kNoSeqNum});
	}
	
	// selected might consist of multiple ranges
	// output per chain

	stable_sort(selected.begin(), selected.end(), [](auto& a, auto& b) -> bool
		{
			int d = a.chainID.compare(b.chainID);
			if (d == 0)
				d = a.seqNr - b.seqNr;
			return d < 0;
		});
	
	CollectResidues(db, selected);
	
	std::vector<std::tuple<std::string,int,int>> result;

	auto i = selected.begin();
	
	while (i != selected.end())
	{
		auto b = std::find_if(i, selected.end(), [](auto s) -> bool { return s.selected; });
		if (b == selected.end())
			break;
		
		auto e = std::find_if(b, selected.end(), [b](auto s) -> bool { return s.asymID != b->asymID or not s.selected; });

		// return ranges with strict increasing sequence numbers.
		// So when there's a gap in the sequence we split the range.
		// Beware of iCodes though
		result.push_back(make_tuple(b->asymID, b->seqID, b->seqID));
		for (auto j = b + 1; j != e; ++j)
		{
			if (j->seqID == std::get<2>(result.back()) + 1)
				std::get<2>(result.back()) = j->seqID;
			else if (j->seqID != std::get<2>(result.back()))		// probably an insertion code
				result.push_back(make_tuple(b->asymID, j->seqID, j->seqID));
		}

		i = e;
	}
	
	return result;
}

struct TLSSelectionNot : public TLSSelection
{
	TLSSelectionNot(TLSSelectionPtr selection)
		: selection(selection.release()) {}
	
	virtual void CollectResidues(Datablock& db, std::vector<TLSResidue>& residues, std::size_t indentLevel) const
	{
		selection->CollectResidues(db, residues, indentLevel + 1);
		
		for (auto& r: residues)
			r.selected = not r.selected;

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "NOT" << std::endl;
			DumpSelection(residues, indentLevel);
		}
	} 

	TLSSelectionPtr selection;
};

struct TLSSelectionAll : public TLSSelection
{
	TLSSelectionAll() {}
	
	virtual void CollectResidues(Datablock& db, std::vector<TLSResidue>& residues, std::size_t indentLevel) const
	{
		for (auto& r: residues)
			r.selected = true;

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "ALL" << std::endl;
			DumpSelection(residues, indentLevel);
		}
	} 
};

struct TLSSelectionChain : public TLSSelectionAll
{
	TLSSelectionChain(const std::string& chainID)
		: m_chain(chainID) {}
	
	virtual void CollectResidues(Datablock& db, std::vector<TLSResidue>& residues, std::size_t indentLevel) const
	{
		bool allChains = m_chain == "*";
		
		for (auto& r: residues)
			r.selected = allChains or r.chainID == m_chain;

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "CHAIN " << m_chain << std::endl;
			DumpSelection(residues, indentLevel);
		}
	} 

	std::string	m_chain;
};

struct TLSSelectionResID : public TLSSelectionAll
{
	TLSSelectionResID(int seqNr, char iCode)
		: m_seq_nr(seqNr), m_icode(iCode) {}
	
	virtual void CollectResidues(Datablock& db, std::vector<TLSResidue>& residues, std::size_t indentLevel) const
	{
		for (auto& r: residues)
			r.selected = r.seqNr == m_seq_nr and r.iCode == m_icode;

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "ResID " << m_seq_nr << (m_icode ? std::string { m_icode} : "") << std::endl;
			DumpSelection(residues, indentLevel);
		}
	} 

	int		m_seq_nr;
	char	m_icode;
};

struct TLSSelectionRangeSeq : public TLSSelectionAll
{
	TLSSelectionRangeSeq(int first, int last)
		: m_first(first), m_last(last) {}
	
	virtual void CollectResidues(Datablock& db, std::vector<TLSResidue>& residues, std::size_t indentLevel) const
	{
		for (auto& r: residues)
		{
			r.selected = ((r.seqNr >= m_first or m_first == kResidueNrWildcard) and
						  (r.seqNr <= m_last or m_last == kResidueNrWildcard));
		}

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Range " << m_first << ':' << m_last << std::endl;
			DumpSelection(residues, indentLevel);
		}
	} 

	int		m_first, m_last;
};

struct TLSSelectionRangeID : public TLSSelectionAll
{
	TLSSelectionRangeID(int first, int last, char icodeFirst = 0, char icodeLast = 0)
		: m_first(first), m_last(last), m_icode_first(icodeFirst), m_icode_last(icodeLast) {}
	
	virtual void CollectResidues(Datablock& db, std::vector<TLSResidue>& residues, std::size_t indentLevel) const
	{
		// need to do this per chain
		std::set<std::string> chains;
		for (auto& r: residues)
			chains.insert(r.chainID);
		
		for (std::string chain: chains)
		{
			auto f = std::find_if(residues.begin(), residues.end(),
				[=,this](auto r) -> bool {
					return r.chainID == chain and r.seqNr == m_first and r.iCode == m_icode_first;
				});
			
			auto l = std::find_if(residues.begin(), residues.end(),
				[=,this](auto r) -> bool {
					return r.chainID == chain and r.seqNr == m_last and r.iCode == m_icode_last;
				});
			
			if (f != residues.end() and l != residues.end() and f <= l)
			{
				++l;
				
				for (; f != l; ++f)
					f->selected = true;
			}
		}

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Through " << m_first << ':' << m_last << std::endl;
			DumpSelection(residues, indentLevel);
		}
	} 

	int		m_first, m_last;
	char	m_icode_first, m_icode_last;
};

struct TLSSelectionUnion : public TLSSelection
{
	TLSSelectionUnion(TLSSelectionPtr& lhs, TLSSelectionPtr& rhs)
		: lhs(lhs.release()), rhs(rhs.release()) {}
	
	TLSSelectionUnion(TLSSelectionPtr& lhs, TLSSelectionPtr&& rhs)
		: lhs(lhs.release()), rhs(rhs.release()) {}
	
	virtual void CollectResidues(Datablock& db, std::vector<TLSResidue>& residues, std::size_t indentLevel) const
	{
		auto a = residues;
		for_each(a.begin(), a.end(), [](auto& r) { r.selected = false; });
		
		auto b = residues;
		for_each(b.begin(), b.end(), [](auto& r) { r.selected = false; });
		
		lhs->CollectResidues(db, a, indentLevel + 1);
		rhs->CollectResidues(db, b, indentLevel + 1);
		
		for (auto ai = a.begin(), bi = b.begin(), ri = residues.begin(); ri != residues.end(); ++ai, ++bi, ++ri)
			ri->selected = ai->selected or bi->selected;

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Union" << std::endl;
			DumpSelection(residues, indentLevel);
		}
	} 
	
	TLSSelectionPtr	lhs;
	TLSSelectionPtr	rhs;
};

struct TLSSelectionIntersection : public TLSSelection
{
	TLSSelectionIntersection(TLSSelectionPtr& lhs, TLSSelectionPtr& rhs)
		: lhs(lhs.release()), rhs(rhs.release()) {}
	
	TLSSelectionIntersection(TLSSelectionPtr& lhs, TLSSelectionPtr&& rhs)
		: lhs(lhs.release()), rhs(rhs.release()) {}
	
	virtual void CollectResidues(Datablock& db, std::vector<TLSResidue>& residues, std::size_t indentLevel) const
	{
		auto a = residues;
		for_each(a.begin(), a.end(), [](auto& r) { r.selected = false; });
		
		auto b = residues;
		for_each(b.begin(), b.end(), [](auto& r) { r.selected = false; });
		
		lhs->CollectResidues(db, a, indentLevel + 1);
		rhs->CollectResidues(db, b, indentLevel + 1);
		
		for (auto ai = a.begin(), bi = b.begin(), ri = residues.begin(); ri != residues.end(); ++ai, ++bi, ++ri)
			ri->selected = ai->selected and bi->selected;

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Intersection" << std::endl;
			DumpSelection(residues, indentLevel);
		}
	} 
	
	TLSSelectionPtr	lhs;
	TLSSelectionPtr	rhs;
};

struct TLSSelectionByName : public TLSSelectionAll
{
  public:
	TLSSelectionByName(const std::string& resname)
		: m_name(resname) {}

	virtual void CollectResidues(Datablock& db, std::vector<TLSResidue>& residues, std::size_t indentLevel) const
	{
		for (auto& r: residues)
			r.selected = r.name == m_name;

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Name " << m_name << std::endl;
			DumpSelection(residues, indentLevel);
		}
	} 

	std::string	m_name;
};

struct TLSSelectionByElement : public TLSSelectionAll
{
  public:
	TLSSelectionByElement(const std::string& element)
		: m_element(element) {}

	virtual void CollectResidues(Datablock& db, std::vector<TLSResidue>& residues, std::size_t indentLevel) const
	{
		// rationale... We want to select residues only. So we select
		// residues that have just a single atom of type m_element.
		// And we assume these have as residue name... m_element.
		// ... Right?
		
		for (auto& r: residues)
			r.selected = iequals(r.name, m_element);

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Element " << m_element << std::endl;
			DumpSelection(residues, indentLevel);
		}
	} 

	std::string	m_element;
};

// --------------------------------------------------------------------

class TLSSelectionParserImpl
{
  public:
	TLSSelectionParserImpl(const std::string& selection)
		: m_selection(selection), m_p(m_selection.begin()), m_end(m_selection.end()) {}
	
	virtual TLSSelectionPtr Parse() = 0;
	
  protected:

	virtual int GetNextToken() = 0;
	virtual void Match(int token);
	virtual std::string ToString(int token) = 0;

	std::string				m_selection;
	std::string::iterator	m_p, m_end;
	int						m_lookahead = 0;
	std::string				m_token;
};


void TLSSelectionParserImpl::Match(int token)
{
	if (m_lookahead == token)
		m_lookahead = GetNextToken();
	else
	{
		std::string expected;
		if (token >= 256)
			expected = ToString(token);
		else
			expected = { char(token) };
		
		std::string found;
		if (m_lookahead >= 256)
			found = ToString(m_lookahead) + " (" + m_token + ')';
		else
			found = { char(m_lookahead) };
	
		throw std::runtime_error("Expected " + expected + " but found " + found);
	}
}

// --------------------------------------------------------------------

class TLSSelectionParserImplPhenix : public TLSSelectionParserImpl
{
  public:
	TLSSelectionParserImplPhenix(const std::string& selection)
		: TLSSelectionParserImpl(selection)
	{
		m_lookahead = GetNextToken();
	}
	
	virtual TLSSelectionPtr Parse();

  private:

	TLSSelectionPtr ParseAtomSelection();
	TLSSelectionPtr	ParseTerm();
	TLSSelectionPtr	ParseFactor();

	enum TOKEN {
		pt_NONE		= 0,
		pt_IDENT	= 256,
		pt_STRING,
		pt_NUMBER,
		pt_RESID,
		pt_EOLN,
		pt_KW_ALL,
		pt_KW_CHAIN,
		pt_KW_RESSEQ,
		pt_KW_RESID,
		pt_KW_ICODE,
		pt_KW_RESNAME,
		pt_KW_ELEMENT,
		pt_KW_AND,
		pt_KW_OR,
		pt_KW_NOT,
		pt_KW_PDB,
		pt_KW_ENTRY,
		pt_KW_THROUGH
	};

	virtual int GetNextToken();
	virtual std::string ToString(int token);
	
	int m_value_i;
	std::string m_value_s;
	char m_icode;
};

int TLSSelectionParserImplPhenix::GetNextToken()
{
	int result = pt_NONE;
	enum STATE {
		st_START,
		st_RESID	= 200,
		st_NUM		= 300,
		st_IDENT	= 400,
		st_QUOTED	= 500,
		st_DQUOTED	= 550,
		st_OTHER	= 600
	};
	int state = st_START;
	
	m_value_i = 0;
	m_icode = 0;
	m_value_s.clear();
	auto s = m_p;
	
	auto start = state;
	m_token.clear();
	
	auto restart = [&]()
	{
		switch (start)
		{
			case st_START:	state = start = st_RESID;	break;
			case st_RESID:	state = start = st_NUM;		break;
			case st_NUM:	state = start = st_IDENT;	break;
			case st_IDENT:	state = start = st_QUOTED;	break;
			case st_QUOTED:	state = start = st_DQUOTED;	break;
			case st_DQUOTED:state = start = st_OTHER;	break;
		}
		m_token.clear();
		m_p = s;
	};
	
	auto retract = [&]()
	{
		--m_p;
		m_token.pop_back();
	};
	
	while (result == pt_NONE)
	{
		char ch = *m_p++;
		if (m_p > m_end)
			ch = 0;
		else
			m_token += ch;
		
		switch (state)
		{
			// start block
			case st_START:
				if (ch == 0)
					result = pt_EOLN;
				else if (isspace(ch))
				{
					m_token.clear();
					++s;
				}
				else
					restart();
				break;

			// RESID block
			case st_RESID:
				if (ch == '-')
					state = st_RESID + 1;
				else if (isdigit(ch))
				{
					m_value_i = (ch - '0');
					state = st_RESID + 2;
				}
				else
					restart();
				break;
			
			case st_RESID + 1:
				if (isdigit(ch))
				{
					m_value_i = -(ch - '0');
					state = st_RESID + 2;
				}
				else
					restart();
				break;
			
			case st_RESID + 2:
				if (isdigit(ch))
					m_value_i = 10 * m_value_i + (m_value_i < 0 ? -1 : 1) * (ch - '0');
				else if (isalpha(ch))
				{
					m_icode = ch;
					state = st_RESID + 3;
				}
				else
					restart();
				break;
			
			case st_RESID + 3:
				if (isalnum(ch))
					restart();
				else
				{
					retract();
					result = pt_RESID;
				}
				break;
			
			// NUM block
			
			case st_NUM:
				if (ch == '-')
					state = st_NUM + 1;
				else if (isdigit(ch))
				{
					m_value_i = ch - '0';
					state = st_NUM + 2;
				}
				else
					restart();
				break;
			
			case st_NUM + 1:
				if (isdigit(ch))
				{
					m_value_i = -(ch - '0');
					state = st_NUM + 2;
				}
				else
					restart();
				break;
			
			case st_NUM + 2:
				if (isdigit(ch))
					m_value_i = 10 * m_value_i + (m_value_i < 0 ? -1 : 1) * (ch - '0');
				else if (not isalpha(ch))
				{
					result = pt_NUMBER;
					retract();
				}
				else
					restart();
				break;
			
			// IDENT block
			
			case st_IDENT:
				if (isalnum(ch))
				{
					m_value_s = { ch };
					state = st_IDENT + 1;
				}
				else
					restart();
				break;
			
			case st_IDENT + 1:
				if (isalnum(ch) or ch == '\'')
					m_value_s += ch;
				else
				{
					--m_p;
					result = pt_IDENT;
				}
				break;
			
			// QUOTED block
			
			case st_QUOTED:
				if (ch == '\'')
				{
					m_value_s.clear();
					state = st_QUOTED + 1;
				}
				else
					restart();
				break;
			
			case st_QUOTED + 1:
				if (ch == '\'')
					result = pt_STRING;
				else if (ch == 0)
					throw std::runtime_error("Unexpected end of selection, missing quote character?");
				else
					m_value_s += ch;
				break;
			
			// QUOTED block
			
			case st_DQUOTED:
				if (ch == '\"')
				{
					m_value_s.clear();
					state = st_DQUOTED + 1;
				}
				else
					restart();
				break;
			
			case st_DQUOTED + 1:
				if (ch == '\"')
					result = pt_STRING;
				else if (ch == 0)
					throw std::runtime_error("Unexpected end of selection, missing quote character?");
				else
					m_value_s += ch;
				break;
			
			// OTHER block
			case st_OTHER:
				result = ch;
				break;
		}
	}
	
	if (result == pt_IDENT)
	{
		if (iequals(m_value_s, "CHAIN"))
			result = pt_KW_CHAIN;
		else if (iequals(m_value_s, "ALL"))
			result = pt_KW_ALL;
		else if (iequals(m_value_s, "AND"))
			result = pt_KW_AND;
		else if (iequals(m_value_s, "OR"))
			result = pt_KW_OR;
		else if (iequals(m_value_s, "NOT"))
			result = pt_KW_NOT;
		else if (iequals(m_value_s, "RESSEQ"))
			result = pt_KW_RESSEQ;
		else if (iequals(m_value_s, "RESID") or iequals(m_value_s, "RESI"))
			result = pt_KW_RESID;
		else if (iequals(m_value_s, "RESNAME"))
			result = pt_KW_RESNAME;
		else if (iequals(m_value_s, "ELEMENT"))
			result = pt_KW_ELEMENT;
		else if (iequals(m_value_s, "PDB"))
			result = pt_KW_PDB;
		else if (iequals(m_value_s, "ENTRY"))
			result = pt_KW_ENTRY;
		else if (iequals(m_value_s, "THROUGH"))
			result = pt_KW_THROUGH;
	}
		
	return result;
}

std::string TLSSelectionParserImplPhenix::ToString(int token)
{
	switch (token)
	{
		case pt_IDENT:		return "identifier";
		case pt_STRING:		return "string";
		case pt_NUMBER:		return "number";
		case pt_RESID:		return "resid";
		case pt_EOLN:		return "end of line";

		case pt_KW_ALL:		return "ALL";
		case pt_KW_CHAIN:	return "CHAIN";
		case pt_KW_RESSEQ:	return "RESSEQ";
		case pt_KW_RESID:	return "RESID";
		case pt_KW_RESNAME:	return "RESNAME";
		case pt_KW_ELEMENT:	return "ELEMENT";
		case pt_KW_AND:		return "AND";
		case pt_KW_OR:		return "OR";
		case pt_KW_NOT:		return "NOT";
		case pt_KW_PDB:		return "PDB";
		case pt_KW_ENTRY:	return "ENTRY";
		case pt_KW_THROUGH:	return "THROUGH";

		default:			return "character";
	}
}

TLSSelectionPtr TLSSelectionParserImplPhenix::Parse()
{
	if (m_lookahead == pt_KW_PDB)
	{
		Match(pt_KW_PDB);
//		Match(pt_KW_ENTRY);
		
		throw std::runtime_error("Unimplemented PDB ENTRY specification");
	}

	TLSSelectionPtr result = ParseAtomSelection();

	bool extraParenthesis = false;

	if (m_lookahead == ')')
	{
		extraParenthesis = true;
		m_lookahead = GetNextToken();
	}

	Match(pt_EOLN);
	
	if (extraParenthesis and cif::VERBOSE > 0)
		std::cerr << "WARNING: too many closing parenthesis in TLS selection statement" << std::endl;
	
	return result;
}

TLSSelectionPtr TLSSelectionParserImplPhenix::ParseAtomSelection()
{
	TLSSelectionPtr result = ParseTerm();
	
	while (m_lookahead == pt_KW_OR)
	{
		Match(pt_KW_OR);
		result.reset(new TLSSelectionUnion(result, ParseTerm()));
	}
	
	return result;
}

TLSSelectionPtr TLSSelectionParserImplPhenix::ParseTerm()
{
	TLSSelectionPtr result = ParseFactor();
	
	while (m_lookahead == pt_KW_AND)
	{
		Match(pt_KW_AND);
		result.reset(new TLSSelectionIntersection(result, ParseFactor()));
	}
	
	return result;
}

TLSSelectionPtr TLSSelectionParserImplPhenix::ParseFactor()
{
	TLSSelectionPtr result;

	switch (m_lookahead)
	{
		case '(':
			Match('(');
			result = ParseAtomSelection();
			if (m_lookahead == pt_EOLN and cif::VERBOSE > 0)
				std::cerr << "WARNING: missing closing parenthesis in TLS selection statement" << std::endl;
			else
				Match(')');
			break;
		
		case pt_KW_NOT:
			Match(pt_KW_NOT);
			result.reset(new TLSSelectionNot(ParseAtomSelection()));
			break;
		
		case pt_KW_CHAIN:
		{
			Match(pt_KW_CHAIN);
			
			std::string chainID = m_value_s;
			if (m_lookahead == pt_NUMBER)	// sigh
			{
				chainID = std::to_string(m_value_i);
				Match(pt_NUMBER);
			}
			else
				Match(m_lookahead == pt_STRING ? pt_STRING : pt_IDENT);
			
			result.reset(new TLSSelectionChain(chainID));
			break;
		}
		
		case pt_KW_RESNAME:
		{
			Match(pt_KW_RESNAME);
			std::string name = m_value_s;
			Match(pt_IDENT);
			result.reset(new TLSSelectionByName(name));
			break;
		}
		
		case pt_KW_ELEMENT:
		{
			Match(pt_KW_ELEMENT);
			std::string element = m_value_s;
			Match(pt_IDENT);
			result.reset(new TLSSelectionByElement(element));
			break;
		}
		
		case pt_KW_RESSEQ:
		{
			Match(pt_KW_RESSEQ);
			
			int from = m_value_i;
			Match(pt_NUMBER);
			
			int to = from;
			if (m_lookahead == ':')
			{
				Match(':');
				to = m_value_i;
				Match(pt_NUMBER);
			}
			
			result.reset(new TLSSelectionRangeSeq(from, to));
			break;
		}
		
		case pt_KW_RESID:
		{
			Match(pt_KW_RESID);
			
			int from, to;
			char icode_from = 0, icode_to = 0;
			bool through = false;
			
			from = to = m_value_i;

			if (m_lookahead == pt_NUMBER)
				Match(pt_NUMBER);
			else
			{
				icode_from = m_icode;
				Match(pt_RESID);
			}

			if (m_lookahead == ':' or m_lookahead == pt_KW_THROUGH or m_lookahead == '-')
			{
				through = m_lookahead == pt_KW_THROUGH;

				Match(m_lookahead);

				to = m_value_i;
				if (m_lookahead == pt_NUMBER)
					Match(pt_NUMBER);
				else
				{
					icode_to = m_icode;
					Match(pt_RESID);
				}

				if (through)
					result.reset(new TLSSelectionRangeID(from, to, icode_from, icode_to));
				else
				{
					if (cif::VERBOSE > 0 and (icode_from or icode_to))
						std::cerr << "Warning, ignoring insertion codes" << std::endl;
					
					result.reset(new TLSSelectionRangeSeq(from, to));
				}
			}
			else
				result.reset(new TLSSelectionResID(from, icode_from));
			
			break;
		}
		
		case pt_KW_ALL:
			Match(pt_KW_ALL);
			result.reset(new TLSSelectionAll());
			break;
		
		default:
			throw std::runtime_error("Unexpected token " + ToString(m_lookahead) + " (" + m_token + ')');
	}
	
	return result;
}

// --------------------------------------------------------------------

class TLSSelectionParserImplBuster : public TLSSelectionParserImpl
{
  public:
	TLSSelectionParserImplBuster(const std::string& selection);
	
	virtual TLSSelectionPtr Parse();
	
  protected:
	
	enum TOKEN {
		bt_NONE	= 0,
		bt_IDENT	= 256,
		bt_NUMBER,
		bt_EOLN,
	};
	
	virtual int GetNextToken();
	virtual std::string ToString(int token);
		
	TLSSelectionPtr ParseGroup();
	std::tuple<std::string,int> ParseAtom();

	TLSSelectionPtr ParseOldGroup();
	
	int m_value_i;
	std::string m_value_s;
	bool m_parsing_old_style = false;
};

TLSSelectionParserImplBuster::TLSSelectionParserImplBuster(const std::string& selection)
	: TLSSelectionParserImpl(selection)
{
	m_lookahead = GetNextToken();
}

int TLSSelectionParserImplBuster::GetNextToken()
{
	int result = bt_NONE;
	enum STATE { st_START, st_NEGATE, st_NUM, st_IDENT } state = st_START;
	
	m_value_i = 0;
	m_value_s.clear();
	bool negative = false;
	
	while (result == bt_NONE)
	{
		char ch = *m_p++;
		if (m_p > m_end)
			ch = 0;
		
		switch (state)
		{
			case st_START:
				if (ch == 0)
					result = bt_EOLN;
				else if (isspace(ch))
					continue;
				else if (isdigit(ch))
				{
					m_value_i = ch - '0';
					state = st_NUM;
				}
				else if (isalpha(ch))
				{
					m_value_s = { ch };
					state = st_IDENT;
				}
				else if (ch == '-')
				{
					state = st_NEGATE;
				}
				else
					result = ch;
				break;
			
			case st_NEGATE:
				if (isdigit(ch))
				{
					m_value_i = ch - '0';
					state = st_NUM;
					negative = true;
				}
				else
				{
					--m_p;
					result = '-';
				}
				break;
			
			case st_NUM:
				if (isdigit(ch))
					m_value_i = 10 * m_value_i + (ch - '0');
				else
				{
					if (negative)
						m_value_i = -m_value_i;

					result = bt_NUMBER;
					--m_p;
				}
				break;
			
			case st_IDENT:
				if (isalnum(ch))
					m_value_s += ch;
				else
				{
					--m_p;
					result = bt_IDENT;
				}
				break;
		}
	}
	
	return result;
}

std::string TLSSelectionParserImplBuster::ToString(int token)
{
	switch (token)
	{
		case bt_IDENT:	return "identifier (" + m_value_s + ')';
		case bt_NUMBER:	return "number (" + std::to_string(m_value_i) + ')';
		case bt_EOLN:	return "end of line";

		default:
			assert(false);
			return "unknown token";
	}
}

TLSSelectionPtr TLSSelectionParserImplBuster::ParseGroup()
{
	TLSSelectionPtr result;
	
	auto add = [&result](const std::string& chainID, int from, int to)
	{
		TLSSelectionPtr sc(new TLSSelectionChain(chainID));
		TLSSelectionPtr sr(new TLSSelectionRangeSeq(from, to));
		TLSSelectionPtr s(new TLSSelectionIntersection(sc, sr));
		
		if (result == nullptr)
			result.reset(s.release());
		else
			result.reset(new TLSSelectionUnion{result, s });
	};
	
	Match('{');

	do
	{
		std::string chain1;
		int seqNr1;
		std::tie(chain1, seqNr1) = ParseAtom();
		
		if (m_lookahead == '-')
		{
			std::string chain2;
			int seqNr2 = seqNr1;

			Match('-');
			
			if (m_lookahead == bt_NUMBER)
			{
				seqNr2 = m_value_i;
				Match(bt_NUMBER);
			}
			else
			{
				std::tie(chain2, seqNr2) = ParseAtom();
				if (chain1 != chain2)
				{
					if (cif::VERBOSE > 0)
						std::cerr << "Warning, ranges over multiple chains detected" << std::endl;
					
					TLSSelectionPtr sc1(new TLSSelectionChain(chain1));
					TLSSelectionPtr sr1(new TLSSelectionRangeSeq(seqNr1, kResidueNrWildcard));
					TLSSelectionPtr s1(new TLSSelectionIntersection(sc1, sr1));

					TLSSelectionPtr sc2(new TLSSelectionChain(chain2));
					TLSSelectionPtr sr2(new TLSSelectionRangeSeq(kResidueNrWildcard, seqNr2));
					TLSSelectionPtr s2(new TLSSelectionIntersection(sc2, sr2));
					
					TLSSelectionPtr s(new TLSSelectionUnion(s1, s2));
					
					if (result == nullptr)
						result.reset(s.release());
					else
						result.reset(new TLSSelectionUnion{result, s });
					
					chain1.clear();
				}
			}

			if (not chain1.empty())
				add(chain1, seqNr1, seqNr2);
		}
		else
			add(chain1, seqNr1, seqNr1);
	}
	while (m_lookahead != '}');
	
	Match('}');
	
	return result;
}

std::tuple<std::string,int> TLSSelectionParserImplBuster::ParseAtom()
{
	std::string chain = m_value_s;
	int seqNr = kResidueNrWildcard;
	
	if (m_lookahead == '*')
		Match('*');
	else
		Match(bt_IDENT);
	
	Match('|');
	
	if (m_lookahead == '*')
		Match('*');
	else
	{
		seqNr = m_value_i;
		Match(bt_NUMBER);
		
		if (m_lookahead == ':')
		{
			Match(':');
			std::string atom = m_value_s;
			
			if (cif::VERBOSE > 0)
				std::cerr << "Warning: ignoring atom ID '" << atom << "' in TLS selection" << std::endl;
			
			Match(bt_IDENT);
		}
	}
	
	return make_tuple(chain, seqNr);
}

TLSSelectionPtr TLSSelectionParserImplBuster::Parse()
{
	TLSSelectionPtr result = ParseGroup();
	Match(bt_EOLN);
	return result;
}

// --------------------------------------------------------------------

class TLSSelectionParserImplBusterOld : public TLSSelectionParserImpl
{
  public:
	TLSSelectionParserImplBusterOld(const std::string& selection)
		: TLSSelectionParserImpl(selection)
	{
		m_lookahead = GetNextToken();
	}
	
	virtual TLSSelectionPtr Parse();

  private:

	TLSSelectionPtr ParseAtomSelection();
	TLSSelectionPtr	ParseTerm();
	TLSSelectionPtr	ParseFactor();
	
	TLSSelectionPtr ParseResid();
	TLSSelectionPtr ParseChainResid();

	enum TOKEN {
		pt_NONE		= 0,
		pt_IDENT	= 256,
		pt_CHAINRESID,
		pt_STRING,
		pt_NUMBER,
		pt_RANGE,
		pt_EOLN,
		
		pt_KW_ALL,
		pt_KW_CHAIN,
		pt_KW_RESSEQ,
		pt_KW_RESID,
		pt_KW_RESNAME,
		pt_KW_ELEMENT,
		pt_KW_AND,
		pt_KW_OR,
		pt_KW_NOT,
		pt_KW_PDB,
		pt_KW_ENTRY,
		pt_KW_THROUGH
	};

	virtual int GetNextToken();
	virtual std::string ToString(int token);
	
	int m_value_i;
	std::string m_value_s;
	int m_value_r[2];
};

int TLSSelectionParserImplBusterOld::GetNextToken()
{
	int result = pt_NONE;
	enum STATE { st_START, st_NEGATE, st_NUM, st_RANGE, st_IDENT_1, st_IDENT, st_CHAINRESID, st_QUOTED_1, st_QUOTED_2 } state = st_START;
	
	m_value_i = 0;
	m_value_s.clear();
	
	bool negative = false;
	
	while (result == pt_NONE)
	{
		char ch = *m_p++;
		if (m_p > m_end)
			ch = 0;
		
		switch (state)
		{
			case st_START:
				if (ch == 0)
					result = pt_EOLN;
				else if (isspace(ch))
					continue;
				else if (isdigit(ch))
				{
					m_value_i = ch - '0';
					state = st_NUM;
				}
				else if (isalpha(ch))
				{
					m_value_s = { ch };
					state = st_IDENT_1;
				}
				else if (ch == '-')
				{
					state = st_NEGATE;
				}
				else if (ch == '\'')
				{
					state = st_QUOTED_1;
				}
				else
					result = ch;
				break;
			
			case st_NEGATE:
				if (isdigit(ch))
				{
					m_value_i = ch - '0';
					state = st_NUM;
					negative = true;
				}
				else
				{
					--m_p;
					result = '-';
				}
				break;
			
			case st_NUM:
				if (isdigit(ch))
					m_value_i = 10 * m_value_i + (ch - '0');
				else if (ch == '-' or ch == ':')
				{
					if (negative)
						m_value_i = -m_value_i;

					m_value_r[0] = m_value_i;
					m_value_r[1] = 0;
					state = st_RANGE;
				}
				else
				{
					if (negative)
						m_value_i = -m_value_i;

					result = pt_NUMBER;
					--m_p;
				}
				break;
			
			case st_RANGE:		// TODO: question, is "-2--1" a valid range? We do not support that, yet
				if (isdigit(ch))
					m_value_r[1] = 10 * m_value_r[1] + (ch - '0');
				else if (m_value_r[1] != 0)
				{
					result = pt_RANGE;
					--m_p;
				}
				else
				{
					--m_p;
					--m_p;
					result = pt_NUMBER;
				}
				break;
			
			case st_IDENT_1:
				if (isalpha(ch))
				{
					m_value_s += ch;
					state = st_IDENT;
				}
				else if (isdigit(ch))
				{
					m_value_i = (ch - '0');
					state = st_CHAINRESID;
				}
				else
				{
					--m_p;
					result = pt_IDENT;
				}
				break;
			
			case st_CHAINRESID:
				if (isalpha(ch))
				{
					m_value_s += std::to_string(m_value_i);
					m_value_s += ch;
					state = st_IDENT;
				}
				else if (isdigit(ch))
					m_value_i = 10 * m_value_i + (ch - '0');
				else
				{
					--m_p;
					result = pt_CHAINRESID;
				}
				break;
			
			case st_IDENT:
				if (isalnum(ch))
					m_value_s += ch;
				else
				{
					--m_p;
					result = pt_IDENT;
				}
				break;
			
			case st_QUOTED_1:
				if (ch == '\'')
				{
					--m_p;
					result = '\'';
				}
				else
				{
					m_value_s = { ch };
					state = st_QUOTED_2;
				}
				break;
			
			case st_QUOTED_2:
				if (ch == '\'')
					result = pt_STRING;
				else if (ch == 0)
					throw std::runtime_error("Unexpected end of selection, missing quote character?");
				else
					m_value_s += ch;
				break;
				
		}
	}
	
	if (result == pt_IDENT)
	{
		if (iequals(m_value_s, "CHAIN"))
			result = pt_KW_CHAIN;
		else if (iequals(m_value_s, "ALL"))
			result = pt_KW_ALL;
		else if (iequals(m_value_s, "AND"))
			result = pt_KW_AND;
		else if (iequals(m_value_s, "OR"))
			result = pt_KW_OR;
		else if (iequals(m_value_s, "NOT"))
			result = pt_KW_NOT;
		else if (iequals(m_value_s, "RESSEQ"))
			result = pt_KW_RESSEQ;
		else if (iequals(m_value_s, "RESID") or iequals(m_value_s, "RESI") or iequals(m_value_s, "RESIDUES"))
			result = pt_KW_RESID;
		else if (iequals(m_value_s, "RESNAME"))
			result = pt_KW_RESNAME;
		else if (iequals(m_value_s, "PDB"))
			result = pt_KW_PDB;
		else if (iequals(m_value_s, "ENTRY"))
			result = pt_KW_ENTRY;
		else if (iequals(m_value_s, "THROUGH"))
			result = pt_KW_THROUGH;
	}
	
	return result;
}

std::string TLSSelectionParserImplBusterOld::ToString(int token)
{
	switch (token)
	{
		case pt_IDENT:	return "identifier (" + m_value_s + ')';
		case pt_STRING:	return "string (" + m_value_s + ')';
		case pt_NUMBER:	return "number (" + std::to_string(m_value_i) + ')';
		case pt_RANGE:	return "range (" + std::to_string(m_value_r[0]) + ':' + std::to_string(m_value_r[1]) + ')'; 
		case pt_EOLN:	return "end of line";

		case pt_KW_ALL:		return "ALL";
		case pt_KW_CHAIN:	return "CHAIN";
		case pt_KW_RESSEQ:	return "RESSEQ";
		case pt_KW_RESID:	return "RESID";
		case pt_KW_RESNAME:	return "RESNAME";
		case pt_KW_ELEMENT:	return "ELEMENT";
		case pt_KW_AND:		return "AND";
		case pt_KW_OR:		return "OR";
		case pt_KW_NOT:		return "NOT";
		case pt_KW_PDB:		return "PDB";
		case pt_KW_ENTRY:	return "ENTRY";
		case pt_KW_THROUGH:	return "THROUGH";
		default:
			assert(false);
			return "unknown token";
	}
}

TLSSelectionPtr TLSSelectionParserImplBusterOld::Parse()
{
	if (m_lookahead == pt_KW_PDB)
	{
		Match(pt_KW_PDB);
//		Match(pt_KW_ENTRY);
		
		throw std::runtime_error("Unimplemented PDB ENTRY specification");
	}

	TLSSelectionPtr result = ParseAtomSelection();

	Match(pt_EOLN);
	
	return result;
}

TLSSelectionPtr TLSSelectionParserImplBusterOld::ParseAtomSelection()
{
	TLSSelectionPtr result = ParseTerm();
	
	while (m_lookahead == pt_KW_OR)
	{
		Match(pt_KW_OR);
		result.reset(new TLSSelectionUnion(result, ParseTerm()));
	}
	
	return result;
}

TLSSelectionPtr TLSSelectionParserImplBusterOld::ParseTerm()
{
	TLSSelectionPtr result = ParseFactor();
	
	while (m_lookahead == pt_KW_AND)
	{
		Match(pt_KW_AND);
		result.reset(new TLSSelectionIntersection(result, ParseFactor()));
	}
	
	return result;
}

TLSSelectionPtr TLSSelectionParserImplBusterOld::ParseFactor()
{
	TLSSelectionPtr result;

	switch (m_lookahead)
	{
		case '(':
			Match('(');
			result = ParseAtomSelection();
			Match(')');
			break;
		
		case pt_KW_NOT:
			Match(pt_KW_NOT);
			result.reset(new TLSSelectionNot(ParseAtomSelection()));
			break;
		
		case pt_KW_CHAIN:
		{
			Match(pt_KW_CHAIN);

			std::string chainID = m_value_s;
			if (m_lookahead == pt_NUMBER)	// sigh
			{
				chainID = std::to_string(m_value_i);
				Match(pt_NUMBER);
			}
			else
				Match(m_lookahead == pt_STRING ? pt_STRING : pt_IDENT);
			
			result.reset(new TLSSelectionChain(chainID));
			break;
		}
		
		case pt_KW_RESNAME:
		{
			Match(pt_KW_RESNAME);
			std::string name = m_value_s;
			Match(pt_IDENT);
			result.reset(new TLSSelectionByName(name));
			break;
		}
		
		case pt_KW_RESSEQ:
			Match(pt_KW_RESSEQ);
			result = ParseResid();
			break;
	
		case pt_KW_RESID:
			Match(pt_KW_RESID);
			result = ParseResid();
			break;
		
		case pt_KW_ALL:
			Match(pt_KW_ALL);
			result.reset(new TLSSelectionAll());
			break;
		
		case pt_CHAINRESID:
			result = ParseChainResid();
			break;
		
		default:
			throw std::runtime_error("Unexpected token " + ToString(m_lookahead));
	}
	
	return result;
}

TLSSelectionPtr TLSSelectionParserImplBusterOld::ParseResid()
{
	TLSSelectionPtr result;
	
	for (;;)
	{
		int from, to;
				
		if (m_lookahead == pt_RANGE)
		{
			from = m_value_r[0];
			to = m_value_r[1];
			Match(pt_RANGE);
		}
		else
		{
			from = m_value_i;
			Match(pt_NUMBER);
			
			to = from;
			if (m_lookahead == ':' or m_lookahead == '-' or m_lookahead == pt_KW_THROUGH)
			{
				Match(m_lookahead);
				to = m_value_i;
				Match(pt_NUMBER);
			}
		}
		
		TLSSelectionPtr range(new TLSSelectionRangeSeq(from, to));
		
		if (result)
			result.reset(new TLSSelectionUnion(result, range));
		else
			result.reset(range.release());
		
		if (m_lookahead == ',')
		{
			Match(',');
			continue;
		}
		
		break;
	}
	
	return result;
}

TLSSelectionPtr TLSSelectionParserImplBusterOld::ParseChainResid()
{
	TLSSelectionPtr result;
	
	for (;;)
	{
		int from, to;
		
		from = to = m_value_i;
		std::string chainID = m_value_s;
		
		Match(pt_CHAINRESID);

		if (m_lookahead == '-')
		{
			Match(m_lookahead);
			to = m_value_i;
			
			if (m_value_s != chainID)
				throw std::runtime_error("Cannot have two different chainIDs in a range selection");
			
			Match(pt_CHAINRESID);
		}
		
		TLSSelectionPtr sc(new TLSSelectionChain(chainID));
		TLSSelectionPtr sr(new TLSSelectionRangeSeq(from, to));
		TLSSelectionPtr range(new TLSSelectionIntersection(sc, sr));
		
		if (result)
			result.reset(new TLSSelectionUnion(result, range));
		else
			result.reset(range.release());
		
		if (m_lookahead == ',')
		{
			Match(',');
			continue;
		}
		
		break;
	}
	
	return result;
}

// --------------------------------------------------------------------

class TLSSelectionParserBase
{
  public:
	virtual TLSSelectionPtr Parse(const std::string& selection) const = 0;
	virtual ~TLSSelectionParserBase() {}
};

template<typename IMPL>
class TLSSelectionParser
{
  public:
	virtual TLSSelectionPtr Parse(const std::string& selection) const
	{
		TLSSelectionPtr result;
		
		try
		{
			IMPL p(selection);
			result = p.Parse();
		}
		catch (const std::exception& ex)
		{
			if (cif::VERBOSE >= 0)
				std::cerr << "ParseError: " << ex.what() << std::endl;
		}
		
		return result;
	}
};

// --------------------------------------------------------------------


TLSSelectionPtr ParseSelectionDetails(const std::string& program, const std::string& selection)
{
	TLSSelectionParser<TLSSelectionParserImplPhenix> phenix;
	TLSSelectionParser<TLSSelectionParserImplBuster> buster;
	TLSSelectionParser<TLSSelectionParserImplBusterOld> busterOld;

	TLSSelectionPtr result;

	if (ba::icontains(program, "buster"))
	{
		result = buster.Parse(selection);

		if (not result)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Falling back to old BUSTER" << std::endl;
			result = busterOld.Parse(selection);
		}
		
		if (not result)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Falling back to PHENIX" << std::endl;
			result = phenix.Parse(selection);
		}
	}
	else if (ba::icontains(program, "phenix"))
	{
		result = phenix.Parse(selection);

		if (not result)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Falling back to BUSTER" << std::endl;
			result = buster.Parse(selection);
		}

		if (not result)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Falling back to old BUSTER" << std::endl;
			result = busterOld.Parse(selection);
		}
	}
	else
	{
		if (cif::VERBOSE > 0)
			std::cerr << "No known program specified, trying PHENIX" << std::endl;

		result = phenix.Parse(selection);

		if (not result)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Falling back to BUSTER" << std::endl;
			result = buster.Parse(selection);
		}

		if (not result)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Falling back to old BUSTER" << std::endl;
			result = busterOld.Parse(selection);
		}
	}
	
	return result;
}

}
