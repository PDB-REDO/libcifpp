/*
    Created by: Maarten L. Hekkelman
    Date: dinsdag 07 november, 2017

    Copyright 2017 NKI AVL

    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to
    permit persons to whom the Software is furnished to do so, subject to
    the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

// #include <sys/ioctl.h>
// #include <termios.h>

#include <iomanip>
#include <iostream>

#include <cif++.hpp>
#include <cif++/pdb/tls.hpp>

namespace cif
{

const int
	kResidueNrWildcard = std::numeric_limits<int>::min(),
	kNoSeqNum = std::numeric_limits<int>::max() - 1;

// --------------------------------------------------------------------
// We parse selection statements and create a selection expression tree
// which is then interpreted by setting the selected flag for the
// residues. After that, the selected ranges are collected and printed.

struct tls_residue
{
	std::string chainID;
	int seqNr = 0;
	char iCode;
	std::string name;
	bool selected;

	std::string asymID;
	int seqID = 0;

	bool operator==(const tls_residue &rhs) const
	{
		return chainID == rhs.chainID and
		       seqNr == rhs.seqNr and
		       iCode == rhs.iCode and
		       iequals(name, rhs.name) and
		       selected == rhs.selected;
	}
};

void dump_selection(const std::vector<tls_residue> &selected, int indentLevel)
{
	std::string indent(indentLevel * 2, ' ');

	auto i = selected.begin();
	bool first = true;

	// First print in PDB space
	while (i != selected.end())
	{
		auto b = find_if(i, selected.end(), [](auto s) -> bool
			{ return s.selected; });
		if (b == selected.end())
			break;

		if (first)
			std::cout << indent << "PDB:" << std::endl;
		first = false;

		auto e = find_if(b, selected.end(), [b](auto s) -> bool
			{ return s.chainID != b->chainID or not s.selected; });

		std::cout << indent << " >> " << b->chainID << ' ' << b->seqNr << ':' << (e - 1)->seqNr << std::endl;
		i = e;
	}

	// Then in mmCIF space

	if (not first)
		std::cout << indent << "mmCIF:" << std::endl;

	i = selected.begin();
	while (i != selected.end())
	{
		auto b = find_if(i, selected.end(), [](auto s) -> bool
			{ return s.selected; });
		if (b == selected.end())
			break;

		auto e = find_if(b, selected.end(), [b](auto s) -> bool
			{ return s.asymID != b->asymID or not s.selected; });

		std::string asymID = b->asymID;
		int from = b->seqID, to = from;

		for (auto j = b + 1; j != e; ++j)
		{
			if (j->seqID == to + 1)
				to = j->seqID;
			else if (j->seqID != to) // probably an insertion code
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
			std::cout << indent << "Empty selection" << std::endl;
	}
}

std::vector<std::tuple<std::string, int, int>> tls_selection::get_ranges(cif::datablock &db, bool pdbNamespace) const
{
	std::vector<tls_residue> selected;

	// Collect the residues from poly seq scheme...
	for (auto r : db["pdbx_poly_seq_scheme"])
	{
		std::string chain, seqNr, iCode, name;

		std::string asymID;
		int seqID = 0;

		if (pdbNamespace)
			cif::tie(chain, seqNr, iCode, name, asymID, seqID) = r.get("pdb_strand_id", "pdb_seq_num", "pdb_ins_code", "pdb_mon_id", "asym_id", "seq_id");
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

		selected.push_back({ chain, stoi(seqNr), iCode[0], name, false, asymID, seqID });
	}

	// ... those from the nonpoly scheme
	for (auto r : db["pdbx_nonpoly_scheme"])
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

		selected.push_back({ chain, stoi(seqNr), iCode[0], name, false, asymID, kNoSeqNum });
	}

	// ... those from the nonpoly scheme
	for (auto r : db["pdbx_branch_scheme"])
	{
		std::string chain, seqNr, iCode, name, asymID;

		if (pdbNamespace)
		{
			cif::tie(chain, seqNr, iCode, name, asymID) = r.get("auth_asym_id", "pdb_seq_num", "pdb_ins_code", "pdb_mon_id", "asym_id");
			if (seqNr.empty())
				continue;
		}
		else
		{
			cif::tie(chain, name) = r.get("asym_id", "mon_id");
			asymID = chain;
			seqNr = "0";
		}

		if (iCode.length() > 1)
			throw std::runtime_error("invalid iCode");

		selected.push_back({ chain, stoi(seqNr), iCode[0], name, false, asymID, kNoSeqNum });
	}

	// selected might consist of multiple ranges
	// output per chain

	stable_sort(selected.begin(), selected.end(), [](auto &a, auto &b) -> bool
		{
			int d = a.chainID.compare(b.chainID);
			if (d == 0)
				d = a.seqNr - b.seqNr;
			return d < 0; });

	collect_residues(db, selected);

	std::vector<std::tuple<std::string, int, int>> result;

	if (pdbNamespace)
	{
		auto i = selected.begin();

		while (i != selected.end())
		{
			auto b = find_if(i, selected.end(), [](auto s) -> bool
				{ return s.selected; });
			if (b == selected.end())
				break;

			auto e = find_if(b, selected.end(), [b](auto s) -> bool
				{ return s.chainID != b->chainID or not s.selected; });

			// return ranges with strict increasing sequence numbers.
			// So when there's a gap in the sequence we split the range.
			// Beware of iCodes though
			result.push_back(std::make_tuple(b->chainID, b->seqNr, b->seqNr));
			for (auto j = b + 1; j != e; ++j)
			{
				if (j->seqNr == std::get<2>(result.back()) + 1)
					std::get<2>(result.back()) = j->seqNr;
				else if (j->seqNr != std::get<2>(result.back())) // probably an insertion code
					result.push_back(std::make_tuple(b->chainID, j->seqNr, j->seqNr));
			}

			i = e;
		}
	}
	else
	{
		auto i = selected.begin();

		while (i != selected.end())
		{
			auto b = find_if(i, selected.end(), [](auto s) -> bool
				{ return s.selected; });
			if (b == selected.end())
				break;

			auto e = find_if(b, selected.end(), [b](auto s) -> bool
				{ return s.asymID != b->asymID or not s.selected; });

			// return ranges with strict increasing sequence numbers.
			// So when there's a gap in the sequence we split the range.
			// Beware of iCodes though
			result.push_back(std::make_tuple(b->asymID, b->seqID, b->seqID));
			for (auto j = b + 1; j != e; ++j)
			{
				if (j->seqID == std::get<2>(result.back()) + 1)
					std::get<2>(result.back()) = j->seqID;
				else if (j->seqID != std::get<2>(result.back())) // probably an insertion code
					result.push_back(std::make_tuple(b->asymID, j->seqID, j->seqID));
			}

			i = e;
		}
	}

	for (auto &&[name, i1, i2] : result)
	{
		if (i1 == kNoSeqNum) i1 = 0;
		if (i2 == kNoSeqNum) i2 = 0;
	}

	return result;
}

struct tls_selection_not : public tls_selection
{
	tls_selection_not(std::unique_ptr<tls_selection> selection)
		: selection(selection.release())
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, size_t indentLevel) const override
	{
		selection->collect_residues(db, residues, indentLevel + 1);

		for (auto &r : residues)
			r.selected = not r.selected;

		if (cif::VERBOSE)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "NOT" << std::endl;
			dump_selection(residues, indentLevel);
		}
	}

	std::unique_ptr<tls_selection> selection;
};

struct tls_selection_all : public tls_selection
{
	tls_selection_all() {}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, size_t indentLevel) const override
	{
		for (auto &r : residues)
			r.selected = true;

		if (cif::VERBOSE)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "ALL" << std::endl;
			dump_selection(residues, indentLevel);
		}
	}
};

struct tls_selection_chain : public tls_selection_all
{
	tls_selection_chain(const std::string &chainID)
		: m_chain(chainID)
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, size_t indentLevel) const override
	{
		bool allChains = m_chain == "*";

		for (auto &r : residues)
			r.selected = allChains or r.chainID == m_chain;

		if (cif::VERBOSE)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "CHAIN " << m_chain << std::endl;
			dump_selection(residues, indentLevel);
		}
	}

	std::string m_chain;
};

struct tls_selection_res_id : public tls_selection_all
{
	tls_selection_res_id(int seqNr, char iCode)
		: m_seq_nr(seqNr)
		, m_icode(iCode)
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, size_t indentLevel) const override
	{
		for (auto &r : residues)
			r.selected = r.seqNr == m_seq_nr and r.iCode == m_icode;

		if (cif::VERBOSE)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "ResID " << m_seq_nr << (m_icode ? std::string{ m_icode } : "") << std::endl;
			dump_selection(residues, indentLevel);
		}
	}

	int m_seq_nr;
	char m_icode;
};

struct tls_selection_range_seq : public tls_selection_all
{
	tls_selection_range_seq(int first, int last)
		: m_first(first)
		, m_last(last)
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, size_t indentLevel) const override
	{
		for (auto &r : residues)
		{
			r.selected = ((r.seqNr >= m_first or m_first == kResidueNrWildcard) and
						  (r.seqNr <= m_last or m_last == kResidueNrWildcard));
		}

		if (cif::VERBOSE)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Range " << m_first << ':' << m_last << std::endl;
			dump_selection(residues, indentLevel);
		}
	}

	int m_first, m_last;
};

struct tls_selection_range_id : public tls_selection_all
{
	tls_selection_range_id(int first, int last, char icodeFirst = 0, char icodeLast = 0)
		: m_first(first)
		, m_last(last)
		, m_icode_first(icodeFirst)
		, m_icode_last(icodeLast)
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, size_t indentLevel) const override
	{
		// need to do this per chain
		std::set<std::string> chains;
		for (auto &r : residues)
			chains.insert(r.chainID);

		for (std::string chain : chains)
		{
			auto f = find_if(residues.begin(), residues.end(),
				[this,chain](auto r) -> bool
				{
					return r.chainID == chain and r.seqNr == m_first and r.iCode == m_icode_first;
				});

			auto l = find_if(residues.begin(), residues.end(),
				[this,chain](auto r) -> bool
				{
					return r.chainID == chain and r.seqNr == m_last and r.iCode == m_icode_last;
				});

			if (f != residues.end() and l != residues.end() and f <= l)
			{
				++l;

				for (; f != l; ++f)
					f->selected = true;
			}
		}

		if (cif::VERBOSE)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Through " << m_first << ':' << m_last << std::endl;
			dump_selection(residues, indentLevel);
		}
	}

	int m_first, m_last;
	char m_icode_first, m_icode_last;
};

struct tls_selection_union : public tls_selection
{
	tls_selection_union(std::unique_ptr<tls_selection> &lhs, std::unique_ptr<tls_selection> &rhs)
		: lhs(lhs.release())
		, rhs(rhs.release())
	{
	}

	tls_selection_union(std::unique_ptr<tls_selection> &lhs, std::unique_ptr<tls_selection> &&rhs)
		: lhs(lhs.release())
		, rhs(rhs.release())
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, size_t indentLevel) const override
	{
		auto a = residues;
		for_each(a.begin(), a.end(), [](auto &r)
			{ r.selected = false; });

		auto b = residues;
		for_each(b.begin(), b.end(), [](auto &r)
			{ r.selected = false; });

		lhs->collect_residues(db, a, indentLevel + 1);
		rhs->collect_residues(db, b, indentLevel + 1);

		for (auto ai = a.begin(), bi = b.begin(), ri = residues.begin(); ri != residues.end(); ++ai, ++bi, ++ri)
			ri->selected = ai->selected or bi->selected;

		if (cif::VERBOSE)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Union" << std::endl;
			dump_selection(residues, indentLevel);
		}
	}

	std::unique_ptr<tls_selection> lhs;
	std::unique_ptr<tls_selection> rhs;
};

struct tls_selection_intersection : public tls_selection
{
	tls_selection_intersection(std::unique_ptr<tls_selection> &lhs, std::unique_ptr<tls_selection> &rhs)
		: lhs(lhs.release())
		, rhs(rhs.release())
	{
	}

	tls_selection_intersection(std::unique_ptr<tls_selection> &lhs, std::unique_ptr<tls_selection> &&rhs)
		: lhs(lhs.release())
		, rhs(rhs.release())
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, size_t indentLevel) const override
	{
		auto a = residues;
		for_each(a.begin(), a.end(), [](auto &r)
			{ r.selected = false; });

		auto b = residues;
		for_each(b.begin(), b.end(), [](auto &r)
			{ r.selected = false; });

		lhs->collect_residues(db, a, indentLevel + 1);
		rhs->collect_residues(db, b, indentLevel + 1);

		for (auto ai = a.begin(), bi = b.begin(), ri = residues.begin(); ri != residues.end(); ++ai, ++bi, ++ri)
			ri->selected = ai->selected and bi->selected;

		if (cif::VERBOSE)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Intersection" << std::endl;
			dump_selection(residues, indentLevel);
		}
	}

	std::unique_ptr<tls_selection> lhs;
	std::unique_ptr<tls_selection> rhs;
};

struct tls_selection_by_name : public tls_selection_all
{
  public:
	tls_selection_by_name(const std::string &resname)
		: m_name(resname)
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, size_t indentLevel) const override
	{
		for (auto &r : residues)
			r.selected = r.name == m_name;

		if (cif::VERBOSE)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Name " << m_name << std::endl;
			dump_selection(residues, indentLevel);
		}
	}

	std::string m_name;
};

struct tls_selection_by_element : public tls_selection_all
{
  public:
	tls_selection_by_element(const std::string &element)
		: m_element(element)
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, size_t indentLevel) const override
	{
		// rationale... We want to select residues only. So we select
		// residues that have just a single atom of type m_element.
		// And we assume these have as residue name... m_element.
		// ... Right?

		for (auto &r : residues)
			r.selected = iequals(r.name, m_element);

		if (cif::VERBOSE)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Element " << m_element << std::endl;
			dump_selection(residues, indentLevel);
		}
	}

	std::string m_element;
};

// --------------------------------------------------------------------

class tls_selection_parser_impl
{
  public:
	tls_selection_parser_impl(const std::string &selection)
		: m_selection(selection)
		, m_p(m_selection.begin())
		, m_end(m_selection.end())
	{
	}

	virtual std::unique_ptr<tls_selection> Parse() = 0;

  protected:
	virtual int get_next_token() = 0;
	virtual void match(int token);
	virtual std::string to_string(int token) = 0;

	std::string m_selection;
	std::string::iterator m_p, m_end;
	int m_lookahead;
	std::string m_token;
};

void tls_selection_parser_impl::match(int token)
{
	if (m_lookahead == token)
		m_lookahead = get_next_token();
	else
	{
		std::string expected;
		if (token >= 256)
			expected = to_string(token);
		else
			expected = { char(token) };

		std::string found;
		if (m_lookahead >= 256)
			found = to_string(m_lookahead) + " (" + m_token + ')';
		else
			found = { char(m_lookahead) };

		throw std::runtime_error("Expected " + expected + " but found " + found);
	}
}

// --------------------------------------------------------------------

class TLSSelectionParserImplPhenix : public tls_selection_parser_impl
{
  public:
	TLSSelectionParserImplPhenix(const std::string &selection)
		: tls_selection_parser_impl(selection)
	{
		m_lookahead = get_next_token();
	}

	virtual std::unique_ptr<tls_selection> Parse();

  private:
	std::unique_ptr<tls_selection> ParseAtomSelection();
	std::unique_ptr<tls_selection> ParseTerm();
	std::unique_ptr<tls_selection> ParseFactor();

	enum TOKEN
	{
		pt_NONE = 0,
		pt_IDENT = 256,
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

	virtual int get_next_token();
	virtual std::string to_string(int token);

	int m_value_i;
	std::string m_value_s;
	char m_icode;
};

int TLSSelectionParserImplPhenix::get_next_token()
{
	int result = pt_NONE;
	enum STATE
	{
		st_START,
		st_RESID = 200,
		st_NUM = 300,
		st_IDENT = 400,
		st_QUOTED = 500,
		st_DQUOTED = 550,
		st_OTHER = 600
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
			case st_START: state = start = st_RESID; break;
			case st_RESID: state = start = st_NUM; break;
			case st_NUM: state = start = st_IDENT; break;
			case st_IDENT: state = start = st_QUOTED; break;
			case st_QUOTED: state = start = st_DQUOTED; break;
			case st_DQUOTED: state = start = st_OTHER; break;
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

std::string TLSSelectionParserImplPhenix::to_string(int token)
{
	switch (token)
	{
		case pt_IDENT: return "identifier";
		case pt_STRING: return "std::string";
		case pt_NUMBER: return "number";
		case pt_RESID: return "resid";
		case pt_EOLN: return "end of line";

		case pt_KW_ALL: return "ALL";
		case pt_KW_CHAIN: return "CHAIN";
		case pt_KW_RESSEQ: return "RESSEQ";
		case pt_KW_RESID: return "RESID";
		case pt_KW_RESNAME: return "RESNAME";
		case pt_KW_ELEMENT: return "ELEMENT";
		case pt_KW_AND: return "AND";
		case pt_KW_OR: return "OR";
		case pt_KW_NOT: return "NOT";
		case pt_KW_PDB: return "PDB";
		case pt_KW_ENTRY: return "ENTRY";
		case pt_KW_THROUGH: return "THROUGH";

		default: return "character";
	}
}

std::unique_ptr<tls_selection> TLSSelectionParserImplPhenix::Parse()
{
	if (m_lookahead == pt_KW_PDB)
	{
		match(pt_KW_PDB);
		//		Match(pt_KW_ENTRY);

		throw std::runtime_error("Unimplemented PDB ENTRY specification");
	}

	std::unique_ptr<tls_selection> result = ParseAtomSelection();

	bool extraParenthesis = false;

	if (m_lookahead == ')')
	{
		extraParenthesis = true;
		m_lookahead = get_next_token();
	}

	match(pt_EOLN);

	if (extraParenthesis)
		std::cerr << "WARNING: too many closing parenthesis in TLS selection statement" << std::endl;

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplPhenix::ParseAtomSelection()
{
	std::unique_ptr<tls_selection> result = ParseTerm();

	while (m_lookahead == pt_KW_OR)
	{
		match(pt_KW_OR);
		result.reset(new tls_selection_union(result, ParseTerm()));
	}

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplPhenix::ParseTerm()
{
	std::unique_ptr<tls_selection> result = ParseFactor();

	while (m_lookahead == pt_KW_AND)
	{
		match(pt_KW_AND);
		result.reset(new tls_selection_intersection(result, ParseFactor()));
	}

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplPhenix::ParseFactor()
{
	std::unique_ptr<tls_selection> result;

	switch (m_lookahead)
	{
		case '(':
			match('(');
			result = ParseAtomSelection();
			if (m_lookahead == pt_EOLN)
				std::cerr << "WARNING: missing closing parenthesis in TLS selection statement" << std::endl;
			else
				match(')');
			break;

		case pt_KW_NOT:
			match(pt_KW_NOT);
			result.reset(new tls_selection_not(ParseAtomSelection()));
			break;

		case pt_KW_CHAIN:
		{
			match(pt_KW_CHAIN);

			std::string chainID = m_value_s;
			if (m_lookahead == pt_NUMBER) // sigh
			{
				chainID = to_string(m_value_i);
				match(pt_NUMBER);
			}
			else
				match(m_lookahead == pt_STRING ? pt_STRING : pt_IDENT);

			result.reset(new tls_selection_chain(chainID));
			break;
		}

		case pt_KW_RESNAME:
		{
			match(pt_KW_RESNAME);
			std::string name = m_value_s;
			match(pt_IDENT);
			result.reset(new tls_selection_by_name(name));
			break;
		}

		case pt_KW_ELEMENT:
		{
			match(pt_KW_ELEMENT);
			std::string element = m_value_s;
			match(pt_IDENT);
			result.reset(new tls_selection_by_element(element));
			break;
		}

		case pt_KW_RESSEQ:
		{
			match(pt_KW_RESSEQ);

			int from = m_value_i;
			match(pt_NUMBER);

			int to = from;
			if (m_lookahead == ':')
			{
				match(':');
				to = m_value_i;
				match(pt_NUMBER);
			}

			result.reset(new tls_selection_range_seq(from, to));
			break;
		}

		case pt_KW_RESID:
		{
			match(pt_KW_RESID);

			int from, to;
			char icode_from = 0, icode_to = 0;
			bool through = false;

			from = to = m_value_i;

			if (m_lookahead == pt_NUMBER)
				match(pt_NUMBER);
			else
			{
				icode_from = m_icode;
				match(pt_RESID);
			}

			if (m_lookahead == ':' or m_lookahead == pt_KW_THROUGH or m_lookahead == '-')
			{
				through = m_lookahead == pt_KW_THROUGH;

				match(m_lookahead);

				to = m_value_i;
				if (m_lookahead == pt_NUMBER)
					match(pt_NUMBER);
				else
				{
					icode_to = m_icode;
					match(pt_RESID);
				}

				if (through)
					result.reset(new tls_selection_range_id(from, to, icode_from, icode_to));
				else
				{
					if (cif::VERBOSE and (icode_from or icode_to))
						std::cerr << "Warning, ignoring insertion codes" << std::endl;

					result.reset(new tls_selection_range_seq(from, to));
				}
			}
			else
				result.reset(new tls_selection_res_id(from, icode_from));

			break;
		}

		case pt_KW_ALL:
			match(pt_KW_ALL);
			result.reset(new tls_selection_all());
			break;

		default:
			throw std::runtime_error("Unexpected token " + to_string(m_lookahead) + " (" + m_token + ')');
	}

	return result;
}

// --------------------------------------------------------------------

class TLSSelectionParserImplBuster : public tls_selection_parser_impl
{
  public:
	TLSSelectionParserImplBuster(const std::string &selection);

	virtual std::unique_ptr<tls_selection> Parse();

  protected:
	enum TOKEN
	{
		bt_NONE = 0,
		bt_IDENT = 256,
		bt_NUMBER,
		bt_EOLN,
	};

	virtual int get_next_token();
	virtual std::string to_string(int token);

	std::unique_ptr<tls_selection> ParseGroup();
	std::tuple<std::string, int> ParseAtom();

	std::unique_ptr<tls_selection> ParseOldGroup();

	int m_value_i;
	std::string m_value_s;
	bool m_parsing_old_style = false;
};

TLSSelectionParserImplBuster::TLSSelectionParserImplBuster(const std::string &selection)
	: tls_selection_parser_impl(selection)
{
	m_lookahead = get_next_token();
}

int TLSSelectionParserImplBuster::get_next_token()
{
	int result = bt_NONE;
	enum STATE
	{
		st_START,
		st_NEGATE,
		st_NUM,
		st_IDENT
	} state = st_START;

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

std::string TLSSelectionParserImplBuster::to_string(int token)
{
	switch (token)
	{
		case bt_IDENT: return "identifier (" + m_value_s + ')';
		case bt_NUMBER: return "number (" + to_string(m_value_i) + ')';
		case bt_EOLN: return "end of line";

		default:
			assert(false);
			return "unknown token";
	}
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBuster::ParseGroup()
{
	std::unique_ptr<tls_selection> result;

	auto add = [&result](const std::string &chainID, int from, int to)
	{
		std::unique_ptr<tls_selection> sc(new tls_selection_chain(chainID));
		std::unique_ptr<tls_selection> sr(new tls_selection_range_seq(from, to));
		std::unique_ptr<tls_selection> s(new tls_selection_intersection(sc, sr));

		if (result == nullptr)
			result.reset(s.release());
		else
			result.reset(new tls_selection_union{ result, s });
	};

	match('{');

	do
	{
		std::string chain1;
		int seqNr1;
		std::tie(chain1, seqNr1) = ParseAtom();

		if (m_lookahead == '-')
		{
			std::string chain2;
			int seqNr2 = seqNr1;

			match('-');

			if (m_lookahead == bt_NUMBER)
			{
				seqNr2 = m_value_i;
				match(bt_NUMBER);
			}
			else
			{
				std::tie(chain2, seqNr2) = ParseAtom();
				if (chain1 != chain2)
				{
					std::cerr << "Warning, ranges over multiple chains detected" << std::endl;

					std::unique_ptr<tls_selection> sc1(new tls_selection_chain(chain1));
					std::unique_ptr<tls_selection> sr1(new tls_selection_range_seq(seqNr1, kResidueNrWildcard));
					std::unique_ptr<tls_selection> s1(new tls_selection_intersection(sc1, sr1));

					std::unique_ptr<tls_selection> sc2(new tls_selection_chain(chain2));
					std::unique_ptr<tls_selection> sr2(new tls_selection_range_seq(kResidueNrWildcard, seqNr2));
					std::unique_ptr<tls_selection> s2(new tls_selection_intersection(sc2, sr2));

					std::unique_ptr<tls_selection> s(new tls_selection_union(s1, s2));

					if (result == nullptr)
						result.reset(s.release());
					else
						result.reset(new tls_selection_union{ result, s });

					chain1.clear();
				}
			}

			if (not chain1.empty())
				add(chain1, seqNr1, seqNr2);
		}
		else
			add(chain1, seqNr1, seqNr1);
	} while (m_lookahead != '}');

	match('}');

	return result;
}

std::tuple<std::string, int> TLSSelectionParserImplBuster::ParseAtom()
{
	std::string chain = m_value_s;
	int seqNr = kResidueNrWildcard;

	if (m_lookahead == '*')
		match('*');
	else
		match(bt_IDENT);

	match('|');

	if (m_lookahead == '*')
		match('*');
	else
	{
		seqNr = m_value_i;
		match(bt_NUMBER);

		if (m_lookahead == ':')
		{
			match(':');
			std::string atom = m_value_s;

			if (cif::VERBOSE)
				std::cerr << "Warning: ignoring atom ID '" << atom << "' in TLS selection" << std::endl;

			match(bt_IDENT);
		}
	}

	return std::make_tuple(chain, seqNr);
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBuster::Parse()
{
	std::unique_ptr<tls_selection> result = ParseGroup();
	match(bt_EOLN);
	return result;
}

// --------------------------------------------------------------------

class TLSSelectionParserImplBusterOld : public tls_selection_parser_impl
{
  public:
	TLSSelectionParserImplBusterOld(const std::string &selection)
		: tls_selection_parser_impl(selection)
	{
		m_lookahead = get_next_token();
	}

	virtual std::unique_ptr<tls_selection> Parse();

  private:
	std::unique_ptr<tls_selection> ParseAtomSelection();
	std::unique_ptr<tls_selection> ParseTerm();
	std::unique_ptr<tls_selection> ParseFactor();

	std::unique_ptr<tls_selection> ParseResid();
	std::unique_ptr<tls_selection> ParseChainResid();

	enum TOKEN
	{
		pt_NONE = 0,
		pt_IDENT = 256,
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

	virtual int get_next_token();
	virtual std::string to_string(int token);

	int m_value_i;
	std::string m_value_s;
	int m_value_r[2];
};

int TLSSelectionParserImplBusterOld::get_next_token()
{
	int result = pt_NONE;
	enum STATE
	{
		st_START,
		st_NEGATE,
		st_NUM,
		st_RANGE,
		st_IDENT_1,
		st_IDENT,
		st_CHAINRESID,
		st_QUOTED_1,
		st_QUOTED_2
	} state = st_START;

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

			case st_RANGE: // TODO: question, is "-2--1" a valid range? We do not support that, yet
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
					m_value_s += to_string(m_value_i);
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

std::string TLSSelectionParserImplBusterOld::to_string(int token)
{
	switch (token)
	{
		case pt_IDENT: return "identifier (" + m_value_s + ')';
		case pt_STRING: return "std::string (" + m_value_s + ')';
		case pt_NUMBER: return "number (" + to_string(m_value_i) + ')';
		case pt_RANGE: return "range (" + to_string(m_value_r[0]) + ':' + to_string(m_value_r[1]) + ')';
		case pt_EOLN: return "end of line";

		case pt_KW_ALL: return "ALL";
		case pt_KW_CHAIN: return "CHAIN";
		case pt_KW_RESSEQ: return "RESSEQ";
		case pt_KW_RESID: return "RESID";
		case pt_KW_RESNAME: return "RESNAME";
		case pt_KW_ELEMENT: return "ELEMENT";
		case pt_KW_AND: return "AND";
		case pt_KW_OR: return "OR";
		case pt_KW_NOT: return "NOT";
		case pt_KW_PDB: return "PDB";
		case pt_KW_ENTRY: return "ENTRY";
		case pt_KW_THROUGH: return "THROUGH";
		default:
			assert(false);
			return "unknown token";
	}
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBusterOld::Parse()
{
	if (m_lookahead == pt_KW_PDB)
	{
		match(pt_KW_PDB);
		//		Match(pt_KW_ENTRY);

		throw std::runtime_error("Unimplemented PDB ENTRY specification");
	}

	std::unique_ptr<tls_selection> result = ParseAtomSelection();

	match(pt_EOLN);

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBusterOld::ParseAtomSelection()
{
	std::unique_ptr<tls_selection> result = ParseTerm();

	while (m_lookahead == pt_KW_OR)
	{
		match(pt_KW_OR);
		result.reset(new tls_selection_union(result, ParseTerm()));
	}

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBusterOld::ParseTerm()
{
	std::unique_ptr<tls_selection> result = ParseFactor();

	while (m_lookahead == pt_KW_AND)
	{
		match(pt_KW_AND);
		result.reset(new tls_selection_intersection(result, ParseFactor()));
	}

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBusterOld::ParseFactor()
{
	std::unique_ptr<tls_selection> result;

	switch (m_lookahead)
	{
		case '(':
			match('(');
			result = ParseAtomSelection();
			match(')');
			break;

		case pt_KW_NOT:
			match(pt_KW_NOT);
			result.reset(new tls_selection_not(ParseAtomSelection()));
			break;

		case pt_KW_CHAIN:
		{
			match(pt_KW_CHAIN);

			std::string chainID = m_value_s;
			if (m_lookahead == pt_NUMBER) // sigh
			{
				chainID = to_string(m_value_i);
				match(pt_NUMBER);
			}
			else
				match(m_lookahead == pt_STRING ? pt_STRING : pt_IDENT);

			result.reset(new tls_selection_chain(chainID));
			break;
		}

		case pt_KW_RESNAME:
		{
			match(pt_KW_RESNAME);
			std::string name = m_value_s;
			match(pt_IDENT);
			result.reset(new tls_selection_by_name(name));
			break;
		}

		case pt_KW_RESSEQ:
			match(pt_KW_RESSEQ);
			result = ParseResid();
			break;

		case pt_KW_RESID:
			match(pt_KW_RESID);
			result = ParseResid();
			break;

		case pt_KW_ALL:
			match(pt_KW_ALL);
			result.reset(new tls_selection_all());
			break;

		case pt_CHAINRESID:
			result = ParseChainResid();
			break;

		default:
			throw std::runtime_error("Unexpected token " + to_string(m_lookahead));
	}

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBusterOld::ParseResid()
{
	std::unique_ptr<tls_selection> result;

	for (;;)
	{
		int from, to;

		if (m_lookahead == pt_RANGE)
		{
			from = m_value_r[0];
			to = m_value_r[1];
			match(pt_RANGE);
		}
		else
		{
			from = m_value_i;
			match(pt_NUMBER);

			to = from;
			if (m_lookahead == ':' or m_lookahead == '-' or m_lookahead == pt_KW_THROUGH)
			{
				match(m_lookahead);
				to = m_value_i;
				match(pt_NUMBER);
			}
		}

		std::unique_ptr<tls_selection> range(new tls_selection_range_seq(from, to));

		if (result)
			result.reset(new tls_selection_union(result, range));
		else
			result.reset(range.release());

		if (m_lookahead == ',')
		{
			match(',');
			continue;
		}

		break;
	}

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBusterOld::ParseChainResid()
{
	std::unique_ptr<tls_selection> result;

	for (;;)
	{
		int from, to;

		from = to = m_value_i;
		std::string chainID = m_value_s;

		match(pt_CHAINRESID);

		if (m_lookahead == '-')
		{
			match(m_lookahead);
			to = m_value_i;

			if (m_value_s != chainID)
				throw std::runtime_error("Cannot have two different chainIDs in a range selection");

			match(pt_CHAINRESID);
		}

		std::unique_ptr<tls_selection> sc(new tls_selection_chain(chainID));
		std::unique_ptr<tls_selection> sr(new tls_selection_range_seq(from, to));
		std::unique_ptr<tls_selection> range(new tls_selection_intersection(sc, sr));

		if (result)
			result.reset(new tls_selection_union(result, range));
		else
			result.reset(range.release());

		if (m_lookahead == ',')
		{
			match(',');
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
	virtual std::unique_ptr<tls_selection> Parse(const std::string &selection) const = 0;
	virtual ~TLSSelectionParserBase() {}
};

template <typename IMPL>
class TLSSelectionParser
{
  public:
	virtual std::unique_ptr<tls_selection> Parse(const std::string &selection) const
	{
		std::unique_ptr<tls_selection> result;

		try
		{
			IMPL p(selection);
			result = p.Parse();
		}
		catch (const std::exception &ex)
		{
			std::cerr << "ParseError: " << ex.what() << std::endl;
		}

		return result;
	}
};

// --------------------------------------------------------------------

std::unique_ptr<tls_selection> parse_tls_selection_details(const std::string &program, const std::string &selection)
{
	TLSSelectionParser<TLSSelectionParserImplPhenix> phenix;
	TLSSelectionParser<TLSSelectionParserImplBuster> buster;
	TLSSelectionParser<TLSSelectionParserImplBusterOld> busterOld;

	std::unique_ptr<tls_selection> result;

	if (cif::icontains(program, "buster"))
	{
		result = buster.Parse(selection);

		if (not result)
		{
			if (cif::VERBOSE)
				std::cerr << "Falling back to old BUSTER" << std::endl;
			result = busterOld.Parse(selection);
		}

		if (not result)
		{
			if (cif::VERBOSE)
				std::cerr << "Falling back to PHENIX" << std::endl;
			result = phenix.Parse(selection);
		}
	}
	else if (cif::icontains(program, "phenix"))
	{
		result = phenix.Parse(selection);

		if (not result)
		{
			if (cif::VERBOSE)
				std::cerr << "Falling back to BUSTER" << std::endl;
			result = buster.Parse(selection);
		}

		if (not result)
		{
			if (cif::VERBOSE)
				std::cerr << "Falling back to old BUSTER" << std::endl;
			result = busterOld.Parse(selection);
		}
	}
	else
	{
		if (cif::VERBOSE)
			std::cerr << "No known program specified, trying PHENIX" << std::endl;

		result = phenix.Parse(selection);

		if (not result)
		{
			if (cif::VERBOSE)
				std::cerr << "Falling back to BUSTER" << std::endl;
			result = buster.Parse(selection);
		}

		if (not result)
		{
			if (cif::VERBOSE)
				std::cerr << "Falling back to old BUSTER" << std::endl;
			result = busterOld.Parse(selection);
		}
	}

	return result;
}

} // namespace cif
