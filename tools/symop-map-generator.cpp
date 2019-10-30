/*
cd ~/projects/pdb-redo/libcif++/tools/
clang++ -I ~/my-clipper/include -L ~/my-clipper/lib -o symop-map-generator symop-map-generator.cpp -lclipper-core
./symop-map-generator

*/

#include <cassert>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <regex>

#include <cstdlib>

#include <clipper/clipper.h>

using namespace std;

std::regex kNameRx(R"(^(\d+) +(\d+) +(\d+) +(\S+) +(\S+) +(\S+) +'([^']+)'( +'([^']+)')?(?: +!.+)?$)");

class SymopParser
{
  public:
	SymopParser() {}

	array<int,15> parse(const string& s)
	{
		m_p = s.begin();
		m_e = s.end();
		m_lookahead = next_token();

		parsepart(0);
		match((Token)',');
		parsepart(1);
		match((Token)',');
		parsepart(2);

		if (m_lookahead != 0 or m_p != m_e)
			throw runtime_error("symmetry expression contains more data than expected");

		return {
			m_rot[0][0], m_rot[0][1], m_rot[0][2], 
			m_rot[1][0], m_rot[1][1], m_rot[1][2], 
			m_rot[2][0], m_rot[2][1], m_rot[2][2], 
			m_trn[0][0], m_trn[0][1],
			m_trn[1][0], m_trn[1][1],
			m_trn[2][0], m_trn[2][1]
		};
	}

  private:

	enum Token : int { Eof = 0, Number = 256, XYZ };

	string to_string(Token t)
	{
		switch (t)
		{
			case Eof:		return "end of expression";
			case Number:	return "number";
			case XYZ:		return "'x', 'y' or 'z'";
			default:
				if (isprint(t))
					return string({'\'', static_cast<char>(t), '\''});
				return "invalid character " + std::to_string(static_cast<int>(t));
		}
	}
	
	Token next_token()
	{
		Token result = Eof;
		while (m_p != m_e)
		{
			char ch = *m_p++;
			if (ch == ' ')
				continue;

			switch (ch)
			{
				case 'x':
				case 'X':
					result = XYZ;
					m_nr = 0;
					break;
				
				case 'y':
				case 'Y':
					result = XYZ;
					m_nr = 1;
					break;
				
				case 'z':
				case 'Z':
					result = XYZ;
					m_nr = 2;
					break;
				
				default:
					if (isdigit(ch))
					{
						m_nr = ch - '0';
						result = Number;
					}
					else
						result = (Token)ch;
					break;
			}
			break;
		}

		return result;
	}

	void match(Token token)
	{
		if (m_lookahead != token)
			throw runtime_error("Unexpected character " + to_string(m_lookahead) + " expected " + to_string(token));
		
		m_lookahead = next_token();
	}

	void parsepart(int row)
	{
		bool first = false;

		do
		{
			int sign = m_lookahead == '-' ? -1 : 1;
			if (m_lookahead == '-' or m_lookahead == '+')
				match(m_lookahead);

			if (m_lookahead == Number)
			{
				m_trn[row][0] = sign * m_nr;
				match(Number);

				match((Token)'/');

				m_trn[row][1] = m_nr;
				match(Number);
			}
			else
			{
				m_rot[row][m_nr] = sign;
				match(XYZ);
			}
		}
		while (m_lookahead == '+' or m_lookahead == '-');
	}

	Token m_lookahead;
	char m_nr;

	string m_s;
	string::const_iterator m_p, m_e;

	int m_rot[3][3] = {};
	int m_trn[3][2] = {};
};

int main()
{
	try
	{
		const char* CLIBD = getenv("CLIBD");
		
		if (CLIBD == nullptr)
			throw runtime_error("CCP4 not sourced");

		ifstream file(CLIBD + "/symop.lib"s);
		if (not file.is_open())
			throw runtime_error("Could not open symop.lib file");
		
		string line;
		clipper::Spacegroup spacegroup;

		enum class State { Initial, InSpacegroup, Error } state = State::Initial;
		int symopnr = 0;

		cout << R"(
// This file was generated from $CLIBD/symop.lib

struct SymopNrTable
{
	uint8_t spacegroupNr;
	uint8_t	rotationalNr:7;
	int8_t rot[3][3]:2;
	uint8_t trn[3][2]:4;
};

kSymopNrTable[] = {
)" << endl;

		while (getline(file, line))
		{
			if (line.empty())
				throw runtime_error("Invalid symop.lib file, contains empty line");
			
			switch (state)
			{
				case State::Error:
				case State::InSpacegroup:
					if (line[0] == ' ')
					{
						if (state == State::Error)
							continue;
						
						try
						{
							SymopParser p;
							auto symop = p.parse(line);

							cout << "    { " << setw(3) << spacegroup.spacegroup_number() << '\t' << symopnr++;
							for (auto i: symop)
								cout << ',' << setw(2) << i;
							cout << " }," << endl;
						}
						catch (const exception& e)
						{
							cerr << line << endl
								 << e.what() << endl;
						}


						continue;
					}
					// fall through
				
				case State::Initial:
				{
					smatch m;
					if (not regex_match(line, m, kNameRx))
					{
						cerr << line << endl;
						throw runtime_error("Name line does not match regular expression");
					}

					int nr = stoi(m[1].str());
					if (nr > 230)
					{
						state = State::Error;
						continue;
					}

					try
					{
						spacegroup = clipper::Spacegroup(clipper::Spgr_descr(stoi(m[1].str())));
						symopnr = 1;
						
						cout << "    // " << spacegroup.symbol_hm() << endl;

						if (spacegroup.num_symops() != stoi(m[2]))
						{
							cerr << line << endl
								 << "Num of symops do not match: " << spacegroup.num_symops() << " <> " << m[2] << endl;
							state = State::Error;
							continue;
						}

						state = State::InSpacegroup;
					}
					catch (const clipper::Message_fatal& e)
					{
						cerr << line << endl
							 << e.text() << endl;
						state = State::Error;
					}
					break;
				}
			}
		}
		cout << "};" << endl;
	}
	catch (const exception& ex)
	{
		cerr << endl
			 << "Program terminated due to error:" << endl
			 << ex.what() << endl;
	}
	
	
	return 0;
}
