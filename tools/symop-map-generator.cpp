/*
cd ~/projects/pdb-redo/libcif++/tools/
clang++ -I ~/my-clipper/include -L ~/my-clipper/lib -o symop-map-generator symop-map-generator.cpp -lclipper-core
./symop-map-generator

*/

#include <cassert>

#include <iostream>
#include <fstream>
#include <regex>

#include <cstdlib>

#include <clipper/clipper.h>

using namespace std;

std::regex kNameRx(R"(^(\d+) +(\d+) +(\d+) +(\S+) +(\S+) +(\S+) +'([^']+)'( +'([^']+)')?(?: +!.+)?$)");

class SymopParser
{
  public:
	SymopParser(const std::string& s)
		: m_s(s)
		, m_p(s.begin())
		, m_e(s.end())
	{
		m_lookahead = next_char();
	}

	void parse()
	{
		parsepart(0);
		match(',');
		parsepart(1);
		match(',');
		parsepart(2);
		if (m_lookahead != 0)
			throw runtime_error("symmetry expression contains more data than expected");
	}

  private:

	char next_char()
	{
		char result = 0;
		while (m_p != m_e)
		{
			char ch = *m_p++;
			if (ch == ' ')
				continue;

			result = ch;
			break;
		}

		return result;
	}

	void retract()
	{
		--m_p;
	}

	void match(char c)
	{
		if (m_lookahead != c)
			throw runtime_error("Unexpected character " + m_lookahead + " expected " + c);
		
		m_lookahead = next_char();
	}

	void parsepart(int xyz)
	{
		if (isdigit(m_lookahead))
		{
			parserational(xyz);

			switch (m_lookahead)
			{
				case '-':
					match('-');
					break;
			}

		}


		}
		if (m_lookahead == '+')
		{
			match('+');
		}
		else if (m_lookahead = '-')
		{

		}


	}

	tuple<int,int> parserational();

	char m_lookahead;

	string m_s;
	string::iterator m_p, m_e;

	int m_rot[3][3] = {};
	int m_trn[3] = {};
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
		
	}
	catch (const exception& ex)
	{
		cerr << endl
			 << "Program terminated due to error:" << endl
			 << ex.what() << endl;
	}
	
	
	return 0;
}
