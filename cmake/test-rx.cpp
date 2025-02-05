// See: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86164

#include <iostream>
#include <regex>

int main()
{
	std::string s(100'000, '*');
	std::smatch m;
	std::regex r("^(.*?)$");

	std::regex_search(s, m, r);

	std::cout << s.substr(0, 10) << '\n';
	std::cout << m.str(1).substr(0, 10) << '\n';

	return 0;
}
