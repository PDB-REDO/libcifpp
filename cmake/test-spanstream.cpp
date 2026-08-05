#include <iostream>
#include <spanstream>

int main()
{
	std::ispanstream is("hello, world!\n");
	std::cout << is.rdbuf();
	return 0;
}