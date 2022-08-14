#include "../include/cif++/cif.hpp"
#include "../include/cif++/pdb/PDB2Cif.hpp""

#include <iostream>
#include <fstream>

// #include "pdb2cif.h"

int main(int argc, char* argv[])
{
	using namespace std::literals;

	if (argc != 2)
	{
		std::cerr << "Usage: pdb2cif-test <input-file>" << std::endl;
		exit(1);
	}
	
	std::ifstream is(argv[1]);
	if (not is.is_open())
		throw std::runtime_error("Could not open file "s + argv[1]);
	
	cif::File f;
	ReadPDBFile(is, f);
	f.save(std::cout);
	
	return 0;	
}
