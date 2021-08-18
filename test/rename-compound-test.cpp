#include "../include/cif++/Cif++.hpp"
#include "../include/cif++/PDB2Cif.hpp"
#include "../include/cif++/Structure.hpp"

#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
	cif::VERBOSE = 3;
	
	mmcif::CompoundFactory::instance().pushDictionary("RXA.cif");

	mmcif::File f("../examples/1cbs.cif.gz");
	mmcif::Structure structure(f);

	auto &res = structure.getResidue("B", "REA");
	structure.changeResidue(res, "RXA", {});

	structure.cleanupEmptyCategories();

	f.file().save(std::cout);
	
	return 0;	
}
