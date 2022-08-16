#include "../include/cif++/Cif++.hpp"
#include "../include/cif++/PDB2Cif.hpp"
#include "../include/cif++/Structure.hpp"

#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{
	cif::VERBOSE = 3;

	try
	{
		std::filesystem::path testdir = std::filesystem::current_path();

		if (argc == 3)
			testdir = argv[2];

		if (std::filesystem::exists(testdir / ".." / "data" / "ccd-subset.cif"))
			cif::add_file_resource("components.cif", testdir / ".." / "data" / "ccd-subset.cif");

		if (std::filesystem::exists(testdir / ".." / "rsrc" / "mmcif_pdbx.dic"))
			cif::add_file_resource("mmcif_pdbx.dic", testdir / ".." / "rsrc" / "mmcif_pdbx.dic");

		mmcif::CompoundFactory::instance().pushDictionary(testdir / "REA.cif");
		mmcif::CompoundFactory::instance().pushDictionary(testdir / "RXA.cif");

		mmcif::File f(testdir / ".."/"examples"/"1cbs.cif.gz");
		mmcif::Structure structure(f);

		auto &res = structure.getResidue("B");
		structure.changeResidue(res, "RXA", {});

		structure.cleanupEmptyCategories();

		f.save(std::cout);
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		exit(1);
	}
	
	return 0;	
}
