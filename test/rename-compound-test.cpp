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

	try
	{
		std::filesystem::path testdir = std::filesystem::current_path();

		if (argc == 3)
			testdir = argv[2];

		if (std::filesystem::exists(testdir / ".." / "data" / "ccd-subset.cif"))
			cif::addFileResource("components.cif", testdir / ".." / "data" / "ccd-subset.cif");

		if (std::filesystem::exists(testdir / ".." / "rsrc" / "mmcif_pdbx_v50.dic"))
			cif::addFileResource("mmcif_pdbx_v50.dic", testdir / ".." / "rsrc" / "mmcif_pdbx_v50.dic");

		mmcif::CompoundFactory::instance().pushDictionary(testdir / "REA.cif");
		mmcif::CompoundFactory::instance().pushDictionary(testdir / "RXA.cif");

		mmcif::File f(testdir / ".."/"examples"/"1cbs.cif.gz");
		mmcif::Structure structure(f);

		auto &res = structure.getResidue("B", "REA");
		structure.changeResidue(res, "RXA", {});

		structure.cleanupEmptyCategories();

		f.file().save(std::cout);
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		exit(1);
	}
	
	return 0;	
}
