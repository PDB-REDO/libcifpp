#if __has_include("../src/Config.hpp")
#include "../src/Config.hpp"
#endif
#include "../include/cif++/Cif++.hpp"
#include "../include/cif++/PDB2Cif.hpp"

#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>

// #include "pdb2cif.h"

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
	using namespace std::literals;

	po::options_description desc("pdb2cif-test "s + PACKAGE_VERSION + " options");
	desc.add_options()
		("input,i",		po::value<std::string>(),	"Input file")
		("help,h",									"Display help message")
		("version",									"Print version")
		("verbose,v",								"Verbose output")
		("debug,d",		po::value<int>(),			"Debug level (for even more verbose output)");

	po::positional_options_description p;
	p.add("input", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		std::cout << argv[0] << " version " PACKAGE_VERSION << std::endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0)
	{
		std::cerr << desc << std::endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();
	
	std::ifstream is(vm["input"].as<std::string>());
	if (not is.is_open())
		throw std::runtime_error("Could not open file " + vm["input"].as<std::string>());
	
	cif::File f;
	ReadPDBFile(is, f);
	f.save(std::cout);
	
	return 0;	
}
