#define CATCH_CONFIG_RUNNER 1

#include "test-main.hpp"

#include <cif++.hpp>

std::filesystem::path gTestDir = std::filesystem::current_path();

int main(int argc, char *argv[])
{
	Catch::Session session; // There must be exactly one instance

	// Build a new parser on top of Catch2's
#if CATCH22
	using namespace Catch::clara;
#else
	// Build a new parser on top of Catch2's
	using namespace Catch::Clara;
#endif

	auto cli = session.cli()                               // Get Catch2's command line parser
	           | Opt(gTestDir, "data-dir")                 // bind variable to a new option, with a hint string
	                 ["-D"]["--data-dir"]                  // the option names it will respond to
	           ("The directory containing the data files") // description string for the help output
	           | Opt(cif::VERBOSE, "verbose")["-v"]["--cif-verbose"]("Flag for cif::VERBOSE");

	// Now pass the new composite back to Catch2 so it uses that
	session.cli(cli);

	// Let Catch2 (using Clara) parse the command line
	int returnCode = session.applyCommandLine(argc, argv);
	if (returnCode != 0) // Indicates a command line error
		return returnCode;

	// do this now, avoids the need for installing
	cif::add_file_resource("mmcif_pdbx.dic", gTestDir / ".." / "rsrc" / "mmcif_pdbx.dic");

	// initialize CCD location
	cif::add_file_resource("components.cif", gTestDir / ".." / "rsrc" / "ccd-subset.cif");

	cif::compound_factory::instance().push_dictionary(gTestDir / "HEM.cif");

	return session.run();
}