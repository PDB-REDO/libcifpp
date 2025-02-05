#include <filesystem>
#include <iostream>

#include <cif++.hpp>

namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		std::cerr << "Usage: example <inputfile>\n";
		exit(1);
	}

	cif::file file(argv[1]);

	if (file.empty())
	{
		std::cerr << "Empty file\n";
		exit(1);
	}

	auto &db = file.front();
	auto &atom_site = db["atom_site"];
	auto n = atom_site.find(cif::key("label_atom_id") == "OXT").size();

	std::cout << "File contains " << atom_site.size() << " atoms of which " << n << (n == 1 ? " is" : " are") << " OXT\n"
			  << "residues with an OXT are:\n";

	for (const auto &[asym, comp, seqnr] : atom_site.find<std::string, std::string, int>(
			 cif::key("label_atom_id") == "OXT", "label_asym_id", "label_comp_id", "label_seq_id"))
	{
		std::cout << asym << ' ' << comp << ' ' << seqnr << '\n';
	}

	return 0;
}
