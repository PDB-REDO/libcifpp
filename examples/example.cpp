#include <filesystem>
#include <iostream>

#include <cif++.hpp>

namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
	if (argc != 2)
		exit(1);

	cif::file file = cif::pdb::read(argv[1]);

	if (file.empty())
	{
		std::cerr << "Empty file" << std::endl;
		exit(1);
	}

	auto &db = file.front();
	auto &atom_site = db["atom_site"];
	auto n = atom_site.find(cif::key("label_atom_id") == "OXT").size();

	std::cout << "File contains " << atom_site.size() << " atoms of which " << n << (n == 1 ? " is" : " are") << " OXT" << std::endl
			  << "residues with an OXT are:" << std::endl;

	for (const auto &[asym, comp, seqnr] : atom_site.find<std::string, std::string, int>(
			 cif::key("label_atom_id") == "OXT", "label_asym_id", "label_comp_id", "label_seq_id"))
	{
		std::cout << asym << ' ' << comp << ' ' << seqnr << std::endl;
	}

	return 0;
}
