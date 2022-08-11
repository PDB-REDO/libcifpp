#include <iostream>
#include <filesystem>

#include <cif++/cif.hpp>

namespace fs = std::filesystem;

int main()
{
	cif::file file;
	file.load_dictionary("mmcif_pdbx_v50");
	file.load("1cbs.cif.gz");

	auto& db = file.front();
	auto &atom_site = db["atom_site"];
	auto n = atom_site.find(cif::key("label_atom_id") == "OXT").size();

	std::cout << "File contains " << atom_site.size() << " atoms of which " << n << (n == 1 ? " is" : " are") << " OXT" << std::endl
		<< "residues with an OXT are:" << std::endl;
	
	for (const auto& [asym, comp, seqnr]: atom_site.find<std::string,std::string,int>(
			cif::key("label_atom_id") == "OXT", "label_asym_id", "label_comp_id", "label_seq_id"))
	{
		std::cout << asym << ' ' << comp << ' ' << seqnr << std::endl;
	}

	return 0;
}
