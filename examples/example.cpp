#include <iostream>
#include <filesystem>

#include <cif++/Cif++.hpp>

namespace fs = std::filesystem;

int main()
{
	fs::path in("1cbs.cif.gz");

	cif::File file;

	file.loadDictionary("mmcif_pdbx_v50");

	file.load("1cbs.cif.gz");

	auto& db = file.firstDatablock()["atom_site"];
	auto n = db.find(cif::Key("label_atom_id") == "OXT").size();

	std::cout << "File contains " << db.size() << " atoms of which " << n << (n == 1 ? " is" : " are") << " OXT" << std::endl
		<< "residues with an OXT are:" << std::endl;
	
	for (const auto& [asym, comp, seqnr]: db.find<std::string,std::string,int>(
			cif::Key("label_atom_id") == "OXT",
			{ "label_asym_id", "label_comp_id", "label_seq_id" }
		))
	{
		std::cout << asym << ' ' << comp << ' ' << seqnr << std::endl;
	}

	return 0;
}
