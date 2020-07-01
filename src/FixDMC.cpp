// --------------------------------------------------------------------

#include <set>

#include <cif++/Structure.h>

namespace mmcif
{

void addC(Monomer& mon)
{

}

void addCA(Monomer& mon)
{

}

void addN(Monomer& mon)
{

}

void addO(Monomer& mon)
{

}

void CreateMissingBackboneAtoms(Structure& structure, bool simplified)
{
	for (auto& poly: structure.polymers())
	{
		for (auto& mon: poly)
		{
			if (mon.isComplete() or mon.hasAlternateBackboneAtoms())
				continue;
			
			std::set<std::string> missing;
			if (not mon.hasAtomWithID("C"))		missing.insert("C");
			if (not mon.hasAtomWithID("CA"))	missing.insert("CA");
			if (not mon.hasAtomWithID("N"))		missing.insert("N");
			if (not mon.hasAtomWithID("O"))		missing.insert("O");

			switch (missing.size())
			{
				case 1:
					if (missing.count("O"))
						addO(mon);
					else if (missing.count("N"))
						addN(mon);
					else if (missing.count("CA"))
						addCA(mon);
					else if (missing.count("C"))
						addC(mon);
					break;
			}
		}
	}
}

}