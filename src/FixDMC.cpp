// --------------------------------------------------------------------

#include <set>

#include <cif++/Structure.hpp>

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
			
			auto atomC = mon.atomByID("C");
			auto atomCA = mon.atomByID("CA");
			auto atomN = mon.atomByID("N");
			auto atomO = mon.atomByID("O");

			int missing = (atomC ? 0 : 1) + (atomCA ? 0 : 1) + (atomN ? 0 : 1) + (atomO ? 0 : 1);

			switch (missing)
			{
				case 1:
					if (not atomO)
						addO(mon);
					else if (not atomN)
						addN(mon);
					else if (not atomCA)
						addCA(mon);
					else if (not atomC)
						addC(mon);
					break;
			}
		}
	}
}

}