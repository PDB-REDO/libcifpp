// copyright

#pragma once

#include <unordered_map>

#include "cif++/Structure.h"

namespace mmcif
{

class BondMap
{
  public:
	BondMap(const Structure& p);
	
	BondMap(const BondMap&) = delete;
	BondMap& operator=(const BondMap&) = delete;

	bool operator()(const Atom& a, const Atom& b) const
	{
		return isBonded(index.at(a.id()), index.at(b.id()));
	}

	bool is1_4(const Atom& a, const Atom& b) const;
	
  private:

	bool isBonded(size_t ai, size_t bi) const;

	size_t dim;
	std::vector<bool> bond;
	std::unordered_map<std::string,size_t> index;
};

}
