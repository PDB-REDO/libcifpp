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

	bool operator()(const Atom& a, const Atom& b) const;
	
  private:

	uint32 dim;
	std::vector<bool> bond;
	std::unordered_map<std::string,size_t> index;
};

}
