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

	bool is1_4(const Atom& a, const Atom& b) const
	{
		uint32_t ixa = index.at(a.id());
		uint32_t ixb = index.at(b.id());
	
		return bond_1_4.count(key(ixa, ixb));
	}
	
	// links coming from the struct_conn records:
	std::vector<std::string> linked(const Atom& a) const;
	
  private:

	bool isBonded(uint32_t ai, uint32_t bi) const
	{
		return bond.count(key(ai, bi)) != 0;
	}

	uint64_t key(uint32_t a, uint32_t b) const
	{
		if (a > b)
			std::swap(a, b);
		return static_cast<uint64_t>(a) | (static_cast<uint64_t>(b) << 32);
	}
	
	std::tuple<uint32_t,uint32_t> dekey(uint64_t k) const
	{
		return std::make_tuple(
			static_cast<uint32_t>(k >> 32),
			static_cast<uint32_t>(k)
		);
	}
	
	uint32_t dim;
	std::unordered_map<std::string,uint32_t> index;
	std::set<uint64_t> bond, bond_1_4;

	std::map<std::string,std::set<std::string>> link;
};

}
