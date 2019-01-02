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
		uint32 ixa = index.at(a.id());
		uint32 ixb = index.at(b.id());
	
		return bond_1_4.count(key(ixa, ixb));
	}
	
	// links coming from the struct_conn records:
	std::vector<std::string> linked(const Atom& a) const;
	
  private:

	bool isBonded(uint32 ai, uint32 bi) const
	{
		return bond.count(key(ai, bi)) != 0;
	}

	uint64 key(uint32 a, uint32 b) const
	{
		if (a > b)
			std::swap(a, b);
		return static_cast<uint64>(a) | (static_cast<uint64>(b) << 32);
	}
	
	std::tuple<uint32,uint32> dekey(uint64 k) const
	{
		return std::make_tuple(
			static_cast<uint32>(k >> 32),
			static_cast<uint32>(k)
		);
	}
	
	uint32 dim;
	std::unordered_map<std::string,uint32> index;
	std::set<uint64> bond, bond_1_4;

	std::map<std::string,std::set<std::string>> link;
};

}
