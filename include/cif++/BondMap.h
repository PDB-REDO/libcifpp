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
		return get(index.at(a.id()), index.at(b.id())) == 1;
	}

	bool operator()(const std::string& id_a, const std::string& id_b) const
	{
		return get(index.at(id_a), index.at(id_b)) == 1;
	}

	bool is1_4(const Atom& a, const Atom& b) const
	{
		return get(index.at(a.id()), index.at(b.id())) == 3;
	}
	
	bool is1_4(const std::string& id_a, const std::string& id_b) const
	{
		return get(index.at(id_a), index.at(id_b)) == 3;
	}

//	bool is1_4(const Atom& a, const Atom& b) const;
//	bool is1_4(const std::string& id_a, const std::string& id_b) const;
	
  private:

	uint8 get(size_t ia, size_t ib) const
	{
		uint8 result = 0;
		if (ia != ib)
		{
			if (ib < ia)
				std::swap(ia, ib);
			
			size_t ix = ib + ia * dim - ia * (ia + 1) / 2;
			assert(ix < bond.size());
			result = bond[ix];
		}
		return result;
	}

	size_t dim;
	std::vector<uint8> bond;
	std::unordered_map<std::string,size_t> index;
};

}
