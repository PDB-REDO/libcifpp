// copyright

#pragma once

#include <unordered_map>

#include <clipper/clipper.h>

#include "cif++/Structure.h"

namespace mmcif
{

class DistanceMap
{
  public:
	DistanceMap(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell);

	// simplified version for subsets of atoms (used in refining e.g.)	
	DistanceMap(const std::vector<Atom>& atoms);
	
	DistanceMap(const DistanceMap&) = delete;
	DistanceMap& operator=(const DistanceMap&) = delete;

	float operator()(const Atom& a, const Atom& b) const;
	std::vector<Atom> near(const Atom& a, float maxDistance = 3.5f) const;

  private:

	size_t dim;
	std::vector<float> dist;
	std::unordered_map<std::string,size_t> index;
};

}
