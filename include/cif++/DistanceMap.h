// copyright

#pragma once

#include <unordered_map>

#include <clipper/clipper.h>

#include "cif++/Structure.h"

namespace libcif
{

class DistanceMap
{
  public:
	DistanceMap(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell);
	
	DistanceMap(const DistanceMap&) = delete;
	DistanceMap& operator=(const DistanceMap&) = delete;

	float operator()(const Atom& a, const Atom& b) const;
	std::vector<Atom> near(const Atom& a, float maxDistance = 3.5f) const;

  private:

	uint32 dim;
	std::vector<float> dist;
	std::unordered_map<std::string,size_t> index;
};

}
