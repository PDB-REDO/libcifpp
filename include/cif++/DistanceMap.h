// copyright

#pragma once

#include <unordered_map>

//#include <boost/hash/hash.hpp>

#include <clipper/clipper.h>

#include "cif++/Structure.h"

namespace mmcif
{

class DistanceMap
{
  public:
	DistanceMap(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell);

	// simplified version for subsets of atoms (used in refining e.g.)	
//	DistanceMap(const Structure& p, const std::vector<Atom>& atoms);
	
	DistanceMap(const DistanceMap&) = delete;
	DistanceMap& operator=(const DistanceMap&) = delete;

	float operator()(const Atom& a, const Atom& b) const;

	std::vector<Atom> near(const Atom& a, float maxDistance = 3.5f) const;
	std::vector<Atom> near(const Point& p, float maxDistance = 3.5f) const;

  private:

	const Structure& structure;
	size_t dim;
	std::unordered_map<std::string,size_t> index;
	std::map<size_t,std::string> rIndex;
	
	struct key_hash
	{
		size_t operator()(std::tuple<size_t,size_t> const& k) const noexcept
		{
			size_t i[2];
			std::tie(i[0], i[1]) = k;
			return boost::hash_range(i, i + 1);
		}
	};
	
	std::vector<std::tuple<float,size_t>>	mA;
	std::vector<size_t>						mIA, mJA;
	Point									mD;			// needed to move atoms to center
	std::vector<clipper::RTop_orth>			mRtOrth;
};

}
