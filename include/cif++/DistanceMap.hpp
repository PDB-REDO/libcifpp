// copyright

#pragma once

#include <unordered_map>

#include <clipper/clipper.h>

#include "cif++/Structure.hpp"

namespace mmcif
{

class DistanceMap
{
  public:
	DistanceMap(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell,
		float maxDistance);

	// simplified version for subsets of atoms (used in refining e.g.)
//	DistanceMap(const Structure& p, const std::vector<Atom>& atoms);
	
	DistanceMap(const DistanceMap&) = delete;
	DistanceMap& operator=(const DistanceMap&) = delete;

	float operator()(const Atom& a, const Atom& b) const;

	std::vector<Atom> near(const Atom& a, float maxDistance = 3.5f) const;
	std::vector<Atom> near(const Point& p, float maxDistance = 3.5f) const;

	static clipper::Coord_orth
		CalculateOffsetForCell(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell);

	static std::vector<clipper::RTop_orth>
		AlternativeSites(const clipper::Spacegroup& spacegroup, const clipper::Cell& cell);

  private:

	typedef std::map<std::tuple<size_t,size_t>,std::tuple<float,int32_t>> DistMap;

	void AddDistancesForAtoms(const Residue& a, const Residue& b, DistMap& dm, int32_t rtix);

	const Structure&						structure;
	size_t									dim;
	std::unordered_map<std::string,size_t>	index;
	std::map<size_t,std::string>			rIndex;
	
	float									mMaxDistance, mMaxDistanceSQ;
	
	std::vector<std::tuple<float,int32_t>>	mA;
	std::vector<size_t>						mIA, mJA;
	Point									mD;			// needed to move atoms to center
	std::vector<clipper::RTop_orth>			mRtOrth;
};

}
