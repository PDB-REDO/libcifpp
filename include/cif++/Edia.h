// copyright

#pragma once

#include <zeep/xml/document.hpp>
#include <clipper/clipper.h>

#include "cif++/Structure.h"
#include "cif++/DistanceMap.h"
#include "cif++/BondMap.h"

namespace libcif
{

// --------------------------------------------------------------------
// Code to calculate EDIA (electron density for individual atoms)

// --------------------------------------------------------------------
// AtomRadius can be used to retrieve the precalculated atom radius
// for an atom with a certain charge at a certain resolution.

class AtomRadius
{
  public:
	AtomRadius();

	float operator()(AtomType a, int charge, float resolution);
	float operator()(Atom a, float resolution)
	{
		return operator()(a.type(), a.charge(), resolution);
	}

  private:

	// key consists of atom_type, charge and resolution

	typedef std::tuple<AtomType,int,float>	Key;
	typedef std::map<Key,float>				Cache;

	zeep::xml::document mCurves;
	Cache mCache;
};

// --------------------------------------------------------------------

float CalculateEDIA(const Atom& atom, const clipper::Xmap<float>& xmap,
	float resolution, float meanDensity, float rmsDensity,
	const DistanceMap& dm, const BondMap& bm);

}
