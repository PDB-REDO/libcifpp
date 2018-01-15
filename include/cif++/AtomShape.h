// AtomShape, analogue to the similarly named code in clipper

#pragma once

#include "cif++/Structure.h"

namespace libcif
{

// --------------------------------------------------------------------
// Class used in calculating radii

class AtomShape
{
  public:
	AtomShape(const Atom& atom, float resHigh, float resLow);

	float radius() const;

  private:
	AtomType	mSymbol;
	int			mCharge;
	float		mUIso, mOccupancy;
	float		mResHigh, mResLow;
};
	
}
