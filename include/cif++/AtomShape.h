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
	AtomShape(const Atom& atom, float resHigh, float resLow,
		float bFactor);

	~AtomShape();
	
	AtomShape(const AtomShape&) = delete;
	AtomShape& operator=(const AtomShape&) = delete;

	float radius() const;
	float calculatedDensity(float r) const;

  private:
	struct AtomShapeImpl*	mImpl;
};
	
}
