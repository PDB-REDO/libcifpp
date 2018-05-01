// Lib for working with structures as contained in mmCIF and PDB files

#pragma once

#include "cif++/ResolutionCalculator.h"

#include <clipper/clipper.h>

#include <cmath>

namespace mmcif
{

// --------------------------------------------------------------------

class ResolutionCalculator
{
  public:
	ResolutionCalculator(double a, double b, double c,
		double alpha, double beta, double gamma);
	ResolutionCalculator(const clipper::Cell& cell);
	
	double operator()(int h, int k, int l) const
	{
		double tmpres = h * h * mCoefs[0] + h * k * mCoefs[1] +
					    h * l * mCoefs[2] + k * k * mCoefs[3] +
					    k * l * mCoefs[4] + l * l * mCoefs[5];
		
		return 1.0 / std::sqrt(tmpres);
	}
	
  private:
	double	mCoefs[6];
};

}
