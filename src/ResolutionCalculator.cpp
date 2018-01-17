// --------------------------------------------------------------------

#include "cif++/ResolutionCalculator.h"

#include "cif++/Point.h"

using namespace std;

namespace libcif
{

ResolutionCalculator::ResolutionCalculator(const clipper::Cell& cell)
	: ResolutionCalculator(cell.a(), cell.b(), cell.c(),
		180 * cell.alpha() / kPI,
		180 * cell.beta() / kPI,
		180 * cell.gamma() / kPI)
{
}

ResolutionCalculator::ResolutionCalculator(double a, double b, double c,
		double alpha, double beta, double gamma)
{
	double deg2rad = atan(1.0) / 45.0;
	
	double ca = cos(deg2rad * alpha);
	double sa = sin(deg2rad * alpha);
	double cb = cos(deg2rad * beta);
	double sb = sin(deg2rad * beta);
	double cg = cos(deg2rad * gamma);
	double sg = sin(deg2rad * gamma);
	
	double cast = (cb * cg - ca) / (sb * sg);
	double cbst = (cg * ca - cb) / (sg * sa);
	double cgst = (ca * cb - cg) / (sa * sb);
	
	double sast = sqrt(1 - cast * cast);
	double sbst = sqrt(1 - cbst * cbst);
	double sgst = sqrt(1 - cgst * cgst);

	double ast = 1 / (a * sb * sgst);
	double bst = 1 / (b * sg * sast);
	double cst = 1 / (c * sa * sbst);
	
	mCoefs[0] = ast * ast;
	mCoefs[1] = 2 * ast * bst * cgst;
	mCoefs[2] = 2 * ast * cst * cbst;
	mCoefs[3] = bst * bst;
	mCoefs[4] = 2 * bst * cst * cast;
	mCoefs[5] = cst * cst;
}

}
