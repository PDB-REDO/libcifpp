/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "cif++/ResolutionCalculator.hpp"

#include "cif++/Point.hpp"

using namespace std;

namespace mmcif
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
