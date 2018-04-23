// Lib for working with structures as contained in mmCIF and PDB files

#include <random>

#include "cif++/Point.h"

using namespace std;

namespace libcif
{

// --------------------------------------------------------------------

quaternion Normalize(quaternion q)
{
	valarray<double> t(4);
	
	t[0] = q.R_component_1();
	t[1] = q.R_component_2();
	t[2] = q.R_component_3();
	t[3] = q.R_component_4();
	
	t *= t;
	
	double length = sqrt(t.sum());

	if (length > 0.001)
		q /= length;
	else
		q = quaternion(1, 0, 0, 0);

	return q;
}

// --------------------------------------------------------------------

float DihedralAngle(const Point& p1, const Point& p2, const Point& p3, const Point& p4)
{
	Point v12 = p1 - p2;	// vector from p2 to p1
	Point v43 = p4 - p3;	// vector from p3 to p4
	
	Point z = p2 - p3;		// vector from p3 to p2
	
	Point p = CrossProduct(z, v12);
	Point x = CrossProduct(z, v43);
	Point y = CrossProduct(z, x);
	
	double u = DotProduct(x, x);
	double v = DotProduct(y, y);
	
	double result = 360;
	if (u > 0 and v > 0)
	{
		u = DotProduct(p, x) / sqrt(u);
		v = DotProduct(p, y) / sqrt(v);
		if (u != 0 or v != 0)
			result = atan2(v, u) * 180 / kPI;
	}
	
	return result;
}

float CosinusAngle(const Point& p1, const Point& p2, const Point& p3, const Point& p4)
{
	Point v12 = p1 - p2;
	Point v34 = p3 - p4;
	
	double result = 0;
	
	double x = DotProduct(v12, v12) * DotProduct(v34, v34);
	if (x > 0)
		result = DotProduct(v12, v34) / sqrt(x);
	
	return result;
}

// --------------------------------------------------------------------

tuple<double,Point> QuaternionToAngleAxis(quaternion q)
{
	if (q.R_component_1() > 1)
		q = Normalize(q);

	// angle:
	double angle = 2 * acos(q.R_component_1());
	angle = angle * 180 / kPI;

	// axis:
	double s = sqrt(1 - q.R_component_1() * q.R_component_1());
	if (s < 0.001)
		s = 1;
	
	Point axis(q.R_component_2() / s, q.R_component_3() / s, q.R_component_4() / s);

	return make_tuple(angle, axis);
}

Point CenterPoints(vector<Point>& Points)
{
	Point t;
	
	for (Point& pt : Points)
	{
		t.getX() += pt.getX();
		t.getY() += pt.getY();
		t.getZ() += pt.getZ();
	}
	
	t.getX() /= Points.size();
	t.getY() /= Points.size();
	t.getZ() /= Points.size();
	
	for (Point& pt : Points)
	{
		pt.getX() -= t.getX();
		pt.getY() -= t.getY();
		pt.getZ() -= t.getZ();
	}
	
	return t;
}

Point Centroid(vector<Point>& Points)
{
	Point result;
	
	for (Point& pt : Points)
		result += pt;
	
	result /= Points.size();
	
	return result;
}

double RMSd(const vector<Point>& a, const vector<Point>& b)
{
	double sum = 0;
	for (uint32 i = 0; i < a.size(); ++i)
	{
		valarray<double> d(3);
		
		d[0] = b[i].getX() - a[i].getX();
		d[1] = b[i].getY() - a[i].getY();
		d[2] = b[i].getZ() - a[i].getZ();

		d *= d;
		
		sum += d.sum();
	}
	
	return sqrt(sum / a.size());
}

// The next function returns the largest solution for a quartic equation
// based on Ferrari's algorithm.
// A depressed quartic is of the form:
//
//   x^4 + ax^2 + bx + c = 0
//
// (since I'm too lazy to find out a better way, I've implemented the
//  routine using complex values to avoid nan's as a result of taking
//  sqrt of a negative number)
double LargestDepressedQuarticSolution(double a, double b, double c)
{
	complex<double> P = - (a * a) / 12 - c;
	complex<double> Q = - (a * a * a) / 108 + (a * c) / 3 - (b * b) / 8;
	complex<double> R = - Q / 2.0 + sqrt((Q * Q) / 4.0 + (P * P * P) / 27.0);
	
	complex<double> U = pow(R, 1 / 3.0);
	
	complex<double> y;
	if (U == 0.0)
		y = -5.0 * a / 6.0 + U - pow(Q, 1.0 / 3.0);
	else
		y = -5.0 * a / 6.0 + U - P / (3.0 * U);

	complex<double> W = sqrt(a + 2.0 * y);
	
	// And to get the final result:
	// result = (±W + sqrt(-(3 * alpha + 2 * y ± 2 * beta / W))) / 2;
	// We want the largest result, so:

	valarray<double> t(4);

	t[0] = (( W + sqrt(-(3.0 * a + 2.0 * y + 2.0 * b / W))) / 2.0).real();
	t[1] = (( W + sqrt(-(3.0 * a + 2.0 * y - 2.0 * b / W))) / 2.0).real();
	t[2] = ((-W + sqrt(-(3.0 * a + 2.0 * y + 2.0 * b / W))) / 2.0).real();
	t[3] = ((-W + sqrt(-(3.0 * a + 2.0 * y - 2.0 * b / W))) / 2.0).real();

	return t.max();
}

//quaternion AlignPoints(const vector<Point>& pa, const vector<Point>& pb)
//{
//	// First calculate M, a 3x3 matrix containing the sums of products of the coordinates of A and B
//	matrix<double> M(3, 3, 0);
//
//	for (uint32 i = 0; i < pa.size(); ++i)
//	{
//		const Point& a = pa[i];
//		const Point& b = pb[i];
//		
//		M(0, 0) += a.getX() * b.getX();	M(0, 1) += a.getX() * b.getY();	M(0, 2) += a.getX() * b.getZ();
//		M(1, 0) += a.getY() * b.getX();	M(1, 1) += a.getY() * b.getY();	M(1, 2) += a.getY() * b.getZ();
//		M(2, 0) += a.getZ() * b.getX();	M(2, 1) += a.getZ() * b.getY();	M(2, 2) += a.getZ() * b.getZ();
//	}
//	
//	// Now calculate N, a symmetric 4x4 matrix
//	symmetric_matrix<double> N(4);
//	
//	N(0, 0) =  M(0, 0) + M(1, 1) + M(2, 2);
//	N(0, 1) =  M(1, 2) - M(2, 1);
//	N(0, 2) =  M(2, 0) - M(0, 2);
//	N(0, 3) =  M(0, 1) - M(1, 0);
//	
//	N(1, 1) =  M(0, 0) - M(1, 1) - M(2, 2);
//	N(1, 2) =  M(0, 1) + M(1, 0);
//	N(1, 3) =  M(0, 2) + M(2, 0);
//	
//	N(2, 2) = -M(0, 0) + M(1, 1) - M(2, 2);
//	N(2, 3) =  M(1, 2) + M(2, 1);
//	
//	N(3, 3) = -M(0, 0) - M(1, 1) + M(2, 2);
//
//	// det(N - λI) = 0
//	// find the largest λ (λm)
//	//
//	// Aλ4 + Bλ3 + Cλ2 + Dλ + E = 0
//	// A = 1
//	// B = 0
//	// and so this is a so-called depressed quartic
//	// solve it using Ferrari's algorithm
//	
//	double C = -2 * (
//		M(0, 0) * M(0, 0) + M(0, 1) * M(0, 1) + M(0, 2) * M(0, 2) +
//		M(1, 0) * M(1, 0) + M(1, 1) * M(1, 1) + M(1, 2) * M(1, 2) +
//		M(2, 0) * M(2, 0) + M(2, 1) * M(2, 1) + M(2, 2) * M(2, 2));
//	
//	double D = 8 * (M(0, 0) * M(1, 2) * M(2, 1) +
//					M(1, 1) * M(2, 0) * M(0, 2) +
//					M(2, 2) * M(0, 1) * M(1, 0)) -
//			   8 * (M(0, 0) * M(1, 1) * M(2, 2) +
//					M(1, 2) * M(2, 0) * M(0, 1) +
//					M(2, 1) * M(1, 0) * M(0, 2));
//	
//	double E = 
//		(N(0,0) * N(1,1) - N(0,1) * N(0,1)) * (N(2,2) * N(3,3) - N(2,3) * N(2,3)) +
//		(N(0,1) * N(0,2) - N(0,0) * N(2,1)) * (N(2,1) * N(3,3) - N(2,3) * N(1,3)) +
//		(N(0,0) * N(1,3) - N(0,1) * N(0,3)) * (N(2,1) * N(2,3) - N(2,2) * N(1,3)) +
//		(N(0,1) * N(2,1) - N(1,1) * N(0,2)) * (N(0,2) * N(3,3) - N(2,3) * N(0,3)) +
//		(N(1,1) * N(0,3) - N(0,1) * N(1,3)) * (N(0,2) * N(2,3) - N(2,2) * N(0,3)) +
//		(N(0,2) * N(1,3) - N(2,1) * N(0,3)) * (N(0,2) * N(1,3) - N(2,1) * N(0,3));
//	
//	// solve quartic
//	double lm = LargestDepressedQuarticSolution(C, D, E);
//	
//	// calculate t = (N - λI)
//	matrix<double> li = identity_matrix<double>(4) * lm;
//	matrix<double> t = N - li;
//	
//	// calculate a matrix of cofactors for t
//	matrix<double> cf(4, 4);
//
//	const uint32 ixs[4][3] =
//	{
//		{ 1, 2, 3 },
//		{ 0, 2, 3 },
//		{ 0, 1, 3 },
//		{ 0, 1, 2 }
//	};
//
//	uint32 maxR = 0;
//	for (uint32 r = 0; r < 4; ++r)
//	{
//		const uint32* ir = ixs[r];
//		
//		for (uint32 c = 0; c < 4; ++c)
//		{
//			const uint32* ic = ixs[c];
//
//			cf(r, c) =
//				t(ir[0], ic[0]) * t(ir[1], ic[1]) * t(ir[2], ic[2]) +
//				t(ir[0], ic[1]) * t(ir[1], ic[2]) * t(ir[2], ic[0]) +
//				t(ir[0], ic[2]) * t(ir[1], ic[0]) * t(ir[2], ic[1]) -
//				t(ir[0], ic[2]) * t(ir[1], ic[1]) * t(ir[2], ic[0]) -
//				t(ir[0], ic[1]) * t(ir[1], ic[0]) * t(ir[2], ic[2]) -
//				t(ir[0], ic[0]) * t(ir[1], ic[2]) * t(ir[2], ic[1]);
//		}
//		
//		if (r > maxR and cf(r, 0) > cf(maxR, 0))
//			maxR = r;
//	}
//	
//	// NOTE the negation of the y here, why? Maybe I swapped r/c above?
//	quaternion q(cf(maxR, 0), cf(maxR, 1), -cf(maxR, 2), cf(maxR, 3));
//	q = Normalize(q);
//	
//	return q;
//}

// --------------------------------------------------------------------

Point Nudge(Point p, float offset)
{
	static std::random_device rd;
	static mt19937_64 rng(rd());

	uniform_real_distribution<> randomAngle(0, 2 * kPI);
	normal_distribution<> randomOffset(0, offset);

	float theta = randomAngle(rng);
	float phi1 = randomAngle(rng) - kPI;
	float phi2 = randomAngle(rng) - kPI;
		
	quaternion q = boost::math::spherical(1.0f, theta, phi1, phi2);

	Point r{ 0, 0, 1 };
	r.rotate(q);
	r *= randomOffset(rng);
	
	return p + r;
}

}
