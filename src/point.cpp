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

#include "cif++/cif++.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <initializer_list>
#include <numbers>
#include <optional>
#include <random>
#include <stdexcept>
#include <tuple>
#include <valarray>
#include <vector>

namespace cif
{

// --------------------------------------------------------------------

point center_points(std::vector<point> &Points)
{
	point t;

	for (point &pt : Points)
	{
		t.x += pt.x;
		t.y += pt.y;
		t.z += pt.z;
	}

	t.x /= static_cast<float>(Points.size());
	t.y /= static_cast<float>(Points.size());
	t.z /= static_cast<float>(Points.size());

	for (point &pt : Points)
	{
		pt.x -= t.x;
		pt.y -= t.y;
		pt.z -= t.z;
	}

	return t;
}

quaternion construct_for_dihedral_angle(point p1, point p2, point p3, point p4,
	float angle, float /*esd*/)
{
	p1 -= p3;
	p2 -= p3;
	p4 -= p3;
	p3 -= p3;

	auto axis = -p2;
	float dh = dihedral_angle(p1, p2, p3, p4);

	return construct_from_angle_axis(angle - dh, axis);
}

point centroid(const std::vector<point> &pts)
{
	point result;

	for (auto &pt : pts)
		result += pt;

	result /= static_cast<float>(pts.size());

	return result;
}

double RMSd(const std::vector<point> &a, const std::vector<point> &b)
{
	double sum = 0;
	for (uint32_t i = 0; i < a.size(); ++i)
	{
		std::valarray<double> d(3);

		d[0] = b[i].x - a[i].x;
		d[1] = b[i].y - a[i].y;
		d[2] = b[i].z - a[i].z;

		d *= d;

		sum += d.sum();
	}

	return std::sqrt(sum / static_cast<double>(a.size()));
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
	std::complex<double> P = -(a * a) / 12 - c;
	std::complex<double> Q = -(a * a * a) / 108 + (a * c) / 3 - (b * b) / 8;
	std::complex<double> R = -Q / 2.0 + std::sqrt((Q * Q) / 4.0 + (P * P * P) / 27.0);

	std::complex<double> U = std::pow(R, 1 / 3.0);

	std::complex<double> y;
	if (U == 0.0)
		y = -5.0 * a / 6.0 + U - std::pow(Q, 1.0 / 3.0);
	else
		y = -5.0 * a / 6.0 + U - P / (3.0 * U);

	std::complex<double> W = std::sqrt(a + 2.0 * y);

	// And to get the final result:
	// result = (±W + std::sqrt(-(3 * alpha + 2 * y ± 2 * beta / W))) / 2;
	// We want the largest result, so:

	std::valarray<double> t(4);

	t[0] = ((W + std::sqrt(-(3.0 * a + 2.0 * y + 2.0 * b / W))) / 2.0).real();
	t[1] = ((W + std::sqrt(-(3.0 * a + 2.0 * y - 2.0 * b / W))) / 2.0).real();
	t[2] = ((-W + std::sqrt(-(3.0 * a + 2.0 * y + 2.0 * b / W))) / 2.0).real();
	t[3] = ((-W + std::sqrt(-(3.0 * a + 2.0 * y - 2.0 * b / W))) / 2.0).real();

	return t.max();
}

quaternion align_points(const std::vector<point> &pa, const std::vector<point> &pb)
{
	// First calculate M, a 3x3 matrix containing the sums of products of the coordinates of A and B
	matrix3x3<double> M;

	for (uint32_t i = 0; i < pa.size(); ++i)
	{
		const point &a = pa[i];
		const point &b = pb[i];

		M(0, 0) += a.x * b.x;
		M(0, 1) += a.x * b.y;
		M(0, 2) += a.x * b.z;
		M(1, 0) += a.y * b.x;
		M(1, 1) += a.y * b.y;
		M(1, 2) += a.y * b.z;
		M(2, 0) += a.z * b.x;
		M(2, 1) += a.z * b.y;
		M(2, 2) += a.z * b.z;
	}

	// Now calculate N, a symmetric 4x4 matrix
	symmetric_matrix4x4<double> N(4);

	N(0, 0) = M(0, 0) + M(1, 1) + M(2, 2);
	N(0, 1) = M(1, 2) - M(2, 1);
	N(0, 2) = M(2, 0) - M(0, 2);
	N(0, 3) = M(0, 1) - M(1, 0);

	N(1, 1) = M(0, 0) - M(1, 1) - M(2, 2);
	N(1, 2) = M(0, 1) + M(1, 0);
	N(1, 3) = M(0, 2) + M(2, 0);

	N(2, 2) = -M(0, 0) + M(1, 1) - M(2, 2);
	N(2, 3) = M(1, 2) + M(2, 1);

	N(3, 3) = -M(0, 0) - M(1, 1) + M(2, 2);

	// det(N - λI) = 0
	// find the largest λ (λm)
	//
	// Aλ4 + Bλ3 + Cλ2 + Dλ + E = 0
	// A = 1
	// B = 0
	// and so this is a so-called depressed quartic
	// solve it using Ferrari's algorithm

	double C = -2 * (M(0, 0) * M(0, 0) + M(0, 1) * M(0, 1) + M(0, 2) * M(0, 2) +
						M(1, 0) * M(1, 0) + M(1, 1) * M(1, 1) + M(1, 2) * M(1, 2) +
						M(2, 0) * M(2, 0) + M(2, 1) * M(2, 1) + M(2, 2) * M(2, 2));

	double D = 8 * (M(0, 0) * M(1, 2) * M(2, 1) +
					   M(1, 1) * M(2, 0) * M(0, 2) +
					   M(2, 2) * M(0, 1) * M(1, 0)) -
	           8 * (M(0, 0) * M(1, 1) * M(2, 2) +
					   M(1, 2) * M(2, 0) * M(0, 1) +
					   M(2, 1) * M(1, 0) * M(0, 2));

	// E is the determinant of N:
	double E =
		(N(0, 0) * N(1, 1) - N(0, 1) * N(0, 1)) * (N(2, 2) * N(3, 3) - N(2, 3) * N(2, 3)) +
		(N(0, 1) * N(0, 2) - N(0, 0) * N(2, 1)) * (N(2, 1) * N(3, 3) - N(2, 3) * N(1, 3)) +
		(N(0, 0) * N(1, 3) - N(0, 1) * N(0, 3)) * (N(2, 1) * N(2, 3) - N(2, 2) * N(1, 3)) +
		(N(0, 1) * N(2, 1) - N(1, 1) * N(0, 2)) * (N(0, 2) * N(3, 3) - N(2, 3) * N(0, 3)) +
		(N(1, 1) * N(0, 3) - N(0, 1) * N(1, 3)) * (N(0, 2) * N(2, 3) - N(2, 2) * N(0, 3)) +
		(N(0, 2) * N(1, 3) - N(2, 1) * N(0, 3)) * (N(0, 2) * N(1, 3) - N(2, 1) * N(0, 3));

	// solve quartic
	double lambda = LargestDepressedQuarticSolution(C, D, E);

	// calculate t = (N - λI)
	matrix<double> t(N - identity_matrix(4) * lambda);

	// calculate a matrix of cofactors for t
	auto cf = matrix_cofactors(t);

	int maxR = 0;
	double maxCF = std::abs(cf(0, 0));

	for (int r = 1; r < 4; ++r)
	{
		auto cfr = std::abs(cf(r, 0));
		if (maxCF < cfr)
		{
			maxCF = cfr;
			maxR = r;
		}
	}

	quaternion q(
		static_cast<float>(cf(maxR, 0)),
		static_cast<float>(cf(maxR, 1)),
		static_cast<float>(cf(maxR, 2)),
		static_cast<float>(cf(maxR, 3)));
	q = normalize(q);

	return q;
}

// --------------------------------------------------------------------

std::tuple<point, float> smallest_sphere_around_2_points(std::array<point, 2> pts)
{
	return { (pts[0] + pts[1]) / 2.0f, distance(pts[0], pts[1]) / 2.0f };
}

std::tuple<point, float> smallest_sphere_around_3_points(std::array<point, 3> pts)
{
	// Find two bisectors
	auto vz = cross(pts[1] - pts[0], pts[2] - pts[0]);

	auto bs1 = cross(vz, pts[1] - pts[0]);
	normalize(bs1);

	auto v1 = (pts[1] - pts[0]);
	normalize(v1);

	auto s1 = pts[0] + (distance(pts[1], pts[0]) / 2) * v1;

	auto bs2 = cross(vz, pts[2] - pts[0]);
	normalize(bs2);

	auto v2 = (pts[2] - pts[0]);
	normalize(v2);

	auto s2 = pts[0] + (distance(pts[2], pts[0]) / 2) * v2;

	auto c = line_line_intersection(s1, s1 + bs1, s2, s2 + bs2);
	if (c)
		return { *c, distance(*c, pts[0]) };

	// Colinear points I guess, try something else
	auto l1 = distance_squared(pts[0], pts[1]);
	auto l2 = distance_squared(pts[0], pts[2]);
	auto l3 = distance_squared(pts[1], pts[2]);

	if (l1 > l2 and l1 > l3)
		return smallest_sphere_around_2_points({ pts[0], pts[1] });
	else if (l2 > l1 and l2 > l3)
		return smallest_sphere_around_2_points({ pts[0], pts[2] });
	else
		return smallest_sphere_around_2_points({ pts[1], pts[2] });
}

std::tuple<point, float> smallest_sphere_around_4_points(std::array<point, 4> pts)
{
	auto t0 = -norm_squared(pts[0]);
	auto t1 = -norm_squared(pts[1]);
	auto t2 = -norm_squared(pts[2]);
	auto t3 = -norm_squared(pts[3]);

	// clang-format off
	matrix4x4<float> Tm({
		pts[0].x, pts[0].y, pts[0].z, 1,
		pts[1].x, pts[1].y, pts[1].z, 1,
		pts[2].x, pts[2].y, pts[2].z, 1,
		pts[3].x, pts[3].y, pts[3].z, 1
	});
	auto T = determinant(Tm);

	if (T != 0)
	{
		matrix4x4<float> Dm({
			t0, pts[0].y, pts[0].z, 1,
			t1, pts[1].y, pts[1].z, 1,
			t2, pts[2].y, pts[2].z, 1,
			t3, pts[3].y, pts[3].z, 1
		});
		auto D = determinant(Dm) / T;
		
		matrix4x4<float> Em({
			pts[0].x, t0, pts[0].z, 1,
			pts[1].x, t1, pts[1].z, 1,
			pts[2].x, t2, pts[2].z, 1,
			pts[3].x, t3, pts[3].z, 1
		});
		auto E = determinant(Em) / T;
		
		matrix4x4<float> Fm({
			pts[0].x, pts[0].y, t0, 1,
			pts[1].x, pts[1].y, t1, 1,
			pts[2].x, pts[2].y, t2, 1,
			pts[3].x, pts[3].y, t3, 1
		});
		
		auto F = determinant(Fm) / T;
		
		matrix4x4<float> Gm({
			pts[0].x, pts[0].y, pts[0].z, t0,
			pts[1].x, pts[1].y, pts[1].z, t1,
			pts[2].x, pts[2].y, pts[2].z, t2,
			pts[3].x, pts[3].y, pts[3].z, t3
		});
		auto G = determinant(Gm) / T;
		
		point center{ -D / 2, -E / 2, -F / 2 };
		float radius = std::sqrt(D * D + E * E + F * F - 4 * G) / 2;

		// clang-format on

		return { center, radius };
	}

	// Perhaps some colinear points, try something else:

	for (auto ix : std::initializer_list<std::array<size_t, 4>>{
			 { 1, 2, 3, 0 },
			 { 0, 2, 3, 1 },
			 { 0, 1, 3, 2 },
			 { 0, 1, 2, 3 },
		 })
	{
		auto [center, radius] =
			smallest_sphere_around_3_points({ pts[ix[0]], pts[ix[1]], pts[ix[2]] });

		if (distance(pts[ix[3]], center) <= radius)
			return { center, radius };
	}

	assert(false);
	exit(1);
}

std::tuple<point, float> smallest_sphere_around_all_points(std::vector<point> P, std::vector<point> R)
{
	if (P.empty() or R.size() == 4)
	{
		switch (R.size())
		{
			case 1:
				return { R[0], 0 };

			case 2:
				return smallest_sphere_around_2_points({ R[0], R[1] });

			case 3:
				return smallest_sphere_around_3_points({ R[0], R[1], R[2] });

			case 4:
				return smallest_sphere_around_4_points({ R[0], R[1], R[2], R[3] });

			default:
				assert(false);
		}
	}

	auto p = P.back();
	P.pop_back();

	auto [c, r] = smallest_sphere_around_all_points(P, R);
	assert(not std::isnan(r));
	if (distance(c, p) <= r)
		return { c, r };

	R.emplace_back(p);
	return smallest_sphere_around_all_points(P, R);
}

bool point_in_circle(point p, std::vector<point> c)
{
	switch (c.size())
	{
		case 0:
			return false;

		case 1:
			return p == c.front();

		case 2:
		{
			auto [center, radius] = smallest_sphere_around_2_points({ c[0], c[1] });
			return distance_squared(p, center) <= radius * radius;
		}

		case 3:
		{
			auto [center, radius] = smallest_sphere_around_3_points({ c[0], c[1], c[2] });
			return distance_squared(p, center) <= radius * radius;
		}

		case 4:
		{
			auto [center, radius] = smallest_sphere_around_4_points({ c[0], c[1], c[2], c[3] });
			return distance_squared(p, center) <= radius * radius;
		}

		default:
			assert(false);
			throw std::runtime_error("Error finding smallest sphere");
	}
}

std::tuple<point, float> smallest_sphere_around_points(std::vector<point> pts)
{
	std::random_device rd;
	std::mt19937 g(rd());

	std::shuffle(pts.begin(), pts.end(), g);

	std::vector<size_t> cix;

	auto cirle_points = [&]()
	{
		std::vector<point> result;
		for (auto ix : cix)
			result.emplace_back(pts[ix]);
		return result;
	};

	size_t i = 0;
	while (i < pts.size())
	{
		if (std::ranges::find(cix, i) != cix.end() or
			point_in_circle(pts[i], cirle_points()))
		{
			++i;
		}
		else
		{
			std::erase_if(cix, [i](size_t j)
				{ return j < i; });
			cix.push_back(i);
			if (cix.size() < 4)
				i = 0;
			else
				++i;
		}
	}

	switch (cix.size())
	{
		case 1:
			return { pts[cix[0]], 0 };
		case 2:
			return smallest_sphere_around_2_points({ pts[cix[0]], pts[cix[1]] });
		case 3:
			return smallest_sphere_around_3_points({ pts[cix[0]], pts[cix[1]], pts[cix[2]] });
		case 4:
			return smallest_sphere_around_4_points({ pts[cix[0]], pts[cix[1]], pts[cix[2]], pts[cix[3]] });
		default:
			assert(false);
			throw std::runtime_error(std::format("Error finding smallest sphere (cix size: {}, pts size: {})",
				cix.size(), pts.size()));
	}
}

} // namespace cif
