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

#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <functional>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <limits>
#include <numbers>
#include <optional>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <utility>
#include <valarray>
#include <vector>

#if __has_include(<clipper/core/coords.h>)
# define HAVE_LIBCLIPPER 1
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wignored-qualifiers"
# include <clipper/core/clipper_types.h>
# include <clipper/core/coords.h>
# pragma GCC diagnostic pop
#endif

/** \file point.hpp
 *
 * This file contains the definition for *point* as well as
 * lots of routines and classes that can manipulate points.
 */

namespace cif
{

// Using glm now

template <typename T>
using quaternion_type = glm::qua<T>;

using quaternion = quaternion_type<float>;

template <typename T>
using point_type = glm::vec<3, T>;

using point = point_type<float>;

// --------------------------------------------------------------------
// several standard 3d operations

/// \brief return the squared distance between points @a a and @a b
template <typename F1, typename F2>
constexpr auto distance_squared(const point_type<F1> &a, const point_type<F2> &b)
{
	return (a.x - b.x) * (a.x - b.x) +
	       (a.y - b.y) * (a.y - b.y) +
	       (a.z - b.z) * (a.z - b.z);
}

/// \brief return the squared norm of point @a p
template <typename F>
constexpr F norm_squared(const point_type<F> &p)
{
	return p.x * p.x + p.y * p.y + p.z * p.z;
}

/// \brief return the norm of point @a p
template <typename F>
constexpr point_type<F> norm(const point_type<F> &p)
{
	return std::sqrt(norm_squared(p));
}

/// \brief return the point where two lines intersect, or an empty value if they don't intersect at all
template <typename F>
std::optional<point> line_line_intersection(const point_type<F> &p1,
	const point_type<F> &p2, const point_type<F> &p3, const point_type<F> &p4)
{
	auto p13 = p1 - p3;
	auto p43 = p4 - p3;
	if (std::abs(p43.x) < std::numeric_limits<F>::epsilon() and std::abs(p43.y) < std::numeric_limits<F>::epsilon() and std::abs(p43.z) < std::numeric_limits<F>::epsilon())
		return std::nullopt;

	auto p21 = p2 - p1;
	if (std::abs(p21.x) < std::numeric_limits<F>::epsilon() and std::abs(p21.y) < std::numeric_limits<F>::epsilon() and std::abs(p21.z) < std::numeric_limits<F>::epsilon())
		return std::nullopt;

	auto d1343 = dot(p43, p13);
	auto d4321 = dot(p43, p21);
	auto d1321 = dot(p13, p21);
	auto d4343 = dot(p43, p43);
	auto d2121 = dot(p21, p21);

	auto denom = d2121 * d4343 - d4321 * d4321;
	if (std::abs(denom) < std::numeric_limits<F>::epsilon())
		return std::nullopt;

	auto numer = d1343 * d4321 - d1321 * d4343;

	auto mua = numer / denom;
	auto mub = (d1343 + d4321 * mua) / d4343;

	auto pa = p1 + mua * p21;
	auto pb = p3 + mub * p43;

	return { (pa + pb) / 2.0f };
}

/// \brief return the angle in degrees between the vectors from point @a p2 to @a p1 and @a p2 to @a p3
template <typename F>
constexpr auto angle(const point_type<F> &p1, const point_type<F> &p2, const point_type<F> &p3)
{
	point_type<F> v1 = p1 - p2;
	point_type<F> v2 = p3 - p2;

	return std::acos(dot(v1, v2) / (v1.length() * v2.length())) * 180 / std::numbers::pi_v<F>;
}

/// \brief return the dihedral angle in degrees for the four points @a p1, @a p2, @a p3 and @a p4
///
/// See https://en.wikipedia.org/wiki/Dihedral_angle for an explanation of what a dihedral angle is
template <typename F>
constexpr auto dihedral_angle(const point_type<F> &p1, const point_type<F> &p2, const point_type<F> &p3, const point_type<F> &p4)
{
	point_type<F> v12 = p1 - p2; // vector from p2 to p1
	point_type<F> v43 = p4 - p3; // vector from p3 to p4

	point_type<F> z = p2 - p3; // vector from p3 to p2

	point_type<F> p = cross(z, v12);
	point_type<F> x = cross(z, v43);
	point_type<F> y = cross(z, x);

	auto u = dot(x, x);
	auto v = dot(y, y);

	F result = 360;
	if (u > 0 and v > 0)
	{
		u = dot(p, x) / std::sqrt(u);
		v = dot(p, y) / std::sqrt(v);
		if (u != 0 or v != 0)
			result = std::atan2(v, u) * static_cast<F>(180 / std::numbers::pi_v<F>);
	}

	return result;
}

/// \brief return the cosinus angle for the four points @a p1, @a p2, @a p3 and @a p4
template <typename F>
constexpr auto cosinus_angle(const point_type<F> &p1, const point_type<F> &p2, const point_type<F> &p3, const point_type<F> &p4)
{
	point_type<F> v12 = p1 - p2;
	point_type<F> v34 = p3 - p4;

	auto x = dot(v12, v12) * dot(v34, v34);

	return x > 0 ? dot(v12, v34) / std::sqrt(x) : 0;
}

/// \brief return the distance from point @a p to the line from @a l1 to @a l2
template <typename F>
constexpr auto distance_point_to_line(const point_type<F> &l1, const point_type<F> &l2, const point_type<F> &p)
{
	auto line = l2 - l1;
	auto p_to_l1 = p - l1;
	auto p_to_l2 = p - l2;
	auto cross = glm::cross(p_to_l1, p_to_l2);
	return cross.length() / line.length();
}

/// \brief return the smallest sphere around the points in @a pts
std::tuple<point, float> smallest_sphere_around_points(std::vector<point> pts);

// --------------------------------------------------------------------

/// \brief Return a quaternion created from angle @a angle and axis @a axis
quaternion construct_from_angle_axis(float angle, point axis);

/// \brief Return a tuple of an angle and an axis for quaternion @a q
std::tuple<float, point> quaternion_to_angle_axis(quaternion q);

/// @brief Given four points and an angle, return the quaternion required to rotate
/// point p4 along the p2-p3 axis and around point p3 to obtain the required within
/// an accuracy of esd
quaternion construct_for_dihedral_angle(point p1, point p2, point p3, point p4,
	float angle, float esd);

/// \brief Return the point that is the centroid of all the points in @a pts
point centroid(const std::vector<point> &pts);

/// \brief Move all the points in @a pts so that their centroid is at the origin
/// (0, 0, 0) and return the offset used (the former centroid)
point center_points(std::vector<point> &pts);

/// \brief Returns how the two sets of points \a a and \b b can be aligned
///
/// \param a	The first set of points
/// \param b    The second set of points
/// \result     The quaternion which should be applied to the points in \a a to
///             obtain the best superposition.
quaternion align_points(const std::vector<point> &a, const std::vector<point> &b);

/// \brief The RMSd for the points in \a a and \a b
double RMSd(const std::vector<point> &a, const std::vector<point> &b);


} // namespace cif
