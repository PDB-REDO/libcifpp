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

#include "cif++/exports.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <functional>
#include <valarray>

#if __has_include(<clipper/core/coords.h>)
#define HAVE_LIBCLIPPER 1
#include <clipper/core/coords.h>
#endif

namespace cif
{

// --------------------------------------------------------------------

const double
	kPI = 3.141592653589793238462643383279502884;

// --------------------------------------------------------------------
// A stripped down quaternion implementation, based on boost::math::quaternion
// We use quaternions to do rotations in 3d space

template <typename T>
class quaternion_type
{
  public:
	using value_type = T;

	constexpr explicit quaternion_type(value_type const &value_a = {}, value_type const &value_b = {}, value_type const &value_c = {}, value_type const &value_d = {})
		: a(value_a)
		, b(value_b)
		, c(value_c)
		, d(value_d)
	{
	}

	constexpr explicit quaternion_type(std::complex<value_type> const &z0, std::complex<value_type> const &z1 = std::complex<value_type>())
		: a(z0.real())
		, b(z0.imag())
		, c(z1.real())
		, d(z1.imag())
	{
	}

	constexpr quaternion_type(quaternion_type const &) = default;
	constexpr quaternion_type(quaternion_type &&) = default;

	template <typename X>
	constexpr explicit quaternion_type(quaternion_type<X> const &rhs)
		: a(static_cast<value_type>(rhs.a))
		, b(static_cast<value_type>(rhs.b))
		, c(static_cast<value_type>(rhs.c))
		, d(static_cast<value_type>(rhs.d))
	{
	}

	// accessors
	//
	// Note:    Like complex number, quaternions do have a meaningful notion of "real part",
	//            but unlike them there is no meaningful notion of "imaginary part".
	//            Instead there is an "unreal part" which itself is a quaternion, and usually
	//            nothing simpler (as opposed to the complex number case).
	//            However, for practicality, there are accessors for the other components
	//            (these are necessary for the templated copy constructor, for instance).

	constexpr value_type real() const
	{
		return a;
	}

	constexpr quaternion_type unreal() const
	{
		return { 0, b, c, d };
	}

	constexpr void swap(quaternion_type &o)
	{
		std::swap(a, o.a);
		std::swap(b, o.b);
		std::swap(c, o.c);
		std::swap(d, o.d);
	}

	// assignment operators

	template <typename X>
	constexpr quaternion_type &operator=(quaternion_type<X> const &rhs)
	{
		a = static_cast<value_type>(rhs.a);
		b = static_cast<value_type>(rhs.b);
		c = static_cast<value_type>(rhs.c);
		d = static_cast<value_type>(rhs.d);

		return *this;
	}

	constexpr quaternion_type &operator=(quaternion_type const &rhs)
	{
		a = rhs.a;
		b = rhs.b;
		c = rhs.c;
		d = rhs.d;

		return *this;
	}

	constexpr quaternion_type &operator=(value_type const &rhs)
	{
		a = rhs;

		b = c = d = static_cast<value_type>(0);

		return *this;
	}

	constexpr quaternion_type &operator=(std::complex<value_type> const &rhs)
	{
		a = rhs.real();
		b = rhs.imag();

		c = d = static_cast<value_type>(0);

		return *this;
	}

	// other assignment-related operators
	//
	// NOTE:    Quaternion multiplication is *NOT* commutative;
	//            symbolically, "q *= rhs;" means "q = q * rhs;"
	//            and "q /= rhs;" means "q = q * inverse_of(rhs);"

	constexpr quaternion_type &operator+=(value_type const &rhs)
	{
		a += rhs;
		return *this;
	}

	constexpr quaternion_type &operator+=(std::complex<value_type> const &rhs)
	{
		a += std::real(rhs);
		b += std::imag(rhs);
		return *this;
	}

	template <class X>
	constexpr quaternion_type &operator+=(quaternion_type<X> const &rhs)
	{
		a += rhs.a;
		b += rhs.b;
		c += rhs.c;
		d += rhs.d;
		return *this;
	}

	constexpr quaternion_type &operator-=(value_type const &rhs)
	{
		a -= rhs;
		return *this;
	}

	constexpr quaternion_type &operator-=(std::complex<value_type> const &rhs)
	{
		a -= std::real(rhs);
		b -= std::imag(rhs);
		return *this;
	}

	template <class X>
	constexpr quaternion_type &operator-=(quaternion_type<X> const &rhs)
	{
		a -= rhs.a;
		b -= rhs.b;
		c -= rhs.c;
		d -= rhs.d;
		return *this;
	}

	constexpr quaternion_type &operator*=(value_type const &rhs)
	{
		a *= rhs;
		b *= rhs;
		c *= rhs;
		d *= rhs;
		return *this;
	}

	constexpr quaternion_type &operator*=(std::complex<value_type> const &rhs)
	{
		value_type ar = rhs.real();
		value_type br = rhs.imag();
		quaternion_type result(a * ar - b * br, a * br + b * ar, c * ar + d * br, -c * br + d * ar);
		swap(result);
		return *this;
	}

	constexpr friend quaternion_type operator*(const quaternion_type &a, const quaternion_type &b)
	{
		auto result = a;
		result *= b;
		return result;
	}

	template <typename X>
	constexpr quaternion_type &operator*=(quaternion_type<X> const &rhs)
	{
		value_type ar = static_cast<value_type>(rhs.a);
		value_type br = static_cast<value_type>(rhs.b);
		value_type cr = static_cast<value_type>(rhs.c);
		value_type dr = static_cast<value_type>(rhs.d);

		quaternion_type result(a * ar - b * br - c * cr - d * dr, a * br + b * ar + c * dr - d * cr, a * cr - b * dr + c * ar + d * br, a * dr + b * cr - c * br + d * ar);
		swap(result);
		return *this;
	}

	constexpr quaternion_type &operator/=(value_type const &rhs)
	{
		a /= rhs;
		b /= rhs;
		c /= rhs;
		d /= rhs;
		return *this;
	}

	constexpr quaternion_type &operator/=(std::complex<value_type> const &rhs)
	{
		value_type ar = rhs.real();
		value_type br = rhs.imag();
		value_type denominator = ar * ar + br * br;
		quaternion_type result((+a * ar + b * br) / denominator, (-a * br + b * ar) / denominator, (+c * ar - d * br) / denominator, (+c * br + d * ar) / denominator);
		swap(result);
		return *this;
	}

	template <typename X>
	constexpr quaternion_type &operator/=(quaternion_type<X> const &rhs)
	{
		value_type ar = static_cast<value_type>(rhs.a);
		value_type br = static_cast<value_type>(rhs.b);
		value_type cr = static_cast<value_type>(rhs.c);
		value_type dr = static_cast<value_type>(rhs.d);

		value_type denominator = ar * ar + br * br + cr * cr + dr * dr;
		quaternion_type result((+a * ar + b * br + c * cr + d * dr) / denominator, (-a * br + b * ar - c * dr + d * cr) / denominator, (-a * cr + b * dr + c * ar - d * br) / denominator, (-a * dr - b * cr + c * br + d * ar) / denominator);
		swap(result);
		return *this;
	}

	constexpr friend quaternion_type normalize(quaternion_type q)
	{
		std::valarray<value_type> t(4);

		t[0] = q.a;
		t[1] = q.b;
		t[2] = q.c;
		t[3] = q.d;

		t *= t;

		value_type length = std::sqrt(t.sum());

		if (length > 0.001)
			q /= static_cast<value_type>(length);
		else
			q = quaternion_type(1, 0, 0, 0);

		return q;
	}

	constexpr friend quaternion_type conj(quaternion_type q)
	{
		return quaternion_type{ +q.a, -q.b, -q.c, -q.d };
	}

	constexpr value_type get_a() const { return a; }
	constexpr value_type get_b() const { return b; }
	constexpr value_type get_c() const { return c; }
	constexpr value_type get_d() const { return d; }

	constexpr bool operator==(const quaternion_type &rhs) const
	{
		return a == rhs.a and b == rhs.b and c == rhs.c and d == rhs.d;
	}

	constexpr bool operator!=(const quaternion_type &rhs) const
	{
		return a != rhs.a or b != rhs.b or c != rhs.c or d != rhs.d;
	}

	constexpr operator bool() const
	{
		return a != 0 or b != 0 or c != 0 or d != 0;
	}

  private:
	value_type a, b, c, d;
};

template <typename T>
inline quaternion_type<T> spherical(T const &rho, T const &theta, T const &phi1, T const &phi2)
{
	T cos_phi1 = std::cos(phi1);
	T cos_phi2 = std::cos(phi2);

	T a = std::cos(theta) * cos_phi1 * cos_phi2;
	T b = std::sin(theta) * cos_phi1 * cos_phi2;
	T c = std::sin(phi1) * cos_phi2;
	T d = std::sin(phi2);

	quaternion_type result(a, b, c, d);
	result *= rho;

	return result;
}

using quaternion = quaternion_type<float>;

// --------------------------------------------------------------------

//	point, a location with x, y and z coordinates as floating point.
//	This one is derived from a tuple<float,float,float> so
//	you can do things like:
//
//	float x, y, z;
//	tie(x, y, z) = atom.loc();

template <typename F>
struct point_type
{
	using value_type = F;

	value_type m_x, m_y, m_z;

	constexpr point_type()
		: m_x(0)
		, m_y(0)
		, m_z(0)
	{
	}

	constexpr point_type(value_type x, value_type y, value_type z)
		: m_x(x)
		, m_y(y)
		, m_z(z)
	{
	}

	template <typename PF>
	constexpr point_type(const point_type<PF> &pt)
		: m_x(static_cast<F>(pt.m_x))
		, m_y(static_cast<F>(pt.m_y))
		, m_z(static_cast<F>(pt.m_z))
	{
	}

	constexpr point_type(const std::tuple<value_type, value_type, value_type> &pt)
		: point_type(std::get<0>(pt), std::get<1>(pt), std::get<2>(pt))
	{
	}

#if HAVE_LIBCLIPPER
	constexpr point_type(const clipper::Coord_orth &pt)
		: m_x(pt[0])
		, m_y(pt[1])
		, m_z(pt[2])
	{
	}

	constexpr point_type &operator=(const clipper::Coord_orth &rhs)
	{
		m_x = rhs[0];
		m_y = rhs[1];
		m_z = rhs[2];
		return *this;
	}
#endif

	template <typename PF>
	constexpr point_type &operator=(const point_type<PF> &rhs)
	{
		m_x = static_cast<F>(rhs.m_x);
		m_y = static_cast<F>(rhs.m_y);
		m_z = static_cast<F>(rhs.m_z);
		return *this;
	}

	constexpr value_type &get_x() { return m_x; }
	constexpr value_type get_x() const { return m_x; }
	constexpr void set_x(value_type x) { m_x = x; }

	constexpr value_type &get_y() { return m_y; }
	constexpr value_type get_y() const { return m_y; }
	constexpr void set_y(value_type y) { m_y = y; }

	constexpr value_type &get_z() { return m_z; }
	constexpr value_type get_z() const { return m_z; }
	constexpr void set_z(value_type z) { m_z = z; }

	constexpr point_type &operator+=(const point_type &rhs)
	{
		m_x += rhs.m_x;
		m_y += rhs.m_y;
		m_z += rhs.m_z;

		return *this;
	}

	constexpr point_type &operator+=(value_type d)
	{
		m_x += d;
		m_y += d;
		m_z += d;

		return *this;
	}

	constexpr point_type &operator-=(const point_type &rhs)
	{
		m_x -= rhs.m_x;
		m_y -= rhs.m_y;
		m_z -= rhs.m_z;

		return *this;
	}

	constexpr point_type &operator-=(value_type d)
	{
		m_x -= d;
		m_y -= d;
		m_z -= d;

		return *this;
	}

	constexpr point_type &operator*=(value_type rhs)
	{
		m_x *= rhs;
		m_y *= rhs;
		m_z *= rhs;
		return *this;
	}

	constexpr point_type &operator/=(value_type rhs)
	{
		m_x /= rhs;
		m_y /= rhs;
		m_z /= rhs;
		return *this;
	}

	constexpr value_type normalize()
	{
		auto length = m_x * m_x + m_y * m_y + m_z * m_z;
		if (length > 0)
		{
			length = std::sqrt(length);
			operator/=(length);
		}
		return length;
	}

	constexpr void rotate(const quaternion &q)
	{
		quaternion_type<value_type> p(0, m_x, m_y, m_z);

		p = q * p * conj(q);

		m_x = p.get_b();
		m_y = p.get_c();
		m_z = p.get_d();
	}

	constexpr void rotate(const quaternion &q, point_type pivot)
	{
		operator-=(pivot);
		rotate(q);
		operator+=(pivot);
	}

#if HAVE_LIBCLIPPER
	operator clipper::Coord_orth() const
	{
		return clipper::Coord_orth(m_x, m_y, m_z);
	}
#endif

	constexpr operator std::tuple<const value_type &, const value_type &, const value_type &>() const
	{
		return std::make_tuple(std::ref(m_x), std::ref(m_y), std::ref(m_z));
	}

	constexpr operator std::tuple<value_type &, value_type &, value_type &>()
	{
		return std::make_tuple(std::ref(m_x), std::ref(m_y), std::ref(m_z));
	}

	constexpr bool operator==(const point_type &rhs) const
	{
		return m_x == rhs.m_x and m_y == rhs.m_y and m_z == rhs.m_z;
	}

	// consider point as a vector... perhaps I should rename point?
	constexpr value_type length_sq() const
	{
		return m_x * m_x + m_y * m_y + m_z * m_z;
	}

	constexpr value_type length() const
	{
		return std::sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
	}
};

using point = point_type<float>;

template <typename F>
inline constexpr std::ostream &operator<<(std::ostream &os, const point_type<F> &pt)
{
	os << '(' << pt.m_x << ',' << pt.m_y << ',' << pt.m_z << ')';
	return os;
}

template <typename F>
inline constexpr point_type<F> operator+(const point_type<F> &lhs, const point_type<F> &rhs)
{
	return point_type<F>(lhs.m_x + rhs.m_x, lhs.m_y + rhs.m_y, lhs.m_z + rhs.m_z);
}

template <typename F>
inline constexpr point_type<F> operator-(const point_type<F> &lhs, const point_type<F> &rhs)
{
	return point_type<F>(lhs.m_x - rhs.m_x, lhs.m_y - rhs.m_y, lhs.m_z - rhs.m_z);
}

template <typename F>
inline constexpr point_type<F> operator-(const point_type<F> &pt)
{
	return point_type<F>(-pt.m_x, -pt.m_y, -pt.m_z);
}

template <typename F>
inline constexpr point_type<F> operator*(const point_type<F> &pt, F f)
{
	return point_type<F>(pt.m_x * f, pt.m_y * f, pt.m_z * f);
}

template <typename F>
inline constexpr point_type<F> operator*(F f, const point_type<F> &pt)
{
	return point_type<F>(pt.m_x * f, pt.m_y * f, pt.m_z * f);
}

template <typename F>
inline constexpr point_type<F> operator/(const point_type<F> &pt, F f)
{
	return point_type<F>(pt.m_x / f, pt.m_y / f, pt.m_z / f);
}

// --------------------------------------------------------------------
// several standard 3d operations

template <typename F>
inline constexpr auto distance_squared(const point_type<F> &a, const point_type<F> &b)
{
	return (a.m_x - b.m_x) * (a.m_x - b.m_x) +
	       (a.m_y - b.m_y) * (a.m_y - b.m_y) +
	       (a.m_z - b.m_z) * (a.m_z - b.m_z);
}

template <typename F>
inline constexpr auto distance(const point_type<F> &a, const point_type<F> &b)
{
	return std::sqrt(
		(a.m_x - b.m_x) * (a.m_x - b.m_x) +
		(a.m_y - b.m_y) * (a.m_y - b.m_y) +
		(a.m_z - b.m_z) * (a.m_z - b.m_z));
}

template <typename F>
inline constexpr auto dot_product(const point_type<F> &a, const point_type<F> &b)
{
	return a.m_x * b.m_x + a.m_y * b.m_y + a.m_z * b.m_z;
}

template <typename F>
inline constexpr point_type<F> cross_product(const point_type<F> &a, const point_type<F> &b)
{
	return point_type<F>(a.m_y * b.m_z - b.m_y * a.m_z,
		a.m_z * b.m_x - b.m_z * a.m_x,
		a.m_x * b.m_y - b.m_x * a.m_y);
}

template <typename F>
constexpr auto angle(const point_type<F> &p1, const point_type<F> &p2, const point_type<F> &p3)
{
	point_type<F> v1 = p1 - p2;
	point_type<F> v2 = p3 - p2;

	return std::acos(dot_product(v1, v2) / (v1.length() * v2.length())) * 180 / kPI;
}

template <typename F>
constexpr auto dihedral_angle(const point_type<F> &p1, const point_type<F> &p2, const point_type<F> &p3, const point_type<F> &p4)
{
	point_type<F> v12 = p1 - p2; // vector from p2 to p1
	point_type<F> v43 = p4 - p3; // vector from p3 to p4

	point_type<F> z = p2 - p3; // vector from p3 to p2

	point_type<F> p = cross_product(z, v12);
	point_type<F> x = cross_product(z, v43);
	point_type<F> y = cross_product(z, x);

	auto u = dot_product(x, x);
	auto v = dot_product(y, y);

	F result = 360;
	if (u > 0 and v > 0)
	{
		u = dot_product(p, x) / std::sqrt(u);
		v = dot_product(p, y) / std::sqrt(v);
		if (u != 0 or v != 0)
			result = std::atan2(v, u) * static_cast<F>(180 / kPI);
	}

	return result;
}

template <typename F>
constexpr auto cosinus_angle(const point_type<F> &p1, const point_type<F> &p2, const point_type<F> &p3, const point_type<F> &p4)
{
	point_type<F> v12 = p1 - p2;
	point_type<F> v34 = p3 - p4;

	auto x = dot_product(v12, v12) * dot_product(v34, v34);

	return x > 0 ? dot_product(v12, v34) / std::sqrt(x) : 0;
}

template <typename F>
constexpr auto distance_point_to_line(const point_type<F> &l1, const point_type<F> &l2, const point_type<F> &p)
{
	auto line = l2 - l1;
	auto p_to_l1 = p - l1;
	auto p_to_l2 = p - l2;
	auto cross = cross_product(p_to_l1, p_to_l2);
	return cross.length() / line.length();
}

// --------------------------------------------------------------------
// For e.g. simulated annealing, returns a new point that is moved in
// a random direction with a distance randomly chosen from a normal
// distribution with a stddev of offset.

point nudge(point p, float offset);

// --------------------------------------------------------------------

quaternion construct_from_angle_axis(float angle, point axis);
std::tuple<double, point> quaternion_to_angle_axis(quaternion q);

/// @brief Given four points and an angle, return the quaternion required to rotate
/// point p4 along the p2-p3 axis and around point p3 to obtain the required within
/// an accuracy of esd
quaternion construct_for_dihedral_angle(point p1, point p2, point p3, point p4,
	float angle, float esd);

point centroid(const std::vector<point> &Points);
point center_points(std::vector<point> &Points);

/// \brief Returns how the two sets of points \a a and \b b can be aligned
///
/// \param a	The first set of points
/// \param b    The second set of points
/// \result     The quaternion which should be applied to the points in \a a to
///             obtain the best superposition.
quaternion align_points(const std::vector<point> &a, const std::vector<point> &b);

/// \brief The RMSd for the points in \a a and \a b
double RMSd(const std::vector<point> &a, const std::vector<point> &b);

// --------------------------------------------------------------------
// Helper class to generate evenly divided points on a sphere
// we use a fibonacci sphere to calculate even distribution of the dots

template <int N>
class spherical_dots
{
  public:
	constexpr static int P = 2 * N * 1;

	using array_type = typename std::array<point, P>;
	using iterator = typename array_type::const_iterator;

	static spherical_dots &instance()
	{
		static spherical_dots sInstance;
		return sInstance;
	}

	size_t size() const { return m_points.size(); }
	const point operator[](uint32_t inIx) const { return m_points[inIx]; }
	iterator begin() const { return m_points.begin(); }
	iterator end() const { return m_points.end(); }

	double weight() const { return m_weight; }

	spherical_dots()
	{
		const double
			kGoldenRatio = (1 + std::sqrt(5.0)) / 2;

		m_weight = (4 * kPI) / P;

		auto p = m_points.begin();

		for (int32_t i = -N; i <= N; ++i)
		{
			double lat = std::asin((2.0 * i) / P);
			double lon = std::fmod(i, kGoldenRatio) * 2 * kPI / kGoldenRatio;

			p->m_x = std::sin(lon) * std::cos(lat);
			p->m_y = std::cos(lon) * std::cos(lat);
			p->m_z = std::sin(lat);

			++p;
		}
	}

  private:
	array_type m_points;
	double m_weight;
};

} // namespace cif
