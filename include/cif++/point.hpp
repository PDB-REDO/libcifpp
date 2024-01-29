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
#include <functional>
#include <valarray>

#if __has_include(<clipper/core/coords.h>)
#define HAVE_LIBCLIPPER 1
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#include <clipper/core/coords.h>
#pragma GCC diagnostic pop
#endif

/** \file point.hpp
 *
 * This file contains the definition for *cif::point* as well as
 * lots of routines and classes that can manipulate points.
 */

namespace cif
{

// --------------------------------------------------------------------

/// \brief Our value for Pi
const double
	kPI = 3.141592653589793238462643383279502884;

// --------------------------------------------------------------------
/**
 * @brief A stripped down quaternion implementation, based on boost::math::quaternion
 *
 * We use quaternions to do rotations in 3d space. Quaternions are faster than
 * matrix calculations and they also suffer less from drift caused by rounding
 * errors.
 *
 * Like complex number, quaternions do have a meaningful notion of "real part",
 * but unlike them there is no meaningful notion of "imaginary part".
 * Instead there is an "unreal part" which itself is a quaternion, and usually
 * nothing simpler (as opposed to the complex number case).
 * However, for practicality, there are accessors for the other components
 * (these are necessary for the templated copy constructor, for instance).
 *
 * @note Quaternion multiplication is *NOT* commutative;
 * symbolically, "q *= rhs;" means "q = q * rhs;"
 * and "q /= rhs;" means "q = q * inverse_of(rhs);"
 */

template <typename T>
class quaternion_type
{
  public:
	/// \brief the value type of the elements, usually this is float
	using value_type = T;

	/// \brief constructor with the four members
	constexpr explicit quaternion_type(value_type const &value_a = {}, value_type const &value_b = {}, value_type const &value_c = {}, value_type const &value_d = {})
		: a(value_a)
		, b(value_b)
		, c(value_c)
		, d(value_d)
	{
	}

	/// \brief constructor taking two complex values as input
	constexpr explicit quaternion_type(std::complex<value_type> const &z0, std::complex<value_type> const &z1 = std::complex<value_type>())
		: a(z0.real())
		, b(z0.imag())
		, c(z1.real())
		, d(z1.imag())
	{
	}

	constexpr quaternion_type(quaternion_type const &) = default; ///< Copy constructor
	constexpr quaternion_type(quaternion_type &&) = default;      ///< Copy constructor

	/// \brief Copy constructor accepting a quaternion with a different value_type
	template <typename X>
	constexpr explicit quaternion_type(quaternion_type<X> const &rhs)
		: a(static_cast<value_type>(rhs.a))
		, b(static_cast<value_type>(rhs.b))
		, c(static_cast<value_type>(rhs.c))
		, d(static_cast<value_type>(rhs.d))
	{
	}

	// accessors

	/// \brief See class description, return the *real* part of the quaternion
	constexpr value_type real() const
	{
		return a;
	}

	/// \brief See class description, return the *unreal* part of the quaternion
	constexpr quaternion_type unreal() const
	{
		return { 0, b, c, d };
	}

	/// \brief swap
	constexpr void swap(quaternion_type &o)
	{
		std::swap(a, o.a);
		std::swap(b, o.b);
		std::swap(c, o.c);
		std::swap(d, o.d);
	}

	// assignment operators

	/// \brief Assignment operator accepting a quaternion with optionally another value_type
	template <typename X>
	constexpr quaternion_type &operator=(quaternion_type<X> const &rhs)
	{
		a = static_cast<value_type>(rhs.a);
		b = static_cast<value_type>(rhs.b);
		c = static_cast<value_type>(rhs.c);
		d = static_cast<value_type>(rhs.d);

		return *this;
	}

	/// \brief Assignment operator
	constexpr quaternion_type &operator=(quaternion_type const &rhs)
	{
		a = rhs.a;
		b = rhs.b;
		c = rhs.c;
		d = rhs.d;

		return *this;
	}

	/// \brief Assignment operator that sets the *real* part to @a rhs and the *unreal* parts to zero
	constexpr quaternion_type &operator=(value_type const &rhs)
	{
		a = rhs;

		b = c = d = static_cast<value_type>(0);

		return *this;
	}

	/// \brief Assignment operator that sets the *real* part to the real part of @a rhs
	/// and the first *unreal* part to the imaginary part of of @a rhs. The other *unreal*
	// parts are set to zero.
	constexpr quaternion_type &operator=(std::complex<value_type> const &rhs)
	{
		a = rhs.real();
		b = rhs.imag();

		c = d = static_cast<value_type>(0);

		return *this;
	}

	// other assignment-related operators

	/// \brief operator += adding value @a rhs to the *real* part
	constexpr quaternion_type &operator+=(value_type const &rhs)
	{
		a += rhs;
		return *this;
	}

	/// \brief operator += adding the real part of @a rhs to the *real* part
	/// and the imaginary part of @a rhs to the first *unreal* part
	constexpr quaternion_type &operator+=(std::complex<value_type> const &rhs)
	{
		a += std::real(rhs);
		b += std::imag(rhs);
		return *this;
	}

	/// \brief operator += adding the parts of @a rhs to the equivalent part of this
	template <class X>
	constexpr quaternion_type &operator+=(quaternion_type<X> const &rhs)
	{
		a += rhs.a;
		b += rhs.b;
		c += rhs.c;
		d += rhs.d;
		return *this;
	}

	/// \brief operator -= subtracting value @a rhs from the *real* part
	constexpr quaternion_type &operator-=(value_type const &rhs)
	{
		a -= rhs;
		return *this;
	}

	/// \brief operator -= subtracting the real part of @a rhs from the *real* part
	/// and the imaginary part of @a rhs from the first *unreal* part
	constexpr quaternion_type &operator-=(std::complex<value_type> const &rhs)
	{
		a -= std::real(rhs);
		b -= std::imag(rhs);
		return *this;
	}

	/// \brief operator -= subtracting the parts of @a rhs from the equivalent part of this
	template <class X>
	constexpr quaternion_type &operator-=(quaternion_type<X> const &rhs)
	{
		a -= rhs.a;
		b -= rhs.b;
		c -= rhs.c;
		d -= rhs.d;
		return *this;
	}

	/// \brief multiply all parts with value @a rhs
	constexpr quaternion_type &operator*=(value_type const &rhs)
	{
		a *= rhs;
		b *= rhs;
		c *= rhs;
		d *= rhs;
		return *this;
	}

	/// \brief multiply with complex number @a rhs
	constexpr quaternion_type &operator*=(std::complex<value_type> const &rhs)
	{
		value_type ar = rhs.real();
		value_type br = rhs.imag();
		quaternion_type result(a * ar - b * br, a * br + b * ar, c * ar + d * br, -c * br + d * ar);
		swap(result);
		return *this;
	}

	/// \brief multiply @a a with @a b and return the result
	friend constexpr quaternion_type operator*(const quaternion_type &a, const quaternion_type &b)
	{
		auto result = a;
		result *= b;
		return result;
	}

	/// \brief multiply with quaternion @a rhs
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

	/// \brief divide all parts by @a rhs
	constexpr quaternion_type &operator/=(value_type const &rhs)
	{
		a /= rhs;
		b /= rhs;
		c /= rhs;
		d /= rhs;
		return *this;
	}

	/// \brief divide by complex number @a rhs
	constexpr quaternion_type &operator/=(std::complex<value_type> const &rhs)
	{
		value_type ar = rhs.real();
		value_type br = rhs.imag();
		value_type denominator = ar * ar + br * br;
		quaternion_type result((+a * ar + b * br) / denominator, (-a * br + b * ar) / denominator, (+c * ar - d * br) / denominator, (+c * br + d * ar) / denominator);
		swap(result);
		return *this;
	}

	/// \brief divide by quaternion @a rhs
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

	/// \brief normalise the values so that the length of the result is exactly 1
	friend constexpr quaternion_type normalize(quaternion_type q)
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

	/// \brief return the conjugate of this
	friend constexpr quaternion_type conj(quaternion_type q)
	{
		return quaternion_type{ +q.a, -q.b, -q.c, -q.d };
	}

	constexpr value_type get_a() const { return a; } ///< Return part a
	constexpr value_type get_b() const { return b; } ///< Return part b
	constexpr value_type get_c() const { return c; } ///< Return part c
	constexpr value_type get_d() const { return d; } ///< Return part d

	/// \brief compare with @a rhs
	constexpr bool operator==(const quaternion_type &rhs) const
	{
		return a == rhs.a and b == rhs.b and c == rhs.c and d == rhs.d;
	}

	/// \brief compare with @a rhs
	constexpr bool operator!=(const quaternion_type &rhs) const
	{
		return a != rhs.a or b != rhs.b or c != rhs.c or d != rhs.d;
	}

	/// \brief test for all zero values
	constexpr operator bool() const
	{
		return a != 0 or b != 0 or c != 0 or d != 0;
	}

  private:
	value_type a, b, c, d;
};

/**
 * @brief This code is similar to the code in boost so I copy the documentation as well:
 *
 * > spherical is a simple transposition of polar, it takes as inputs a (positive)
 * > magnitude and a point on the hypersphere, given by three angles. The first of
 * > these, theta has a natural range of -pi to +pi, and the other two have natural
 * > ranges of -pi/2 to +pi/2 (as is the case with the usual spherical coordinates in
 * > **R**<sup>3</sup>). Due to the many symmetries and periodicities, nothing untoward happens if
 * > the magnitude is negative or the angles are outside their natural ranges. The
 * > expected degeneracies (a magnitude of zero ignores the angles settings...) do
 * > happen however.
 */

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

/// \brief By default we use the float version of a quaternion
using quaternion = quaternion_type<float>;

// --------------------------------------------------------------------

/**
 * @brief 3D point: a location with x, y and z coordinates as floating point.
 *
 * Note that you can simply use structured binding to get access to the
 * individual parts like so:
 *
 * @code{.cpp}
 * float x, y, z;
 * tie(x, y, z) = atom.get_location();
 * @endcode
 */

template <typename F>
struct point_type
{
	/// \brief the value type of the x, y and z members
	using value_type = F;

	value_type m_x, ///< The x part of the location
		m_y,        ///< The y part of the location
		m_z;        ///< The z part of the location

	/// \brief default constructor, initialises the values to zero
	constexpr point_type()
		: m_x(0)
		, m_y(0)
		, m_z(0)
	{
	}

	/// \brief constructor taking three values
	constexpr point_type(value_type x, value_type y, value_type z)
		: m_x(x)
		, m_y(y)
		, m_z(z)
	{
	}

	/// \brief Copy constructor
	template <typename PF>
	constexpr point_type(const point_type<PF> &pt)
		: m_x(static_cast<F>(pt.m_x))
		, m_y(static_cast<F>(pt.m_y))
		, m_z(static_cast<F>(pt.m_z))
	{
	}

	/// \brief constructor taking a tuple of three values
	constexpr point_type(const std::tuple<value_type, value_type, value_type> &pt)
		: point_type(std::get<0>(pt), std::get<1>(pt), std::get<2>(pt))
	{
	}

#if HAVE_LIBCLIPPER
	/// \brief Construct a point using the values in clipper coordinate @a pt
	constexpr point_type(const clipper::Coord_orth &pt)
		: m_x(pt[0])
		, m_y(pt[1])
		, m_z(pt[2])
	{
	}

	/// \brief Assign a point using the values in clipper coordinate @a rhs
	constexpr point_type &operator=(const clipper::Coord_orth &rhs)
	{
		m_x = rhs[0];
		m_y = rhs[1];
		m_z = rhs[2];
		return *this;
	}
#endif

	/// \brief Assignment operator
	template <typename PF>
	constexpr point_type &operator=(const point_type<PF> &rhs)
	{
		m_x = static_cast<F>(rhs.m_x);
		m_y = static_cast<F>(rhs.m_y);
		m_z = static_cast<F>(rhs.m_z);
		return *this;
	}

	constexpr value_type &get_x() { return m_x; }      ///< Get a reference to x
	constexpr value_type get_x() const { return m_x; } ///< Get the value of x
	constexpr void set_x(value_type x) { m_x = x; }    ///< Set the value of x to @a x

	constexpr value_type &get_y() { return m_y; }      ///< Get a reference to y
	constexpr value_type get_y() const { return m_y; } ///< Get the value of y
	constexpr void set_y(value_type y) { m_y = y; }    ///< Set the value of y to @a y

	constexpr value_type &get_z() { return m_z; }      ///< Get a reference to z
	constexpr value_type get_z() const { return m_z; } ///< Get the value of z
	constexpr void set_z(value_type z) { m_z = z; }    ///< Set the value of z to @a z

	/// \brief add @a rhs
	constexpr point_type &operator+=(const point_type &rhs)
	{
		m_x += rhs.m_x;
		m_y += rhs.m_y;
		m_z += rhs.m_z;

		return *this;
	}

	/// \brief add @a d to all members
	constexpr point_type &operator+=(value_type d)
	{
		m_x += d;
		m_y += d;
		m_z += d;

		return *this;
	}

	/// \brief Add the points @a lhs and @a rhs and return the result
	template <typename F2>
	friend constexpr auto operator+(const point_type &lhs, const point_type<F2> &rhs)
	{
		return point_type<std::common_type_t<value_type, F2>>(lhs.m_x + rhs.m_x, lhs.m_y + rhs.m_y, lhs.m_z + rhs.m_z);
	}

	/// \brief subtract @a rhs
	constexpr point_type &operator-=(const point_type &rhs)
	{
		m_x -= rhs.m_x;
		m_y -= rhs.m_y;
		m_z -= rhs.m_z;

		return *this;
	}

	/// \brief subtract @a d from all members
	constexpr point_type &operator-=(value_type d)
	{
		m_x -= d;
		m_y -= d;
		m_z -= d;

		return *this;
	}

	/// \brief Subtract the points @a lhs and @a rhs and return the result
	template <typename F2>
	friend constexpr auto operator-(const point_type &lhs, const point_type<F2> &rhs)
	{
		return point_type<std::common_type_t<value_type, F2>>(lhs.m_x - rhs.m_x, lhs.m_y - rhs.m_y, lhs.m_z - rhs.m_z);
	}

	/// \brief Return the negative copy of @a pt
	friend constexpr point_type operator-(const point_type &pt)
	{
		return point_type(-pt.m_x, -pt.m_y, -pt.m_z);
	}

	/// \brief multiply all members with @a rhs
	constexpr point_type &operator*=(value_type rhs)
	{
		m_x *= rhs;
		m_y *= rhs;
		m_z *= rhs;
		return *this;
	}

	/// \brief multiply point @a pt with value @a f and return the result
	template <typename F2>
	friend constexpr auto operator*(const point_type &pt, F2 f)
	{
		return point_type<std::common_type_t<value_type, F2>>(pt.m_x * f, pt.m_y * f, pt.m_z * f);
	}

	/// \brief multiply point @a pt with value @a f and return the result
	template <typename F2>
	friend constexpr auto operator*(F2 f, const point_type &pt)
	{
		return point_type<std::common_type_t<value_type, F2>>(pt.m_x * f, pt.m_y * f, pt.m_z * f);
	}

	/// \brief divide all members by @a rhs
	constexpr point_type &operator/=(value_type rhs)
	{
		m_x /= rhs;
		m_y /= rhs;
		m_z /= rhs;
		return *this;
	}

	/// \brief divide point @a pt by value @a f and return the result
	template <typename F2>
	friend constexpr auto operator/(const point_type &pt, F2 f)
	{
		return point_type<std::common_type_t<value_type, F2>>(pt.m_x / f, pt.m_y / f, pt.m_z / f);
	}

	/**
	 * @brief looking at this point as a vector, normalise it which
	 * means dividing all members by the length making the length
	 * effectively 1.
	 *
	 * @return The previous length of this vector
	 */
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

	/// \brief Rotate this point using the quaterion @a q
	constexpr void rotate(const quaternion &q)
	{
		quaternion_type<value_type> p(0, m_x, m_y, m_z);

		p = q * p * conj(q);

		m_x = p.get_b();
		m_y = p.get_c();
		m_z = p.get_d();
	}

	/// \brief Rotate this point using the quaterion @a q by first
	/// moving the point to @a pivot and after rotating moving it
	/// back
	constexpr void rotate(const quaternion &q, point_type pivot)
	{
		operator-=(pivot);
		rotate(q);
		operator+=(pivot);
	}

#if HAVE_LIBCLIPPER
	/// \brief Make it possible to pass a point to clipper functions expecting a clipper coordinate
	operator clipper::Coord_orth() const
	{
		return clipper::Coord_orth(m_x, m_y, m_z);
	}
#endif

	/// \brief Allow access to this point as if it is a tuple of three const value_type's
	constexpr operator std::tuple<const value_type &, const value_type &, const value_type &>() const
	{
		return std::make_tuple(std::ref(m_x), std::ref(m_y), std::ref(m_z));
	}

	/// \brief Allow access to this point as if it is a tuple of three value_type's
	constexpr operator std::tuple<value_type &, value_type &, value_type &>()
	{
		return std::make_tuple(std::ref(m_x), std::ref(m_y), std::ref(m_z));
	}

	/// \brief Compare with @a rhs
	constexpr bool operator==(const point_type &rhs) const
	{
		return m_x == rhs.m_x and m_y == rhs.m_y and m_z == rhs.m_z;
	}

	// consider point as a vector... perhaps I should rename point?

	/// \brief looking at the point as if it is a vector, return the squared length
	constexpr value_type length_sq() const
	{
		return m_x * m_x + m_y * m_y + m_z * m_z;
	}

	/// \brief looking at the point as if it is a vector, return the length
	constexpr value_type length() const
	{
		return std::sqrt(length_sq());
	}

	/// \brief Print out the point @a pt to @a os
	friend std::ostream &operator<<(std::ostream &os, const point_type &pt)
	{
		os << '(' << pt.m_x << ',' << pt.m_y << ',' << pt.m_z << ')';
		return os;
	}
};

/// \brief By default we use points with float value_type
using point = point_type<float>;

// --------------------------------------------------------------------
// several standard 3d operations

/// \brief return the squared distance between points @a a and @a b
template <typename F1, typename F2>
constexpr auto distance_squared(const point_type<F1> &a, const point_type<F2> &b)
{
	return (a.m_x - b.m_x) * (a.m_x - b.m_x) +
	       (a.m_y - b.m_y) * (a.m_y - b.m_y) +
	       (a.m_z - b.m_z) * (a.m_z - b.m_z);
}

/// \brief return the distance between points @a a and @a b
template <typename F1, typename F2>
constexpr auto distance(const point_type<F1> &a, const point_type<F2> &b)
{
	return std::sqrt(
		(a.m_x - b.m_x) * (a.m_x - b.m_x) +
		(a.m_y - b.m_y) * (a.m_y - b.m_y) +
		(a.m_z - b.m_z) * (a.m_z - b.m_z));
}

/// \brief return the dot product between the vectors @a a and @a b
template <typename F1, typename F2>
inline constexpr auto dot_product(const point_type<F1> &a, const point_type<F2> &b)
{
	return a.m_x * b.m_x + a.m_y * b.m_y + a.m_z * b.m_z;
}

/// \brief return the cross product between the vectors @a a and @a b
template <typename F1, typename F2>
inline constexpr auto cross_product(const point_type<F1> &a, const point_type<F2> &b)
{
	return point_type<std::common_type_t<F1, F2>>(
		a.m_y * b.m_z - b.m_y * a.m_z,
		a.m_z * b.m_x - b.m_z * a.m_x,
		a.m_x * b.m_y - b.m_x * a.m_y);
}

/// \brief return the angle in degrees between the vectors from point @a p2 to @a p1 and @a p2 to @a p3
template <typename F>
constexpr auto angle(const point_type<F> &p1, const point_type<F> &p2, const point_type<F> &p3)
{
	point_type<F> v1 = p1 - p2;
	point_type<F> v2 = p3 - p2;

	return std::acos(dot_product(v1, v2) / (v1.length() * v2.length())) * 180 / kPI;
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

/// \brief return the cosinus angle for the four points @a p1, @a p2, @a p3 and @a p4
template <typename F>
constexpr auto cosinus_angle(const point_type<F> &p1, const point_type<F> &p2, const point_type<F> &p3, const point_type<F> &p4)
{
	point_type<F> v12 = p1 - p2;
	point_type<F> v34 = p3 - p4;

	auto x = dot_product(v12, v12) * dot_product(v34, v34);

	return x > 0 ? dot_product(v12, v34) / std::sqrt(x) : 0;
}

/// \brief return the distance from point @a p to the line from @a l1 to @a l2
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
/**
 * @brief For e.g. simulated annealing, returns a new point that is moved in
 * a random direction with a distance randomly chosen from a normal
 * distribution with a stddev of offset.
 */
point nudge(point p, float offset);

// --------------------------------------------------------------------

/// \brief Return a quaternion created from angle @a angle and axis @a axis
quaternion construct_from_angle_axis(float angle, point axis);

/// \brief Return a tuple of an angle and an axis for quaternion @a q
std::tuple<double, point> quaternion_to_angle_axis(quaternion q);

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

// --------------------------------------------------------------------
/**
 * @brief Helper class to generate evenly divided points on a sphere
 *
 * We use a fibonacci sphere to calculate even distribution of the dots
 *
 * @tparam N The number of points on the sphere is 2 * N + 1
 */
template <int N>
class spherical_dots
{
  public:
	/// \brief the number of points
	constexpr static int P = 2 * N * 1;

	/// \brief the *weight* of the fibonacci sphere
	constexpr static double W = (4 * kPI) / P;

	/// \brief the internal storage type
	using array_type = typename std::array<point, P>;

	/// \brief iterator type
	using iterator = typename array_type::const_iterator;

	/// \brief singleton instance
	static spherical_dots &instance()
	{
		static spherical_dots sInstance;
		return sInstance;
	}

	/// \brief The number of points
	size_t size() const { return P; }

	/// \brief Access a point by index
	const point operator[](uint32_t inIx) const { return m_points[inIx]; }

	/// \brief iterator pointing to the first point
	iterator begin() const { return m_points.begin(); }

	/// \brief iterator pointing past the last point
	iterator end() const { return m_points.end(); }

	/// \brief return the *weight*,
	double weight() const { return W; }

	spherical_dots()
	{
		const double
			kGoldenRatio = (1 + std::sqrt(5.0)) / 2;

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
};

} // namespace cif
