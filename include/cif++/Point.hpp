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

#include <functional>

#if HAVE_LIBCLIPPER
#include <clipper/core/coords.h>
#endif

#include <boost/math/quaternion.hpp>

namespace mmcif
{

typedef boost::math::quaternion<float> Quaternion;

const double
	kPI = 3.141592653589793238462643383279502884;

// --------------------------------------------------------------------

//	Point, a location with x, y and z coordinates as floating point.
//	This one is derived from a tuple<float,float,float> so
//	you can do things like:
//
//	float x, y, z;
//	tie(x, y, z) = atom.loc();

template <typename F>
struct PointF
{
	typedef F FType;

	FType mX, mY, mZ;

	PointF()
		: mX(0)
		, mY(0)
		, mZ(0)
	{
	}
	PointF(FType x, FType y, FType z)
		: mX(x)
		, mY(y)
		, mZ(z)
	{
	}

	template <typename PF>
	PointF(const PointF<PF> &pt)
		: mX(static_cast<F>(pt.mX))
		, mY(static_cast<F>(pt.mY))
		, mZ(static_cast<F>(pt.mZ))
	{
	}

#if HAVE_LIBCLIPPER
	PointF(const clipper::Coord_orth &pt)
		: mX(pt[0])
		, mY(pt[1])
		, mZ(pt[2])
	{
	}

	PointF &operator=(const clipper::Coord_orth &rhs)
	{
		mX = rhs[0];
		mY = rhs[1];
		mZ = rhs[2];
		return *this;
	}
#endif

	template <typename PF>
	PointF &operator=(const PointF<PF> &rhs)
	{
		mX = static_cast<F>(rhs.mX);
		mY = static_cast<F>(rhs.mY);
		mZ = static_cast<F>(rhs.mZ);
		return *this;
	}

	FType &getX() { return mX; }
	FType getX() const { return mX; }
	void setX(FType x) { mX = x; }

	FType &getY() { return mY; }
	FType getY() const { return mY; }
	void setY(FType y) { mY = y; }

	FType &getZ() { return mZ; }
	FType getZ() const { return mZ; }
	void setZ(FType z) { mZ = z; }

	PointF &operator+=(const PointF &rhs)
	{
		mX += rhs.mX;
		mY += rhs.mY;
		mZ += rhs.mZ;

		return *this;
	}

	PointF &operator+=(FType d)
	{
		mX += d;
		mY += d;
		mZ += d;

		return *this;
	}

	PointF &operator-=(const PointF &rhs)
	{
		mX -= rhs.mX;
		mY -= rhs.mY;
		mZ -= rhs.mZ;

		return *this;
	}

	PointF &operator-=(FType d)
	{
		mX -= d;
		mY -= d;
		mZ -= d;

		return *this;
	}

	PointF &operator*=(FType rhs)
	{
		mX *= rhs;
		mY *= rhs;
		mZ *= rhs;
		return *this;
	}

	PointF &operator/=(FType rhs)
	{
		mX /= rhs;
		mY /= rhs;
		mZ /= rhs;
		return *this;
	}

	FType normalize()
	{
		auto length = mX * mX + mY * mY + mZ * mZ;
		if (length > 0)
		{
			length = std::sqrt(length);
			operator/=(length);
		}
		return length;
	}

	void rotate(const boost::math::quaternion<FType> &q)
	{
		boost::math::quaternion<FType> p(0, mX, mY, mZ);

		p = q * p * boost::math::conj(q);

		mX = p.R_component_2();
		mY = p.R_component_3();
		mZ = p.R_component_4();
	}

#if HAVE_LIBCLIPPER
	operator clipper::Coord_orth() const
	{
		return clipper::Coord_orth(mX, mY, mZ);
	}
#endif

	operator std::tuple<const FType &, const FType &, const FType &>() const
	{
		return std::make_tuple(std::ref(mX), std::ref(mY), std::ref(mZ));
	}

	operator std::tuple<FType &, FType &, FType &>()
	{
		return std::make_tuple(std::ref(mX), std::ref(mY), std::ref(mZ));
	}

	bool operator==(const PointF &rhs) const
	{
		return mX == rhs.mX and mY == rhs.mY and mZ == rhs.mZ;
	}

	// consider point as a vector... perhaps I should rename Point?
	FType lengthsq() const
	{
		return mX * mX + mY * mY + mZ * mZ;
	}

	FType length() const
	{
		return sqrt(mX * mX + mY * mY + mZ * mZ);
	}
};

typedef PointF<float> Point;
typedef PointF<double> DPoint;

template <typename F>
inline std::ostream &operator<<(std::ostream &os, const PointF<F> &pt)
{
	os << '(' << pt.mX << ',' << pt.mY << ',' << pt.mZ << ')';
	return os;
}

template <typename F>
inline PointF<F> operator+(const PointF<F> &lhs, const PointF<F> &rhs)
{
	return PointF<F>(lhs.mX + rhs.mX, lhs.mY + rhs.mY, lhs.mZ + rhs.mZ);
}

template <typename F>
inline PointF<F> operator-(const PointF<F> &lhs, const PointF<F> &rhs)
{
	return PointF<F>(lhs.mX - rhs.mX, lhs.mY - rhs.mY, lhs.mZ - rhs.mZ);
}

template <typename F>
inline PointF<F> operator-(const PointF<F> &pt)
{
	return PointF<F>(-pt.mX, -pt.mY, -pt.mZ);
}

template <typename F>
inline PointF<F> operator*(const PointF<F> &pt, F f)
{
	return PointF<F>(pt.mX * f, pt.mY * f, pt.mZ * f);
}

template <typename F>
inline PointF<F> operator*(F f, const PointF<F> &pt)
{
	return PointF<F>(pt.mX * f, pt.mY * f, pt.mZ * f);
}

template <typename F>
inline PointF<F> operator/(const PointF<F> &pt, F f)
{
	return PointF<F>(pt.mX / f, pt.mY / f, pt.mZ / f);
}

// --------------------------------------------------------------------
// several standard 3d operations

template <typename F>
inline double DistanceSquared(const PointF<F> &a, const PointF<F> &b)
{
	return (a.mX - b.mX) * (a.mX - b.mX) +
	       (a.mY - b.mY) * (a.mY - b.mY) +
	       (a.mZ - b.mZ) * (a.mZ - b.mZ);
}

template <typename F>
inline double Distance(const PointF<F> &a, const PointF<F> &b)
{
	return sqrt(
		(a.mX - b.mX) * (a.mX - b.mX) +
		(a.mY - b.mY) * (a.mY - b.mY) +
		(a.mZ - b.mZ) * (a.mZ - b.mZ));
}

template <typename F>
inline F DotProduct(const PointF<F> &a, const PointF<F> &b)
{
	return a.mX * b.mX + a.mY * b.mY + a.mZ * b.mZ;
}

template <typename F>
inline PointF<F> CrossProduct(const PointF<F> &a, const PointF<F> &b)
{
	return PointF<F>(a.mY * b.mZ - b.mY * a.mZ,
		a.mZ * b.mX - b.mZ * a.mX,
		a.mX * b.mY - b.mX * a.mY);
}

template <typename F>
double Angle(const PointF<F> &p1, const PointF<F> &p2, const PointF<F> &p3)
{
	PointF<F> v1 = p1 - p2;
	PointF<F> v2 = p3 - p2;

	return std::acos(DotProduct(v1, v2) / (v1.length() * v2.length())) * 180 / kPI;
}

template <typename F>
double DihedralAngle(const PointF<F> &p1, const PointF<F> &p2, const PointF<F> &p3, const PointF<F> &p4)
{
	PointF<F> v12 = p1 - p2; // vector from p2 to p1
	PointF<F> v43 = p4 - p3; // vector from p3 to p4

	PointF<F> z = p2 - p3; // vector from p3 to p2

	PointF<F> p = CrossProduct(z, v12);
	PointF<F> x = CrossProduct(z, v43);
	PointF<F> y = CrossProduct(z, x);

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

template <typename F>
double CosinusAngle(const PointF<F> &p1, const PointF<F> &p2, const PointF<F> &p3, const PointF<F> &p4)
{
	PointF<F> v12 = p1 - p2;
	PointF<F> v34 = p3 - p4;

	double result = 0;

	double x = DotProduct(v12, v12) * DotProduct(v34, v34);
	if (x > 0)
		result = DotProduct(v12, v34) / sqrt(x);

	return result;
}

template <typename F>
auto DistancePointToLine(const PointF<F> &l1, const PointF<F> &l2, const PointF<F> &p)
{
	auto line = l2 - l1;
	auto p_to_l1 = p - l1;
	auto p_to_l2 = p - l2;
	auto cross = CrossProduct(p_to_l1, p_to_l2);
	return cross.length() / line.length();
}

// --------------------------------------------------------------------
// For e.g. simulated annealing, returns a new point that is moved in
// a random direction with a distance randomly chosen from a normal
// distribution with a stddev of offset.

template <typename F>
PointF<F> Nudge(PointF<F> p, F offset);

// --------------------------------------------------------------------
// We use quaternions to do rotations in 3d space

Quaternion Normalize(Quaternion q);

std::tuple<double, Point> QuaternionToAngleAxis(Quaternion q);
Point Centroid(std::vector<Point> &Points);
Point CenterPoints(std::vector<Point> &Points);

/// \brief Returns how the two sets of points \a a and \b b can be aligned
///        
/// \param a	The first set of points
/// \param b    The second set of points
/// \result     The quaternion which should be applied to the points in \a a to
///             obtain the best superposition.
Quaternion AlignPoints(const std::vector<Point> &a, const std::vector<Point> &b);

/// \brief The RMSd for the points in \a a and \a b
double RMSd(const std::vector<Point> &a, const std::vector<Point> &b);

// --------------------------------------------------------------------
// Helper class to generate evenly divided Points on a sphere
// we use a fibonacci sphere to calculate even distribution of the dots

template <int N>
class SphericalDots
{
  public:
	enum
	{
		P = 2 * N + 1
	};
	typedef typename std::array<Point, P> array_type;
	typedef typename array_type::const_iterator iterator;

	static SphericalDots &instance()
	{
		static SphericalDots sInstance;
		return sInstance;
	}

	size_t size() const { return mPoints.size(); }
	const Point operator[](uint32_t inIx) const { return mPoints[inIx]; }
	iterator begin() const { return mPoints.begin(); }
	iterator end() const { return mPoints.end(); }

	double weight() const { return mWeight; }

	SphericalDots()
	{

		const double
			kGoldenRatio = (1 + std::sqrt(5.0)) / 2;

		mWeight = (4 * kPI) / P;

		auto p = mPoints.begin();

		for (int32_t i = -N; i <= N; ++i)
		{
			double lat = std::asin((2.0 * i) / P);
			double lon = std::fmod(i, kGoldenRatio) * 2 * kPI / kGoldenRatio;

			p->mX = sin(lon) * cos(lat);
			p->mY = cos(lon) * cos(lat);
			p->mZ = sin(lat);

			++p;
		}
	}

  private:
	array_type mPoints;
	double mWeight;
};

typedef SphericalDots<50> SphericalDots_50;

} // namespace mmcif
