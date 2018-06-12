// Lib for working with structures as contained in mmCIF and PDB files

#pragma once

#include "cif++/Config.h"

#include <boost/filesystem/operations.hpp>
#include <boost/math/quaternion.hpp>

#include "clipper/core/coords.h"

namespace mmcif
{

typedef boost::math::quaternion<float>	quaternion;

const long double
	kPI = 3.141592653589793238462643383279502884L;

// --------------------------------------------------------------------

//	Point, a location with x, y and z coordinates as float.
//	This one is derived from a tuple<float,float,float> so
//	you can do things like:
//
//	float x, y, z;
//	tie(x, y, z) = atom.loc();

struct Point
{
	float mX, mY, mZ;
	
	Point()								: mX(0), mY(0), mZ(0) {}
	Point(float x, float y, float z)	: mX(x), mY(y), mZ(z) {}
	Point(const clipper::Coord_orth& pt): mX(pt[0]), mY(pt[1]), mZ(pt[2]) {}
	
	Point& operator=(const clipper::Coord_orth& rhs)
	{
		mX = rhs[0];
		mY = rhs[1];
		mZ = rhs[2];
		return *this;
	}
	
	float& getX()			{ return mX; }
	float getX() const		{ return mX; }
	void setX(float x)		{ mX = x; }

	float& getY()			{ return mY; }
	float getY() const		{ return mY; }
	void setY(float y)		{ mY = y; }

	float& getZ()			{ return mZ; }
	float getZ() const		{ return mZ; }
	void setZ(float z)		{ mZ = z; }
	
	Point& operator+=(const Point& rhs)
	{
		mX += rhs.mX;
		mY += rhs.mY;
		mZ += rhs.mZ;
		
		return *this;
	}
	
	Point& operator+=(float d)
	{
		mX += d;
		mY += d;
		mZ += d;
		
		return *this;
	}

	Point& operator-=(const Point& rhs)
	{
		mX -= rhs.mX;
		mY -= rhs.mY;
		mZ -= rhs.mZ;
		
		return *this;
	}

	Point& operator-=(float d)
	{
		mX -= d;
		mY -= d;
		mZ -= d;
		
		return *this;
	}

	Point& operator*=(float rhs)
	{
		mX *= rhs;
		mY *= rhs;
		mZ *= rhs;
		return *this;
	}
	
	Point& operator/=(float rhs)
	{
		mX /= rhs;
		mY /= rhs;
		mZ /= rhs;
		return *this;
	}

	float normalize()
	{
		auto length = mX * mX + mY * mY + mZ * mZ;
		if (length > 0)
		{
			length = std::sqrt(length);
			operator/=(length);
		}
		return length;
	}
	
	void rotate(const boost::math::quaternion<float>& q)
	{
		boost::math::quaternion<float> p(0, mX, mY, mZ);
		
		p = q * p * boost::math::conj(q);
	
		mX = p.R_component_2();
		mY = p.R_component_3();
		mZ = p.R_component_4();
	}
	
	operator clipper::Coord_orth() const
	{
		return clipper::Coord_orth(mX, mY, mZ);
	}
	
	operator std::tuple<const float&, const float&, const float&>() const
	{
		return std::make_tuple(std::ref(mX), std::ref(mY), std::ref(mZ));
	}

	operator std::tuple<float&,float&,float&>()
	{
		return std::make_tuple(std::ref(mX), std::ref(mY), std::ref(mZ));
	}
	
	bool operator==(const Point& rhs) const
	{
		return mX == rhs.mX and mY == rhs.mY and mZ == rhs.mZ;
	}
	
	// consider point as a vector... perhaps I should rename Point?
	float lengthsq() const
	{
		return mX * mX + mY * mY + mZ * mZ;
	}

	float length() const
	{
		return sqrt(mX * mX + mY * mY + mZ * mZ);
	}
};

inline std::ostream& operator<<(std::ostream& os, const Point& pt)
{
	os << '(' << pt.mX << ',' << pt.mY << ',' << pt.mZ << ')';
	return os; 
}

inline Point operator+(const Point& lhs, const Point& rhs)
{
	return Point(lhs.mX + rhs.mX, lhs.mY + rhs.mY, lhs.mZ + rhs.mZ);
}

inline Point operator-(const Point& lhs, const Point& rhs)
{
	return Point(lhs.mX - rhs.mX, lhs.mY - rhs.mY, lhs.mZ - rhs.mZ);
}

inline Point operator-(const Point& pt)
{
	return Point(-pt.mX, -pt.mY, -pt.mZ);
}

inline Point operator*(const Point& pt, float f)
{
	return Point(pt.mX * f, pt.mY * f, pt.mZ * f);
}

inline Point operator*(float f, const Point& pt)
{
	return Point(pt.mX * f, pt.mY * f, pt.mZ * f);
}

inline Point operator/(const Point& pt, float f)
{
	return Point(pt.mX / f, pt.mY / f, pt.mZ / f);
}

// --------------------------------------------------------------------
// several standard 3d operations

inline double DistanceSquared(const Point& a, const Point& b)
{
	return
		(a.mX - b.mX) * (a.mX - b.mX) +
		(a.mY - b.mY) * (a.mY - b.mY) +
		(a.mZ - b.mZ) * (a.mZ - b.mZ);
}

inline double Distance(const Point& a, const Point& b)
{
	return sqrt(
		(a.mX - b.mX) * (a.mX - b.mX) +
		(a.mY - b.mY) * (a.mY - b.mY) +
		(a.mZ - b.mZ) * (a.mZ - b.mZ));
}

inline float DotProduct(const Point& a, const Point& b)
{
	return a.mX * b.mX + a.mY * b.mY + a.mZ * b.mZ;
}

inline Point CrossProduct(const Point& a, const Point& b)
{
	return Point(a.mY * b.mZ - b.mY * a.mZ,
				  a.mZ * b.mX - b.mZ * a.mX,
				  a.mX * b.mY - b.mX * a.mY);
}

float DihedralAngle(const Point& p1, const Point& p2, const Point& p3, const Point& p4);
float CosinusAngle(const Point& p1, const Point& p2, const Point& p3, const Point& p4);

// --------------------------------------------------------------------
// For e.g. simulated annealing, returns a new point that is moved in
// a random direction with a distance randomly chosen from a normal
// distribution with a stddev of offset.

Point Nudge(Point p, float offset);

// --------------------------------------------------------------------
// We use quaternions to do rotations in 3d space

quaternion Normalize(quaternion q);

//std::tuple<double,Point> QuaternionToAngleAxis(quaternion q);
Point Centroid(std::vector<Point>& Points);
Point CenterPoints(std::vector<Point>& Points);
quaternion AlignPoints(const std::vector<Point>& a, const std::vector<Point>& b);
double RMSd(const std::vector<Point>& a, const std::vector<Point>& b);

// --------------------------------------------------------------------
// Helper class to generate evenly divided Points on a sphere
// we use a fibonacci sphere to calculate even distribution of the dots

template<int N>
class SphericalDots
{
  public:
	enum { P = 2 * N + 1 };
	typedef typename std::array<Point,P>	array_type;
	typedef typename array_type::const_iterator	iterator;

	static SphericalDots& instance()
	{
		static SphericalDots sInstance;
		return sInstance;
	}
	
	size_t size() const							{ return mPoints.size(); }
	const Point operator[](uint32 inIx) const	{ return mPoints[inIx]; }
	iterator begin() const						{ return mPoints.begin(); }
	iterator end() const						{ return mPoints.end(); }

	double weight() const						{ return mWeight; }

	SphericalDots()
	{
		using namespace std;
		
		const double
			kGoldenRatio = (1 + std::sqrt(5.0)) / 2;
		
		mWeight = (4 * kPI) / P;
		
		auto p = mPoints.begin();
		
		for (int32 i = -N; i <= N; ++i)
		{
			double lat = std::asin((2.0 * i) / P);
			double lon = std::fmod(i, kGoldenRatio) * 2 * kPI / kGoldenRatio;
			
			p->mX = sin(lon) * cos(lat);
			p->mY = cos(lon) * cos(lat);
			p->mZ =            sin(lat);

			++p;
		}
	}

  private:

	array_type				mPoints;
	double					mWeight;
};

typedef SphericalDots<50> SphericalDots_50;

}
