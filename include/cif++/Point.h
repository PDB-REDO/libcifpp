// Lib for working with structures as contained in mmCIF and PDB files

#pragma once

#include "cif++/Config.h"

#include <boost/filesystem/operations.hpp>
#include <boost/math/quaternion.hpp>

#include "clipper/core/coords.h"

namespace libcif
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

struct Point : public std::tuple<float,float,float>
{
	typedef std::tuple<float,float,float>	base_type;
	
	Point()								: base_type(0.f, 0.f, 0.f) {}
	Point(float x, float y, float z)	: base_type(x, y, z) {}
	Point(const clipper::Coord_orth& pt): base_type(pt[0], pt[1], pt[2]) {}
	
	Point& operator=(const clipper::Coord_orth& rhs)
	{
		setX(rhs[0]);
		setY(rhs[1]);
		setZ(rhs[2]);
		return *this;
	}
	
	float& getX()			{ return std::get<0>(*this); }
	float getX() const		{ return std::get<0>(*this); }
	void setX(float x)		{ std::get<0>(*this) = x; }

	float& getY()			{ return std::get<1>(*this); }
	float getY() const		{ return std::get<1>(*this); }
	void setY(float y)		{ std::get<1>(*this) = y; }

	float& getZ()			{ return std::get<2>(*this); }
	float getZ() const		{ return std::get<2>(*this); }
	void setZ(float z)		{ std::get<2>(*this) = z; }
	
	Point& operator+=(const Point& rhs)
	{
		getX() += rhs.getX();
		getY() += rhs.getY();
		getZ() += rhs.getZ();
		return *this;
	}
	
	Point& operator-=(const Point& rhs)
	{
		getX() -= rhs.getX();
		getY() -= rhs.getY();
		getZ() -= rhs.getZ();
		return *this;
	}

	Point& operator*=(float rhs)
	{
		getX() *= rhs;
		getY() *= rhs;
		getZ() *= rhs;
		return *this;
	}
	
	Point& operator/=(float rhs)
	{
		getX() *= rhs;
		getY() *= rhs;
		getZ() *= rhs;
		return *this;
	}

	float normalize()
	{
		auto length = getX() * getX() + getY() * getY() + getZ() * getZ();
		if (length > 0)
		{
			length = std::sqrt(length);
			operator/=(length);
		}
		return length;
	}
	
	void rotate(const boost::math::quaternion<float>& q)
	{
		boost::math::quaternion<float> p(0, getX(), getY(), getZ());
		
		p = q * p * boost::math::conj(q);
	
		getX() = p.R_component_2();
		getY() = p.R_component_3();
		getZ() = p.R_component_4();
	}
	
	operator clipper::Coord_orth() const
	{
		return clipper::Coord_orth(getX(), getY(), getZ());
	}
};

inline std::ostream& operator<<(std::ostream& os, const Point& pt)
{
	os << '(' << pt.getX() << ',' << pt.getY() << ',' << pt.getZ() << ')';
	return os; 
}

inline Point operator+(const Point& lhs, const Point& rhs)
{
	return Point(lhs.getX() + rhs.getX(), lhs.getY() + rhs.getY(), lhs.getZ() + rhs.getZ());
}

inline Point operator-(const Point& lhs, const Point& rhs)
{
	return Point(lhs.getX() - rhs.getX(), lhs.getY() - rhs.getY(), lhs.getZ() - rhs.getZ());
}

inline Point operator-(const Point& pt)
{
	return Point(-pt.getX(), -pt.getY(), -pt.getZ());
}

inline Point operator*(const Point& pt, float f)
{
	return Point(pt.getX() * f, pt.getY() * f, pt.getZ() * f);
}

inline Point operator/(const Point& pt, float f)
{
	return Point(pt.getX() / f, pt.getY() / f, pt.getZ() / f);
}

// --------------------------------------------------------------------
// several standard 3d operations

inline double DistanceSquared(const Point& a, const Point& b)
{
	return
		(a.getX() - b.getX()) * (a.getX() - b.getX()) +
		(a.getY() - b.getY()) * (a.getY() - b.getY()) +
		(a.getZ() - b.getZ()) * (a.getZ() - b.getZ());
}

inline double Distance(const Point& a, const Point& b)
{
	return sqrt(
		(a.getX() - b.getX()) * (a.getX() - b.getX()) +
		(a.getY() - b.getY()) * (a.getY() - b.getY()) +
		(a.getZ() - b.getZ()) * (a.getZ() - b.getZ()));
}

inline float DotProduct(const Point& a, const Point& b)
{
	return a.getX() * b.getX() + a.getY() * b.getY() + a.getZ() * b.getZ();
}

inline Point CrossProduct(const Point& a, const Point& b)
{
	return Point(a.getY() * b.getZ() - b.getY() * a.getZ(),
				  a.getZ() * b.getX() - b.getZ() * a.getX(),
				  a.getX() * b.getY() - b.getX() * a.getY());
}

float DihedralAngle(const Point& p1, const Point& p2, const Point& p3, const Point& p4);
float CosinusAngle(const Point& p1, const Point& p2, const Point& p3, const Point& p4);

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
			
			p->setX(sin(lon) * cos(lat));
			p->setY(cos(lon) * cos(lat));
			p->setZ(           sin(lat));

			++p;
		}
	}

  private:

	array_type				mPoints;
	double					mWeight;
};

typedef SphericalDots<50> SphericalDots_50;


}
