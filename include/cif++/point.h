// Lib for working with structures as contained in mmCIF and PDB files

#pragma once

#include <libcif/config.h>

#include <boost/filesystem/operations.hpp>
#include <boost/math/quaternion.hpp>

#include "clipper/core/coords.h"

namespace libcif
{

typedef boost::math::quaternion<float>	quaternion;

const long double
	kPI = 3.141592653589793238462643383279502884L;

// --------------------------------------------------------------------

//	point, a location with x, y and z coordinates as float.
//	This one is derived from a tuple<float,float,float> so
//	you can do things like:
//
//	float x, y, z;
//	tie(x, y, z) = atom.loc();

struct point : public std::tuple<float,float,float>
{
	typedef std::tuple<float,float,float>	base_type;
	
	point()								: base_type(0.f, 0.f, 0.f) {}
	point(float x, float y, float z)	: base_type(x, y, z) {}
	point(const clipper::Coord_orth& pt): base_type(pt[0], pt[1], pt[2]) {}
	
	point& operator=(const clipper::Coord_orth& rhs)
	{
		x(rhs[0]);
		y(rhs[1]);
		z(rhs[2]);
		return *this;
	}
	
	float& x()				{ return std::get<0>(*this); }
	float x() const			{ return std::get<0>(*this); }
	void x(float x)			{ std::get<0>(*this) = x; }

	float& y()				{ return std::get<1>(*this); }
	float y() const			{ return std::get<1>(*this); }
	void y(float y)			{ std::get<1>(*this) = y; }

	float& z()				{ return std::get<2>(*this); }
	float z() const			{ return std::get<2>(*this); }
	void z(float z)			{ std::get<2>(*this) = z; }
	
	point& operator+=(const point& rhs)
	{
		x() += rhs.x();
		y() += rhs.y();
		z() += rhs.z();
		return *this;
	}
	
	point& operator-=(const point& rhs)
	{
		x() -= rhs.x();
		y() -= rhs.y();
		z() -= rhs.z();
		return *this;
	}

	point& operator*=(float rhs)
	{
		x() *= rhs;
		y() *= rhs;
		z() *= rhs;
		return *this;
	}
	
	point& operator/=(float rhs)
	{
		x() *= rhs;
		y() *= rhs;
		z() *= rhs;
		return *this;
	}

	float normalize()
	{
		auto length = x() * x() + y() * y() + z() * z();
		if (length > 0)
		{
			length = std::sqrt(length);
			operator/=(length);
		}
		return length;
	}
	
	void rotate(const boost::math::quaternion<float>& q)
	{
		boost::math::quaternion<float> p(0, x(), y(), z());
		
		p = q * p * boost::math::conj(q);
	
		x() = p.R_component_2();
		y() = p.R_component_3();
		z() = p.R_component_4();
	}
	
	operator clipper::Coord_orth() const
	{
		return clipper::Coord_orth(x(), y(), z());
	}
};


inline std::ostream& operator<<(std::ostream& os, const point& pt)
{
	os << '(' << pt.x() << ',' << pt.y() << ',' << pt.z() << ')';
	return os; 
}

inline point operator+(const point& lhs, const point& rhs)
{
	return point(lhs.x() + rhs.x(), lhs.y() + rhs.y(), lhs.z() + rhs.z());
}

inline point operator-(const point& lhs, const point& rhs)
{
	return point(lhs.x() - rhs.x(), lhs.y() - rhs.y(), lhs.z() - rhs.z());
}

inline point operator-(const point& pt)
{
	return point(-pt.x(), -pt.y(), -pt.z());
}

inline point operator*(const point& pt, float f)
{
	return point(pt.x() * f, pt.y() * f, pt.z() * f);
}

inline point operator/(const point& pt, float f)
{
	return point(pt.x() / f, pt.y() / f, pt.z() / f);
}

// --------------------------------------------------------------------
// several standard 3d operations

inline double DistanceSquared(const point& a, const point& b)
{
	return
		(a.x() - b.x()) * (a.x() - b.x()) +
		(a.y() - b.y()) * (a.y() - b.y()) +
		(a.z() - b.z()) * (a.z() - b.z());
}

inline double Distance(const point& a, const point& b)
{
	return sqrt(
		(a.x() - b.x()) * (a.x() - b.x()) +
		(a.y() - b.y()) * (a.y() - b.y()) +
		(a.z() - b.z()) * (a.z() - b.z()));
}

inline float DotProduct(const point& a, const point& b)
{
	return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
}

inline point CrossProduct(const point& a, const point& b)
{
	return point(a.y() * b.z() - b.y() * a.z(),
				  a.z() * b.x() - b.z() * a.x(),
				  a.x() * b.y() - b.x() * a.y());
}

float DihedralAngle(const point& p1, const point& p2, const point& p3, const point& p4);
float CosinusAngle(const point& p1, const point& p2, const point& p3, const point& p4);

// --------------------------------------------------------------------
// We use quaternions to do rotations in 3d space

quaternion Normalize(quaternion q);

//std::tuple<double,point> QuaternionToAngleAxis(quaternion q);
point Centroid(std::vector<point>& points);
point CenterPoints(std::vector<point>& points);
quaternion AlignPoints(const std::vector<point>& a, const std::vector<point>& b);
double RMSd(const std::vector<point>& a, const std::vector<point>& b);

// --------------------------------------------------------------------
// Helper class to generate evenly divided points on a sphere
// we use a fibonacci sphere to calculate even distribution of the dots

template<int N>
class spherical_dots
{
  public:
	enum { P = 2 * N + 1 };
	typedef typename std::array<point,P>	array_type;
	typedef typename array_type::const_iterator	iterator;

	static spherical_dots& instance()
	{
		static spherical_dots s_instance;
		return s_instance;
	}
	
	size_t size() const							{ return m_points.size(); }
	const point operator[](uint32 inIx) const	{ return m_points[inIx]; }
	iterator begin() const						{ return m_points.begin(); }
	iterator end() const						{ return m_points.end(); }

	double weight() const						{ return m_weight; }

	spherical_dots()
	{
		using namespace std;
		
		const double
			kGoldenRatio = (1 + std::sqrt(5.0)) / 2;
		
		m_weight = (4 * kPI) / P;
		
		auto p = m_points.begin();
		
		for (int32 i = -N; i <= N; ++i)
		{
			double lat = std::asin((2.0 * i) / P);
			double lon = std::fmod(i, kGoldenRatio) * 2 * kPI / kGoldenRatio;
			
			p->x(sin(lon) * cos(lat));
			p->y(cos(lon) * cos(lat));
			p->z(           sin(lat));

			++p;
		}
	}

  private:

	array_type				m_points;
	double					m_weight;
};

typedef spherical_dots<50> spherical_dots_50;


}
