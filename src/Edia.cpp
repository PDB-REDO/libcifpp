// copyright

#include "cif++/Config.h"

#include <boost/format.hpp>

#include "cif++/mrsrc.h"

#include "cif++/Edia.h"
#include "cif++/CifUtils.h"

using namespace std;

extern int VERBOSE;

namespace libcif
{

// --------------------------------------------------------------------
// Code to calculate EDIA (electron density for individual atoms)

AtomRadius::AtomRadius()
{
	mrsrc::rsrc curves("curves.xml");
	if (not curves)
		throw runtime_error("Missing curves.xml resource");
	
	mCurves.read(string(curves.data(), curves.size()));
}

float AtomRadius::operator()(AtomType a, int charge, float resolution)
{
	Key key(a, charge, resolution);

	try
	{
		if (mCache.count(key) == 0)
		{
			string symbol = AtomTypeTraits(a).symbol();
			
			auto fmt = boost::format("/curves/curve[@symbol='%1%']/cc[@charge='%2%']/radius[@resolution=%3%]");
	
			if (resolution <= 0.5 or resolution >= 3.0)
			{
				auto e = mCurves.find_first((fmt % symbol % charge % (resolution <= 0.5 ? 0.5 : 3.0)).str());
				if (not e)
					throw runtime_error("Missing radius value in curves.xml");
				
				mCache[key] = stof(e->str());
			}
			else
			{
				float r[2];
				
				     if (resolution < 1.0)		{	r[0] = 0.5;	r[1] = 1.0;	}
				else if (resolution < 1.5)		{	r[0] = 1.0;	r[1] = 1.5;	}
				else if (resolution < 2.0)		{	r[0] = 1.5;	r[1] = 2.0;	}
				else if (resolution < 2.5)		{	r[0] = 2.0;	r[1] = 2.5;	}
				else							{	r[0] = 2.5;	r[1] = 3.0;	}
				
				float c[2];
				
				auto e = mCurves.find_first((fmt % symbol % charge % r[0]).str());
				if (not e)
					throw runtime_error("Missing radius value in curves.xml");
				c[0] = stof(e->str());

				e = mCurves.find_first((fmt % symbol % charge % r[1]).str());
				if (not e)
					throw runtime_error("Missing radius value in curves.xml");
				c[1] = stof(e->str());
				
				mCache[key] = c[0] + ((c[1] - c[0]) * (resolution - r[0])) / (r[1] - r[0]);
			}
		}
		
		return mCache.at(key);
	}
	catch (const exception& ex)
	{
		stringstream s;
		s << "Error getting atom radius for " << AtomTypeTraits(a).symbol()
		  << " (charge " << charge << ") with resolution: " << resolution << endl
		  << " exception: " << ex.what() << endl;
		throw runtime_error(s.str());
	}
}

// --------------------------------------------------------------------

class PointWeightFunction
{
  public:
	PointWeightFunction(Point center, float atomRadius)
		: m_Center(center), m_Radius(atomRadius)
	{
		m_P[0] = P{ -1.0f, 0, 1.0f, 1.0822f };
		m_P[1] = P{ 5.1177f, 1.29366f, -0.4f, 1.4043f };
		m_P[2] = P{ -0.9507f, 2, 0, 2 };
	}

	
	float operator()(Point p) const
	{
		float d = Distance(m_Center, p);
		d /= m_Radius;
		
		float result = 0;

		for (int i = 0; i < 3; ++i)
		{
			auto& p = m_P[i];
			
			if (d > p.x)
				continue;
			
			result = p.m * (d - p.c) * (d - p.c) + p.b;
			
			assert(result != 0);
			if (result == 0)
				result = numeric_limits<float>::epsilon();
			
			break;
		}

		return result;
	}

  private:
	
	struct P
	{
		float	m, c, b, x;
	};
	
	Point		m_Center;
	float		m_Radius;
	P			m_P[3];
};

// --------------------------------------------------------------------

float CalculateEDIA(const Atom& atom, const clipper::Xmap<float>& xmap,
	float resolution, float meanDensity, float rmsDensity,
	const DistanceMap& dm, const BondMap& bm)
{
	// internal types for this function
	struct lessAtom
	{
		bool operator()(const Atom& a, const Atom& b) const { return a.id().compare(b.id()) < 0; }
	};
	
	typedef set<Atom, lessAtom> atomSet;

	
	static AtomRadius atomRadius;
	float radius = atomRadius(atom.type(), 0, resolution);
	
	float x, y, z;
	tie(x, y, z) = atom.location();

	// calculate min and max orthogonal coordinates first, based on atom position, radius and bfactor
	clipper::Coord_orth oMin = { x - radius * 2, y - radius * 2, z - radius * 2 },
						oMax = { x + radius * 2, y + radius * 2, z + radius * 2 };

	clipper::Coord_frac fMin = oMin.coord_frac(xmap.cell());
	clipper::Coord_frac fMax = oMax.coord_frac(xmap.cell());

	clipper::Coord_map mMin = fMin.coord_map(xmap.grid_sampling());
	clipper::Coord_map mMax = fMax.coord_map(xmap.grid_sampling());

	clipper::Coord_grid gridMin = mMin.floor();
	clipper::Coord_grid gridMax = mMax.ceil();

	PointWeightFunction w(atom.location(), radius);
	
	vector<Atom> atomsNearBy = dm.near(atom, 3.5f);

	vector<PointWeightFunction> wn;
	for (auto a: atomsNearBy)
		wn.emplace_back(a.location(), atomRadius(a, resolution));

	double sum1 = 0, sum2 = 0;

	auto i0 = clipper::Xmap_base::Map_reference_coord(xmap, gridMin);
	for (auto iu = i0; iu.coord().u() <= gridMax[0]; iu.next_u())
		for (auto iv = iu; iv.coord().v() <= gridMax[1]; iv.next_v())
			for (auto iw = iv; iw.coord().w() <= gridMax[2]; iw.next_w())
			{
				auto p = iw.coord_orth();

				float z = (xmap[iw] - meanDensity) / rmsDensity;
				
				if (z < 0)
					z = 0;
				if (z > 1.2)
					z = 1.2;
				
				float wp = w(p);
				
				// And divide the ownership
				
				atomSet S, D, I;
				
				if (wp != 0)
				{
					if (wp < 0)
						D.insert(atom);
					else
					{
						S.insert(atom);
						I.insert(atom);
					}
				}
				
				for (size_t i = 0; i < atomsNearBy.size(); ++i)
				{
					float w = wn[i](p);
					if (w == 0)
						continue;
					
					if (w < 0)
						D.insert(atomsNearBy[i]);
					else if (w > 0)
					{
						S.insert(atomsNearBy[i]);
						
						if (not bm(atomsNearBy[i], atom))
							I.insert(atomsNearBy[i]);
					}
				}
				
				float o = 0;
				if (wp > 0)
				{
					if (I.size() == 1)
						o = 1;
					else
					{
						float sumpb = accumulate(I.begin(), I.end(), 0.f,
							[p](float s, const Atom& b) -> float
							{
								return s + Distance(p, b.location());
							});

						o = 1 - Distance(atom.location(), p) / sumpb;
					}
				}
				else if (D.count(atom) and S.empty())
				{
					if (D.size() == 1)
						o = 1;
					else
					{
						float sumpb = accumulate(D.begin(), D.end(), 0.f,
							[p](float s, const Atom& b) -> float
							{
								return s + Distance(p, b.location());
							});
	
						o = 1 - Distance(atom.location(), p) / sumpb;
					}
				}

				if (VERBOSE > 2)
					cout << Point(p) << ":\td: " << xmap[iw] << "\tz: " << z << "\to: " << o << "\tzraw: " << ((xmap[iw] - meanDensity) / rmsDensity) << "\twp: " << wp << endl;
				
				sum1 += z * wp * o;
				if (wp > 0)
					sum2 += wp;
			}
	
	float result = sum1 / sum2;
	
	if (VERBOSE > 1)
		cout << string(cif::get_terminal_width(), '-') << endl
			 << "Processing atom " << atom.id() << endl
			 << "Symbol: " << AtomTypeTraits(atom.type()).symbol() << endl
			 << "Location: " << atom.location() << endl
			 << "Radius: " << radius << endl
			 << "EDIA: " << result << endl;
	
	return result;
}

}
