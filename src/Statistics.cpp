/* 
   Created by: Maarten L. Hekkelman
   Date: maandag 04 juni, 2018
*/

#include "cif++/Config.h"

#include <fstream>
#include <numeric>

#include "cif++/Structure.h"
#include "cif++/AtomShape.h"
#include "cif++/Statistics.h"

using namespace std;

namespace mmcif
{

using clipper::Coord_grid;
using clipper::Coord_orth;
using clipper::Xmap;

// --------------------------------------------------------------------

ostream& operator<<(ostream& os, const ResidueStatistics& st)
{
	if (st.compID == "HOH")
		os << st.asymID << '_' << st.authSeqID << '_' << st.compID << '\t';
	else
		os << st.asymID << '_' << st.seqID << '_' << st.compID << '\t';
	 
	os << st.RSR << '\t'
	   << st.SRSR << '\t'
	   << st.RSCCS << '\t'
	   << st.ngrid << '\t'
	   << st.EDIAm << '\t'
	   << st.OPIA;
	 
	 return os;
}

// --------------------------------------------------------------------

double anorm(double x)
{
	return 0.5 * erfc(-x * sqrt(0.5));
}

double phinvs(double p)
{
//
// ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3.
//
// Produces the normal deviate Z corresponding to a given lower tail
// area of P; Z is accurate to about 1 part in 10**16.

	// Coefficients for P close to 0.5
	const double A[8] = {
		3.3871328727963666080,		1.3314166789178437745e+2,	1.9715909503065514427e+3,	1.3731693765509461125e+4,
		4.5921953931549871457e+4,	6.7265770927008700853e+4,	3.3430575583588128105e+4,	2.5090809287301226727e+3
	}, B[8] = {
		0,							4.2313330701600911252e+1,	6.8718700749205790830e+2,	5.3941960214247511077e+3,
		2.1213794301586595867e+4,	3.9307895800092710610e+4,	2.8729085735721942674e+4,	5.2264952788528545610e+3
	};
	
	// Coefficients for P not close to 0, 0.5 or 1.
	const double C[8] = {
		1.42343711074968357734e0,	4.63033784615654529590e0,	5.76949722146069140550e0,	3.64784832476320460504e0,
		1.27045825245236838258e0,	2.41780725177450611770e-1,	2.27238449892691845833e-2,	7.74545014278341407640e-4,
	}, D[8] = {
		0,							2.05319162663775882187e0,	1.67638483018380384940e0,	6.89767334985100004550e-1,
		1.48103976427480074590e-1,	1.51986665636164571966e-2,	5.47593808499534494600e-4,	1.05075007164441684324e-9
	};
	
	// Coefficients for P near 0 or 1.
	const double E[8] = {
		6.65790464350110377720e0,	5.46378491116411436990e0,	1.78482653991729133580e0,	2.96560571828504891230e-1,
		2.65321895265761230930e-2,	1.24266094738807843860e-3,	2.71155556874348757815e-5,	2.01033439929228813265e-7,
	}, F[8] = {
		0,							5.99832206555887937690e-1,	1.36929880922735805310e-1,	1.48753612908506148525e-2,
		7.86869131145613259100e-4,	1.84631831751005468180e-5,	1.42151175831644588870e-7,	2.04426310338993978564e-15
	};

	if (p < 0 or p > 1)
		throw runtime_error("P should be >=0 and <=1");

	double q = p - 0.5;
	double result;
	
	if (abs(q) < 0.425)
	{
        double r = 0.180625e0 - q * q;
        result =
        	q * (((((((A[7] * r + A[6]) * r + A[5]) * r + A[4]) * r + A[3]) * r + A[2]) * r + A[1]) * r + A[0])
        		/ (((((((B[7] * r + B[6]) * r + B[5]) * r + B[4]) * r + B[3]) * r + B[2]) * r + B[1]) * r + 1);
	}
	else
	{
		double r;
		
		if (q < 0)
			r = p;
		else
			r = 1 - p;
        
        r = sqrt(-log(r));

        if (r <= 5)
        {
        	r -= 1.6;
        	result = (((((((C[7] * r + C[6]) * r + C[5]) * r + C[4]) * r + C[3]) * r + C[2]) * r + C[1]) * r + C[0])
        		/ (((((((D[7] * r + D[6]) * r + D[5]) * r + D[4]) * r + D[3]) * r + D[2]) * r + D[1]) * r + 1);
        }
        else
        {
        	r -= 0.5;
        	result = (((((((E[7] * r + E[6]) * r + E[5]) * r + E[4]) * r + E[3]) * r + E[2]) * r + E[1]) * r + E[0])
        		/ (((((((F[7] * r + F[6]) * r + F[5]) * r + F[4]) * r + F[3]) * r + F[2]) * r + F[1]) * r + 1);
        }
        
        if (q < 0)
        	result = -result;
	}
	return result;
}


double errsol(double a)
{
	auto c = sqrt(2.0 / kPI);
	auto b = abs(a);
	
	double result = 0;
	if (b > 3 / c)
	{
		auto x = abs(pow(b, 1/3.0) - 2 * pow(kPI / b, 2));
		if (a < 0)
			x = -x;
		
		for (;;)
		{
			auto xx = x * x;
			auto y = c * exp(-0.5 * xx);
			auto d = (b * (2 * anorm(x) - 1 - x * y) / xx - x) / (b * y - 3);
			x -= d;
			
			if (abs(d) <= 1e-4)
				break;
		}
		
		result = x;
	}
	
	return result;
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

		for (auto& p: m_P)
		{
			if (d > p.x)
				continue;
			
			result = p.m * (d - p.c) * (d - p.c) + p.b;
			
//			assert(result != 0);
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

template<typename F>
void iterateGrid(const Coord_orth& p, float r, const Xmap<float>& m, F&& func)
{
	using namespace clipper;
	
	Coord_frac fp = p.coord_frac(m.cell());

	Coord_frac o = Coord_orth(r, r, r).coord_frac(m.cell());
	o[0] = abs(o[0]);
	o[1] = abs(o[1]);
	o[2] = abs(o[2]);

	Coord_frac fMin = fp - o, fMax = fp + o;
	Coord_map mMin = fMin.coord_map(m.grid_sampling()), mMax = fMax.coord_map(m.grid_sampling());
	Coord_grid gMin = mMin.floor(), gMax = mMax.ceil();

	auto i0 = Xmap_base::Map_reference_coord(m, gMin);
	for (auto iu = i0; iu.coord().u() <= gMax[0]; iu.next_u())
		for (auto iv = iu; iv.coord().v() <= gMax[1]; iv.next_v())
			for (auto iw = iv; iw.coord().w() <= gMax[2]; iw.next_w())
				func(iw);
}

// --------------------------------------------------------------------

struct AtomGridData
{
	AtomGridData(const Coord_grid& gp, double density)
		: p(gp), density(density) {}

	Coord_grid		p;
	double			density;
};

struct AtomDataSums
{
	size_t					ngrid = 0;
	double					rfSums[2] = {};		// sums for R-Factor
	double					edSums[2] = {};		// Sums for ED1 and ED3
	double					ccSums[3] = {};		// Sums for CC calculation
	double					rgSums[2] = {};
	double					swSums[3] = {};		// Sums used for sample CC calculation
	
	AtomDataSums& operator+=(const AtomDataSums& rhs)
	{
		ngrid += rhs.ngrid;
		rfSums[0] += rhs.rfSums[0];
		rfSums[1] += rhs.rfSums[1];
		edSums[0] += rhs.edSums[0];
		edSums[1] += rhs.edSums[1];
		ccSums[0] += rhs.ccSums[0];
		ccSums[1] += rhs.ccSums[1];
		ccSums[2] += rhs.ccSums[2];
		rgSums[0] += rhs.rgSums[0];
		rgSums[1] += rhs.rgSums[1];
		swSums[0] += rhs.swSums[0];
		swSums[1] += rhs.swSums[1];
		swSums[2] += rhs.swSums[2];
		return *this;
	}
	
	double cc() const
	{
		double s = (ccSums[1] - (edSums[0] * edSums[0]) / ngrid) * (ccSums[2] - (edSums[1] * edSums[1]) / ngrid);
		return (ccSums[0] - edSums[0] * edSums[1] / ngrid) / sqrt(s);
	}

	double srg() const
	{
		double rg = sqrt(rgSums[0] / rgSums[1]);
		
		return sqrt(swSums[0] - rg * rg * swSums[1] + 0.5 * rg * rg * rg * rg * swSums[2]) / (rg * rgSums[1]);
	}
};

struct AtomData
{
	AtomData(Atom atom, float radius)
		: atom(atom)
//		, asymID(atom.authAsymID())
//		, seqID(atom.property<string>("auth_seq_id"))
		, asymID(atom.labelAsymID())
		, seqID(atom.labelSeqID())
		, radius(radius) {}
	
	Atom					atom;
	string					asymID;
	int						seqID;
	float					radius;
	vector<AtomGridData>	points;
	double					averageDensity = 0;
	double					edia = 0;
	AtomDataSums			sums;
};

// --------------------------------------------------------------------

tuple<float,float> CalculateMapStatistics(const Xmap<float>& f)
{
	double sum = 0, sum2 = 0;
	int count = 0;
	
	for (auto ix = f.first(); not ix.last(); ix.next())
	{
		auto v = f[ix];
		
		if (isnan(v))
			throw runtime_error("map contains NaN values");
		
		++count;
		sum += v;
		sum2 += v * v;
	}
	
	float meanDensity = static_cast<float>(sum / count);
	float rmsDensity = static_cast<float>(sqrt((sum2 / count) - (meanDensity * meanDensity)));
	
	return make_tuple(meanDensity, rmsDensity);
}

// --------------------------------------------------------------------

class BoundingBox
{
  public:
//	BoundingBox(const Structure& structure, const vector<tuple<string,int,string,string>>& residues, float margin)
//	{
//		mXMin = mYMin = mZMin = numeric_limits<float>::max();
//		mXMax = mYMax = mZMax = numeric_limits<float>::min();
//		
//		for (auto& r: residues)
//		{
//			int seqID;
//			string asymID, compID, pdbID;
//			tie(asymID, seqID, compID, pdbID) = r;
//			
//			Residue res(structure, compID, asymID, seqID);
//			for (auto& atom: res.atoms())
//			{
//				auto l = atom.location();
//				if (mXMin > l.mX)
//					mXMin = l.mX;
//				if (mXMax < l.mX)
//					mXMax = l.mX;
//				if (mYMin > l.mY)
//					mYMin = l.mY;
//				if (mYMax < l.mY)
//					mYMax = l.mY;
//				if (mZMin > l.mZ)
//					mZMin = l.mZ;
//				if (mZMax < l.mZ)
//					mZMax = l.mZ;
//			}
//		}
//		
//		mXMin -= margin;
//		mXMax += margin;
//		mYMin -= margin;
//		mYMax += margin;
//		mZMin -= margin;
//		mZMax += margin;
//	}

	template<class List>
	BoundingBox(const Structure& structure, List atoms, float margin)
	{
		mXMin = mYMin = mZMin = numeric_limits<float>::max();
		mXMax = mYMax = mZMax = numeric_limits<float>::min();
		
		for (auto& atom: atoms)
		{
			auto l = atom.location();
			if (mXMin > l.mX)
				mXMin = l.mX;
			if (mXMax < l.mX)
				mXMax = l.mX;
			if (mYMin > l.mY)
				mYMin = l.mY;
			if (mYMax < l.mY)
				mYMax = l.mY;
			if (mZMin > l.mZ)
				mZMin = l.mZ;
			if (mZMax < l.mZ)
				mZMax = l.mZ;
		}
		
		mXMin -= margin;
		mXMax += margin;
		mYMin -= margin;
		mYMax += margin;
		mZMin -= margin;
		mZMax += margin;
	}
	
	bool contains(const Point& p) const
	{
		return  p.mX >= mXMin and p.mX <= mXMax
			and p.mY >= mYMin and p.mY <= mYMax
			and p.mZ >= mZMin and p.mZ <= mZMax;
	}

  private:
	float mXMin, mXMax, mYMin, mYMax, mZMin, mZMax;
};
// --------------------------------------------------------------------

StatsCollector::StatsCollector(const MapMaker<float>& mm, Structure& structure, bool electronScattering)
	: mStructure(structure), mMapMaker(mm), mElectronScattering(electronScattering)
{
	mSpacegroup = mm.spacegroup();
	mCell = mm.cell();
	mGrid = mm.gridSampling();
	mResHigh = mm.resHigh();
	mResLow = mm.resLow();
	
	initialize();
}

void StatsCollector::initialize()
{
	mMeanDensityFb = mMapMaker.fb().meanDensity();
	mRMSDensityFb = mMapMaker.fb().rmsDensity();
	mRMSDensityFd = mMapMaker.fd().rmsDensity();

	// calculate degrees of freedom
	auto omcd = mCell.matrix_orth();

	mVF = 1;
	mVC = 1;
	
	for (int i = 0; i < 3; ++i)
	{
		mVC *= omcd(i, i);
		mVF *= omcd(i, i) / mGrid[i];
	}	

	mVF *= pow(2 / mResHigh, 3);

	mSZ = 0;
//	double so = 0;
//	const double C = sqrt(2.0 / kPI);
	
	for (auto& a: mStructure.atoms())
	{
		auto t = a.type();
		if (t <= He)
			continue;
		
		float w = a.occupancy() * t;
		
		if (w <= 0)
			continue;
		
		mSZ += w;

//		float bIso = Util::u2b(a.uIso());
//		if (bIso < 4)
//			bIso = 4;
//		float x = sqrt(bIso) / mResHigh;
//		x = w * (2 * anorm(x) - 1 - C * x * exp(-0.5 * pow(x, 2))) / pow(x, 3);
//		
//		so += x;
	}

//	auto bo = mSZ;

	mSZ = mSZ * mSpacegroup.num_symops() / mVC;
//	mMeanBIso = pow(mResHigh * errsol(bo / so), 2);

	// Calculate overall rms data
	vector<AtomData> atomData;

	for (auto atom: mStructure.atoms())
	{
		AtomShape shape(atom, mResHigh, mResLow, mElectronScattering);
		
		float radius = shape.radius();
		
		if (cif::VERBOSE > 2)
			cerr << (atomData.size() + 1) << '\t'
				 << AtomTypeTraits(atom.type()).symbol() << '\t'
				 << radius << endl;

		atomData.emplace_back(atom, radius);
	}
	
	GridPtDataMap gridPointDensity;
	map<string,vector<double>> zScoresPerAsym;
	sumDensity(atomData, gridPointDensity, zScoresPerAsym);

	// Now that we have the density data, we can calculate the correction/rescale factors
	for (auto zsc: zScoresPerAsym)
	{
		// collect array of z-scores
		vector<double>& zdca0 = zsc.second;

//		double qa, qb;
//		tie(qa, qb) = interpolateCumulativeProbabilities(zdca0, mVF);
		auto& z = zdca0;
		auto vf = mVF;

		sort(z.begin(), z.end());
		
		double qa = 0, qb = 1;
		
		size_t nd = z.size();
		size_t n = static_cast<size_t>(round(vf * nd));
		
		if (n > 100)
		{
			size_t i1 = static_cast<size_t>((n + 1) * anorm(-1.5)) + 1;
			size_t i2 = static_cast<size_t>((n + 1) * anorm(1.5));
		
			size_t ns = i2 - i1 + 1;
		
			double vr = (nd - 1) / (n - 1.0);
			double sw = 0, swx = 0, swxs = 0, swy = 0, swxy = 0, swys = 0;
			
			for (auto i = i1; i <= i2; ++i)
			{
				double qx = phinvs(static_cast<double>(i) / (n + 1));
				double x = vr * i;
				size_t j = static_cast<size_t>(x);
				x -= j;
				
		//		assert(j < z.size());
				if (j < 1 or j >= z.size())
					continue;
		
				auto qyd = (1.0 - x) * z[j - 1] + x * z[j] - qx;
				
				auto wx = exp(-0.5 * qx * qx);
				sw += wx;
				swx += wx * qx;
				swxs += wx * qx * qx;
				swy += wx * qyd;
				swxy += wx * qx * qyd;
				swys += wx * qyd * qyd;
			}
			
			double dd = 1.0 / (sw * swxs - swx * swx);
			double qa = dd * (swxs * swy - swx * swxy);
			double qb = dd * (sw * swxy - swx * swy);
			
			if (cif::VERBOSE > 1)
			{
				swys = dd * (swys - (qa * swy + qb * swxy)) / (ns - 2);
				cerr << endl
					 << "Intercept & gradient before LS: " << qa << " (" << sqrt(swys * swxs) << ") " << qb << " (" << sqrt(swys * sw) << ')' << endl;
			}

			qb += 1.0;
		
			if (cif::VERBOSE > 1)
			{
				cerr << endl
					 << "Rescale SD(delta-rho) using Q-Q plot for asym " << zsc.first << ':' << endl
					 << string(54, '=') << endl
					 << "Input & updated SD(delta-rho): " << mRMSDensityFd << " ; " << qb * mRMSDensityFd << endl
					 << endl;
			}
		}
				
		mRmsScaled[zsc.first] = make_pair(qa * mRMSDensityFd, qb * mRMSDensityFd);
	}
}

vector<ResidueStatistics> StatsCollector::collect() const
{
	vector<tuple<string,int,string,string>> residues;
	vector<Atom> atoms;
	
	for (auto atom: mStructure.atoms())
	{
		if (atom.isWater())
			continue;
		
		auto k = make_tuple(atom.labelAsymID(), atom.labelSeqID(), atom.labelCompID(), atom.authSeqID());
//		auto k = make_tuple(atom.authAsymID(), atom.property<string>("auth_seq_id"), atom.authCompID());

		if (residues.empty() or residues.back() != k)
		{
			residues.emplace_back(move(k));
			atoms.emplace_back(move(atom));
		}
	}	
	
	BoundingBox bbox(mStructure, atoms, 5.0f);
	return collect(residues, bbox, true);
}

vector<ResidueStatistics> StatsCollector::collect(const string& asymID, int resFirst, int resLast, bool authNameSpace) const
{
	vector<tuple<string,int,string,string>> residues;
	vector<Atom> atoms;
	
	for (auto atom: mStructure.atoms())
	{
		if (atom.isWater())	
			continue;
		
		if (authNameSpace)
		{
			int authSeqID = stoi(atom.authSeqID());
			
			if (atom.authAsymID() != asymID or authSeqID < resFirst or authSeqID > resLast)
				continue;
		}
		else
		{
			if (atom.labelAsymID() != asymID or atom.labelSeqID() < resFirst or atom.labelSeqID() > resLast)
				continue;
		}

		auto k = make_tuple(atom.labelAsymID(), atom.labelSeqID(), atom.labelCompID(), atom.authSeqID());
//		auto k = make_tuple(atom.authAsymID(), atom.property<string>("auth_seq_id"), atom.authCompID());

		if (residues.empty() or residues.back() != k)
		{
			residues.emplace_back(move(k));
			atoms.emplace_back(move(atom));
		}
	}	
	
	BoundingBox bbox(mStructure, atoms, 5.0f);
	return collect(residues, bbox, false);
}

vector<ResidueStatistics> StatsCollector::collect(const vector<tuple<string,int,string,string>>& residues,
	BoundingBox& bbox, bool addWaters) const
{
	vector<AtomData> atomData;

//	BoundingBox bb(mStructure, residues, 5.0f);

	for (auto atom: mStructure.atoms())
	{
		if (atom.isWater())
		{
			if (not addWaters)
				continue;
		}
		else if (not bbox.contains(atom.location()))
			continue;
		
		AtomShape shape(atom, mResHigh, mResLow, mElectronScattering);
		
		float radius = shape.radius();
		
		if (cif::VERBOSE > 2)
			cerr << (atomData.size() + 1) << '\t'
				 << AtomTypeTraits(atom.type()).symbol() << '\t'
				 << radius << endl;

		atomData.emplace_back(atom, radius);
	}	
	
	calculate(atomData);
	
	set<string> missing;
	vector<ResidueStatistics> result;

	// And now collect the per residue information
	for (auto r: residues)
	{
		int seqID;
		string asymID, compID, authSeqID;
		tie(asymID, seqID, compID, authSeqID) = r;
		
		AtomDataSums sums;
		double ediaSum = 0;
		size_t n = 0, m = 0;

		auto comp = Compound::create(compID);
		if (comp == nullptr)
		{
			if (not missing.count(compID) and compID != "HOH")
				cerr << "Missing information for compound '" << compID << '\'' << endl;
			missing.insert(compID);
			
			for (const auto& d: atomData)
			{
				if (d.asymID != asymID or d.seqID != seqID)
					continue;
				
				sums += d.sums;
				ediaSum += pow(d.edia + 0.1, -2);
				++n;
				
				if (d.edia >= 0.8)
					++m;
			}
		}
		else
		{
			for (auto& compAtom: comp->atoms())
			{
				if (compAtom.typeSymbol == H)
					continue;

				++n;
				
				auto ci = find_if(atomData.begin(), atomData.end(),
					[=](auto& d) { return d.asymID == asymID and d.seqID == seqID and d.atom.labelAtomID() == compAtom.id; });

				if (ci == atomData.end())
				{
					if (compAtom.id == "OXT")
						--n;
					else if (cif::VERBOSE > 1)
						cerr << "Missing atom '" << compAtom.id << "' in residue " << asymID << ':' << seqID << endl;
					continue;
				}
				
				sums += ci->sums;
				ediaSum += pow(ci->edia + 0.1, -2);
				
				if (ci->edia >= 0.8)
					++m;
			}
		}
		
		result.emplace_back(ResidueStatistics{asymID, seqID, compID,
			authSeqID,
			(sums.rfSums[0] / sums.rfSums[1]),				// rsr
			sums.srg(),										// srsr
			sums.cc(),										// rsccs
			1 / sqrt(ediaSum / n) - 0.1,					// ediam
			100. * m / n,									// opia
			static_cast<int>(round(mVF * sums.ngrid))});	// ngrid
	}
	
	if (addWaters)
	{
		for (const auto& d: atomData)
		{
			const Atom& atom = d.atom;
			
			if (not atom.isWater())
				continue;
	
			result.emplace_back(ResidueStatistics{d.asymID, d.seqID, "HOH", 
				atom.authSeqID(),
				(d.sums.rfSums[0] / d.sums.rfSums[1]),			// rsr
				d.sums.srg(),									// srsr
				d.sums.cc(),									// rsccs
				d.edia,											// ediam
				(d.edia > 0.8 ? 100. : 0.),						// opia
				static_cast<int>(round(mVF * d.sums.ngrid))});	// ngrid
	
		}
	}
	
	return result;
}

ResidueStatistics StatsCollector::collect(initializer_list<const Residue*> residues) const
{
	vector<Atom> atoms;
	for (auto& r: residues)
		for (auto a: r->atoms())
			atoms.push_back(a);
	
	return collect(atoms);
}

ResidueStatistics StatsCollector::collect(initializer_list<Atom> atoms) const
{
	vector<Atom> v(atoms);
	return collect(v);
}

ResidueStatistics StatsCollector::collect(const vector<Atom>& atoms) const
{
	vector<AtomData> atomData;

	BoundingBox bb(mStructure, atoms, 4.f);

	for (auto atom: mStructure.atoms())
	{
		if (not bb.contains(atom.location()))
			continue;
		
		AtomShape shape(atom, mResHigh, mResLow, mElectronScattering);
		
		float radius = shape.radius();
		
		if (cif::VERBOSE > 2)
			cerr << (atomData.size() + 1) << '\t'
				 << AtomTypeTraits(atom.type()).symbol() << '\t'
				 << radius << endl;

		atomData.emplace_back(atom, radius);
	}	
	
	calculate(atomData);
	
	AtomDataSums sums;
	size_t n = 0, m = 0;
	double ediaSum = 0;

	for (auto& atom: atoms)
	{
		++n;
		
		auto ci = find_if(atomData.begin(), atomData.end(),
			[=](auto& d) { return d.asymID == atom.labelAsymID() and d.seqID == atom.labelSeqID() and d.atom.labelAtomID() == atom.labelAtomID(); });
		
		if (ci == atomData.end())
			continue;

		sums += ci->sums;
		ediaSum += pow(ci->edia + 0.1, -2);
		
		if (ci->edia >= 0.8)
			++m;
	}

	ResidueStatistics result{"", 0, "", "",
			(sums.rfSums[0] / sums.rfSums[1]),				// rsr
			sums.srg(),										// srsr
			sums.cc(),										// rsccs
			1 / sqrt(ediaSum / n) - 0.1,					// ediam
			100. * m / n,									// opia
			static_cast<int>(round(mVF * sums.ngrid))		// ngrid		
	};
	
	return result;
}

void StatsCollector::sumDensity(vector<AtomData>& atomData,
	GridPtDataMap& gridPointDensity, map<string,vector<double>>& zScoresPerAsym) const
{
	using namespace clipper;
	
	const Xmap<float>& Fb = mMapMaker.fb();
	const Xmap<float>& Fd = mMapMaker.fd();
	
	// First step, iterate over atoms, then over grid points covered by this atom
	// collecting per gridpoint statistics

	for (auto& data: atomData)
	{
		auto& atom = data.atom;
		
		AtomShape shape(atom, mResHigh, mResLow, mElectronScattering);

		string asymID = data.asymID;
		if (atom.isWater())
			asymID = "0";

		auto radius = data.radius;
		double sumDensity = 0;

		iterateGrid(atom.location(), radius, Fb, [&](Xmap_base::Map_reference_coord& iw)
		{
			auto p = Point(iw.coord_orth());
			
			double d = Distance(p, atom.location());

			if (d <= radius)
			{
				double density = shape.calculatedDensity(p);
				
				assert(not isnan(density));
	
				gridPointDensity[iw.coord()] += density;
				data.points.emplace_back(iw.coord(), density);
				
				sumDensity += density;
				
				zScoresPerAsym[data.asymID].push_back(Fd[iw] / (Fd.multiplicity(iw.coord()) * mRMSDensityFd));
			}
		});

		data.averageDensity = sumDensity / data.points.size();
	}
}

void StatsCollector::collectSums(vector<AtomData>& atomData, GridPtDataMap& gridPointDensity) const
{
	using namespace clipper;
	
	const Xmap<float>& Fb = mMapMaker.fb();
	const Xmap<float>& Fd = mMapMaker.fd();

	cif::Progress progress(atomData.size(), "Stats calculation");

	// Iterate over the atom data to collect the sums
	for (auto& d: atomData)
	{
		auto rmsScaledF = mRmsScaled.at(d.asymID);

		for (auto gp: d.points)
		{
			++d.sums.ngrid;
			
			double e = gp.density / gridPointDensity[gp.p];
			double t = e * mSZ / rmsScaledF.second;
			
			Xmap_base::Map_reference_coord ix(Fb, gp.p);
			
			double fb = Fb[ix];
			double fd = Fd[ix];

			double ed1 = e * (fb - rmsScaledF.first) / rmsScaledF.second + t;
			double ed2 = e * (fd - rmsScaledF.first) / rmsScaledF.second;
			double ed3 = ed1 - ed2;

			d.sums.rfSums[0] += abs(ed2);
			d.sums.rfSums[1] += abs(ed1 + ed3);
			
			double w = gp.density / d.averageDensity;
			if (w < 0)
				w = 0;
			if (w > 1)
				w = 1;
			
			d.sums.rgSums[0] += w * ed2 * ed2;
			d.sums.rgSums[1] += w * ed1 * ed1;
			
			d.sums.swSums[0] += (w * ed2) * (w * ed2);
			d.sums.swSums[1] += (w * ed1) * (w * ed2);
			d.sums.swSums[2] += (w * ed1) * (w * ed1);

			ed1 -= t;
			ed3 -= t;
			
			d.sums.ccSums[0] += ed1 * ed3;
			d.sums.ccSums[1] += ed1 * ed1;
			d.sums.ccSums[2] += ed3 * ed3;
			
			d.sums.edSums[0] += ed1;
			d.sums.edSums[1] += ed3;
		}
		
		progress.consumed(1);
	}
}

void StatsCollector::calculate(vector<AtomData>& atomData) const
{
	GridPtDataMap gridPointDensity;
	map<string,vector<double>> zScoresPerAsym;

	sumDensity(atomData, gridPointDensity, zScoresPerAsym);
	collectSums(atomData, gridPointDensity);
}

// --------------------------------------------------------------------

EDIAStatsCollector::EDIAStatsCollector(MapMaker<float>& mm,
	Structure& structure, bool electronScattering, const BondMap& bondMap)
	: StatsCollector(mm, structure, electronScattering)
	, mDistanceMap(structure, mm.spacegroup(), mm.cell(), 3.5f), mBondMap(bondMap)
{
	
	// create a atom radius map, for EDIA

	const double kResolutions[] =
	{
		0.5, 1.0, 1.5, 2.0, 2.5
	};

	// The following numbers were harvested with the application collect-b-factors
	const double kAverageBFactors[] =
	{
		6.31912,		// 0.5
		14.4939,		// 1.0
		20.8827,		// 1.5
		27.7075,		// 2.0
		55.6378			// 2.5
	};
	const int kAverageBFactorCount = sizeof(kAverageBFactors) / sizeof(double);
	
	int i = static_cast<int>(floor(mResHigh / 0.5)) - 1;
	if (i > kAverageBFactorCount - 1)
		i = kAverageBFactorCount - 1;
	if (i < 0)
		i = 0;

	float ediaBFactor;	
	if (i < kAverageBFactorCount - 1)
		ediaBFactor = kAverageBFactors[i] +
			((kAverageBFactors[i + 1] - kAverageBFactors[i]) * (mResHigh - kResolutions[i]) / (kResolutions[i + 1]  - kResolutions[i]));
	else
		ediaBFactor = kAverageBFactors[i];

	if (cif::VERBOSE)
		cerr << "Calculating radii with B Factor " << ediaBFactor << endl;
	
	for (auto atom: mStructure.atoms())
	{
		if (mRadii.count(atom.type()))
			continue;
		
		AtomShape shape(atom, mResHigh, mResLow, mElectronScattering, ediaBFactor);
		mRadii[atom.type()] = shape.radius();
		
		if (cif::VERBOSE)
			cerr << "Radius for atom with type " << AtomTypeTraits(atom.type()).symbol() << " is " << mRadii[atom.type()] << endl;
	}
}

void EDIAStatsCollector::calculate(vector<AtomData>& atomData) const
{
	StatsCollector::calculate(atomData);
	
	const Xmap<float>& Fb = mMapMaker.fb();
//	Xmap<float>& fd = mMapMaker.fd();
		
	struct lessAtom
	{
		bool operator()(const Atom& a, const Atom& b) const { return a.id().compare(b.id()) < 0; }
	};
	
	typedef set<Atom, lessAtom> atomSet;

	// Calculate EDIA scores
	
	cif::Progress progress(atomData.size(), "EDIA calculation");

	for (auto& data: atomData)
	{
		auto& atom = data.atom;
		float radius = mRadii.at(atom.type());
		
//		if (cif::VERBOSE > 2)
//			cerr << (atomData.size() + 1) << '\t'
//				 << AtomTypeTraits(atom.type()).symbol() << '\t'
//				 << radius << endl;
//
		PointWeightFunction w(atom.location(), radius);
		
		vector<Atom> atomsNearBy = mDistanceMap.near(atom, 3.5f);
	
		vector<PointWeightFunction> wn;
		for (auto a: atomsNearBy)
			wn.emplace_back(a.location(), mRadii.at(a.type()));

		float ediaSum[2] = {};

		iterateGrid(atom.location(), radius, Fb, [&](auto iw)
		{
			Point p = iw.coord_orth();
			
			// EDIA calculations
			float z = (Fb[iw] - mMeanDensityFb) / mRMSDensityFb;
			
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
					
					if (not mBondMap(atomsNearBy[i], atom))
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

//					if (cif::VERBOSE > 2)
//						cout << Point(p) << ":\td: " << xmap[iw] << "\tz: " << z << "\to: " << o << "\tzraw: " << ((xmap[iw] - meanDensity) / rmsDensity) << "\twp: " << wp << endl;
			
			ediaSum[0] += z * wp * o;
			if (wp > 0)
				ediaSum[1] += wp;

		});

		data.edia = ediaSum[0] / ediaSum[1];
		
		progress.consumed(1);
	}
}

}
