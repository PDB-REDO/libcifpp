// copyright

#include "cif++/Config.h"

#include <atomic>
#include <mutex>

#include <boost/thread.hpp>

#include "cif++/DistanceMap.h"
#include "cif++/CifUtils.h"

using namespace std;

//#define DEBUG_VOOR_BART

namespace mmcif
{

// --------------------------------------------------------------------

inline ostream& operator<<(ostream& os, const Atom& a)
{
	os << a.labelAsymId() << ':' << a.labelSeqId() << '/' << a.labelAtomId();
	
	return os;
}

// --------------------------------------------------------------------

vector<clipper::RTop_orth> AlternativeSites(const clipper::Spacegroup& spacegroup,
	const clipper::Cell& cell)
{
	vector<clipper::RTop_orth> result;
	
	for (int i = 0; i < spacegroup.num_symops(); ++i)
	{
		const auto& symop = spacegroup.symop(i);
		
		for (int u: { -1, 0, 1})
			for (int v: { -1, 0, 1})
				for (int w: { -1, 0, 1})
				{
					result.push_back(
						clipper::RTop_frac(
							symop.rot(), symop.trn() + clipper::Vec3<>(u, v, w)
						).rtop_orth(cell));
				}
	}
	
	return result;
}

// --------------------------------------------------------------------

DistanceMap::DistanceMap(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell)
	: structure(p), dim(0)
{
	const float kMaxDistance = 5, kMaxDistanceSQ = kMaxDistance * kMaxDistance;
	
	auto atoms = p.atoms();
	dim = atoms.size();
	
	vector<clipper::Coord_orth> locations(dim);
	
	// bounding box
	Point pMin(numeric_limits<float>::max(), numeric_limits<float>::max(), numeric_limits<float>::max()),
		  pMax(numeric_limits<float>::min(), numeric_limits<float>::min(), numeric_limits<float>::min());
	
	for (auto& atom: atoms)
	{
		size_t ix = index.size();
		index[atom.id()] = ix;
		
		locations[ix] = atom.location();
		
		auto p = atom.location();

		if (pMin.mX > p.mX)
			pMin.mX = p.mX;
		if (pMin.mY > p.mY)
			pMin.mY = p.mY;
		if (pMin.mZ > p.mZ)
			pMin.mZ = p.mZ;

		if (pMax.mX < p.mX)
			pMax.mX = p.mX;
		if (pMax.mY < p.mY)
			pMax.mY = p.mY;
		if (pMax.mZ < p.mZ)
			pMax.mZ = p.mZ;
	};
	
	// correct locations so that the median of x, y and z are inside the cell
	vector<float> c(locations.size());
	auto median = [&]()
	{
		return dim % 1 == 0
			? c[dim / 2]
			: (c[dim / 2 - 1] + c[dim / 2]) / 2;
	};
	
	transform(locations.begin(), locations.end(), c.begin(), [](auto& l) { return l[0]; });
	sort(c.begin(), c.end());
	float mx = median();
	
	transform(locations.begin(), locations.end(), c.begin(), [](auto& l) { return l[1]; });
	sort(c.begin(), c.end());
	float my = median();
	
	transform(locations.begin(), locations.end(), c.begin(), [](auto& l) { return l[2]; });
	sort(c.begin(), c.end());
	float mz = median();

	if (VERBOSE > 1)
		cerr << "median position of atoms: " << Point(mx, my, mz) << endl;
	
	auto calculateD = [&](float m, float c)
	{
		float d = 0;
		while (m + d < -(c / 2))
			d += c;
		while (m + d > (c / 2))
			d -= c;
		return d;
	};

	mD.mX = calculateD(mx, cell.a());
	mD.mY = calculateD(my, cell.b());
	mD.mZ = calculateD(mz, cell.c());
	
	if (VERBOSE)
		cerr << "Calculated d: " << mD.mX << ", " << mD.mY << " and " << mD.mZ << endl;
	
	if (mD.mX != 0 or mD.mY != 0 or mD.mZ != 0)
	{
		if (VERBOSE)
			cerr << "moving coorinates by " << mD.mX << ", " << mD.mY << " and " << mD.mZ << endl;
		
		for_each(locations.begin(), locations.end(), [&](auto& p) { p += mD; });
	}
	
	pMin -= kMaxDistance;	// extend bounding box
	pMax += kMaxDistance;

	mRtOrth = AlternativeSites(spacegroup, cell);
	
	cif::Progress progress(locations.size() - 1, "Creating distance map");
	
	boost::thread_group t;
	size_t N = boost::thread::hardware_concurrency();
	atomic<size_t> next(0);
	mutex m;
	
	for (size_t i = 0; i < N; ++i)
		t.create_thread([&]()
		{
			for (;;)
			{
				size_t i = next++;
				
				if (i >= locations.size())
					break;
				
				clipper::Coord_orth pi = locations[i];

				for (size_t j = i + 1; j < locations.size(); ++j)
				{
					// find nearest location based on spacegroup/cell
					double minR2 = numeric_limits<double>::max();
					
					size_t kb = 0;
					for (size_t k = 0; k < mRtOrth.size(); ++k)
					{
						auto& rt = mRtOrth[k];
						
						auto pj = locations[j];
						
						pj = pj.transform(rt);
						double r2 = (pi - pj).lengthsq();

						if (minR2 > r2)
						{
							minR2 = r2;
							kb = k;
						}

#if defined(DEBUG_VOOR_BART)
						if (r2 < 3.5 * 3.5 and not rt.equals(clipper::RTop<>::identity(), 0.1))
							cout << "symmetry contact between "
								 << atoms[i] << " at " << pi << " and "
								 << atoms[j] << " at " << pj << endl;
#endif
					}
					
					if (minR2 < kMaxDistanceSQ)
					{
						float d = sqrt(minR2);
						auto key = make_tuple(i, j);
						
						lock_guard<mutex> lock(m);
						dist[key] = make_tuple(d, kb);
					}
				}

				progress.consumed(1);
			}
		});
	
	t.join_all();
}

//DistanceMap::DistanceMap(const Structure& p, const vector<Atom>& atoms)
//	: structure(p), dim(0)
//{
//	dim = atoms.size();
//	
//	for (auto& atom: atoms)
//	{
//		size_t ix = index.size();
//		index[atom.id()] = ix;
//	};
//
//	for (size_t i = 0; i < dim; ++i)
//	{
//		for (size_t j = i + 1; j < dim; ++j)
//		{
//			dist[make_tuple(i, j)] = Distance(atoms[i].location(), atoms[j].location());
//		}
//	}
//}

float DistanceMap::operator()(const Atom& a, const Atom& b) const
{
	size_t ixa, ixb;
	
	try
	{
		ixa = index.at(a.id());
	}
	catch (const out_of_range& ex)
	{
		throw runtime_error("atom " + a.id() + " not found in distance map");
	}
		
	try
	{
		ixb = index.at(b.id());
	}
	catch (const out_of_range& ex)
	{
		throw runtime_error("atom " + b.id() + " not found in distance map");
	}
	
	if (ixb < ixa)
		swap(ixa, ixb);
	
	tuple<size_t,size_t> k{ ixa, ixb };

	auto ii = dist.find(k);

	float result = 100;
	
	if (ii != dist.end())
		result = get<0>(ii->second);
	
	return result;
}

vector<Atom> DistanceMap::near(const Atom& a, float maxDistance) const
{
	vector<Atom> result;
	
	size_t ixa;
	try
	{
		ixa = index.at(a.id());
	}
	catch (const out_of_range& ex)
	{
		throw runtime_error("atom " + a.id() + " not found in distance map");
	}

	for (auto& i: index)
	{
		size_t ixb = i.second;
		
		if (ixb == ixa)
			continue;
		
		auto ii = 
			ixa < ixb ?
				dist.find(make_tuple(ixa, ixb)) :
				dist.find(make_tuple(ixb, ixa));
		
		if (ii == dist.end())
			continue;
		
		float distance;
		size_t rti;
		tie(distance, rti) = ii->second;
		
		if (distance > maxDistance)
			continue;
		
		Atom a = structure.getAtomById(i.first);
		
		result.push_back(a.symmetryCopy(mD, mRtOrth.at(rti)));
	}
	
	return result;
}

}
