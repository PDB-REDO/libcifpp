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
	
	// to make index 0 equal to identity, no 
	result.push_back(clipper::RTop_orth::identity());
	
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

DistanceMap::DistanceMap(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell,
		float maxDistance)
	: structure(p), dim(0), mMaxDistance(maxDistance), mMaxDistanceSQ(maxDistance * maxDistance)
{
	auto& atoms = p.atoms();
	dim = atoms.size();
	
	vector<clipper::Coord_orth> locations(dim);
	
	// bounding box
	Point pMin(numeric_limits<float>::max(), numeric_limits<float>::max(), numeric_limits<float>::max()),
		  pMax(numeric_limits<float>::min(), numeric_limits<float>::min(), numeric_limits<float>::min());
	
	for (auto& atom: atoms)
	{
		size_t ix = index.size();
		index[atom.id()] = ix;
		rIndex[ix] = atom.id();
		
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
	
	clipper::Coord_orth D = mD;
	
	if (mD.mX != 0 or mD.mY != 0 or mD.mZ != 0)
	{
		if (VERBOSE)
			cerr << "moving coorinates by " << mD.mX << ", " << mD.mY << " and " << mD.mZ << endl;
		
		for_each(locations.begin(), locations.end(), [&](auto& p) { p += mD; });
	}
	
	pMin -= mMaxDistance;	// extend bounding box
	pMax += mMaxDistance;

	mRtOrth = AlternativeSites(spacegroup, cell);
	
	DistMap dist;
	
	vector<const Residue*> residues;
	residues.reserve(p.residues().size());
	for (auto& r: p.residues())
	{
		residues.emplace_back(&r);

		// Add distances for atoms in this residue		
		AddDistancesForAtoms(r, r, dist, 0);
	}
	
	cif::Progress progress(residues.size(), "Creating distance map");
	
	for (size_t i = 0; i + 1 < residues.size(); ++i)
	{
		progress.consumed(1);
		
		auto& ri = *residues[i];
		
		Point centerI;
		float radiusI;
		tie(centerI, radiusI) = ri.centerAndRadius();
		
		for (size_t j = i + 1; j < residues.size(); ++j)
		{
			auto& rj = *residues[j];
			
			// first case, no symmetry operations
			
			Point centerJ;
			float radiusJ;
			tie(centerJ, radiusJ) = rj.centerAndRadius();
			
			auto d = Distance(centerI, centerJ) - radiusI - radiusJ;
			if (d < mMaxDistance)
			{
				AddDistancesForAtoms(ri, rj, dist, 0);
				continue;
			}
			
			// now try all symmetry operations to see if we can move rj close to ri
			
			clipper::Coord_orth cI = centerI;
			clipper::Coord_orth cJ = centerJ;
			
			auto minR2 = d;
			
			int32 kbest = 0;
			for (int32 k = 1; k < static_cast<int32>(mRtOrth.size()); ++k)
			{
				auto& rt = mRtOrth[k];
				
				auto pJ = (cJ + D).transform(rt) - D;
				double r2 = sqrt((cI - pJ).lengthsq()) - radiusI - radiusJ;

				if (minR2 > r2)
				{
					minR2 = r2;
					kbest = k;
				}
			}
			
			if (minR2 < mMaxDistance)
			{
//cout << ri.labelID() << " en " << rj.labelID() << " liggen dicht bij elkaar na symmetrie operatie: " << kbest << endl;
				AddDistancesForAtoms(ri, rj, dist, kbest);
			}
		}
	}

//	cif::Progress progress(locations.size() - 1, "Creating distance map");
//	
//	boost::thread_group t;
//	size_t N = boost::thread::hardware_concurrency();
//	atomic<size_t> next(0);
//	mutex m;
//	
//	for (size_t i = 0; i < N; ++i)
//		t.create_thread([&]()
//		{
//			for (;;)
//			{
//				size_t i = next++;
//				
//				if (i >= locations.size())
//					break;
//				
//				clipper::Coord_orth pi = locations[i];
//
//				for (size_t j = i + 1; j < locations.size(); ++j)
//				{
//					// find nearest location based on spacegroup/cell
//					double minR2 = numeric_limits<double>::max();
//					
//					size_t kb = 0;
//					for (size_t k = 0; k < mRtOrth.size(); ++k)
//					{
//						auto& rt = mRtOrth[k];
//						
//						auto pj = locations[j];
//						
//						pj = pj.transform(rt);
//						double r2 = (pi - pj).lengthsq();
//
//						if (minR2 > r2)
//						{
//							minR2 = r2;
//							kb = k;
//						}
//					}
//					
//					if (minR2 < mMaxDistanceSQ)
//					{
//						float d = sqrt(minR2);
//						auto key = make_tuple(i, j);
//						
//						lock_guard<mutex> lock(m);
//						dist[key] = make_tuple(d, kb);
//					}
//				}
//
//				progress.consumed(1);
//			}
//		});
//	
//	t.join_all();
	
	// Store as a sparse CSR compressed matrix
	
	size_t nnz = dist.size();
	mA.reserve(nnz);
	mIA.reserve(dim + 1);
	mJA.reserve(nnz);
	
	size_t lastR = 0;
	mIA.push_back(0);
	
	for (auto& di: dist)
	{
		size_t c, r;
		tie(r, c) = di.first;

		if (r != lastR)	// new row
		{
			for (size_t ri = lastR; ri < r; ++ri)
				mIA.push_back(mA.size());
			lastR = r;
		}

		mA.push_back(di.second);
		mJA.push_back(c);
	}

	for (size_t ri = lastR; ri < dim; ++ri)
		mIA.push_back(mA.size());

//// debug code
//	
//	assert(mIA.size() == dim + 1);
//	assert(mA.size() == nnz);
//	assert(mJA.size() == nnz);
//	
//cerr << "nnz: " << nnz << endl;
//	
//	auto get = [&](size_t i, size_t j)
//	{
//		if (i > j)
//			std::swap(i, j);
//		
//		for (size_t cix = mIA[i]; cix < mIA[i + 1]; ++cix)
//		{
//			if (mJA[cix] == j)
//				return std::get<0>(mA[cix]);
//		}
//		
//		return 100.f;
//	};
//	
//	auto get2 = [&](size_t ixa, size_t ixb)
//	{
//		if (ixb < ixa)
//			std::swap(ixa, ixb);
//	
//		int32 L = mIA[ixa];
//		int32 R = mIA[ixa + 1] - 1;
//		
//		while (L <= R)
//		{
//			int32 i = (L + R) / 2;
//			
//			if (mJA[i] == ixb)
//				return std::get<0>(mA[i]);
//	
//			if (mJA[i] < ixb)
//				L = i + 1;
//			else
//				R = i - 1;
//		}
//
//		return 100.f;
//	};
//	
//	// test all values
//	for (size_t i = 0; i < dim; ++i)
//		for (size_t j = 0; j < dim; ++j)
//		{
//			float a = get(i, j);
//			
//			auto ixa = i, ixb = j;
//			if (ixb < ixa)
//				std::swap(ixa, ixb);
//			
//			tuple<size_t,size_t> k{ ixa, ixb };
//		
//			auto ii = dist.find(k);
//		
//			float b = 100;
//			
//			if (ii != dist.end())
//				b = std::get<0>(ii->second);
//			
//			assert(a == b);
//			
//			float c = get2(i, j);
//			assert(a == c);
//		}
}

// --------------------------------------------------------------------

void DistanceMap::AddDistancesForAtoms(const Residue& a, const Residue& b, DistMap& dm, int32 rtix)
{
	for (auto& aa: a.atoms())
	{
		clipper::Coord_orth pa = aa.location();
		size_t ixa = index[aa.id()];
		
		for (auto& bb: b.atoms())
		{
			if (aa.id() == bb.id())
				continue;
			
			clipper::Coord_orth pb = bb.location();
			
			if (rtix)
				pb = pb.transform(mRtOrth[rtix]);
			
			auto d = (pa - pb).lengthsq();
			if (d > mMaxDistanceSQ)
				continue;

			d = sqrt(d);

			size_t ixb = index[bb.id()];

			dm[make_tuple(ixa, ixb)] = make_tuple(d, rtix);
			dm[make_tuple(ixb, ixa)] = make_tuple(d, -rtix);
		}
	}
}

//float DistanceMap::operator()(const Atom& a, const Atom& b) const
//{
//	size_t ixa, ixb;
//	
//	try
//	{
//		ixa = index.at(a.id());
//	}
//	catch (const out_of_range& ex)
//	{
//		throw runtime_error("atom " + a.id() + " not found in distance map");
//	}
//		
//	try
//	{
//		ixb = index.at(b.id());
//	}
//	catch (const out_of_range& ex)
//	{
//		throw runtime_error("atom " + b.id() + " not found in distance map");
//	}
//	
////	if (ixb < ixa)
////		std::swap(ixa, ixb);
//
//	size_t L = mIA[ixa];
//	size_t R = mIA[ixa + 1] - 1;
//	
//	while (L <= R)
//	{
//		size_t i = (L + R) / 2;
//
//		if (mJA[i] == ixb)
//			return get<0>(mA[i]);
//
//		if (mJA[i] < ixb)
//			L = i + 1;
//		else
//			R = i - 1;
//	}
//	
//	return 100.f;
//}

vector<Atom> DistanceMap::near(const Atom& a, float maxDistance) const
{
	assert(maxDistance <= mMaxDistance);
	if (maxDistance > mMaxDistance)
		throw runtime_error("Invalid max distance in DistanceMap::near");
	
	size_t ixa;
	try
	{
		ixa = index.at(a.id());
	}
	catch (const out_of_range& ex)
	{
		throw runtime_error("atom " + a.id() + " not found in distance map");
	}

	vector<Atom> result;
	
	for (size_t i = mIA[ixa]; i < mIA[ixa + 1]; ++i)
	{
		float d;
		int32 rti;
		tie(d, rti) = mA[i];

		if (d > maxDistance)
			continue;

		size_t ixb = mJA[i];
		Atom b = structure.getAtomById(rIndex.at(ixb));
		
		clipper::RTop_orth rt = clipper::RTop_orth::identity();
		
		if (rti > 0)
		{
			rt = mRtOrth.at(rti);
			result.emplace_back(b.symmetryCopy(mD, rt));
		}
		else if (rti < 0)
		{
			rt = mRtOrth.at(-rti).inverse();
			result.emplace_back(b.symmetryCopy(mD, rt));
		}
		else
			result.emplace_back(b);

#if 1 //DEBUG
//		if (rti != 0)
//			cerr << "symmetrie contact " << a.labelID() << " en " << result.back().labelID()
//				 << " d: " << d
//				 << " rti: " << rti
//				 << endl;
		
		auto d2 = Distance(a, result.back());
		if (abs(d2 - d) > 0.01)
		{
			cerr << "Voor a: " << a.location() << " en b: " << b.location() << " => " << result.back().location() << endl;
			
			cerr << "Afstand " << a.labelID() << " en " << result.back().labelID()
				 << " is niet gelijk aan verwachtte waarde:"
				 << "d: " << d
				 << " d2: " << d2
				 << " rti: " << rti
				 << endl;
			
			rt = rt.inverse();
			result.back() = b.symmetryCopy(mD, rt);
			
			d2 = Distance(a, result.back());

			cerr << "inverse b: " << result.back().location() << endl;

			
			if (abs(d2 - d) < 0.01)
				cerr << "==> But the inverse is correct" << endl;
		}
#endif
	}
	
	return result;
}

}
