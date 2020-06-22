// copyright

#include "cif++/Config.h"

#include <atomic>
#include <mutex>

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

clipper::Coord_orth DistanceMap::CalculateOffsetForCell(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell)
{
	auto& atoms = p.atoms();
	size_t dim = atoms.size();
	
	vector<clipper::Coord_orth> locations;
	locations.reserve(dim);
	
	// bounding box
	Point pMin(numeric_limits<float>::max(), numeric_limits<float>::max(), numeric_limits<float>::max()),
		  pMax(numeric_limits<float>::min(), numeric_limits<float>::min(), numeric_limits<float>::min());
	
	for (auto& atom: atoms)
	{
		auto p = atom.location();
		locations.push_back(p);

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
	vector<float> c(dim);
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

	if (cif::VERBOSE > 1)
		cerr << "median position of atoms: " << Point(mx, my, mz) << endl;
	
	auto calculateD = [&](float m, float c)
	{
		float d = 0;
		if (c != 0)
		{
			while (m + d < -(c / 2))
				d += c;
			while (m + d > (c / 2))
				d -= c;
		}
		return d;
	};

	Point D;

	if (cell.a() == 0 or cell.b() == 0 or cell.c() == 0)
		throw runtime_error("Invalid cell, contains a dimension that is zero");

	D.mX = calculateD(mx, cell.a());
	D.mY = calculateD(my, cell.b());
	D.mZ = calculateD(mz, cell.c());
	
	if (D.mX != 0 or D.mY != 0 or D.mZ != 0)
	{
		if (cif::VERBOSE)
			cerr << "moving coorinates by " << D.mX << ", " << D.mY << " and " << D.mZ << endl;
	}

	return D;	
}

// --------------------------------------------------------------------

vector<clipper::RTop_orth> DistanceMap::AlternativeSites(const clipper::Spacegroup& spacegroup,
	const clipper::Cell& cell)
{
	vector<clipper::RTop_orth> result;
	
	// to make the operation at index 0 equal to identity
	result.push_back(clipper::RTop_orth::identity());

	for (int i = 0; i < spacegroup.num_symops(); ++i)
	{
		const auto& symop = spacegroup.symop(i);
		
		for (int u: { -1, 0, 1})
			for (int v: { -1, 0, 1})
				for (int w: { -1, 0, 1})
				{
					if (i == 0 and u == 0 and v == 0 and w == 0)
						continue;

					auto rtop = clipper::RTop_frac(
							symop.rot(), symop.trn() + clipper::Vec3<>(u, v, w)
						).rtop_orth(cell);

					result.push_back(move(rtop));
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

	if (cif::VERBOSE > 1)
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
		if (cif::VERBOSE)
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
	
	cif::Progress progress(residues.size() * residues.size(), "Creating distance map");
	
	for (size_t i = 0; i + 1 < residues.size(); ++i)
	{
		auto& ri = *residues[i];
		
		Point centerI;
		float radiusI;
		tie(centerI, radiusI) = ri.centerAndRadius();
		
		for (size_t j = i + 1; j < residues.size(); ++j)
		{
			progress.consumed(1);

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
			
			int32_t kbest = 0;
			for (int32_t k = 1; k < static_cast<int32_t>(mRtOrth.size()); ++k)
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
				AddDistancesForAtoms(ri, rj, dist, kbest);
		}
	}
	
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
}

// --------------------------------------------------------------------

void DistanceMap::AddDistancesForAtoms(const Residue& a, const Residue& b, DistMap& dm, int32_t rtix)
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
	
//	if (ixb < ixa)
//		std::swap(ixa, ixb);

	size_t L = mIA[ixa];
	size_t R = mIA[ixa + 1] - 1;
	
	while (L <= R)
	{
		size_t i = (L + R) / 2;

		if (mJA[i] == ixb)
			return get<0>(mA[i]);

		if (mJA[i] < ixb)
			L = i + 1;
		else
			R = i - 1;
	}
	
	return 100.f;
}

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
		int32_t rti;
		tie(d, rti) = mA[i];

		if (d > maxDistance)
			continue;

		size_t ixb = mJA[i];
		Atom b = structure.getAtomById(rIndex.at(ixb));
		
		if (rti > 0)
			result.emplace_back(b.symmetryCopy(mD, mRtOrth.at(rti)));
		else if (rti < 0)
			result.emplace_back(b.symmetryCopy(mD, mRtOrth.at(-rti).inverse()));
		else
			result.emplace_back(b);
	}
	
	return result;
}

}
