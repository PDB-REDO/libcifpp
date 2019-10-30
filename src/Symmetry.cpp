// copyright

#include "cif++/Config.h"

#include <atomic>
#include <mutex>

#include <boost/thread.hpp>

#include "cif++/Symmetry.h"
#include "cif++/CifUtils.h"

using namespace std;

namespace mmcif
{

// --------------------------------------------------------------------

clipper::Coord_orth CalculateOffsetForCell(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell)
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
		while (m + d < -(c / 2))
			d += c;
		while (m + d > (c / 2))
			d -= c;
		return d;
	};

	Point D;

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

vector<clipper::RTop_orth> AlternativeSites(const clipper::Spacegroup& spacegroup,
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

// -----------------------------------------------------------------------

string SymmetryAtomIteratorFactory::symop_mmcif(const Atom& a) const
{
	string result;

	if (not a.isSymmetryCopy())
		result = "1_555";
	else
	{
		auto rtop_o = a.symop();

		for (int i = 0; i < mSpacegroup.num_symops(); ++i)
		{
			const auto& symop = mSpacegroup.symop(i);

			for (int u: { -1, 0, 1})
				for (int v: { -1, 0, 1})
					for (int w: { -1, 0, 1})
					{
						// if (i == 0 and u == 0 and v == 0 and w == 0)
						// 	continue;

						auto rtop = clipper::RTop_frac(
								symop.rot(), symop.trn() + clipper::Vec3<>(u, v, w)
							).rtop_orth(mCell);
						
						if (rtop.rot().equals(rtop_o.rot(), 0.00001) and rtop.trn().equals(rtop_o.trn(), 0.000001))
						{
							// gotcha

							auto rtop_f = rtop.rtop_frac(mCell);

							uint32_t t[3] = {
								static_cast<uint32_t>(5 + static_cast<int>(rint(rtop_f.trn()[0]))),
								static_cast<uint32_t>(5 + static_cast<int>(rint(rtop_f.trn()[1]))),
								static_cast<uint32_t>(5 + static_cast<int>(rint(rtop_f.trn()[2])))
							};

							if (t[0] > 9 or t[1] > 9 or t[2] > 9)
								throw runtime_error("Symmetry operation has an out-of-range translation.");

							result += to_string(1 + i) + "_"
								   + to_string(t[0])
								   + to_string(t[1])
								   + to_string(t[2]);

							return result;
						}
					}
		}
	}

	return result;
}

}
