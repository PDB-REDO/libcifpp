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

#include "cif++/Config.hpp"

#include <atomic>
#include <mutex>

#include "cif++/Symmetry.hpp"
#include "cif++/CifUtils.hpp"

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
		assert(c != 0);
		if (c != 0)
		{
			while (m + d < -(c / 2))
				d += c;
			while (m + d > (c / 2))
				d -= c;
		}
		return d;
	};

	if (cell.a() == 0 or cell.b() == 0 or cell.c() == 0)
		throw runtime_error("Invalid cell, contains a dimension that is zero");

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

// --------------------------------------------------------------------
// Unfortunately, clipper has a different numbering scheme than PDB
// for rotation numbers. So we created a table to map those.
// Perhaps a bit over the top, but hey....

#include "SymOpTable_data.cpp"

int32_t GetRotationalIndexNumber(int spacegroup, const clipper::RTop_frac& rt)
{
	auto& rot = rt.rot();
	auto& trn = rt.trn();

	auto rte = [&rot](int i, int j) { return static_cast<int8_t>(lrint(rot(i, j))); };

	SymopData k
	{
		{
			rte(0, 0), rte(0, 1), rte(0, 2),
			rte(1, 0), rte(1, 1), rte(1, 2),
			rte(2, 0), rte(2, 1), rte(2, 2)
		}
	};

	for (int i = 0; i < 3; ++i)
	{
		int n = lrint(trn[i] * 24);
		int d = 24;

		if (n == 0 or abs(n) == 24)
			continue;		// is 0, 0 in our table

		for (int i = 5; i > 1; --i)
			if (n % i == 0 and d % i == 0)
			{
				n /= i;
				d /= i;
			}
		
		n = (n + d) % d;

		switch (i)
		{
			case 0: k.trn_0_0 = n; k.trn_0_1 = d; break;
			case 1: k.trn_1_0 = n; k.trn_1_1 = d; break;
			case 2: k.trn_2_0 = n; k.trn_2_1 = d; break;
		}
	}

	const size_t N = sizeof(kSymopNrTable) / sizeof(SymopDataBlock);
	int32_t L = 0, R = static_cast<int32_t>(N - 1);
	while (L <= R)
	{
		int32_t i = (L + R) / 2;
		if (kSymopNrTable[i].spacegroupNr < spacegroup)
			L = i + 1;
		else
			R = i - 1;
	}

	for (size_t i = L; i < N and kSymopNrTable[i].spacegroupNr == spacegroup; ++i)
	{
		if (kSymopNrTable[i].rt.iv == k.iv)
			return kSymopNrTable[i].rotationalNr;
	}

	throw runtime_error("Symmetry operation was not found in table, cannot find rotational number");
}

// --------------------------------------------------------------------
// And unfortunately, clipper does not know all the spacegroup names mentioned in symop.lib

int GetSpacegroupNumber(std::string spacegroup)
{
	if (spacegroup == "P 21 21 2 A")
		spacegroup = "P 21 21 2 (a)";
	else if (spacegroup.empty())
		throw runtime_error("No spacegroup, cannot continue");

	int result = 0;

	const size_t N = sizeof(kSpaceGroups) / sizeof(Spacegroup);
	int32_t L = 0, R = static_cast<int32_t>(N - 1);
	while (L <= R)
	{
		int32_t i = (L + R) / 2;

		int d = spacegroup.compare(kSpaceGroups[i].name);

		if (d > 0)
			L = i + 1;
		else if (d < 0)
			R = i - 1;
		else
		{
			result = kSpaceGroups[i].nr;
			break;
		}
	}

	// not found, see if we can find a match based on xHM name
	if (result == 0)
	{
		for (auto& sp: kSpaceGroups)
		{
			if (sp.xHM == spacegroup)
			{
				result = sp.nr;
				break;
			}
		}
	}

	if (result == 0)
		throw runtime_error("Spacegroup name " + spacegroup + " was not found in table");
	
	return result;
}

// -----------------------------------------------------------------------

std::string SpacegroupToHall(std::string spacegroup)
{
	int nr = GetSpacegroupNumber(spacegroup);

	// yeah, sucks, I know, might be looping three times this way
	string result;
	for (auto& sp: kSpaceGroups)
	{
		if (sp.nr == nr)
		{
			result = sp.Hall;
			break;
		}
	}

	if (result.empty())
		throw runtime_error("Spacegroup name " + spacegroup + " was not found in table");
	
	return result;
}

// clipper::Spgr_descr GetCCP4SpacegroupDescr(int nr)
// {
// 	const size_t N = sizeof(kSymopNrTable) / sizeof(SymopDataBlock);
// 	int32_t L = 0, R = static_cast<int32_t>(N - 1);
// 	while (L <= R)
// 	{
// 		int32_t i = (L + R) / 2;
// 		if (kSymopNrTable[i].spacegroupNr < nr)
// 			L = i + 1;
// 		else
// 			R = i - 1;
// 	}

// 	if (L < 0 or L >= N)
// 		throw runtime_error("Invalid spacegroup number");

// 	clipper::Spgr_descr::Symop_codes symop_codes;

// 	for (size_t i = L; i < N and kSymopNrTable[i].spacegroupNr == nr; ++i)
// 	{
// 		auto& symop = kSymopNrTable[i];

// 		clipper::ftype m[4][4] = {
// 			{ symop.rt.rot_0_0, symop.rt.rot_0_1, symop.rt.rot_0_2, static_cast<float>(symop.rt.trn_0_0) / symop.rt.trn_0_1 },
// 			{ symop.rt.rot_1_0, symop.rt.rot_1_1, symop.rt.rot_1_2, static_cast<float>(symop.rt.trn_1_0) / symop.rt.trn_1_1 },
// 			{ symop.rt.rot_2_0, symop.rt.rot_2_1, symop.rt.rot_2_2, static_cast<float>(symop.rt.trn_2_0) / symop.rt.trn_2_1 },
// 			{ 0, 0, 0, 1 }
// 		};

// 		symop_codes.emplace_back(clipper::Symop(m));
// 	}

// 	return clipper::Spgr_descr(symop_codes);
// }

clipper::Spgr_descr GetCCP4SpacegroupDescr(int nr)
{
	for (auto& sg: kSpaceGroups)
	{
		if (sg.nr == nr)
			return clipper::Spgr_descr(sg.Hall, clipper::Spgr_descr::Hall);
	}

	throw runtime_error("Invalid spacegroup number: " + to_string(nr));
}

// -----------------------------------------------------------------------

SymmetryAtomIteratorFactory::SymmetryAtomIteratorFactory(const Structure& p, int spacegroupNr, const clipper::Cell& cell)
	: mSpacegroupNr(spacegroupNr)
	, mSpacegroup(GetCCP4SpacegroupDescr(spacegroupNr))
	, mD(CalculateOffsetForCell(p, mSpacegroup, cell))
	, mRtOrth(AlternativeSites(mSpacegroup, cell))
	, mCell(cell)
{

}

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

							int rnr = GetRotationalIndexNumber(mSpacegroupNr, rtop_f);

							uint32_t t[3] = {
								static_cast<uint32_t>(5 + static_cast<int>(rint(rtop_f.trn()[0]))),
								static_cast<uint32_t>(5 + static_cast<int>(rint(rtop_f.trn()[1]))),
								static_cast<uint32_t>(5 + static_cast<int>(rint(rtop_f.trn()[2])))
							};

							if (t[0] > 9 or t[1] > 9 or t[2] > 9)
								throw runtime_error("Symmetry operation has an out-of-range translation.");

							result += to_string(rnr) + "_"
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
