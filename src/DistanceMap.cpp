// copyright

#include "cif++/Config.h"

#include <atomic>

#include <boost/thread.hpp>

#include "cif++/DistanceMap.h"
#include "cif++/CifUtils.h"

using namespace std;

namespace mmcif
{

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
	: dim(0)
{
	auto atoms = p.atoms();
	dim = atoms.size();
	
	dist = vector<float>(dim * (dim - 1), 0.f);
	vector<clipper::Coord_orth> locations(dim);
	vector<bool> isWater(dim, false);
	
	for (auto& atom: atoms)
	{
		size_t ix = index.size();
		index[atom.id()] = ix;
		
		locations[ix] = atom.location();
		isWater[ix] = atom.isWater();
	};

	vector<clipper::RTop_orth> rtOrth = AlternativeSites(spacegroup, cell);
	
	cif::Progress progress(locations.size() - 1, "Creating distance map");
	
	boost::thread_group t;
	size_t N = boost::thread::hardware_concurrency();
	atomic<size_t> next(0);

	for (size_t i = 0; i < N; ++i)
		t.create_thread([&]()
		{
			for (;;)
			{
				size_t i = next++;
				
				if (i >= locations.size())
					break;
				
				for (size_t j = i + 1; j < locations.size(); ++j)
				{
//					if (not (isWater[i] or isWater[j]))
//						continue;
					
					// find nearest location based on spacegroup/cell
					double minR2 = numeric_limits<double>::max();
					
					for (const auto& rt: rtOrth)
					{
						auto p = locations[j].transform(rt);
						
						double r2 = (locations[i] - p).lengthsq();
						if (minR2 > r2)
							minR2 = r2;
					}
		
					size_t ix = j + i * dim - i * (i + 1) / 2;
					
					assert(ix < dist.size());
					dist[ix] = sqrt(minR2);
				}

				progress.consumed(1);
			}
		});
	
	t.join_all();
}

DistanceMap::DistanceMap(const vector<Atom>& atoms)
	: dim(0)
{
	dim = atoms.size();
	
	dist = vector<float>(dim * (dim - 1), 0.f);
	
	for (auto& atom: atoms)
	{
		size_t ix = index.size();
		index[atom.id()] = ix;
	};

	for (size_t i = 0; i < dim; ++i)
	{
		for (size_t j = i + 1; j < dim; ++j)
		{
			size_t ix = j + i * dim - i * (i + 1) / 2;
			assert(ix < dist.size());
			dist[ix] = Distance(atoms[i].location(), atoms[j].location());
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
	
	if (ixb < ixa)
		swap(ixa, ixb);
	
	size_t ix = ixb + ixa * dim - ixa * (ixa + 1) / 2;
	
	assert(ix < dist.size());
	return dist[ix];
}

vector<Atom> DistanceMap::near(const Atom& a, float maxDistance) const
{
	vector<Atom> result;
	const File& f = a.getFile();
	
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
		
		size_t ix;
		
		if (ixa < ixb)
			ix = ixb + ixa * dim - ixa * (ixa + 1) / 2;
		else
			ix = ixa + ixb * dim - ixb * (ixb + 1) / 2;
		
		assert(ix < dist.size());
		if (dist[ix] != 0 and dist[ix] <= maxDistance)
			result.emplace_back(f, i.first);
	}
	
	return result;
}

}
