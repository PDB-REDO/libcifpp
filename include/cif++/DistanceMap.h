// copyright

#pragma once

#include <unordered_map>

#include <clipper/clipper.h>

#include "cif++/Structure.h"

namespace mmcif
{

class DistanceMap
{
  public:
	DistanceMap(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell,
		float maxDistance);

	// simplified version for subsets of atoms (used in refining e.g.)
//	DistanceMap(const Structure& p, const std::vector<Atom>& atoms);
	
	DistanceMap(const DistanceMap&) = delete;
	DistanceMap& operator=(const DistanceMap&) = delete;

	float operator()(const Atom& a, const Atom& b) const;

	std::vector<Atom> near(const Atom& a, float maxDistance = 3.5f) const;
	std::vector<Atom> near(const Point& p, float maxDistance = 3.5f) const;

	static clipper::Coord_orth
		CalculateOffsetForCell(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell);

	static std::vector<clipper::RTop_orth>
		AlternativeSites(const clipper::Spacegroup& spacegroup, const clipper::Cell& cell);

  private:

	typedef std::map<std::tuple<size_t,size_t>,std::tuple<float,int32_t>> DistMap;

	void AddDistancesForAtoms(const Residue& a, const Residue& b, DistMap& dm, int32_t rtix);

	const Structure&						structure;
	size_t									dim;
	std::unordered_map<std::string,size_t>	index;
	std::map<size_t,std::string>			rIndex;
	
	float									mMaxDistance, mMaxDistanceSQ;
	
	std::vector<std::tuple<float,int32_t>>	mA;
	std::vector<size_t>						mIA, mJA;
	Point									mD;			// needed to move atoms to center
	std::vector<clipper::RTop_orth>			mRtOrth;
};

// --------------------------------------------------------------------

class SymmetryAtomIteratorFactory
{
  public:
	// SymmetryAtomIteratorFactory(const Structure& p);
	SymmetryAtomIteratorFactory(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell)
		: mD(DistanceMap::CalculateOffsetForCell(p, spacegroup, cell))
		, mRtOrth(DistanceMap::AlternativeSites(spacegroup, cell)) {}

	SymmetryAtomIteratorFactory(const SymmetryAtomIteratorFactory&) = delete;
	SymmetryAtomIteratorFactory& operator=(const SymmetryAtomIteratorFactory&) = delete;

	class SymmetryAtomIterator : public std::iterator<std::forward_iterator_tag, const Atom>
	{
	  public:
		typedef std::iterator<std::forward_iterator_tag, const Atom>	baseType;
		typedef typename baseType::pointer								pointer;
		typedef typename baseType::reference							reference;

		SymmetryAtomIterator(const SymmetryAtomIteratorFactory& factory, const Atom& atom)
			: m_f(&factory), m_i(0), m_a(atom), m_c(atom) {}

		SymmetryAtomIterator(const SymmetryAtomIteratorFactory& factory, const Atom& atom, int)
			: SymmetryAtomIterator(factory, atom)
		{
			m_i = m_f->mRtOrth.size();
		}

		SymmetryAtomIterator(const SymmetryAtomIterator& iter)
			: m_f(iter.m_f), m_i(iter.m_i), m_a(iter.m_a), m_c(iter.m_c) {}

		SymmetryAtomIterator& operator=(const SymmetryAtomIterator& iter)
		{
			if (this != &iter)
			{
				m_f = iter.m_f;
				m_i = iter.m_i;
				m_a = iter.m_a;
				m_c = iter.m_c;
			}
			return *this;
		}

		reference operator*()		{ return m_c; }
		pointer operator->()		{ return &m_c; }

		SymmetryAtomIterator operator++()
		{
			if (++m_i < m_f->mRtOrth.size())
				m_c = m_a.symmetryCopy(m_f->mD, m_f->mRtOrth[m_i]);
			return *this;
		}

		SymmetryAtomIterator operator++(int)
		{
			SymmetryAtomIterator result(*this);
			this->operator++();
			return result;
		}

		bool operator==(const SymmetryAtomIterator& iter) const
		{
			return m_f == iter.m_f and m_i == iter.m_i;
		}

		bool operator!=(const SymmetryAtomIterator& iter) const
		{
			return m_f != iter.m_f or m_i != iter.m_i;
		}

	  private:
		const SymmetryAtomIteratorFactory* m_f;
		size_t m_i;
		Atom m_a, m_c;
	};

	class SymmetryAtomIteratorRange
	{
	  public:
		SymmetryAtomIteratorRange(const SymmetryAtomIteratorFactory& f, const Atom& a)
			: m_f(f), m_a(a) {}

		SymmetryAtomIterator begin()
		{
			return SymmetryAtomIterator(m_f, m_a);
		}

		SymmetryAtomIterator end()
		{
			return SymmetryAtomIterator(m_f, m_a, 1);
		}

	  private:
		const SymmetryAtomIteratorFactory& m_f;
		Atom m_a;
	};

	SymmetryAtomIteratorRange operator()(const Atom& a) const
	{
		return SymmetryAtomIteratorRange(*this, a);
	}

  private:
	Point									mD;			// needed to move atoms to center
	std::vector<clipper::RTop_orth>			mRtOrth;
};

}
