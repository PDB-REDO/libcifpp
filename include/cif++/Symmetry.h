// copyright

#pragma once

#include "cif++/Structure.h"

namespace mmcif
{

// --------------------------------------------------------------------
// Functions to use when working with symmetry stuff

clipper::Coord_orth CalculateOffsetForCell(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell);
std::vector<clipper::RTop_orth> AlternativeSites(const clipper::Spacegroup& spacegroup, const clipper::Cell& cell);
int32_t GetRotationalIndexNumber(const clipper::Spacegroup& spacegroup, const clipper::Cell& cell,
	const clipper::RTop_orth& rt);

// --------------------------------------------------------------------
// To iterate over all near symmetry copies of an atom

class SymmetryAtomIteratorFactory
{
  public:
	// SymmetryAtomIteratorFactory(const Structure& p);
	SymmetryAtomIteratorFactory(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell)
		: mD(CalculateOffsetForCell(p, spacegroup, cell))
		, mRtOrth(AlternativeSites(spacegroup, cell))
		, mSpacegroup(spacegroup)
		, mCell(cell) {}

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

	std::string symop_mmcif(const Atom& a) const;

  private:
	Point									mD;			// needed to move atoms to center
	std::vector<clipper::RTop_orth>			mRtOrth;
	clipper::Spacegroup						mSpacegroup;
	clipper::Cell							mCell;
};

}
