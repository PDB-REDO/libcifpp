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

#pragma once

#include "cif++/Structure.hpp"

namespace mmcif
{

// --------------------------------------------------------------------
// Functions to use when working with symmetry stuff

clipper::Coord_orth CalculateOffsetForCell(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell);
std::vector<clipper::RTop_orth> AlternativeSites(const clipper::Spacegroup& spacegroup, const clipper::Cell& cell);
int GetSpacegroupNumber(std::string spacegroup);	// alternative for clipper's parsing code
// std::string SpacegroupToHall(std::string spacegroup);

// --------------------------------------------------------------------
// To iterate over all symmetry copies of an atom

class SymmetryAtomIteratorFactory
{
  public:
	// SymmetryAtomIteratorFactory(const Structure& p);

	// SymmetryAtomIteratorFactory(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell)
	// 	: mSpacegroupNr(spacegroup.spacegroup_number())
	// 	, mSpacegroup(spacegroup)
	// 	, mD(CalculateOffsetForCell(p, spacegroup, cell))
	// 	, mRtOrth(AlternativeSites(spacegroup, cell))
	// 	, mCell(cell) {}

	SymmetryAtomIteratorFactory(const Structure& p, int spacegroupNr, const clipper::Cell& cell);

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
	int										mSpacegroupNr;
	clipper::Spacegroup						mSpacegroup;
	Point									mD;			// needed to move atoms to center
	std::vector<clipper::RTop_orth>			mRtOrth;
	clipper::Cell							mCell;
};

}
