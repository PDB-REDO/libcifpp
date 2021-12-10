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

// Calculate DSSP-like secondary structure information

#pragma once

namespace mmcif
{
	
class Structure;
class Monomer;

struct Res;

extern const float
	kCouplingConstant, kMinHBondEnergy, kMaxHBondEnergy;

enum SecondaryStructureType : char
{
	ssLoop			= ' ',
	ssAlphahelix	= 'H',
	ssBetabridge	= 'B',
	ssStrand		= 'E',
	ssHelix_3		= 'G',
	ssHelix_5		= 'I',
	ssHelix_PPII	= 'P',
	ssTurn			= 'T',
	ssBend			= 'S'
};

enum class HelixType
{
	rh_3_10, rh_alpha, rh_pi, rh_pp
};

enum class Helix
{
	None, Start, End, StartAndEnd, Middle
};

//struct HBond
//{
//	std::string 				labelAsymID;
//	int							labelSeqID;
//	double						energy;
//};
//
//struct BridgePartner
//{
//	std::string					labelAsymID;
//	int							labelSeqID;
//	int							ladder;
//	bool						parallel;
//};

struct SecondaryStructure
{
	SecondaryStructureType		type;
//	HBond						donor[2], acceptor[2];
//	BridgePartner				beta[2];
//	int							sheet;
//	bool						bend;
};

//void CalculateSecondaryStructure(Structure& s);

const size_t
	kHistogramSize = 30;

struct DSSP_Statistics
{
	uint32_t nrOfResidues, nrOfChains, nrOfSSBridges, nrOfIntraChainSSBridges, nrOfHBonds;
	uint32_t nrOfHBondsInAntiparallelBridges, nrOfHBondsInParallelBridges;
	uint32_t nrOfHBondsPerDistance[11] = {};
	double accessibleSurface = 0;

	uint32_t residuesPerAlphaHelixHistogram[kHistogramSize] = {};
	uint32_t parallelBridgesPerLadderHistogram[kHistogramSize] = {};
	uint32_t antiparallelBridgesPerLadderHistogram[kHistogramSize] = {};
	uint32_t laddersPerSheetHistogram[kHistogramSize] = {};
};

enum class ChainBreak
{
	None, NewChain, Gap
};

class DSSP
{
  public:
	DSSP(const Structure& s, int min_poly_proline_stretch_length, bool calculateSurfaceAccessibility);
	~DSSP();
	
	DSSP(const DSSP&) = delete;
	DSSP& operator=(const DSSP&) = delete;
	
	SecondaryStructureType operator()(const std::string& inAsymID, int inSeqID) const;
	SecondaryStructureType operator()(const Monomer& m) const;
	
	double accessibility(const std::string& inAsymID, int inSeqID) const;
	double accessibility(const Monomer& m) const;

	bool isAlphaHelixEndBeforeStart(const Monomer& m) const;
	bool isAlphaHelixEndBeforeStart(const std::string& inAsymID, int inSeqID) const;

	DSSP_Statistics GetStatistics() const;

	class iterator;
	using res_iterator = typename std::vector<Res>::iterator;

	class ResidueInfo
	{
	  public:
		friend class iterator;

		ResidueInfo()
			: mImpl(nullptr)
		{
		}

		ResidueInfo(const ResidueInfo &rhs)
			: mImpl(rhs.mImpl)
		{
		}

		ResidueInfo& operator=(const ResidueInfo &rhs)
		{
			mImpl = rhs.mImpl;
			return *this;
		}

		explicit operator bool() const		{ return not empty(); }
		bool empty() const					{ return mImpl == nullptr; }

		const Monomer& residue() const;
		std::string alt_id() const;

		/// \brief return 0 if not a break, ' ' in case of a new chain and '*' in case of a broken chain
		ChainBreak chainBreak() const;

		/// \brief the internal number in DSSP
		int nr() const;

		SecondaryStructureType ss() const;
		
		int ssBridgeNr() const;

		Helix helix(HelixType helixType) const;

		bool bend() const;

		double accessibility() const;

		/// \brief returns resinfo, ladder and parallel
		std::tuple<ResidueInfo,int,bool> bridgePartner(int i) const;

		int sheet() const;

		/// \brief return resinfo and the energy of the bond
		std::tuple<ResidueInfo,double> acceptor(int i) const;
		std::tuple<ResidueInfo,double> donor(int i) const;

	  private:
		ResidueInfo(Res* res) : mImpl(res) {}

		Res* mImpl;
	};

	class iterator
	{
	  public:
		using iterator_category = std::input_iterator_tag;
		using value_type = ResidueInfo;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type*;
		using reference = value_type&;

		iterator(const iterator& i);
		iterator(Res* res);
		iterator& operator=(const iterator& i);

		reference operator*()		{ return mCurrent; }
		pointer operator->()		{ return &mCurrent; }

		iterator& operator++();
		iterator operator++(int)
		{
			auto tmp(*this);
			this->operator++();
			return tmp;
		}

		bool operator==(const iterator& rhs) const		{ return mCurrent.mImpl == rhs.mCurrent.mImpl; }
		bool operator!=(const iterator& rhs) const		{ return mCurrent.mImpl != rhs.mCurrent.mImpl; }

	  private:
		ResidueInfo	mCurrent;
	};

	iterator begin() const;
	iterator end() const;

	bool empty() const		{ return begin() == end(); }

  private:
	struct DSSPImpl* mImpl;
};


}
