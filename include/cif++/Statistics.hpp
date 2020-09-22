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

#include "cif++/MapMaker.hpp"
#include "cif++/DistanceMap.hpp"
#include "cif++/BondMap.hpp"

namespace mmcif
{

// --------------------------------------------------------------------

struct AtomData;
class BoundingBox;

struct ResidueStatistics
{
	std::string	asymID;
	int			seqID;
	std::string compID;
	std::string authSeqID;
	
	double RSR, SRSR, RSCCS, EDIAm, OPIA;
	int ngrid;
};

std::ostream& operator<<(std::ostream& os, const ResidueStatistics& st);

// --------------------------------------------------------------------

class StatsCollector
{
  public:
	StatsCollector(const StatsCollector&) = delete;
	StatsCollector& operator=(const StatsCollector&) = delete;

	StatsCollector(const mmcif::MapMaker<float>& mm,
		mmcif::Structure& structure, bool electronScattering);

	virtual std::vector<ResidueStatistics> collect() const;

	virtual std::vector<ResidueStatistics> collect(const std::string& asymID,
		int resFirst, int resLast, bool authNameSpace = true) const;

	virtual ResidueStatistics collect(std::initializer_list<const mmcif::Residue*> residues) const;

	virtual ResidueStatistics collect(std::initializer_list<mmcif::Atom> atoms) const;

	virtual ResidueStatistics collect(const std::vector<mmcif::Atom>& atoms) const;
	
  protected:

	// asym-seqid-compid
	std::vector<ResidueStatistics> collect(
		const std::vector<std::tuple<std::string,int,std::string,std::string>>& residues,
		BoundingBox& bbox, bool addWaters) const;
	
	void initialize();

	virtual void calculate(std::vector<AtomData>& atomData) const;

	struct cmpGPt
	{
		bool operator()(const clipper::Coord_grid& a, const clipper::Coord_grid& b) const
		{
			int d = a.u() - b.u();
			if (d == 0)
				d = a.v() - b.v();
			if (d == 0)
				d = a.w() - b.w();
			return d < 0;
		}
	};

	typedef std::map<clipper::Coord_grid,double,cmpGPt> GridPtDataMap;


	mmcif::Structure& mStructure;
	const mmcif::MapMaker<float>& mMapMaker;

	clipper::Spacegroup mSpacegroup;
	clipper::Cell mCell;
	clipper::Grid_sampling mGrid;
	float mResHigh, mResLow;
	bool mElectronScattering;

	std::map<std::string,std::pair<double,double>> mRmsScaled;

	void collectSums(std::vector<AtomData>& atomData, GridPtDataMap& gridPointDensity) const;
	void sumDensity(std::vector<AtomData>& atomData,
		GridPtDataMap& gridPointDensity, std::map<std::string,std::vector<double>>& zScoresPerAsym) const;
	
	// Other variables we cache
	
	double mMeanDensityFb, mRMSDensityFb, mRMSDensityFd;
	double mSZ;		// average electron density in cell
	double mVF;		// degrees of freedom
	double mVC;		// cell volume?
};

// --------------------------------------------------------------------

class EDIAStatsCollector : public StatsCollector
{
  public:
	EDIAStatsCollector(mmcif::MapMaker<float>& mm,
		mmcif::Structure& structure, bool electronScattering,
		const mmcif::BondMap& bondMap);

  protected:

	virtual void calculate(std::vector<AtomData>& atomData) const;

	mmcif::DistanceMap mDistanceMap;
	const mmcif::BondMap& mBondMap;
	std::map<mmcif::AtomType,float>	mRadii;
};

}
