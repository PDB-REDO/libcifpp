#pragma once

#include <clipper/clipper.h>
#include <boost/filesystem/path.hpp>

#include "cif++/Structure.h"

namespace libcif
{
	
template<typename FTYPE>
class MapMaker
{
  public:
	typedef FTYPE										ftype;
	typedef typename clipper::Xmap<ftype>				Xmap;
	typedef clipper::HKL_data<clipper::data32::F_phi>	FPdata;
	typedef clipper::HKL_data<clipper::data32::F_sigF>	FOdata;
	
	
	typedef clipper::Spacegroup							Spacegroup;
	typedef clipper::Cell								Cell;
	typedef clipper::Grid_sampling						Grid_sampling;
	
	enum AnisoScalingFlag {
		as_None, as_Observed, as_Calculated
	};
	
	MapMaker(Xmap& fb, Xmap& fd);
	~MapMaker();

	void loadFromMTZ(const boost::filesystem::path& mtzFile,
		float samplingRate = 4.5,
		std::initializer_list<std::string> fbLabels = { "FWT", "PHWT" },
		std::initializer_list<std::string> fdLabels = { "DELFWT", "PHDELWT" },
		std::initializer_list<std::string> foLabels = { "FP", "SIGFP" },
		std::initializer_list<std::string> fcLabels = { "FC_ALL", "PHIC_ALL" });
	void recalculateFromMTZ(const boost::filesystem::path& mtzFile,
		const Structure& structure,
		bool noBulk, AnisoScalingFlag anisoScaling,
		float samplingRate = 4.5,
		std::initializer_list<std::string> foLabels = { "FP", "SIGFP" },
		std::initializer_list<std::string> freeLabels = { "FREE" });
	void loadFromMapFiles(const boost::filesystem::path& fbMapFile,
		const boost::filesystem::path& fdMapFile);

	double rmsDensityFb() const			{ return mRMSDensityFb; }
	double meanDensityFb() const		{ return mMeanDensityFb; }
	double rmsDensityFd() const			{ return mRMSDensityFd; }
	double meanDensityFd() const		{ return mMeanDensityFd; }
	
	double resLow() const				{ return mResLow; }
	double resHigh() const				{ return mResHigh; }

	const Spacegroup& spacegroup() const		{ return mSpacegroup; }
	const Cell& cell() const					{ return mCell; }
	const Grid_sampling& gridSampling() const	{ return mGrid; }

  private:

	void fixMTZ(FPdata& fb, FPdata& fd, FOdata& fo, FPdata& fc);

	Xmap&		mFb;
	Xmap&		mFd;

	Spacegroup		mSpacegroup;
	Cell			mCell;
	Grid_sampling	mGrid;

	float		mSamplingRate;

	double		mRMSDensityFb, mRMSDensityFd;
	double		mMeanDensityFb, mMeanDensityFd;
	double		mResLow, mResHigh;
};

}
