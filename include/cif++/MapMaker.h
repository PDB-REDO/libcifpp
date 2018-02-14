#pragma once

#include <clipper/clipper.h>
#include <boost/filesystem/path.hpp>

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
	
	MapMaker(Xmap& fb, Xmap& fd, bool fixMTZ, float samplingRate = 4.5);
	~MapMaker();

	void LoadFromMTZ(const boost::filesystem::path& mtzFile,
		std::initializer_list<std::string> fbLabels = { "FWT", "PHWT" },
		std::initializer_list<std::string> fdLabels = { "DELFWT", "PHDELWT" },
		std::initializer_list<std::string> foLabels = { "FP", "SIGFP" },
		std::initializer_list<std::string> fcLabels = { "FC_ALL", "PHIC_ALL" });
	void RecalculateFromMTZ(
		std::initializer_list<std::string> foLabels = { "FP", "SIGFP" },
		std::initializer_list<std::string> freeLabels = { "FP", "SIGFP" });
	void LoadFromMapFiles();

	double RMSDensityFb() const			{ return mRMSDensityFb; }
	double MeanDensityFb() const		{ return mMeanDensityFb; }
	double RMSDensityFd() const			{ return mRMSDensityFd; }
	double MeanDensityFd() const		{ return mMeanDensityFd; }
	
	double ResLow() const				{ return mResLow; }
	double ResHight() const				{ return mResHigh; }

  private:

	void FixMTZ(FPdata& fb, FPdata& fd, FOdata& fo, FPdata& fc);

	Xmap&		mFb;
	Xmap&		mFd;

	Spacegroup		mSpacegroup;
	Cell			mCell;
	Grid_sampling	mGrid;

	bool		mFixMTZ;
	float		mSamplingRate;

	double		mRMSDensityFb, mRMSDensityFd;
	double		mMeanDensityFb, mMeanDensityFd;
	double		mResLow, mResHigh;
};

}
