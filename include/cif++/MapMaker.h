#pragma once

#include <clipper/clipper.h>

#include "cif++/Structure.h"

namespace mmcif
{

template<typename FTYPE>
class Map
{
  public:
	typedef FTYPE										ftype;
	typedef typename clipper::Xmap<ftype>				Xmap;
	
	Map();
	~Map();

	void calculateStats();

	double rmsDensity() const							{ return mRMSDensity; }
	double meanDensity() const							{ return mMeanDensity; }
	
	operator Xmap& ()									{ return mMap; }
	operator const Xmap& () const						{ return mMap; }
	Xmap& get()											{ return mMap; }
	const Xmap& get() const								{ return mMap; }
	
	// These routines work with CCP4 map files
	void read(const std::string& f);
	void write(const std::string& f);

	void write_masked(std::ostream& os, clipper::Grid_range range);
	void write_masked(const std::string& f, 
		clipper::Grid_range range);
	
	clipper::Spacegroup spacegroup() const				{ return mMap.spacegroup(); }
	clipper::Cell cell() const							{ return mMap.cell(); }

  private:

	Xmap mMap;
	double mMinDensity, mMaxDensity;
	double mRMSDensity, mMeanDensity;
};

using clipper::HKL_info;
using clipper::HKL_data;
using clipper::data32::F_phi;
using clipper::data32::F_sigF;
using clipper::data32::Phi_fom;
using clipper::data32::Flag;

using clipper::Spacegroup;
using clipper::Cell;
using clipper::Grid_sampling;

// --------------------------------------------------------------------

bool IsMTZFile(const std::string& p);

// --------------------------------------------------------------------
	
template<typename FTYPE>
class MapMaker
{
  public:
	typedef Map<FTYPE> MapType;
	typedef typename MapType::Xmap Xmap;

	enum AnisoScalingFlag {
		as_None, as_Observed, as_Calculated
	};
	
	MapMaker();
	~MapMaker();
	
	MapMaker(const MapMaker&) = delete;
	MapMaker& operator=(const MapMaker&) = delete;
	
	void loadMTZ(const std::string& mtzFile,
		float samplingRate,
		std::initializer_list<std::string> fbLabels = { "FWT", "PHWT" },
		std::initializer_list<std::string> fdLabels = { "DELFWT", "PHDELWT" },
		std::initializer_list<std::string> foLabels = { "FP", "SIGFP" },
		std::initializer_list<std::string> fcLabels = { "FC_ALL", "PHIC_ALL" },
		std::initializer_list<std::string> faLabels = { "FAN", "PHAN" });

	void loadMaps(
		const std::string& fbMapFile,
		const std::string& fdMapFile,
		float reshi, float reslo);

	// following works on both mtz files and structure factor files in CIF format
	void calculate(const std::string& hklin,
		const Structure& structure,
		bool noBulk, AnisoScalingFlag anisoScaling,
		float samplingRate, bool electronScattering = false,
		std::initializer_list<std::string> foLabels = { "FP", "SIGFP" },
		std::initializer_list<std::string> freeLabels = { "FREE" });

	void recalc(const Structure& structure,
		bool noBulk, AnisoScalingFlag anisoScaling,
		float samplingRate, bool electronScattering = false);

	void printStats();

	void writeMTZ(const std::string& file,
		const std::string& project, const std::string& crystal);

	MapType& fb()								{ return mFb; }
	MapType& fd()								{ return mFd; }
	MapType& fa()								{ return mFa; }

	const MapType& fb() const					{ return mFb; }
	const MapType& fd() const					{ return mFd; }
	const MapType& fa() const					{ return mFa; }
	
	double resLow() const						{ return mResLow; }
	double resHigh() const						{ return mResHigh; }

	const Spacegroup& spacegroup() const		{ return mHKLInfo.spacegroup(); }
	const Cell& cell() const					{ return mHKLInfo.cell(); }
	const Grid_sampling& gridSampling() const	{ return mGrid; }

  private:

	void loadFoFreeFromReflectionsFile(const std::string& hklin);
	void loadFoFreeFromMTZFile(const std::string& hklin,
		std::initializer_list<std::string> foLabels,
		std::initializer_list<std::string> freeLabels);
	
	void fixMTZ();
	
	MapType				mFb, mFd, mFa;
	Grid_sampling		mGrid;
	float				mSamplingRate;
	double				mResLow, mResHigh;
	int					mNumRefln = 1000, mNumParam = 20;
	
	// Cached raw data
	HKL_info			mHKLInfo;
	HKL_data<F_sigF>	mFoData;
	HKL_data<Flag>		mFreeData;
	HKL_data<F_phi>		mFcData, mFbData, mFdData, mFaData;
	HKL_data<Phi_fom>	mPhiFomData;
};

}
