#pragma once

#include <clipper/clipper.h>
#include <boost/filesystem/path.hpp>

#include "cif++/Structure.h"

namespace libcif
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
	
	// These routines work with CCP4 map files
	void read(const boost::filesystem::path& f);
	void write(const boost::filesystem::path& f);
	
	clipper::Spacegroup spacegroup() const				{ return mMap.spacegroup(); }
	clipper::Cell cell() const							{ return mMap.cell(); }

  private:

	Xmap mMap;
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

bool IsMTZFile(const boost::filesystem::path& p);

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
	
	void loadMTZ(const boost::filesystem::path& mtzFile,
		float samplingRate = 4.5,
		std::initializer_list<std::string> fbLabels = { "FWT", "PHWT" },
		std::initializer_list<std::string> fdLabels = { "DELFWT", "PHDELWT" },
		std::initializer_list<std::string> foLabels = { "FP", "SIGFP" },
		std::initializer_list<std::string> fcLabels = { "FC_ALL", "PHIC_ALL" });

	void loadMaps(
		const boost::filesystem::path& fbMapFile,
		const boost::filesystem::path& fdMapFile,
		float reshi, float reslo);

	// following works on both mtz files and structure factor files in CIF format
	void calculate(const boost::filesystem::path& hklin,
		const Structure& structure,
		bool noBulk, AnisoScalingFlag anisoScaling,
		float samplingRate = 4.5, bool electronScattering = false,
		std::initializer_list<std::string> foLabels = { "FP", "SIGFP" },
		std::initializer_list<std::string> freeLabels = { "FREE" });

	void writeMTZ(const boost::filesystem::path& file,
		const std::string& project, const std::string& crystal);

	MapType& fb()								{ return mFb; }
	MapType& fd()								{ return mFd; }

	const MapType& fb() const					{ return mFb; }
	const MapType& fd() const					{ return mFd; }
	
	double resLow() const						{ return mResLow; }
	double resHigh() const						{ return mResHigh; }

	const Spacegroup& spacegroup() const		{ return mHKLInfo.spacegroup(); }
	const Cell& cell() const					{ return mHKLInfo.cell(); }
	const Grid_sampling& gridSampling() const	{ return mGrid; }

  private:

	void loadFoFreeFromReflectionsFile(const boost::filesystem::path& hklin);
	void loadFoFreeFromMTZFile(const boost::filesystem::path& hklin,
		std::initializer_list<std::string> foLabels,
		std::initializer_list<std::string> freeLabels);
	
	void recalc(const Structure& structure,
		bool noBulk, AnisoScalingFlag anisoScaling,
		float samplingRate = 4.5, bool electronScattering = false);

	void fixMTZ();
	void printStats();
	
	MapType				mFb, mFd;
	Grid_sampling		mGrid;
	float				mSamplingRate;
	double				mResLow, mResHigh;
	int					mNumRefln = 1000, mNumParam = 20;
	
	// Cached raw data
	HKL_info			mHKLInfo;
	HKL_data<F_sigF>	mFoData;
	HKL_data<Flag>		mFreeData;
	HKL_data<F_phi>		mFcData, mFbData, mFdData;
	HKL_data<Phi_fom>	mPhiFomData;
};

}
