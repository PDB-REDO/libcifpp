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
	
template<typename FTYPE>
class MapMaker
{
  public:
	typedef Map<FTYPE> MapType;
	typedef typename MapType::Xmap Xmap;

	typedef clipper::HKL_data<clipper::data32::F_phi>	FPdata;
	typedef clipper::HKL_data<clipper::data32::F_sigF>	FOdata;
	typedef clipper::HKL_data<clipper::data32::Phi_fom>	WData;
	
	
	typedef clipper::Spacegroup							Spacegroup;
	typedef clipper::Cell								Cell;
	typedef clipper::Grid_sampling						Grid_sampling;
	
	enum AnisoScalingFlag {
		as_None, as_Observed, as_Calculated
	};
	
	MapMaker();
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
		float samplingRate = 4.5, bool electronScattering = false,
		std::initializer_list<std::string> foLabels = { "FP", "SIGFP" },
		std::initializer_list<std::string> freeLabels = { "FREE" });

	void loadFromMapFiles(
		const boost::filesystem::path& fbMapFile,
		const boost::filesystem::path& fdMapFile,
		float reshi, float reslo);

	MapType& fb()								{ return mFb; }
	MapType& fd()								{ return mFd; }

	const MapType& fb() const					{ return mFb; }
	const MapType& fd() const					{ return mFd; }
	
	double resLow() const						{ return mResLow; }
	double resHigh() const						{ return mResHigh; }

	const Spacegroup& spacegroup() const		{ return mSpacegroup; }
	const Cell& cell() const					{ return mCell; }
	const Grid_sampling& gridSampling() const	{ return mGrid; }

  private:

	void fixMTZ(FPdata& fb, FPdata& fd, FOdata& fo, FPdata& fc, WData& fom);

	MapType mFb, mFd;
	Spacegroup mSpacegroup;
	Cell mCell;
	Grid_sampling mGrid;
	float mSamplingRate;
	double mResLow, mResHigh;
};

}
