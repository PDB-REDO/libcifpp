#include "cif++/Config.h"

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include <clipper/clipper-ccp4.h>

#include "cif++/MapMaker.h"
#include "cif++/ResolutionCalculator.h"

using namespace std;
namespace fs = boost::filesystem;
namespace ba = boost::algorithm;

extern int VERBOSE;

namespace libcif
{

template<typename FTYPE>
MapMaker<FTYPE>::MapMaker(Xmap& fb, Xmap& fd, bool fixMTZ, float samplingRate)
	: mFb(fb), mFd(fd)
	, mFixMTZ(fixMTZ), mSamplingRate(samplingRate)
{
}

template<typename FTYPE>
MapMaker<FTYPE>::~MapMaker()
{
}

template<typename FTYPE>
void MapMaker<FTYPE>::LoadFromMTZ(const fs::path& mtzFile,
	initializer_list<string> fbLabels, initializer_list<string> fdLabels,
	initializer_list<string> foLabels, initializer_list<string> fcLabels)
{
	if (VERBOSE)
		cerr << "Reading map from " << mtzFile << endl;
	
	const string kBasePath("/%1%/%2%/[%3%]");

	clipper::CCP4MTZfile mtzin;
	mtzin.open_read(mtzFile.string());
	
	clipper::HKL_info hklInfo;
	mtzin.import_hkl_info(hklInfo);
	
	clipper::HKL_data<clipper::data32::F_phi> fbData, fdData, fcData;
	clipper::HKL_data<clipper::data32::F_sigF> foData;

	mtzin.import_hkl_data(fbData,
		(boost::format(kBasePath) % "*" % "*" % ba::join(fbLabels, ",")).str());
	mtzin.import_hkl_data(fdData,
		(boost::format(kBasePath) % "*" % "*" % ba::join(fdLabels, ",")).str());
	mtzin.import_hkl_data(foData,
		(boost::format(kBasePath) % "*" % "*" % ba::join(foLabels, ",")).str());
	mtzin.import_hkl_data(fcData,
		(boost::format(kBasePath) % "*" % "*" % ba::join(fcLabels, ",")).str());

	mtzin.close_read();
	
	mCell = hklInfo.cell();
	mSpacegroup = hklInfo.spacegroup();

	ResolutionCalculator rc(mCell);
	mResHigh = 99;
	mResLow = 0;
	
	for (auto hi = foData.first_data(); not hi.last(); hi = foData.next_data(hi))
	{
		float res = rc(hi.hkl().h(), hi.hkl().k(), hi.hkl().l());
		
		if (mResHigh > res)
			mResHigh = res;

		if (mResLow < res)
			mResLow = res;
	}

	if (VERBOSE > 1)
		cerr << "calculated reshi = " << mResHigh << " reslo = " << mResLow << endl;
		
	if (mFixMTZ)
		FixMTZ(fbData, fdData, foData, fcData);

	float samplingRate = mSamplingRate / 2;
	
	mGrid = Grid_sampling(mSpacegroup, mCell,
		hklInfo.resolution(), samplingRate);	// define grid
	
	mFb = Xmap(mSpacegroup, mCell, mGrid);	// define map
	mFb.fft_from(fbData);														// generate map
	
	mFd = Xmap(mSpacegroup, mCell, mGrid);	// define map
	mFd.fft_from(fdData);														// generate map

	if (VERBOSE)
		cerr << "Read Xmaps with sampling rate: " << samplingRate << endl
			 << "  resolution: " << mResHigh
			 << endl;
}

template<typename FTYPE>
void MapMaker<FTYPE>::FixMTZ(FPdata& fb, FPdata& fd, FOdata& fo, FPdata& fc)
{
	const clipper::data32::F_phi fzero(0, 0);
	
	// mtzfix...
	for (auto ih = fb.first(); not ih.last(); ih.next())
	{
		if (fo[ih].missing())
		{
			fb[ih] = fc[ih];
			fd[ih] = fzero;
			continue;
		}
		
		if (fb[ih].missing() or fd[ih].missing())
			continue;

		if (abs(fmod(abs(fb[ih].phi() - fc[ih].phi()) + 180, 360) - 180) > 90)
			fb[ih] = -fb[ih];

		if (abs(fmod(abs(fd[ih].phi() - fc[ih].phi()) + 180, 360) - 180) > 90)
			fd[ih] = -fd[ih]; 

		clipper::HKL_class cls(mSpacegroup, ih.hkl());

		auto mFo = fb[ih] - fd[ih];
		auto DFc = mFo - fd[ih];

		if (cls.centric())
			fb[ih] = mFo;
		else
			fd[ih] = fd[ih] + fd[ih];
		
		fc[ih] = DFc;
	}
}

//void MapMaker::RecalculateFromMTZ(bool fixMtz)
//{
//	
//}
//
//void MapMaker::LoadFromMapFiles(bool fixMtz)
//{
//	
//}

template class MapMaker<float>;
template class MapMaker<double>;

}
