#include "cif++/Config.h"

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include <clipper/clipper-contrib.h>
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
MapMaker<FTYPE>::MapMaker(Xmap& fb, Xmap& fd)
	: mFb(fb), mFd(fd)
{
}

template<typename FTYPE>
MapMaker<FTYPE>::~MapMaker()
{
}

template<typename FTYPE>
void MapMaker<FTYPE>::loadFromMTZ(const fs::path& mtzFile, float samplingRate,
	initializer_list<string> fbLabels, initializer_list<string> fdLabels,
	initializer_list<string> foLabels, initializer_list<string> fcLabels)
{
	if (VERBOSE)
		cerr << "Reading map from " << mtzFile << endl;
	
	const string kBasePath("/%1%/%2%/[%3%]");

	using clipper::HKL_info;
	using clipper::CCP4MTZfile;
	using clipper::HKL_data;
	using clipper::data32::F_phi;
	using clipper::data32::F_sigF;

	CCP4MTZfile mtzin;
	mtzin.open_read(mtzFile.string());
	
	HKL_info hklInfo;
	mtzin.import_hkl_info(hklInfo);
	
	HKL_data<F_phi> fbData, fdData, fcData;
	HKL_data<F_sigF> foData;

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
		
	fixMTZ(fbData, fdData, foData, fcData);

	samplingRate /= 2;	// clipper's way of interpreting?
	
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
void MapMaker<FTYPE>::fixMTZ(FPdata& fb, FPdata& fd, FOdata& fo, FPdata& fc)
{
#warning("WARNING: Need the check first to see if fix is necessary!")

	using clipper::HKL_class;
	using clipper::data32::F_phi;

	const F_phi fzero(0, 0);
	
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

		HKL_class cls(mSpacegroup, ih.hkl());

		auto mFo = fb[ih] - fd[ih];
		auto DFc = mFo - fd[ih];

		if (cls.centric())
			fb[ih] = mFo;
		else
			fd[ih] = fd[ih] + fd[ih];
		
		fc[ih] = DFc;
	}
}

template<typename FTYPE>
void MapMaker<FTYPE>::recalculateFromMTZ(const fs::path& mtzFile,
		const Structure& structure, bool noBulk, AnisoScalingFlag anisoScaling, float samplingRate,
		bool electronScattering, initializer_list<string> foLabels, initializer_list<string> freeLabels)
{
	if (VERBOSE)
		cerr << "Recalculating maps from " << mtzFile << endl;
	
	const string kBasePath("/%1%/%2%/[%3%]");

	using clipper::HKL_info;
	using clipper::CCP4MTZfile;
	using clipper::HKL_data;
	using clipper::data32::F_phi;
	using clipper::data32::F_sigF;
	using clipper::data32::Phi_fom;
	using clipper::data32::Flag;

	CCP4MTZfile mtzin;
	mtzin.open_read(mtzFile.string());
	
	HKL_info hklInfo;
//	MTZcrystal crystal;
	mtzin.import_hkl_info(hklInfo);
	
	HKL_data<F_sigF> fo;
	HKL_data<Flag> free;

//	mtzin.import_crystal(crystal, "/*/*/*");
	mtzin.import_hkl_data(fo,
		(boost::format(kBasePath) % "*" % "*" % ba::join(foLabels, ",")).str());
	mtzin.import_hkl_data(free,
		(boost::format(kBasePath) % "*" % "*" % ba::join(freeLabels, ",")).str());

	mtzin.close_read();
	
	mCell = hklInfo.cell();
	mSpacegroup = hklInfo.spacegroup();

	// The calculation work
	vector<clipper::Atom> atoms;

	for (auto a: structure.atoms())
		atoms.push_back(a.toClipper());

	HKL_data<F_phi> fc(hklInfo, mCell);


	if (not electronScattering)
	{
		auto& exptl = structure.getFile().data()["exptl"];
		electronScattering = not exptl.empty() and exptl.front()["method"] == "ELECTRON CRYSTALLOGRAPHY";
	}

	clipper::ScatteringFactors::selectScattteringFactorsType(
		electronScattering ? clipper::SF_ELECTRON : clipper::SF_WAASMAIER_KIRFEL);
		
	if (noBulk)
	{
		clipper::SFcalc_aniso_fft<float> sfc;
		sfc(fc, atoms);
	}
	else
	{
		clipper::SFcalc_obs_bulk<float> sfcb;
		sfcb(fc, fo, atoms);
		
		if (VERBOSE)
			cerr << "Bulk correction volume: " << sfcb.bulk_frac() << endl
				 << "Bulk correction factor: " << sfcb.bulk_scale() << endl;
	}
	
	if (anisoScaling != as_None)
	{
		clipper::SFscale_aniso<float>::TYPE F = clipper::SFscale_aniso<float>::F;
		clipper::SFscale_aniso<float> sfscl;
		if (anisoScaling == as_Observed)
			sfscl(fo, fc);  // scale Fobs
		else
			sfscl(fc, fo);  // scale Fcal
			
		if (VERBOSE)
			cerr << "Anisotropic scaling:" << endl
				 << sfscl.u_aniso_orth(F).format() << endl;
	}

	// now do sigmaa calc
	HKL_data<F_phi>   fb(hklInfo, mCell), fd(hklInfo, mCell);
	HKL_data<Phi_fom> phiw(hklInfo, mCell);
	HKL_data<Flag>    flag(hklInfo, mCell);

	const int freeflag = 0;
	for (auto ih = flag.first(); not ih.last(); ih.next())
	{
		if (not fo[ih].missing() and (free[ih].missing() or free[ih].flag() == freeflag))
			flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
		else
			flag[ih].flag() = clipper::SFweight_spline<float>::NONE;
	}
	
	int nRefln = 1000;
//	if (vm.count("num-reflns"))
//		nRefln = vm["num-reflns"].as<int>();
	
	int nParam = 20;
//	if (vm.count("num-params"))
//		nParam = vm["num-param"].as<int>();

	// do sigmaa calc
	clipper::SFweight_spline<float> sfw(nRefln, nParam);
	sfw(fb, fd, phiw, fo, fc, flag);

	// fb now contains 2mFo - DFc
	// fd now contains  mFo - DFc

	fixMTZ(fb, fd, fo, fc);

	ResolutionCalculator rc(mCell);
	mResHigh = 99; mResLow = 0;
	
	for (auto hi = fo.first_data(); not hi.last(); hi = fo.next_data(hi))
	{
		float res = rc(hi.hkl().h(), hi.hkl().k(), hi.hkl().l());
		
		if (mResHigh > res)
			mResHigh = res;

		if (mResLow < res)
			mResLow = res;
	}

	if (VERBOSE > 1)
		cerr << "calculated reshi = " << mResHigh << " reslo = " << mResLow << endl;

	samplingRate /= 2;
	
	mGrid = Grid_sampling(mSpacegroup, mCell,
		hklInfo.resolution(), samplingRate);		// define grid
	
	mFb = Xmap(mSpacegroup, mCell, mGrid);			// define map
	mFb.fft_from(fb);								// generate map
	
	mFd = Xmap(mSpacegroup, mCell, mGrid);			// define map
	mFd.fft_from(fd);								// generate map

	if (VERBOSE)
		cerr << "Read Xmaps with sampling rate: " << samplingRate << endl
			 << "  resolution: " << mResHigh
			 << endl;
//
//#if DEBUG
//	char tmpFoFileName[] = "/tmp/fo-XXXXXX.map";
//	if (mkstemps(tmpFoFileName, 4) < 0)
//		throw runtime_error(string("Could not create temp file for map: ") + strerror(errno));
//	
//	CCP4MAPfile foFile;
//	foFile.open_write(tmpFoFileName);
//	foFile.export_xmap(mFb);
//	foFile.close_write();
//	
//	char tmpFcFileName[] = "/tmp/df-XXXXXX.map";
//	if (mkstemps(tmpFcFileName, 4) < 0)
//		throw runtime_error(string("Could not create temp file for map: ") + strerror(errno));
//	
//	CCP4MAPfile fcFile;
//	fcFile.open_write(tmpFcFileName);
//	fcFile.export_xmap(mFd);
//	fcFile.close_write();
//#endif
}

template<typename FTYPE>
void MapMaker<FTYPE>::loadFromMapFiles(const fs::path& fbMapFile, const fs::path& fdMapFile)
{
	using clipper::CCP4MAPfile;

	CCP4MAPfile fbFile;
	fbFile.open_read(fbMapFile.string());
	fbFile.import_xmap(mFb);
	mGrid = fbFile.grid_sampling();
	mCell = fbFile.cell();
	mSpacegroup = fbFile.spacegroup();
	fbFile.close_read();
	
	CCP4MAPfile fdFile;
	fdFile.open_read(fdMapFile.string());
	fdFile.import_xmap(mFd);
	fdFile.close_read();
}

template class MapMaker<float>;
template class MapMaker<double>;

}
