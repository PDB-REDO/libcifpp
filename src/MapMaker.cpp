#include "cif++/Config.h"

#include <iomanip>

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
		cerr << "Reading map from " << mtzFile << endl
			 << "  with labels: FB: " << ba::join(fbLabels, ",") << endl
			 << "  with labels: FD: " << ba::join(fdLabels, ",") << endl
			 << "  with labels: FO: " << ba::join(foLabels, ",") << endl
			 << "  with labels: FC: " << ba::join(fcLabels, ",") << endl;
	
	const string kBasePath("/%1%/%2%/[%3%]");

	using clipper::HKL_info;
	using clipper::CCP4MTZfile;
	using clipper::HKL_data;
	using clipper::data32::F_phi;
	using clipper::data32::F_sigF;
	using clipper::data32::Phi_fom;

	CCP4MTZfile mtzin;
	mtzin.open_read(mtzFile.string());
	
	HKL_info hklInfo;
	mtzin.import_hkl_info(hklInfo);
	
	HKL_data<F_phi> fbData, fdData, fcData;
	HKL_data<F_sigF> foData;
	HKL_data<Phi_fom> fomData;

	mtzin.import_hkl_data(fbData,
		(boost::format(kBasePath) % "*" % "*" % ba::join(fbLabels, ",")).str());
	mtzin.import_hkl_data(fdData,
		(boost::format(kBasePath) % "*" % "*" % ba::join(fdLabels, ",")).str());
	mtzin.import_hkl_data(foData,
		(boost::format(kBasePath) % "*" % "*" % ba::join(foLabels, ",")).str());
	mtzin.import_hkl_data(fcData,
		(boost::format(kBasePath) % "*" % "*" % ba::join(fcLabels, ",")).str());
	mtzin.import_hkl_data(fomData,
		(boost::format(kBasePath) % "*" % "*" % "PHWT,FOM").str());

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

	fixMTZ(fbData, fdData, foData, fcData, fomData);

	samplingRate /= 2;	// clipper's way of interpreting?
	
	mGrid.init(mSpacegroup, mCell,
		hklInfo.resolution(), samplingRate);	// define grid
	
	mFb.init(mSpacegroup, mCell, mGrid);	// define map
	mFb.fft_from(fbData);					// generate map
	
	mFd.init(mSpacegroup, mCell, mGrid);	// define map
	mFd.fft_from(fdData);					// generate map

	if (VERBOSE)
		cerr << "Read Xmaps with sampling rate: " << samplingRate << endl
		     << "  stored resolution: " << hklInfo.resolution().limit() << endl
			 << "  calculated reshi = " << mResHigh << " reslo = " << mResLow << endl
		     << "  spacegroup: " << mSpacegroup.symbol_hm() << endl
		     << "  cell: " << mCell.format() << endl;

#if DEBUG
	char tmpFoFileName[] = "/tmp/fo-XXXXXX.map";
	if (mkstemps(tmpFoFileName, 4) < 0)
		throw runtime_error(string("Could not create temp file for map: ") + strerror(errno));
	
	using clipper::CCP4MAPfile;

	CCP4MAPfile foFile;
	foFile.open_write(tmpFoFileName);
	foFile.export_xmap(mFb);
	foFile.close_write();
	
	char tmpFcFileName[] = "/tmp/df-XXXXXX.map";
	if (mkstemps(tmpFcFileName, 4) < 0)
		throw runtime_error(string("Could not create temp file for map: ") + strerror(errno));
	
	CCP4MAPfile fcFile;
	fcFile.open_write(tmpFcFileName);
	fcFile.export_xmap(mFd);
	fcFile.close_write();
	
	cout << "Wrote fo map to: " << tmpFoFileName << endl
		 << "  and df map to: " << tmpFcFileName << endl;
#endif
}

ostream& operator<<(ostream& os, const clipper::HKL& hkl)
{
	os << "h: " << hkl.h() << ", "
	   << "k: " << hkl.k() << ", "
	   << "l: " << hkl.l();
	
	return os;
};


template<typename FTYPE>
void MapMaker<FTYPE>::fixMTZ(FPdata& fb, FPdata& fd, FOdata& fo, FPdata& fc, WData& fom)
{
	enum {
		A1,  // A1:  FC = 2mFo - FM              
		A2,  // A2:  FC >= 2mFo - FM             
		A3,  // A3:  FD = FM - mFo               
		A4,  // A4:  FD = 2(FM - mFo)            
		C5,  // C5:  FC = 2mFo - FM              
		C6,  // C6:  FM = mFo                    
		C7,  // C7:  FD = mFo - FC               
		C8,  // C8:  FD = 2(mFo - FC)            
		C9,  // C9:  FD <= mFo - FC              
		T10, // 10:  FM = FC (unobserved only)   
		T11, // 11:  FD = 0 (unobserved only)    
		TestCount
	};
	vector<bool> tests(TestCount, true);
	
	// first run the tests to see if we need to fix anything
	
	if (VERBOSE)
		cerr << "Testing MTZ file" << endl;
	
	for (auto ih = fb.first(); not ih.last(); ih.next())
	{
		clipper::HKL_class cls(mSpacegroup, ih.hkl());

		auto W = fom[ih].fom();

		auto FM = fb[ih].f();
		auto PM = fb[ih].phi() * 180 / kPI;
		auto FD = fd[ih].f();
		auto PD = fd[ih].phi() * 180 / kPI;
		auto FO = fo[ih].f();
		auto FC = fc[ih].f();
		auto PC = fc[ih].phi() * 180 / kPI;

		auto WFO = W * FO;

		if (abs(fmod(abs(PM - PC) + 180, 360) - 180) > 90)
			FM = -FM;

		if (abs(fmod(abs(PD - PC) + 180, 360) - 180) > 90)
			FD = -FD; 
			
		if (fo[ih].missing() or W == 0)
		{
			if (tests[T10] and abs(FM - FC) > 0.05)
			{
				tests[T10] = false;
				if (VERBOSE) cerr << "Test 10 failed at " << ih.hkl() << endl;
			}
			
			if (tests[T11] and abs(FD) > 0.05)
			{
				tests[T11] = false;
				if (VERBOSE) cerr << "Test 11 failed at " << ih.hkl() << endl;
			}
		}
		else if (cls.centric())
		{
			if (tests[C5] and abs(FC + FM - 2 * WFO) > 0.05)
			{
				tests[C5] = false;
				if (VERBOSE) cerr << "Test C5 failed at " << ih.hkl() << endl;
			}
			
			if (tests[C6] and abs(FM - WFO) > 0.05)
			{
				tests[C6] = false;
				if (VERBOSE) cerr << "Test C6 failed at " << ih.hkl() << endl;
			}
			
			if (tests[C7] and abs(FC + FD - WFO) > 0.05)
			{
				tests[C7] = false;
				if (VERBOSE) cerr << "Test C7 failed at " << ih.hkl() << endl;
			}
			
			if (tests[C8] and abs(FC + 0.5 * FD - WFO) > 0.05)
			{
				tests[C8] = false;
				if (VERBOSE) cerr << "Test C8 failed at " << ih.hkl() << endl;
			}
			
			if (tests[C9] and (1.01 * FC + Gd - WFO) < -0.05)
			{
				tests[C9] = false;
				if (VERBOSE) cerr << "Test C9 failed at " << ih.hkl() << endl;
			}
			
		}
		else
		{
			if (tests[A1] and abs(FC + FM - 2 * WFO) > 0.05)
			{
				tests[A1] = false;
				if (VERBOSE) cerr << "Test A1 failed at " << ih.hkl() << endl;
			}
			
			if (tests[A2] and 1.01 * FC + FM - 2 * WFO < -0.05)
			{
				tests[A2] = false;
				if (VERBOSE) cerr << "Test A2 failed at " << ih.hkl() << endl;
			}
			
			if (tests[A3] and abs(FM - FD - WFO) > 0.05)
			{
				tests[A3] = false;
				if (VERBOSE) cerr << "Test A3 failed at " << ih.hkl() << endl;
			}
			
			if (tests[A4] and abs(FM - 0.5 * FD - WFO) > 0.05)
			{
				tests[A4] = false;
				if (VERBOSE) cerr << "Test A4 failed at " << ih.hkl() << endl;
			}
		}
	}	

	using clipper::HKL_class;
	using clipper::data32::F_phi;

	const F_phi fzero(0, 0);
	
	// mtzfix...
	for (auto ih = fb.first(); not ih.last(); ih.next())
	{
		if (fb[ih].missing() or fd[ih].missing())
			continue;

		auto PM = fb[ih].phi() * 180 / kPI;
		auto PD = fd[ih].phi() * 180 / kPI;
		auto PC = fc[ih].phi() * 180 / kPI;

		if (abs(fmod(abs(PM - PC) + 180, 360) - 180) > 90)
		{
			fb[ih].f() = -fb[ih].f();
			fb[ih].phi() = fc[ih].phi();
		}

		if (abs(fmod(abs(PD - PC) + 180, 360) - 180) > 90)
		{
			fd[ih].f() = -fd[ih].f();
			fd[ih].phi() = fc[ih].phi();
		}
		
		auto mFo = fb[ih] - fd[ih];

		HKL_class cls(mSpacegroup, ih.hkl());

		if (not fo[ih].missing() and fom[ih].fom() > 0)
		{
			if (cls.centric())
			{
				if (not tests[C6])
					fb[ih] = mFo;
				if (not tests[C7] and tests[C8])
					fd[ih].f() = fd[ih].f() / 2;
			}
			else
			{
				if (tests[A3] and not tests[A4])
					fd[ih] = fd[ih] + fd[ih];
			}
		}
		else
		{
			if (not tests[T10])
			{
				if ((not cls.centric() and tests[A1]) or
					(cls.centric() and (tests[C5] or tests[C7] or tests[C8])))
				{
					fb[ih] = fc[ih];
				}
			}
			
			if (not tests[T11])
				fd[ih] = fzero;
		}
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

	fixMTZ(fb, fd, fo, fc, phiw);

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
	
	mGrid.init(mSpacegroup, mCell,
		hklInfo.resolution(), samplingRate);		// define grid
	
	mFb.init(mSpacegroup, mCell, mGrid);			// define map
	mFb.fft_from(fb);								// generate map
	
	mFd.init(mSpacegroup, mCell, mGrid);			// define map
	mFd.fft_from(fd);								// generate map

	if (VERBOSE)
		cerr << "Read Xmaps with sampling rate: " << samplingRate << endl
			 << "  resolution: " << mResHigh
			 << endl;

#if DEBUG
	char tmpFoFileName[] = "/tmp/fo-XXXXXX.map";
	if (mkstemps(tmpFoFileName, 4) < 0)
		throw runtime_error(string("Could not create temp file for map: ") + strerror(errno));
	
	using clipper::CCP4MAPfile;

	CCP4MAPfile foFile;
	foFile.open_write(tmpFoFileName);
	foFile.export_xmap(mFb);
	foFile.close_write();
	
	char tmpFcFileName[] = "/tmp/df-XXXXXX.map";
	if (mkstemps(tmpFcFileName, 4) < 0)
		throw runtime_error(string("Could not create temp file for map: ") + strerror(errno));
	
	CCP4MAPfile fcFile;
	fcFile.open_write(tmpFcFileName);
	fcFile.export_xmap(mFd);
	fcFile.close_write();
	
	cout << "Wrote fo map to: " << tmpFoFileName << endl
		 << "  and df map to: " << tmpFcFileName << endl;
#endif
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
