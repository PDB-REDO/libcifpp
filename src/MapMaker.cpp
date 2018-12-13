#include "cif++/Config.h"

#include <iomanip>

#include <boost/format.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>

#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>

#include "cif++/MapMaker.h"
#include "cif++/ResolutionCalculator.h"

using namespace std;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

extern int VERBOSE;

namespace mmcif
{
	
// --------------------------------------------------------------------

bool IsMTZFile(const fs::path& p)
{
	bool result = false;
	
	fs::ifstream f(p);
	if (f.is_open())
	{
		char sig[5] = {};
		f.read(sig, 4);
		result = sig == string("MTZ ");
	}

	return result;
}

// --------------------------------------------------------------------

template<typename FTYPE>
Map<FTYPE>::Map()
{
}

template<typename FTYPE>
Map<FTYPE>::~Map()
{
}

template<typename FTYPE>
void Map<FTYPE>::calculateStats()
{
	double sum = 0, sum2 = 0;
	int count = 0;
	
	for (auto ix = mMap.first(); not ix.last(); ix.next())
	{
		auto v = mMap[ix];
		
		if (isnan(v))
			throw runtime_error("map contains NaN values");
		
		++count;
		sum += v;
		sum2 += v * v;
	}
	
	mMeanDensity = sum / count;
	mRMSDensity = sqrt((sum2 / count) - (mMeanDensity * mMeanDensity));
}

template<typename FTYPE>
void Map<FTYPE>::read(const fs::path& mapFile)
{
	fs::path dataFile = mapFile;

	if (VERBOSE)
		cout << "Reading map from " << mapFile << endl;
	
	if (mapFile.extension() == ".gz" or mapFile.extension() == ".bz2")
	{
		// file is compressed
		
		fs::path p = mapFile.parent_path();
		string s = mapFile.filename().string();
		
		io::filtering_stream<io::input> in;
		
		fs::ifstream fi(mapFile);
		
		if (mapFile.extension() == ".gz")
			in.push(io::gzip_decompressor());
		else
			in.push(io::bzip2_decompressor());
			
		in.push(fi);

		char tmpFileName[] = "/tmp/map-tmp-XXXXXX";
		if (mkstemp(tmpFileName) < 0)
			throw runtime_error(string("Could not create temp file for map: ") + strerror(errno));
		
		dataFile = fs::path(tmpFileName);
		fs::ofstream out(dataFile);
		io::copy(in, out);
	}
	
	if (not fs::exists(dataFile))
		throw runtime_error("Could not open map file " + mapFile.string());
	
	using namespace clipper;

	CCP4MAPfile mapin;
	mapin.open_read(dataFile.string());
	mapin.import_xmap(mMap);
	mapin.close_read();
	
	if (dataFile != mapFile)
		fs::remove(dataFile);

	calculateStats();
}

template<typename FTYPE>
void Map<FTYPE>::write(const fs::path& f)
{
	assert(false);
}

template class Map<float>;
template class Map<double>;

// --------------------------------------------------------------------

template<typename FTYPE>
MapMaker<FTYPE>::MapMaker()
{
}

template<typename FTYPE>
MapMaker<FTYPE>::~MapMaker()
{
}

template<typename FTYPE>
void MapMaker<FTYPE>::loadMTZ(const fs::path& hklin, float samplingRate,
	initializer_list<string> fbLabels, initializer_list<string> fdLabels,
	initializer_list<string> foLabels, initializer_list<string> fcLabels)
{
	if (VERBOSE)
		cerr << "Reading map from " << hklin << endl
			 << "  with labels: FB: " << ba::join(fbLabels, ",") << endl
			 << "  with labels: FD: " << ba::join(fdLabels, ",") << endl
			 << "  with labels: FO: " << ba::join(foLabels, ",") << endl
			 << "  with labels: FC: " << ba::join(fcLabels, ",") << endl;

	fs::path dataFile = hklin;
	
	if (hklin.extension() == ".gz" or hklin.extension() == ".bz2")
	{
		// file is compressed
		
		fs::path p = hklin.parent_path();
		string s = hklin.filename().string();
		
		io::filtering_stream<io::input> in;
		
		fs::ifstream fi(hklin);
		
		if (hklin.extension() == ".gz")
			in.push(io::gzip_decompressor());
		else
			in.push(io::bzip2_decompressor());
			
		in.push(fi);

		char tmpFileName[] = "/tmp/mtz-tmp-XXXXXX";
		if (mkstemp(tmpFileName) < 0)
			throw runtime_error(string("Could not create temp file for mtz: ") + strerror(errno));
		
		dataFile = fs::path(tmpFileName);
		fs::ofstream out(dataFile);
		io::copy(in, out);
	}
	
	if (not fs::exists(dataFile))
		throw runtime_error("Could not open mtz file " + hklin.string());
	
	const string kBasePath("/%1%/%2%/[%3%]");

	using clipper::CCP4MTZfile;

	CCP4MTZfile mtzin;
	mtzin.open_read(dataFile.string());
	
	mtzin.import_hkl_info(mHKLInfo);
	
	mtzin.import_hkl_data(mFbData,
		(boost::format(kBasePath) % "*" % "*" % ba::join(fbLabels, ",")).str());
	mtzin.import_hkl_data(mFdData,
		(boost::format(kBasePath) % "*" % "*" % ba::join(fdLabels, ",")).str());
	mtzin.import_hkl_data(mFoData,
		(boost::format(kBasePath) % "*" % "*" % ba::join(foLabels, ",")).str());
	mtzin.import_hkl_data(mFcData,
		(boost::format(kBasePath) % "*" % "*" % ba::join(fcLabels, ",")).str());
	mtzin.import_hkl_data(mFreeData,
		(boost::format(kBasePath) % "*" % "*" % "FREE").str());
	mtzin.import_hkl_data(mPhiFomData,
		(boost::format(kBasePath) % "*" % "*" % "PHWT,FOM").str());

	mtzin.close_read();

	if (dataFile != hklin)
		fs::remove(dataFile);
	
	Cell cell = mHKLInfo.cell();
	Spacegroup spacegroup = mHKLInfo.spacegroup();

	ResolutionCalculator rc(cell);
	mResHigh = 99;
	mResLow = 0;
	
	for (auto hi = mFoData.first_data(); not hi.last(); hi = mFoData.next_data(hi))
	{
		float res = rc(hi.hkl().h(), hi.hkl().k(), hi.hkl().l());
		
		if (mResHigh > res)
			mResHigh = res;

		if (mResLow < res)
			mResLow = res;
	}

//	fixMTZ();

	mGrid.init(spacegroup, cell,
		mHKLInfo.resolution(), samplingRate);	// define grid
	
	clipper::Xmap<FTYPE>& fbMap = mFb;
	clipper::Xmap<FTYPE>& fdMap = mFd;
	
	fbMap.init(spacegroup, cell, mGrid);	// define map
	fbMap.fft_from(mFbData);				// generate map
	
	fdMap.init(spacegroup, cell, mGrid);	// define map
	fdMap.fft_from(mFdData);				// generate map

	if (VERBOSE)
	{
		cerr << "Read Xmaps with sampling rate: " << samplingRate << endl
		     << "  stored resolution: " << mHKLInfo.resolution().limit() << endl
			 << "  calculated reshi = " << mResHigh << " reslo = " << mResLow << endl
		     << "  spacegroup: " << spacegroup.symbol_hm() << endl
		     << "  cell: " << cell.format() << endl
		     << "  grid: " << mGrid.format() << endl;

		printStats();
	}

	mFb.calculateStats();
	mFd.calculateStats();
}

// --------------------------------------------------------------------

template<typename FTYPE>
void MapMaker<FTYPE>::loadMaps(
	const fs::path& fbMapFile, const fs::path& fdMapFile, float reshi, float reslo)
{
	mResHigh = reshi;
	mResLow = reslo;
	
	mFb.read(fbMapFile);
	mFd.read(fdMapFile);
	
	if (not mFb.cell().equals(mFd.cell()))
		throw runtime_error("Fb and Fd map do not contain the same cell");
}

// --------------------------------------------------------------------

ostream& operator<<(ostream& os, const clipper::HKL& hkl)
{
	os << "h: " << hkl.h() << ", "
	   << "k: " << hkl.k() << ", "
	   << "l: " << hkl.l();
	
	return os;
};

// --------------------------------------------------------------------

template<typename FTYPE>
void MapMaker<FTYPE>::calculate(const fs::path& hklin,
	const Structure& structure, bool noBulk, AnisoScalingFlag anisoScaling,
	float samplingRate, bool electronScattering,
	initializer_list<std::string> foLabels, initializer_list<std::string> freeLabels)
{
	if (IsMTZFile(hklin))
		loadFoFreeFromMTZFile(hklin, foLabels, freeLabels);
	else
		loadFoFreeFromReflectionsFile(hklin);
	
	recalc(structure, noBulk, anisoScaling, samplingRate, electronScattering);
}

// --------------------------------------------------------------------
	
template<typename FTYPE>
void MapMaker<FTYPE>::loadFoFreeFromReflectionsFile(const fs::path& hklin)
{
	using clipper::HKL;
	
	cif::File reflnsFile(hklin);
	auto& reflns = reflnsFile.firstDatablock();
	
//	m_xname = reflns["exptl_crystal"].front()["id"].as<string>();
//	m_pname = reflns["entry"].front()["id"].as<string>();

	float a, b, c, alpha, beta, gamma;
	cif::tie(a, b, c, alpha, beta, gamma) = reflns["cell"].front().get(
		"length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma");
	
	using clipper::Cell_descr;
	Cell cell = Cell(Cell_descr{a, b, c, alpha, beta, gamma});

//	if (not cell2.equals(m_cell))
//		throw runtime_error("Reflections file and coordinates file do not agree upon the cell parameters");

	// --------------------------------------------------------------------

	// Read reflections file to calculate resolution low and high
	ResolutionCalculator rc(a, b, c, alpha, beta, gamma);
	double hires = 99;
	
	for (auto r: reflns["refln"])
	{
		int h, k, l;
		
		cif::tie(h, k, l) = r.get("index_h", "index_k", "index_l");
		
		double res = rc(h, k, l);
		
		if (hires > res)
			hires = res;
	}
	
	string spacegroupDescr = reflns["symmetry"].front()["space_group_name_H-M"].as<string>();
	auto spacegroup = Spacegroup(clipper::Spgr_descr{spacegroupDescr});
	mHKLInfo = HKL_info(spacegroup, cell, clipper::Resolution{hires}, true);
	
//	m_crystal = MTZcrystal(m_xname, m_pname, m_cell);

	mFoData.init(mHKLInfo, mHKLInfo.cell());
	mFreeData.init(mHKLInfo, mHKLInfo.cell());

	for (auto ih = mFreeData.first(); not ih.last(); ih.next())
		mFreeData[ih].set_null();
	
	// --------------------------------------------------------------------

	enum FreeRConvention { frXPLO, frCCP4 } freeRConvention = frXPLO;
	int freeRefl = 1, workRefl = 0;
	
	if (false /*m_statusXPLO*/)
	{
		freeRConvention = frCCP4;
		freeRefl = 0;
		workRefl = 1;
	}

	bool first = false;
	for (auto r: reflns["refln"])
	{
		int h, k, l;
		char flag;
		float F, sigF;
		
		cif::tie(h, k, l, flag, F, sigF) = r.get("index_h", "index_k", "index_l", "status", "F_meas_au", "F_meas_sigma_au");

		int ix = mHKLInfo.index_of(HKL{h, k, l});
		
		if (ix < 0)
		{
			if (VERBOSE)
				cerr << "Ignoring hkl(" << h << ", " << k << ", " << l << ")" << endl;
			continue;
		}
		
		if (first and (flag == freeRefl or flag == workRefl))
		{
			cerr << "Non-standard _refln.status column detected" << endl
				 << "Assuming " << (freeRConvention == frXPLO ? "XPLOR" : "CCP4") << " convention for free R flag" << endl;
			first = false;
		}
		
		mFoData[ix] = F_sigF(F, sigF);
		
		switch (flag)
		{
			case 'o':
			case 'h':
			case 'l':
				mFreeData[ix] = Flag(1);
				break;
			
			case 'f':
				mFreeData[ix] = Flag(0);
				break;
			
			case '0':
			case '1':
				mFreeData[ix] = Flag(workRefl == flag ? 1 : 0);
				break;

			default:
				if (VERBOSE > 1)
					cerr << "Unexpected value in status: '" << flag << "' for hkl(" << h << ", " << k << ", " << l << ")" << endl;
				break;
		}
	}
}

// --------------------------------------------------------------------

template<typename FTYPE>
void MapMaker<FTYPE>::loadFoFreeFromMTZFile(const fs::path& hklin,
	initializer_list<std::string> foLabels, initializer_list<std::string> freeLabels)
{
	if (VERBOSE)
		cerr << "Recalculating maps from " << hklin << endl;
	
	const string kBasePath("/%1%/%2%/[%3%]");

	using clipper::CCP4MTZfile;

	CCP4MTZfile mtzin;
	mtzin.open_read(hklin.string());
	
	mtzin.import_hkl_info(mHKLInfo);
	mtzin.import_hkl_data(mFoData,
		(boost::format(kBasePath) % "*" % "*" % ba::join(foLabels, ",")).str());
	mtzin.import_hkl_data(mFreeData,
		(boost::format(kBasePath) % "*" % "*" % ba::join(freeLabels, ",")).str());

	mtzin.close_read();
}

// --------------------------------------------------------------------

template<typename FTYPE>
void MapMaker<FTYPE>::recalc(const Structure& structure,
		bool noBulk, AnisoScalingFlag anisoScaling,
		float samplingRate, bool electronScattering)
{
	Cell cell = mHKLInfo.cell();
	Spacegroup spacegroup = mHKLInfo.spacegroup();

	// The calculation work
	vector<clipper::Atom> atoms;

	for (auto a: structure.atoms())
		atoms.push_back(a.toClipper());

	mFcData.init(mHKLInfo, cell);

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
		sfc(mFcData, atoms);
	}
	else
	{
		clipper::SFcalc_obs_bulk<float> sfcb;
		sfcb(mFcData, mFoData, atoms);
		
		if (VERBOSE)
			cerr << "Bulk correction volume: " << sfcb.bulk_frac() << endl
				 << "Bulk correction factor: " << sfcb.bulk_scale() << endl;
	}
	
	if (anisoScaling != as_None)
	{
		clipper::SFscale_aniso<float>::TYPE F = clipper::SFscale_aniso<float>::F;
		clipper::SFscale_aniso<float> sfscl;
		if (anisoScaling == as_Observed)
			sfscl(mFoData, mFcData);  // scale Fobs
		else
			sfscl(mFcData, mFoData);  // scale Fcal
			
		if (VERBOSE)
			cerr << "Anisotropic scaling:" << endl
				 << sfscl.u_aniso_orth(F).format() << endl;
	}

	// now do sigmaa calc
	mFbData.init(mHKLInfo, cell);
	mFdData.init(mHKLInfo, cell);
	mPhiFomData.init(mHKLInfo, cell);
	
	HKL_data<Flag> flag(mHKLInfo, cell);

	const int freeflag = 0;
	for (auto ih = mFreeData.first(); not ih.last(); ih.next())
	{
		if (not mFoData[ih].missing() and (mFreeData[ih].missing() or mFreeData[ih].flag() == freeflag))
			flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
		else
			flag[ih].flag() = clipper::SFweight_spline<float>::NONE;
	}
	
	// do sigmaa calc
	clipper::SFweight_spline<float> sfw(mNumRefln, mNumParam);
	sfw(mFbData, mFdData, mPhiFomData, mFoData, mFcData, flag);

	// mFbData now contains 2mFo - DFc
	// mFdData now contains  mFo - DFc

	fixMTZ();

	ResolutionCalculator rc(cell);
	mResHigh = 99; mResLow = 0;
	
	for (auto hi = mFoData.first_data(); not hi.last(); hi = mFoData.next_data(hi))
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
	
	mGrid.init(spacegroup, cell,
		mHKLInfo.resolution(), samplingRate);		// define grid
	
	clipper::Xmap<FTYPE>& fbMap = mFb;
	clipper::Xmap<FTYPE>& fdMap = mFd;

	fbMap.init(spacegroup, cell, mGrid);			// define map
	fbMap.fft_from(mFbData);								// generate map
	
	fdMap.init(spacegroup, cell, mGrid);			// define map
	fdMap.fft_from(mFdData);								// generate map

	if (VERBOSE)
	{
		cerr << "Read Xmaps with sampling rate: " << samplingRate << endl
			 << "  resolution: " << mResHigh
			 << endl;

		printStats();
	}

	mFb.calculateStats();
	mFd.calculateStats();
}

template<typename FTYPE>
void MapMaker<FTYPE>::fixMTZ()
{
	Spacegroup spacegroup = mHKLInfo.spacegroup();

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
	
	for (auto ih = mFbData.first(); not ih.last(); ih.next())
	{
		clipper::HKL_class cls(spacegroup, ih.hkl());

		auto W = mPhiFomData[ih].fom();

		auto FM = mFbData[ih].f();
		auto PM = mFbData[ih].phi() * 180 / kPI;
		auto FD = mFdData[ih].f();
		auto PD = mFdData[ih].phi() * 180 / kPI;
		auto FO = mFoData[ih].f();
		auto FC = mFcData[ih].f();
		auto PC = mFcData[ih].phi() * 180 / kPI;

		auto WFO = W * FO;

		if (abs(fmod(abs(PM - PC) + 180, 360) - 180) > 90)
			FM = -FM;

		if (abs(fmod(abs(PD - PC) + 180, 360) - 180) > 90)
			FD = -FD; 
			
		if (mFoData[ih].missing() or W == 0)
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
	for (auto ih = mFbData.first(); not ih.last(); ih.next())
	{
		if (mFbData[ih].missing() or mFdData[ih].missing())
			continue;

		auto PM = mFbData[ih].phi() * 180 / kPI;
		auto PD = mFdData[ih].phi() * 180 / kPI;
		auto PC = mFcData[ih].phi() * 180 / kPI;

		if (abs(fmod(abs(PM - PC) + 180, 360) - 180) > 90)
		{
			mFbData[ih].f() = -mFbData[ih].f();
			mFbData[ih].phi() = mFcData[ih].phi();
		}

		if (abs(fmod(abs(PD - PC) + 180, 360) - 180) > 90)
		{
			mFdData[ih].f() = -mFdData[ih].f();
			mFdData[ih].phi() = mFcData[ih].phi();
		}
		
		auto mFo = mFbData[ih] - mFdData[ih];

		HKL_class cls(spacegroup, ih.hkl());

		if (not mFoData[ih].missing() and mPhiFomData[ih].fom() > 0)
		{
			if (cls.centric())
			{
				if (not tests[C6])
					mFbData[ih] = mFo;
				if (not tests[C7] and tests[C8])
					mFdData[ih].f() = mFdData[ih].f() / 2;
			}
			else
			{
				if (tests[A3] and not tests[A4])
					mFdData[ih] = mFdData[ih] + mFdData[ih];
			}
		}
		else
		{
			if (not tests[T10])
			{
				if ((not cls.centric() and tests[A1]) or
					(cls.centric() and (tests[C5] or tests[C7] or tests[C8])))
				{
					mFbData[ih] = mFcData[ih];
				}
			}
			
			if (not tests[T11])
				mFdData[ih] = fzero;
		}
	}
}

template<typename FTYPE>
void MapMaker<FTYPE>::printStats()
{
	// calc R and R-free
	vector<double> params(mNumParam, 1.0);

	clipper::BasisFn_spline basisfn(mFoData, mNumParam, 1.0);
	clipper::TargetFn_scaleF1F2<clipper::data32::F_phi,clipper::data32::F_sigF> targetfn(mFcData, mFoData);
	clipper::ResolutionFn rfn(mHKLInfo, basisfn, targetfn, params);

	double r1w = 0, f1w = 0, r1f = 0, f1f = 0;
	const int freeflag = 0;

	for (auto ih = mFoData.first_data(); not ih.last(); ih = mFoData.next_data(ih))
	{
		if (mFcData[ih].missing())
			continue;
//			throw runtime_error("missing Fc");
		
		double Fo = mFoData[ih].f();
		double Fc = sqrt(rfn.f(ih)) * mFcData[ih].f();

		if (mFreeData[ih].flag() == freeflag)
		{
			r1f += fabs(Fo - Fc);
			f1f += Fo;
		}
		else
		{
			r1w += fabs(Fo - Fc);
			f1w += Fo;
		}
	}
	
	if (f1f < 0.1)
		f1f = 0.1;
	r1f /= f1f;
	
	if (f1w < 0.1)
		f1w = 0.1;
	r1w /= f1w;

	cerr << "R-factor      : " << r1w << endl
		 << "Free R-factor : " << r1f << endl;
}

template<typename FTYPE>
void MapMaker<FTYPE>::writeMTZ(const fs::path& file, const string& pname, const string& cname)
{
	if (mHKLInfo.is_null())
		throw runtime_error("HKL info not initialized");
	
	clipper::CCP4MTZfile mtz;
	clipper::MTZdataset dataset(pname, 0);
	clipper::MTZcrystal crystal(cname, pname, mHKLInfo.cell());

	const string col = "/" + pname + "/" + cname + "/";
	
	mtz.open_write(file.string());
	mtz.export_hkl_info(mHKLInfo);
	mtz.export_crystal(crystal, col);
	mtz.export_dataset(dataset, col);
	if (not mFreeData.is_null())		mtz.export_hkl_data(mFreeData, col + "[FREE]");
	if (not mFoData.is_null())			mtz.export_hkl_data(mFoData, col + "[FP,SIGFP]");
	if (not mFcData.is_null())			mtz.export_hkl_data(mFcData, col + "[FC_ALL,PHIC_ALL]");
	if (not mFbData.is_null())			mtz.export_hkl_data(mFbData, col + "[FWT,PHWT]");
	if (not mFdData.is_null())			mtz.export_hkl_data(mFdData, col + "[DELFWT,PHDELWT]");
	if (not mPhiFomData.is_null())		mtz.export_hkl_data(mPhiFomData, col + "[PHI,FOM]");
	mtz.close_write();
}

template class MapMaker<float>;
template class MapMaker<double>;

}
