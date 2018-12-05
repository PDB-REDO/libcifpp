// Lib for working with structures as contained in file and PDB files

#include "cif++/Structure.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/format.hpp>

#include "cif++/PDB2Cif.h"
#include "cif++/CifParser.h"
#include "cif++/Cif2PDB.h"
#include "cif++/AtomShape.h"

using namespace std;

namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;

extern int VERBOSE;

namespace mmcif
{
	
// --------------------------------------------------------------------
// FileImpl
	
struct FileImpl
{
	cif::File			mData;
	cif::Datablock*		mDb = nullptr;
	
	void load(fs::path p);
	void save(fs::path p);
};

void FileImpl::load(fs::path p)
{
	fs::ifstream inFile(p, ios_base::in | ios_base::binary);
	if (not inFile.is_open())
		throw runtime_error("No such file: " + p.string());
	
	io::filtering_stream<io::input> in;
	string ext = p.extension().string();
	
	if (p.extension() == ".bz2")
	{
		in.push(io::bzip2_decompressor());
		ext = p.stem().extension().string();
	}
	else if (p.extension() == ".gz")
	{
		in.push(io::gzip_decompressor());
		ext = p.stem().extension().string();
	}
	
	in.push(inFile);

	try
	{
		// OK, we've got the file, now create a protein
		if (ext == ".cif")
			mData.load(in);
		else if (ext == ".pdb" or ext == ".ent")
			ReadPDBFile(in, mData);
		else
		{
			try
			{
				if (VERBOSE)
					cerr << "unrecognized file extension, trying cif" << endl;
	
				mData.load(in);
			}
			catch (const cif::CifParserError& e)
			{
				if (VERBOSE)
					cerr << "Not cif, trying plain old PDB" << endl;
	
				// pffft...
				in.reset();
	
				if (inFile.is_open())
					inFile.seekg(0);
				else
					inFile.open(p, ios_base::in | ios::binary);
	
				if (p.extension() == ".bz2")
					in.push(io::bzip2_decompressor());
				else if (p.extension() == ".gz")
					in.push(io::gzip_decompressor());
				
				in.push(inFile);
	
				ReadPDBFile(in, mData);
			}
		}
	}
	catch (const exception& ex)
	{
		cerr << "Failed trying to load file " << p << endl;
		throw;
	}
	
	// Yes, we've parsed the data. Now locate the datablock.
	mDb = &mData.firstDatablock();
	
	// And validate, otherwise lots of functionality won't work
//	if (mData.getValidator() == nullptr)
	mData.loadDictionary("mmcif_pdbx");
	if (not mData.isValid())
		cerr << "Invalid mmCIF file" << (VERBOSE ? "." : " use --verbose option to see errors") << endl;
}

void FileImpl::save(fs::path p)
{
	fs::ofstream outFile(p, ios_base::out | ios_base::binary);
	io::filtering_stream<io::output> out;
	
	if (p.extension() == ".gz")
	{
		out.push(io::gzip_compressor());
		p = p.stem();
	}
	else if (p.extension() == ".bz2")
	{
		out.push(io::bzip2_compressor());
		p = p.stem();
	}
	
	out.push(outFile);
	
	if (p.extension() == ".pdb")
		WritePDBFile(out, mData);
	else
		mData.save(out);
}


// --------------------------------------------------------------------
// Atom

struct AtomImpl
{
	AtomImpl(const File& f, const string& id)
		: mFile(f), mId(id), mRefcount(1), mCompound(nullptr)
	{
		auto& db = *mFile.impl().mDb;
		auto& cat = db["atom_site"];
		
		mRow = cat[cif::Key("id") == mId];

		prefetch();
	}
	
	AtomImpl(const File& f, const string& id, cif::Row row)
		: mFile(f), mId(id), mRefcount(1), mRow(row), mCompound(nullptr)
	{
		prefetch();
	}

	AtomImpl(const AtomImpl& impl, const Point& d, const clipper::RTop_orth& rt)
		: mFile(impl.mFile), mId(impl.mId), mType(impl.mType), mAtomID(impl.mAtomID)
		, mCompID(impl.mCompID), mAsymID(impl.mAsymID), mSeqID(impl.mSeqID)
		, mAltID(impl.mAltID), mLocation(impl.mLocation), mRefcount(impl.mRefcount)
		, mRow(impl.mRow), mCompound(impl.mCompound), mRadius(impl.mRadius)
		, mCachedProperties(impl.mCachedProperties)
		, mSymmetryCopy(true)
	{
		mLocation += d;
		mLocation = ((clipper::Coord_orth)mLocation).transform(rt);
		mLocation -= d;
	}

	void prefetch()
	{
		// Prefetch some data
		string symbol;
		cif::tie(symbol, mAtomID, mCompID, mAsymID, mSeqID, mAltID) =
			mRow.get("type_symbol", "label_atom_id", "label_comp_id", "label_asym_id", "label_seq_id", "label_alt_id");
		
		mType = AtomTypeTraits(symbol).type();

		float x, y, z;
		cif::tie(x, y, z) = mRow.get("Cartn_x", "Cartn_y", "Cartn_z");
		
		mLocation = Point(x, y, z);

		string compId;
		cif::tie(compId) = mRow.get("label_comp_id");
		
		mCompound = Compound::create(compId);
	}
	
	clipper::Atom toClipper() const
	{
		clipper::Atom result;
		result.set_coord_orth(mLocation);
		
		if (mRow["occupancy"].empty())
			result.set_occupancy(1.0);
		else
			result.set_occupancy(mRow["occupancy"].as<float>());
		
		string element = mRow["type_symbol"].as<string>();
		if (not mRow["pdbx_formal_charge"].empty())
		{
			int charge = mRow["pdbx_formal_charge"].as<int>();
			if (abs(charge > 1))
				element += to_string(charge);
			if (charge < 0)
				element += '-';
			else
				element += '+';
		}
		result.set_element(element);
		
		if (not mRow["U_iso_or_equiv"].empty())
			result.set_u_iso(mRow["U_iso_or_equiv"].as<float>());
		else if (not mRow["B_iso_or_equiv"].empty())
			result.set_u_iso(mRow["B_iso_or_equiv"].as<float>() / (8 * kPI * kPI));
		else
			throw runtime_error("Missing B_iso or U_iso");	
		
		auto& db = *mFile.impl().mDb;
		auto& cat = db["atom_site_anisotrop"];
		auto r = cat[cif::Key("id") == mId];
		if (r.empty())
			result.set_u_aniso_orth(clipper::U_aniso_orth(nan("0"), 0, 0, 0, 0, 0));
		else
		{
			float u11, u12, u13, u22, u23, u33;
			cif::tie(u11, u12, u13, u22, u23, u33) =
				r.get("U[1][1]", "U[1][2]", "U[1][3]", "U[2][2]", "U[2][3]", "U[3][3]");
			
			result.set_u_aniso_orth(clipper::U_aniso_orth(u11, u22, u33, u12, u13, u23));
		}
		
		return result;
	}

	void reference()
	{
		++mRefcount;
	}
	
	void release()
	{
		if (--mRefcount < 0)
			delete this;
	}
	
	bool getAnisoU(float anisou[6]) const
	{
		auto& db = *mFile.impl().mDb;
		auto& cat = db["atom_site_anisotrop"];
		auto r = cat[cif::Key("id") == mId];
		bool result = false;

		if (not r.empty())
		{
			result = true;
			cif::tie(anisou[0], anisou[1], anisou[2], anisou[3], anisou[4], anisou[5]) =
				r.get("U[1][1]", "U[1][2]", "U[1][3]", "U[2][2]", "U[2][3]", "U[3][3]");
		}
		
		return result;
	}
	
	void moveTo(const Point& p)
	{
		assert(not mSymmetryCopy);
		if (mSymmetryCopy)
			throw runtime_error("Moving symmetry copy");
		
		boost::format kPosFmt("%.3f");

		mRow["Cartn_x"] = (kPosFmt % p.getX()).str();
		mRow["Cartn_y"] = (kPosFmt % p.getY()).str();
		mRow["Cartn_z"] = (kPosFmt % p.getZ()).str();

		mLocation = p;
	}
	
	const Compound& comp()
	{
		if (mCompound == nullptr)
		{
			string compId;
			cif::tie(compId) = mRow.get("label_comp_id");
			
			mCompound = Compound::create(compId);
			
			if (VERBOSE and mCompound == nullptr)
				cerr << "Compound not found: '" << compId << '\'' << endl;
		}
		
		if (mCompound == nullptr)
			throw runtime_error("no compound");
	
		return *mCompound;
	}
	
	bool isWater() const
	{
		return mCompound != nullptr and mCompound->isWater();
	}

	float radius() const
	{
		return mRadius;
	}

	const string& property(const string& name)
	{
		static string kEmptyString;
		
		auto i = mCachedProperties.find(name);
		if (i == mCachedProperties.end())
		{
			auto v = mRow[name];
			if (v.empty())
				return kEmptyString;
			
			return mCachedProperties[name] = v.as<string>();
		}
		else
			return i->second;
	}
	
	const File&			mFile;
	string				mId;
	AtomType			mType;

	string				mAtomID;
	string				mCompID;
	string				mAsymID;
	int					mSeqID;
	string				mAltID;
	
	Point				mLocation;
	int					mRefcount;
	cif::Row			mRow;
	const Compound*		mCompound;
	float				mRadius = nan("4");
	map<string,string>	mCachedProperties;
	
	bool				mSymmetryCopy = false;
};

//Atom::Atom(const File& f, const string& id)
//	: mImpl(new AtomImpl(f, id))
//{
//}
//
Atom::Atom(AtomImpl* impl)
	: mImpl(impl)
{
}

Atom::Atom(const Atom& rhs)
	: mImpl(rhs.mImpl)
{
	if (mImpl)
		mImpl->reference();
}

Atom::~Atom()
{
	if (mImpl)
		mImpl->release();
}

Atom& Atom::operator=(const Atom& rhs)
{
	if (this != &rhs)
	{
		if (mImpl)
			mImpl->release();
		mImpl = rhs.mImpl;
		if (mImpl)
			mImpl->reference();
	}

	return *this;
}

template<>
string Atom::property<string>(const string& name) const
{
	return mImpl->property(name);
}

template<>
int Atom::property<int>(const string& name) const
{
	auto v = mImpl->property(name);
	return v.empty() ? 0 : stoi(v);
}

template<>
float Atom::property<float>(const string& name) const
{
	return stof(mImpl->property(name));
}

const string& Atom::id() const
{
	return mImpl->mId;
}

AtomType Atom::type() const
{
	return mImpl->mType;
}

int Atom::charge() const
{
	return property<int>("pdbx_formal_charge");
}

string Atom::energyType() const
{
	string result;
	
	if (mImpl and mImpl->mCompound)
		result = mImpl->mCompound->getAtomById(mImpl->mAtomID).typeEnergy;
	
	return result;
}

float Atom::uIso() const
{
	float result;
	
	if (not property<string>("U_iso_or_equiv").empty())
		result =  property<float>("U_iso_or_equiv");
	else if (not property<string>("B_iso_or_equiv").empty())
		result = property<float>("B_iso_or_equiv") / (8 * kPI * kPI);
	else
		throw runtime_error("Missing B_iso or U_iso");
	
	return result;
}

bool Atom::getAnisoU(float anisou[6]) const
{
	return mImpl->getAnisoU(anisou);
}

float Atom::occupancy() const
{
	return property<float>("occupancy");
}

string Atom::labelAtomId() const
{
	return mImpl->mAtomID;
}

string Atom::labelCompId() const
{
	return mImpl->mCompID;
}

string Atom::labelAsymId() const
{
	return mImpl->mAsymID;
}

int Atom::labelSeqId() const
{
	return mImpl->mSeqID;
}

string Atom::authAsymId() const
{
	return property<string>("auth_asym_id");
}

string Atom::authAtomId() const
{
	return property<string>("auth_atom_id");
}

string Atom::authCompId() const
{
	return property<string>("auth_comp_id");
}

int Atom::authSeqId() const
{
	return property<int>("auth_seq_id");
}

string Atom::labelID() const
{
	return property<string>("label_comp_id") + '_' + mImpl->mAsymID + '_' + to_string(mImpl->mSeqID) + ':' + mImpl->mAtomID;
}

string Atom::pdbID() const
{
	return
		property<string>("auth_comp_id") + '_' +
		property<string>("auth_asym_id") + '_' +
		property<string>("auth_seq_id") + 
		property<string>("pdbx_PDB_ins_code");
}

Point Atom::location() const
{
	return mImpl->mLocation;
}

void Atom::location(Point p)
{
	mImpl->moveTo(p);
}

Atom Atom::symmetryCopy(const Point& d, const clipper::RTop_orth& rt)
{
	return Atom(new AtomImpl(*mImpl, d, rt));
}

const Compound& Atom::comp() const
{
	return mImpl->comp();
}

bool Atom::isWater() const
{
	return mImpl->isWater();
}

bool Atom::operator==(const Atom& rhs) const
{
	return mImpl == rhs.mImpl or
		(&mImpl->mFile == &rhs.mImpl->mFile and mImpl->mId == rhs.mImpl->mId); 	
}

const File& Atom::getFile() const
{
	assert(mImpl);
	return mImpl->mFile;
}

clipper::Atom Atom::toClipper() const
{
	return mImpl->toClipper();
}

void Atom::calculateRadius(float resHigh, float resLow, float perc)
{
	AtomShape shape(*this, resHigh, resLow, false);
	mImpl->mRadius = shape.radius();

	// verbose
	if (VERBOSE > 1)
		cout << "Calculated radius for " << AtomTypeTraits(mImpl->mType).name() << " with charge " << charge() << " is " << mImpl->mRadius << endl;
}

float Atom::radius() const
{
	return mImpl->mRadius;
}

// --------------------------------------------------------------------
// residue

Residue::Residue()
	: mStructure(nullptr), mSeqID(0)
{
}

Residue::Residue(const Residue& rhs)
	: mStructure(rhs.mStructure)
	, mCompoundID(rhs.mCompoundID), mAsymID(rhs.mAsymID), mAltID(rhs.mAltID), mSeqID(rhs.mSeqID)
{
}

Residue& Residue::operator=(const Residue& rhs)
{
	if (this != &rhs)
	{
		mStructure = rhs.mStructure;
		mCompoundID = rhs.mCompoundID;
		mAsymID = rhs.mAsymID;
		mSeqID = rhs.mSeqID;
		mAltID = rhs.mAltID;
		
		mAtoms = rhs.mAtoms;
	}
	
	return *this;
}

const Compound& Residue::compound() const
{
	auto result = Compound::create(mCompoundID);
	if (result == nullptr)
		throw runtime_error("Failed to create compound " + mCompoundID);
	return *result;
}

AtomView Residue::atoms() const
{
	if (mStructure == nullptr)
		throw runtime_error("Invalid Residue object");

	if (mAtoms.empty())
	{
		for (auto& a: mStructure->atoms())
		{
			if (mSeqID > 0 and a.labelSeqId() != mSeqID)
				continue;
			
			if (a.labelAsymId() != mAsymID or
				a.labelCompId() != mCompoundID)
					continue;
			
			if (not mAltID.empty() and a.property<string>("label_alt_id") != mAltID)
				continue;
			
			mAtoms.push_back(a);
		}
//
//		auto& atomSites = mStructure->category("atom_site");
//		
//		auto query = cif::Key("label_asym_id") == mAsymID and cif::Key("label_comp_id") == mCompoundID;
//		
//		if (not mSeqID.empty())
//			query = move(query) and cif::Key("label_seq_id") == mSeqID;
//		
//		if (not mAltID.empty())
//			query = move(query) and (cif::Key("label_alt_id") == cif::Empty() or cif::Key("label_alt_id") == mAltID);
//	
//		auto cifAtoms = atomSites.find(move(query));
//
//		set<string> ids;
//		for (auto cifAtom: cifAtoms)
//			ids.insert(cifAtom["id"].as<string>());
//		
//		for (auto& a: mStructure->atoms())
//		{
//			if (ids.count(a.id()))
//				mAtoms.push_back(a);
//		}
	}
	
	return mAtoms;
	
//	auto& atomSites = mStructure->category("atom_site");
//	
//	auto query = cif::Key("label_asym_id") == mAsymID and cif::Key("label_comp_id") == mCompoundID;
//	
//	if (not mSeqID.empty())
//		query = move(query) and cif::Key("label_seq_id") == mSeqID;
//	
//	if (not mAltID.empty())
//		query = move(query) and (cif::Key("label_alt_id") == cif::Empty() or cif::Key("label_alt_id") == mAltID);
//
//	return atomSites.find(move(query));
}

Atom Residue::atomByID(const string& atomID) const
{
	for (auto& a: atoms())
	{
		if (a.labelAtomId() == atomID)
			return a;
	}
	
	throw runtime_error("Atom with atom_id " + atomID + " not found in residue " + mAsymID + ':' + to_string(mSeqID));
}

// Residue is a single entity if the atoms for the asym with mAsymID is equal
// to the number of atoms in this residue...  hope this is correct....
bool Residue::isEntity() const
{
	auto& db = mStructure->datablock();
	
	auto a1 = db["atom_site"].find(cif::Key("label_asym_id") == mAsymID);
	auto a2 = atoms();
	
	return a1.size() == a2.size();
}

string Residue::authID() const
{
	string result;
	
	try
	{
		char chainID, iCode;
		int seqNum;
		
		tie(chainID, seqNum, iCode) = mStructure->MapLabelToAuth(mAsymID, mSeqID);
		
		result = chainID + to_string(seqNum);
		if (iCode != ' ' and iCode != 0)
			result += iCode;
	}
	catch (...)
	{
		result = mAsymID + to_string(mSeqID);
	}
	
	return result;
}

// --------------------------------------------------------------------
// monomer

Monomer::Monomer()
	: mPolymer(nullptr), mIndex(0)
{
}

Monomer::Monomer(const Monomer& rhs)
	: Residue(rhs), mPolymer(rhs.mPolymer), mIndex(rhs.mIndex)
{
}

Monomer& Monomer::operator=(const Monomer& rhs)
{
	if (this != &rhs)
	{
		Residue::operator=(rhs);
		mPolymer = rhs.mPolymer;
		mIndex = rhs.mIndex;
	}
	
	return *this;
}

Monomer::Monomer(const Polymer& polymer, uint32 index)
	: Residue(*polymer.structure())
	, mPolymer(&polymer)
	, mIndex(index)
{
	
}

Monomer::Monomer(const Polymer& polymer, uint32 index, int seqID, const string& compoundID, const string& altID)
	: Residue(*polymer.structure(), compoundID, polymer.asymID(), seqID, altID)
	, mPolymer(&polymer)
	, mIndex(index)
{
}

float Monomer::phi() const
{
	float result = 360;
	if (mIndex > 0)
	{
		Monomer prev = mPolymer->operator[](mIndex - 1);
		if (prev.mSeqID + 1 == mSeqID)
			result = DihedralAngle(prev.C().location(), N().location(), CAlpha().location(), C().location()); 
	}

	return result;
}

float Monomer::psi() const
{
	float result = 360;
	if (mIndex + 1 < mPolymer->size())
	{
		Monomer next = mPolymer->operator[](mIndex + 1);
		if (mSeqID + 1 == next.mSeqID)
			result = DihedralAngle(N().location(), CAlpha().location(), C().location(), next.N().location()); 
	}

	return result;
}

float Monomer::alpha() const
{
	float result = 360;
	
	if (mIndex > 1 and mIndex + 2 < mPolymer->size())
	{
		Monomer prev = mPolymer->operator[](mIndex - 1);
		Monomer next = mPolymer->operator[](mIndex + 1);
		Monomer nextNext = mPolymer->operator[](mIndex + 2);
		
		result = DihedralAngle(prev.CAlpha().location(), CAlpha().location(), next.CAlpha().location(), nextNext.CAlpha().location());
	}
	
	return result;
}

float Monomer::kappa() const
{
	double result = 360;
	
	try
	{
		if (mIndex > 2 and mIndex + 2 < mPolymer->size())
		{
			Monomer prevPrev = mPolymer->operator[](mIndex - 2);
			Monomer nextNext = mPolymer->operator[](mIndex + 2);
			
			if (prevPrev.mSeqID + 4 == nextNext.mSeqID)
			{
				double ckap = CosinusAngle(CAlpha().location(), prevPrev.CAlpha().location(), nextNext.CAlpha().location(), CAlpha().location());
				double skap = sqrt(1 - ckap * ckap);
				result = atan2(skap, ckap) * 180 / kPI;
			}
		}
	}
	catch (const exception& ex)
	{
		if (VERBOSE)
			cerr << "When trying to calculate kappa for " << asymID() << ':' << seqID() << ": "
				 << ex.what() << endl;
	}

	return result;
}

bool Monomer::areBonded(const Monomer& a, const Monomer& b, float errorMargin)
{
	bool result = false;
	
	try
	{
		Point atoms[4] = {
			a.atomByID("CA").location(),
			a.atomByID("C").location(),
			b.atomByID("N").location(),
			b.atomByID("CA").location()
		};
		
		auto distanceCACA = Distance(atoms[0], atoms[3]);
		double omega = DihedralAngle(atoms[0], atoms[1], atoms[2], atoms[3]);

		bool cis = abs(omega) <= 30.0;
		float maxCACADistance = cis ? 3.0 : 3.8;
		
		result = abs(distanceCACA - maxCACADistance) < errorMargin;
	}
	catch (...) {}

	return result;
}

bool Monomer::isCis(const mmcif::Monomer& a, const mmcif::Monomer& b)
{
	bool result = false;
	
	try
	{
		double omega = DihedralAngle(
			a.atomByID("CA").location(),
			a.atomByID("C").location(),
			b.atomByID("N").location(),
			b.atomByID("CA").location());
		result = abs(omega) <= 30.0;
	}
	catch (...) {}

	return result;
}

// --------------------------------------------------------------------
// polymer

Polymer::iterator::iterator(const Polymer& p, uint32 index)
	: mPolymer(&p), mIndex(index), mCurrent(p, index)
{
	auto& polySeq = mPolymer->mPolySeq;

	if (index < polySeq.size())
	{
		int seqID;
		string asymID, monID;
		cif::tie(asymID, seqID, monID) =
			polySeq[mIndex].get("asym_id", "seq_id", "mon_id");
		
		mCurrent = Monomer(*mPolymer, index, seqID, monID, "");
	}
}

Monomer Polymer::operator[](size_t index) const
{
	if (index >= mPolySeq.size())
		throw out_of_range("Invalid index for residue in polymer");

	string compoundID;
	int seqID;
	
	auto r = mPolySeq[index];
	
	cif::tie(seqID, compoundID) =
		r.get("seq_id", "mon_id");
	
	return Monomer(const_cast<Polymer&>(*this), index, seqID, compoundID, "");
}

Polymer::iterator::iterator(const iterator& rhs)
	: mPolymer(rhs.mPolymer), mIndex(rhs.mIndex), mCurrent(rhs.mCurrent)
{
}

Polymer::iterator& Polymer::iterator::operator=(const iterator& rhs)
{
	if (this != &rhs)
	{
		mPolymer = rhs.mPolymer;
		mIndex = rhs.mIndex;
		mCurrent = rhs.mCurrent;
	}
	
	return *this;
}

Polymer::iterator& Polymer::iterator::operator++()
{
	auto& polySeq = mPolymer->mPolySeq;

	if (mIndex < polySeq.size())
		++mIndex;

	if (mIndex < polySeq.size())
	{
		int seqID;
		string asymID, monID;
		cif::tie(asymID, seqID, monID) =
			polySeq[mIndex].get("asym_id", "seq_id", "mon_id");
		
		mCurrent = Monomer(*mPolymer, mIndex, seqID, monID, "");
	}
	
	return *this;
}

Polymer::Polymer(const Polymer& rhs)
	: mStructure(rhs.mStructure), mEntityID(rhs.mEntityID), mAsymID(rhs.mAsymID), mPolySeq(rhs.mPolySeq)
{
	
}

Polymer::Polymer(const Structure& s, const string& asymID)
	: mStructure(const_cast<Structure*>(&s)), mAsymID(asymID)
	, mPolySeq(s.category("pdbx_poly_seq_scheme").find(cif::Key("asym_id") == mAsymID))
{
	mEntityID = mPolySeq.front()["entity_id"].as<string>();

#if DEBUG
	for (auto r: mPolySeq)
		assert(r["entity_id"] == mEntityID);
#endif
	
}

Polymer::Polymer(const Structure& s, const string& entityID, const string& asymID)
	: mStructure(const_cast<Structure*>(&s)), mEntityID(entityID), mAsymID(asymID)
	, mPolySeq(s.category("pdbx_poly_seq_scheme").find(cif::Key("asym_id") == mAsymID and cif::Key("entity_id") == mEntityID))
{
}

Polymer::iterator Polymer::begin() const
{
	return iterator(*this, 0);
}

Polymer::iterator Polymer::end() const
{
	return iterator(*this, mPolySeq.size());
}

Monomer Polymer::getBySeqID(int seqID) const
{
	for (auto m: *this)
		if (m.seqID() == seqID)
			return m;

	if (VERBOSE)
		cerr << "Monomer with seqID " << seqID << " not found in polymer " << mAsymID << endl;

	return Monomer();
}

int Polymer::Distance(const Monomer& a, const Monomer& b) const
{
	int result = numeric_limits<int>::max();
	
	if (a.asymID() == b.asymID())
	{
		int ixa = numeric_limits<int>::max(), ixb = numeric_limits<int>::max();
	
		int ix = 0, f = 0;
		for (auto m: *this)
		{
			if (m.seqID() == a.seqID())
				ixa = ix, ++f;
			if (m.seqID() == b.seqID())
				ixb = ix, ++f;
			if (f == 2)
			{
				result = abs(ixa - ixb);
				break;
			}
		}
	}

	return result;
}

// --------------------------------------------------------------------
// File

File::File()
	: mImpl(new FileImpl)
{
}

File::File(fs::path File)
	: mImpl(new FileImpl)
{
	load(File);
}

File::~File()
{
	delete mImpl;
}

void File::load(fs::path p)
{
	mImpl->load(p);
	
//	// all data is now in mFile, construct atoms and others
//	
//	auto& db = mFile.firstDatablock();
//	
//	// the entities
//	
//	struct entity
//	{
//		string				id;
//		string				type;
//	};
//	vector<entity> entities;
//	
//	for (auto& _e: db["entity"])
//	{
//		string type = _e["type"];
//		ba::to_lower(type);
//		entities.push_back({ _e["id"], type });
//	}
//
//	auto& atomSites = db["atom_site"];
//	for (auto& atomSite: atomSites)
//	{
//		AtomPtr ap(new Atom(this, atom_site));
//
//		string entity_id = atom_site["entity_id"];
//		
//		auto e = find_if(entities.begin(), entities.end(), [=](entity& e) -> bool { return e.id == entity_id; });
//		if (e == entities.end())
//			throw runtime_error("Entity " + entity_id + " is not defined");
//
//		string comp_id, asym_id, seq_id;
//		cif::tie(comp_id, seq_id) = atom_site.get("label_comp_id", "label_asym_id", "label_seq_id");
//
//		auto r = find_if(m_residues.begin(), m_residues.end(), [=](residue_ptr& res) -> bool
//		{
////			return res.entities
//			return false;
//		});
//		
//		if (e->type == "water")
//			;
//		else if (e->type == "polymer")
//			;
//		else	
//			;
//		
//		m_atoms.push_back(ap);
//	}
	
}

void File::save(boost::filesystem::path File)
{
	mImpl->save(File);
}

cif::Datablock& File::data()
{
	assert(mImpl);
	assert(mImpl->mDb);
	
	if (mImpl == nullptr or mImpl->mDb == nullptr)
		throw runtime_error("No data loaded");
	
	return *mImpl->mDb;
}

cif::File& File::file()
{
	assert(mImpl);
	
	if (mImpl == nullptr)
		throw runtime_error("No data loaded");
	
	return mImpl->mData;
}

// --------------------------------------------------------------------
//	Structure

struct StructureImpl
{
	StructureImpl(Structure& s, File& f, uint32 modelNr)
		: mFile(&f), mModelNr(modelNr)
	{
		auto& db = *mFile->impl().mDb;
		auto& atomCat = db["atom_site"];
		
		for (auto& a: atomCat)
		{
			auto modelNr = a["pdbx_PDB_model_num"];
			
			if (modelNr.empty() or modelNr.as<uint32>() == mModelNr)
				mAtoms.emplace_back(new AtomImpl(f, a["id"].as<string>(), a));
		}
	}
	
	void removeAtom(Atom& a);
	void swapAtoms(Atom& a1, Atom& a2);
	void moveAtom(Atom& a, Point p);
	void changeResidue(Residue& res, const string& newCompound,
		const vector<tuple<string,string>>& remappedAtoms);

	void insertCompound(const string& compoundID, bool isEntity);

	File*			mFile;
	uint32			mModelNr;
	AtomView		mAtoms;
};

// --------------------------------------------------------------------

void StructureImpl::insertCompound(const string& compoundID, bool isEntity)
{
	auto compound = Compound::create(compoundID);
	if (compound == nullptr)
		throw runtime_error("Trying to insert unknown compound " + compoundID + " (not found in CCP4 monomers lib)");

	cif::Datablock& db = *mFile->impl().mDb;
	
	auto& chemComp = db["chem_comp"];
	auto r = chemComp.find(cif::Key("id") == compoundID);
	if (r.empty())
	{
		chemComp.emplace({
			{ "id", compoundID },
			{ "name", compound->name() },
			{ "formula", compound->formula() },
			{ "type", compound->type() }
		});
	}
	
	if (isEntity)
	{
		auto& pdbxEntityNonpoly = db["pdbx_entity_nonpoly"];
		if (pdbxEntityNonpoly.find(cif::Key("comp_id") == compoundID).empty())
		{
			auto& entity = db["entity"];
			string entityID = to_string(entity.size() + 1);
			
			entity.emplace({
				{ "id", entityID },
				{ "type", "non-polymer" },
				{ "pdbx_description", compound->name() },
				{ "formula_weight", compound->formulaWeight() }
			});
			
			pdbxEntityNonpoly.emplace({
				{ "entity_id", entityID  },
				{ "name", compound->name() },
				{ "comp_id", compoundID }
			});
		}
	}
}

void StructureImpl::removeAtom(Atom& a)
{
	cif::Datablock& db = *mFile->impl().mDb;
	
	auto& atomSites = db["atom_site"];
	
	for (auto i = atomSites.begin(); i != atomSites.end(); ++i)
	{
		string id;
		cif::tie(id) = i->get("id");
		
		if (id == a.id())
		{
			atomSites.erase(i);
			break;
		}
	}
	
	mAtoms.erase(remove(mAtoms.begin(), mAtoms.end(), a), mAtoms.end());
}

void StructureImpl::swapAtoms(Atom& a1, Atom& a2)
{
	cif::Datablock& db = *mFile->impl().mDb;
	auto& atomSites = db["atom_site"];

	auto r1 = atomSites.find(cif::Key("id") == a1.id());
	auto r2 = atomSites.find(cif::Key("id") == a2.id());
	
	if (r1.size() != 1)
		throw runtime_error("Cannot swap atoms since the number of atoms with id " + a1.id() + " is " + to_string(r1.size()));
	
	if (r2.size() != 1)
		throw runtime_error("Cannot swap atoms since the number of atoms with id " + a2.id() + " is " + to_string(r2.size()));
	
	auto l1 = r1.front()["label_atom_id"];
	auto l2 = r2.front()["label_atom_id"];
	swap(l1, l2);
	
	auto l3 = r1.front()["auth_atom_id"];
	auto l4 = r2.front()["auth_atom_id"];
	swap(l3, l4);
}

void StructureImpl::moveAtom(Atom& a, Point p)
{
	a.location(p);
}

void StructureImpl::changeResidue(Residue& res, const string& newCompound,
		const vector<tuple<string,string>>& remappedAtoms)
{
	cif::Datablock& db = *mFile->impl().mDb;

	string entityID;

	// First make sure the compound is already known or insert it.
	// And if the residue is an entity, we must make sure it exists
	insertCompound(newCompound, res.isEntity());
	
	auto& atomSites = db["atom_site"];
	auto atoms = res.atoms();

	for (auto& a: remappedAtoms)
	{
		string a1, a2;
		tie(a1, a2) = a;
		
		auto i = find_if(atoms.begin(), atoms.end(), [&](const Atom& a) { return a.labelAtomId() == a1; });
		if (i == atoms.end())
		{
			cerr << "Missing atom for atom ID " << a1 << endl;
			continue;
		}

		auto r = atomSites.find(cif::Key("id") == i->id());

		if (r.size() != 1)
			continue;
		
		if (a1 != a2)
			r.front()["label_atom_id"] = a2;
	}
	
	for (auto a: atoms)
	{
		auto r = atomSites.find(cif::Key("id") == a.id());
		assert(r.size() == 1);

		if (r.size() != 1)
			continue;

		r.front()["label_comp_id"] = newCompound;
		if (not entityID.empty())
			r.front()["label_entity_id"] = entityID;
	}
}

// --------------------------------------------------------------------

Structure::Structure(File& f, uint32 modelNr)
	: mImpl(new StructureImpl(*this, f, modelNr))
{
}

Structure::~Structure()
{
	delete mImpl;
}

AtomView Structure::atoms() const
{
	return mImpl->mAtoms;
}

AtomView Structure::waters() const
{
	AtomView result;
	
	auto& db = datablock();
	
	// Get the entity id for water
	auto& entityCat = db["entity"];
	string waterEntityId;
	for (auto& e: entityCat)
	{
		string id, type;
		cif::tie(id, type) = e.get("id", "type");
		if (ba::iequals(type, "water"))
		{
			waterEntityId = id;
			break;
		}
	}

	for (auto& a: mImpl->mAtoms)
	{
		if (a.property<string>("label_entity_id") == waterEntityId)
			result.push_back(a);
	}
	
	return result;
}

vector<Polymer> Structure::polymers() const
{
	vector<Polymer> result;
	
	auto& polySeqScheme = category("pdbx_poly_seq_scheme");
	
	for (auto& r: polySeqScheme)
	{
		string asymID, entityID, seqID, monID;
		cif::tie(asymID, entityID, seqID, monID) =
			r.get("asym_id", "entity_id", "seq_id", "mon_id");
		
		if (result.empty() or result.back().asymID() != asymID or result.back().entityID() != entityID)
			result.emplace_back(*this, entityID, asymID);
	}
	
	return result;
}

vector<Residue> Structure::nonPolymers() const
{
	vector<Residue> result;

	auto& nonPolyScheme = category("pdbx_nonpoly_scheme");
	
	for (auto& r: nonPolyScheme)
	{
		string asymID, monID;
		cif::tie(asymID, monID) =
			r.get("asym_id", "mon_id");
		
		if (result.empty() or result.back().asymID() != asymID)
			result.emplace_back(*this, monID, asymID);
	}
	
	return result;
}

Atom Structure::getAtomById(string id) const
{
	for (auto& a: mImpl->mAtoms)
	{
		if (a.id() == id)
			return a;
	}
	
	throw out_of_range("Could not find atom with id " + id);
}

File& Structure::getFile() const
{
	return *mImpl->mFile;
}

cif::Category& Structure::category(const char* name) const
{
	auto& db = datablock();
	return db[name];
}

tuple<char,int,char> Structure::MapLabelToAuth(
	const string& asymId, int seqId) const
{
	auto& db = *getFile().impl().mDb;
	
	tuple<char,int,char> result;
	bool found = false;
	
	for (auto r: db["pdbx_poly_seq_scheme"].find(
						cif::Key("asym_id") == asymId and
						cif::Key("seq_id") == seqId))
	{
		string auth_asym_id, pdb_ins_code;
		int pdb_seq_num;
		
		cif::tie(auth_asym_id, pdb_seq_num, pdb_ins_code) =
			r.get("pdb_strand_id", "pdb_seq_num", "pdb_ins_code");

		result = make_tuple(auth_asym_id.front(), pdb_seq_num,
			pdb_ins_code.empty() ? ' ' : pdb_ins_code.front());

		found = true;
		break;
	}
						
	if (not found)
	{
		for (auto r: db["pdbx_nonpoly_scheme"].find(
							cif::Key("asym_id") == asymId and
							cif::Key("seq_id") == seqId))
		{
			string pdb_strand_id, pdb_ins_code;
			int pdb_seq_num;
			
			cif::tie(pdb_strand_id, pdb_seq_num, pdb_ins_code) =
				r.get("pdb_strand_id", "pdb_seq_num", "pdb_ins_code");
	
			result = make_tuple(pdb_strand_id.front(), pdb_seq_num,
				pdb_ins_code.empty() ? ' ' : pdb_ins_code.front());
	
			found = true;
			break;
		}
	}

	return result;
}

tuple<string,int,string,string> Structure::MapLabelToPDB(
	const string& asymId, int seqId, const string& monId) const
{
	auto& db = datablock();
	
	tuple<string,int,string,string> result;
	
	for (auto r: db["pdbx_poly_seq_scheme"].find(
						cif::Key("asym_id") == asymId and
						cif::Key("seq_id") == seqId and
						cif::Key("mon_id") == monId))
	{
		string pdbStrandId, pdbMonId, pdbInsCode;
		int pdbSeqNum;
		
		cif::tie(pdbStrandId, pdbSeqNum, pdbMonId, pdbInsCode) =
			r.get("pdb_strand_id", "pdb_seq_num", "pdb_mon_id", "pdb_ins_code");

		result = make_tuple(pdbStrandId, pdbSeqNum, pdbMonId, pdbInsCode);

		break;
	}
						
	for (auto r: db["pdbx_nonpoly_scheme"].find(
						cif::Key("asym_id") == asymId and
						cif::Key("seq_id") == seqId and
						cif::Key("mon_id") == monId))
	{
		string pdbStrandId, pdbMonId, pdbInsCode;
		int pdbSeqNum;
		
		cif::tie(pdbStrandId, pdbSeqNum, pdbMonId, pdbInsCode) =
			r.get("pdb_strand_id", "pdb_seq_num", "pdb_mon_id", "pdb_ins_code");

		result = make_tuple(pdbStrandId, pdbSeqNum, pdbMonId, pdbInsCode);

		break;
	}

	return result;
}

cif::Datablock& Structure::datablock() const
{
	return *mImpl->mFile->impl().mDb;
}

// --------------------------------------------------------------------
// actions

void Structure::removeAtom(Atom& a)
{
	mImpl->removeAtom(a);
}

void Structure::swapAtoms(Atom& a1, Atom& a2)
{
	mImpl->swapAtoms(a1, a2);
}

void Structure::moveAtom(Atom& a, Point p)
{
	mImpl->moveAtom(a, p);
}

void Structure::changeResidue(Residue& res, const string& newCompound,
		const vector<tuple<string,string>>& remappedAtoms)
{
	mImpl->changeResidue(res, newCompound, remappedAtoms);
}

}
