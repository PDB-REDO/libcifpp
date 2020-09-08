// Lib for working with structures as contained in file and PDB files

#include "cif++/Structure.hpp"

#include <numeric>
#include <fstream>

#if __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#elif __has_include(<boost/filesystem.hpp>)
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#endif

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/format.hpp>

#include "cif++/Config.hpp"
#include "cif++/PDB2Cif.hpp"
#include "cif++/CifParser.hpp"
#include "cif++/Cif2PDB.hpp"
#include "cif++/AtomShape.hpp"

using namespace std;

namespace ba = boost::algorithm;
namespace io = boost::iostreams;

extern int cif::VERBOSE;

namespace mmcif
{
	
// --------------------------------------------------------------------
// FileImpl
	
struct FileImpl
{
	cif::File			mData;
	cif::Datablock*		mDb = nullptr;
	
	void load(const std::string& p);
	void save(const std::string& p);
};

void FileImpl::load(const std::string& p)
{
	fs::path path(p);

	std::ifstream inFile(path, ios_base::in | ios_base::binary);
	if (not inFile.is_open())
		throw runtime_error("No such file: " + path.string());
	
	io::filtering_stream<io::input> in;
	string ext = path.extension().string();
	
	if (path.extension() == ".bz2")
	{
		in.push(io::bzip2_decompressor());
		ext = path.stem().extension().string();
	}
	else if (path.extension() == ".gz")
	{
		in.push(io::gzip_decompressor());
		ext = path.stem().extension().string();
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
				if (cif::VERBOSE)
					cerr << "unrecognized file extension, trying cif" << endl;
	
				mData.load(in);
			}
			catch (const cif::CifParserError& e)
			{
				if (cif::VERBOSE)
					cerr << "Not cif, trying plain old PDB" << endl;
	
				// pffft...
				in.reset();
	
				if (inFile.is_open())
					inFile.seekg(0);
				else
					inFile.open(path, ios_base::in | ios::binary);
	
				if (path.extension() == ".bz2")
					in.push(io::bzip2_decompressor());
				else if (path.extension() == ".gz")
					in.push(io::gzip_decompressor());
				
				in.push(inFile);
	
				ReadPDBFile(in, mData);
			}
		}
	}
	catch (const exception& ex)
	{
		cerr << "Error trying to load file " << path << endl;
		throw;
	}
	
	// Yes, we've parsed the data. Now locate the datablock.
	mDb = &mData.firstDatablock();
	
	// And validate, otherwise lots of functionality won't work
//	if (mData.getValidator() == nullptr)
	mData.loadDictionary("mmcif_pdbx");
	if (not mData.isValid())
		cerr << "Invalid mmCIF file" << (cif::VERBOSE ? "." : " use --verbose option to see errors") << endl;
}

void FileImpl::save(const std::string& p)
{
	fs::path path(p);

	std::ofstream outFile(p, ios_base::out | ios_base::binary);
	io::filtering_stream<io::output> out;
	
	if (path.extension() == ".gz")
	{
		out.push(io::gzip_compressor());
		path = path.stem();
	}
	else if (path.extension() == ".bz2")
	{
		out.push(io::bzip2_compressor());
		path = path.stem();
	}
	
	out.push(outFile);
	
	if (path.extension() == ".pdb")
		WritePDBFile(out, mData);
	else
		mData.save(out);
}


// --------------------------------------------------------------------
// Atom

struct AtomImpl
{
	AtomImpl(const AtomImpl& i)
		: mFile(i.mFile), mID(i.mID), mType(i.mType)
		, mAtomID(i.mAtomID), mCompID(i.mCompID), mAsymID(i.mAsymID)
		, mSeqID(i.mSeqID), mAltID(i.mAltID), mLocation(i.mLocation)
		, mRefcount(1), mRow(i.mRow), mCompound(i.mCompound)
		, mRadius(i.mRadius), mCachedProperties(i.mCachedProperties)
		, mSymmetryCopy(i.mSymmetryCopy), mClone(true)
		, mRTop(i.mRTop), mD(i.mD)
	{
	}
	
	AtomImpl(const File& f, const string& id)
		: mFile(f), mID(id), mRefcount(1), mCompound(nullptr)
	{
		auto& db = *mFile.impl().mDb;
		auto& cat = db["atom_site"];
		
		mRow = cat[cif::Key("id") == mID];

		prefetch();
	}
	
	AtomImpl(const File& f, const string& id, cif::Row row)
		: mFile(f), mID(id), mRefcount(1), mRow(row), mCompound(nullptr)
	{
		prefetch();
	}

	AtomImpl(const AtomImpl& impl, const Point& d, const clipper::RTop_orth& rt)
		: mFile(impl.mFile), mID(impl.mID), mType(impl.mType), mAtomID(impl.mAtomID)
		, mCompID(impl.mCompID), mAsymID(impl.mAsymID), mSeqID(impl.mSeqID)
		, mAltID(impl.mAltID), mLocation(impl.mLocation), mRefcount(1)
		, mRow(impl.mRow), mCompound(impl.mCompound), mRadius(impl.mRadius)
		, mCachedProperties(impl.mCachedProperties)
		, mSymmetryCopy(true), mRTop(rt), mD(d)
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
		
		if (symbol != "X")
			mType = AtomTypeTraits(symbol).type();

		float x, y, z;
		cif::tie(x, y, z) = mRow.get("Cartn_x", "Cartn_y", "Cartn_z");
		
		mLocation = Point(x, y, z);

		string compID;
		cif::tie(compID) = mRow.get("label_comp_id");
		
		mCompound = Compound::create(compID);
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
			if (abs(charge) > 1)
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
		auto r = cat[cif::Key("id") == mID];
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
		if (--mRefcount <= 0)
			delete this;
	}
	
	bool getAnisoU(float anisou[6]) const
	{
		auto& db = *mFile.impl().mDb;
		auto& cat = db["atom_site_anisotrop"];
		auto r = cat[cif::Key("id") == mID];
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
		
		if (not mClone)
		{
			mRow["Cartn_x"] = p.getX();
			mRow["Cartn_y"] = p.getY();
			mRow["Cartn_z"] = p.getZ();
		}
		
//		boost::format kPosFmt("%.3f");
//
//		mRow["Cartn_x"] = (kPosFmt % p.getX()).str();
//		mRow["Cartn_y"] = (kPosFmt % p.getY()).str();
//		mRow["Cartn_z"] = (kPosFmt % p.getZ()).str();

		mLocation = p;
	}
	
	const Compound& comp() const
	{
		if (mCompound == nullptr)
		{
			string compID;
			cif::tie(compID) = mRow.get("label_comp_id");
			
			mCompound = Compound::create(compID);
			
			if (cif::VERBOSE and mCompound == nullptr)
				cerr << "Compound not found: '" << compID << '\'' << endl;
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

	const string& property(const string& name) const
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

	int compare(const AtomImpl& b) const
	{
		int d = mAsymID.compare(b.mAsymID);
		if (d == 0)
			d = mSeqID - b.mSeqID;
		if (d == 0)
			d = mAtomID.compare(b.mAtomID);
		return d;
	}

	void swapAtomLabels(AtomImpl& b)
	{
		std::swap(mAtomID, b.mAtomID);
	}

	const File&			mFile;
	string				mID;
	AtomType			mType;

	string				mAtomID;
	string				mCompID;
	string				mAsymID;
	int					mSeqID;
	string				mAltID;
	
	Point				mLocation;
	int					mRefcount;
	cif::Row			mRow;
	mutable const Compound*		mCompound;
	float				mRadius = nan("4");
	mutable map<string,string>	mCachedProperties;
	
	bool				mSymmetryCopy = false;
	bool				mClone = false;

	clipper::RTop_orth	mRTop;
	Point				mD;
	int32_t				mRTix;
};

//Atom::Atom(const File& f, const string& id)
//	: mImpl(new AtomImpl(f, id))
//{
//}
//

Atom::Atom()
	: mImpl_(nullptr)
{
}

Atom::Atom(AtomImpl* impl)
	: mImpl_(impl)
{
}

AtomImpl* Atom::impl()
{
	if (mImpl_ == nullptr)
		throw std::runtime_error("atom is not set");
	return mImpl_;
}

const AtomImpl* Atom::impl() const
{
	if (mImpl_ == nullptr)
		throw std::runtime_error("atom is not set");
	return mImpl_;
}

Atom Atom::clone() const
{
	return Atom(mImpl_ ? new AtomImpl(*mImpl_) : nullptr);
}

Atom::Atom(const Atom& rhs)
	: mImpl_(rhs.mImpl_)
{
	if (mImpl_)
		mImpl_->reference();
}

Atom::~Atom()
{
	if (mImpl_)
		mImpl_->release();
}

Atom& Atom::operator=(const Atom& rhs)
{
	if (this != &rhs)
	{
		if (mImpl_)
			mImpl_->release();
		mImpl_ = rhs.mImpl_;
		if (mImpl_)
			mImpl_->reference();
	}

	return *this;
}

template<>
string Atom::property<string>(const string& name) const
{
	return impl()->property(name);
}

template<>
int Atom::property<int>(const string& name) const
{
	auto v = impl()->property(name);
	return v.empty() ? 0 : stoi(v);
}

template<>
float Atom::property<float>(const string& name) const
{
	return stof(impl()->property(name));
}

const string& Atom::id() const
{
	return impl()->mID;
}

AtomType Atom::type() const
{
	return impl()->mType;
}

int Atom::charge() const
{
	return property<int>("pdbx_formal_charge");
}

string Atom::energyType() const
{
	string result;
	
	if (impl() and impl()->mCompound)
		result = impl()->mCompound->getAtomByID(impl()->mAtomID).typeEnergy;
	
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
	return impl()->getAnisoU(anisou);
}

float Atom::occupancy() const
{
	return property<float>("occupancy");
}

string Atom::labelAtomID() const
{
	return impl()->mAtomID;
}

string Atom::labelCompID() const
{
	return impl()->mCompID;
}

string Atom::labelAsymID() const
{
	return impl()->mAsymID;
}

string Atom::labelAltID() const
{
	return impl()->mAltID;
}

bool Atom::isAlternate() const
{
	return not impl()->mAltID.empty();
}

int Atom::labelSeqID() const
{
	return impl()->mSeqID;
}

string Atom::authAsymID() const
{
	return property<string>("auth_asym_id");
}

string Atom::authAtomID() const
{
	return property<string>("auth_atom_id");
}

string Atom::pdbxAuthAltID() const
{
	return property<string>("pdbx_auth_alt_id");
}

string Atom::pdbxAuthInsCode() const
{
	return property<string>("pdbx_PDB_ins_code");
}

string Atom::authCompID() const
{
	return property<string>("auth_comp_id");
}

string Atom::authSeqID() const
{
	return property<string>("auth_seq_id");
}

string Atom::labelID() const
{
	return property<string>("label_comp_id") + '_' + impl()->mAsymID + '_' + to_string(impl()->mSeqID) + ':' + impl()->mAtomID;
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
	return impl()->mLocation;
}

void Atom::location(Point p)
{
	impl()->moveTo(p);
}

Atom Atom::symmetryCopy(const Point& d, const clipper::RTop_orth& rt)
{
	return Atom(new AtomImpl(*impl(), d, rt));
}

bool Atom::isSymmetryCopy() const
{
	return impl()->mSymmetryCopy;
}

string Atom::symmetry() const
{
	return clipper::Symop(impl()->mRTop).format() + "\n" + impl()->mRTop.format();
}

const clipper::RTop_orth& Atom::symop() const
{
	return impl()->mRTop;
}

const Compound& Atom::comp() const
{
	return impl()->comp();
}

bool Atom::isWater() const
{
	return impl()->isWater();
}

bool Atom::operator==(const Atom& rhs) const
{
	return impl() == rhs.impl() or
		(&impl()->mFile == &rhs.impl()->mFile and impl()->mID == rhs.impl()->mID); 	
}

clipper::Atom Atom::toClipper() const
{
	return impl()->toClipper();
}

void Atom::calculateRadius(float resHigh, float resLow, float perc)
{
	AtomShape shape(*this, resHigh, resLow, false);
	impl()->mRadius = shape.radius();

	// verbose
	if (cif::VERBOSE > 1)
		cout << "Calculated radius for " << AtomTypeTraits(impl()->mType).name() << " with charge " << charge() << " is " << impl()->mRadius << endl;
}

float Atom::radius() const
{
	return impl()->mRadius;
}

int Atom::compare(const Atom& b) const
{
	return impl() == b.impl() ? 0 : impl()->compare(*b.impl());
}

void Atom::setID(int id)
{
	impl()->mID = to_string(id);
}

std::ostream& operator<<(std::ostream& os, const Atom& atom)
{
	os << atom.labelCompID() << ' ' << atom.labelAsymID() << ':' << atom.labelSeqID();

	if (atom.isAlternate())
		os << '(' << atom.labelAltID() << ')';
	if (atom.authAsymID() != atom.labelAsymID() or atom.authSeqID() != std::to_string(atom.labelSeqID()) or atom.pdbxAuthInsCode().empty() == false)
		os << " [" << atom.authAsymID() << ':' << atom.authSeqID() << atom.pdbxAuthInsCode() << ']';

	return os;
}

// --------------------------------------------------------------------
// residue

Residue::Residue(const Structure& structure, const string& compoundID,
	const string& asymID, const string& authSeqID)
	: mStructure(&structure), mCompoundID(compoundID)
	, mAsymID(asymID), mAuthSeqID(authSeqID)
{
	assert(mCompoundID == "HOH");
	
	for (auto& a: mStructure->atoms())
	{
		if (a.labelAsymID() != mAsymID or
			a.labelCompID() != mCompoundID)
				continue;

		if (not mAuthSeqID.empty() and a.authSeqID() != mAuthSeqID)		// water!
			continue;
		
		mAtoms.push_back(a);
	}
	
	assert(not mAtoms.empty());
}

Residue::Residue(const Structure& structure, const string& compoundID,
	const string& asymID, int seqID)
	: mStructure(&structure), mCompoundID(compoundID)
	, mAsymID(asymID), mSeqID(seqID)
{
	assert(mCompoundID != "HOH");
	
	for (auto& a: mStructure->atoms())
	{
		if (mSeqID > 0 and a.labelSeqID() != mSeqID)
			continue;
		
		if (a.labelAsymID() != mAsymID or
			a.labelCompID() != mCompoundID)
				continue;
		
		mAtoms.push_back(a);
	}
}

Residue::Residue(Residue&& rhs)
	: mStructure(rhs.mStructure), mCompoundID(move(rhs.mCompoundID)), mAsymID(move(rhs.mAsymID))
	, mSeqID(rhs.mSeqID), mAuthSeqID(rhs.mAuthSeqID), mAtoms(move(rhs.mAtoms))
{
//cerr << "move constructor residue" << endl;
	rhs.mStructure = nullptr;
}

Residue& Residue::operator=(Residue&& rhs)
{
//cerr << "move assignment residue" << endl;
	mStructure = rhs.mStructure;		rhs.mStructure = nullptr;
	mCompoundID = move(rhs.mCompoundID);
	mAsymID = move(rhs.mAsymID);
	mSeqID = rhs.mSeqID;
	mAuthSeqID = rhs.mAuthSeqID;
	mAtoms = move(rhs.mAtoms);

	return *this;
}

Residue::~Residue()
{
//cerr << "~Residue" << endl;	
}

string Residue::authInsCode() const
{
	assert(mStructure);

	string result;
	
	try
	{
		char iCode;
		tie(ignore, ignore, iCode) = mStructure->MapLabelToAuth(mAsymID, mSeqID);
		
		result = string{ iCode };
		ba::trim(result);
	}
	catch (...)
	{
	}
	
	return result;
}

string Residue::authAsymID() const
{
	assert(mStructure);

	string result;
	
	try
	{
		tie(result, ignore, ignore) = mStructure->MapLabelToAuth(mAsymID, mSeqID);
	}
	catch (...)
	{
		result = mAsymID;
	}
	
	return result;
}

string Residue::authSeqID() const
{
	assert(mStructure);

	string result;
	
	try
	{
		int seqID;
		tie(ignore, seqID, ignore) = mStructure->MapLabelToAuth(mAsymID, mSeqID);
		result = to_string(seqID);
	}
	catch (...)
	{
	}
	
	return result;
}

const Compound& Residue::compound() const
{
	auto result = Compound::create(mCompoundID);
	if (result == nullptr)
		throw runtime_error("Failed to create compound " + mCompoundID);
	return *result;
}

const AtomView& Residue::atoms() const
{
	if (mStructure == nullptr)
		throw runtime_error("Invalid Residue object");

	return mAtoms;
}

std::string Residue::unique_alt_id() const
{
	if (mStructure == nullptr)
		throw runtime_error("Invalid Residue object");

	auto firstAlt = std::find_if(mAtoms.begin(), mAtoms.end(), [](auto& a) { return not a.labelAltID().empty(); });

	return firstAlt != mAtoms.end() ? firstAlt->labelAltID() : "";
}

AtomView Residue::unique_atoms() const
{
	if (mStructure == nullptr)
		throw runtime_error("Invalid Residue object");

	AtomView result;
	std::string firstAlt;

	for (auto& atom: mAtoms)
	{
		auto alt = atom.labelAltID();
		if (alt.empty())
		{
			result.push_back(atom);
			continue;
		}

		if (firstAlt.empty())
			firstAlt = alt;
		else if (alt != firstAlt)
		{
			if (cif::VERBOSE)
				std::cerr << "skipping alternate atom " << atom << std::endl;
			continue;
		}

		result.push_back(atom);
	}
	
	return result;
}

Atom Residue::atomByID(const string& atomID) const
{
	Atom result;

	for (auto& a: mAtoms)
	{
		if (a.labelAtomID() == atomID)
		{
			result = a;
			break;
		}
	}

	if (not result and cif::VERBOSE)
		std::cerr << "Atom with atom_id " << atomID << " not found in residue " << mAsymID << ':' << mSeqID << std::endl;

	return result;
}

// Residue is a single entity if the atoms for the asym with mAsymID is equal
// to the number of atoms in this residue...  hope this is correct....
bool Residue::isEntity() const
{
	auto& db = mStructure->datablock();
	
	auto a1 = db["atom_site"].find(cif::Key("label_asym_id") == mAsymID);
//	auto a2 = atoms();
	auto& a2 = mAtoms;
	
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

string Residue::labelID() const
{
	if (mCompoundID == "HOH")
		return mAsymID + mAuthSeqID;
	else
		return mAsymID + to_string(mSeqID);
}

tuple<Point,float> Residue::centerAndRadius() const
{
	vector<Point> pts;
	for (auto& a: mAtoms)
		pts.push_back(a.location());
	
	auto center = Centroid(pts);
	float radius = 0;
	
	for (auto& pt: pts)
	{
		auto d = Distance(pt, center);
		if (radius < d)
			radius = d;
	}
	
	return make_tuple(center, radius);
}

bool Residue::hasAlternateAtoms() const
{
	return std::find_if(mAtoms.begin(), mAtoms.end(), [](const Atom& atom) { return atom.isAlternate(); }) != mAtoms.end();
}

std::ostream& operator<<(std::ostream& os, const Residue& res)
{
	os << res.compoundID() << ' ' << res.asymID() << ':' << res.seqID();

	if (res.authAsymID() != res.asymID() or res.authSeqID() != std::to_string(res.seqID()))
		os << " [" << res.authAsymID() << ':' << res.authSeqID() << ']';

	return os;
}

// --------------------------------------------------------------------
// monomer

//Monomer::Monomer(Monomer&& rhs)
//	: Residue(move(rhs)), mPolymer(rhs.mPolymer), mIndex(rhs.mIndex)
//{
//}

Monomer::Monomer(const Polymer& polymer, uint32_t index, int seqID, const string& compoundID)
	: Residue(*polymer.structure(), compoundID, polymer.asymID(), seqID)
	, mPolymer(&polymer)
	, mIndex(index)
{
}

Monomer::Monomer(Monomer&& rhs)
	: Residue(move(rhs)), mPolymer(rhs.mPolymer), mIndex(rhs.mIndex)
{
cerr << "move constructor monomer" << endl;

//	mStructure = rhs.mStructure;			rhs.mStructure = nullptr;
//	mCompoundID = move(rhs.mCompoundID);
//	mAsymID = move(rhs.mAsymID);
//	mSeqID = rhs.mSeqID;
//	mAtoms = move(rhs.mAtoms);
//	
//	mPolymer = rhs.mPolymer; rhs.mPolymer = nullptr;
//	mIndex = rhs.mIndex;
	rhs.mPolymer = nullptr;
}

Monomer& Monomer::operator=(Monomer&& rhs)
{
cerr << "move assignment monomer" << endl;

	Residue::operator=(move(rhs));
	mPolymer = rhs.mPolymer;		rhs.mPolymer = nullptr;
	mIndex = rhs.mIndex;
	
	return *this;
}

bool Monomer::is_first_in_chain() const
{
	return mIndex == 0;
}

bool Monomer::is_last_in_chain() const
{
	return mIndex + 1 == mPolymer->size();
}

bool Monomer::has_alpha() const
{
	return
		mIndex >= 1 and mIndex + 2 < mPolymer->size();
}

bool Monomer::has_kappa() const
{
	return mIndex >= 2 and mIndex + 2 < mPolymer->size();
}

float Monomer::phi() const
{
	float result = 360;

	try
	{
		if (mIndex > 0)
		{
			auto& prev = mPolymer->operator[](mIndex - 1);
			if (prev.mSeqID + 1 == mSeqID)
				result = DihedralAngle(prev.C().location(), N().location(), CAlpha().location(), C().location()); 
		}
	}
	catch (const exception& ex)
	{
		if (cif::VERBOSE)
			cerr << ex.what() << endl;
	}

	return result;
}

float Monomer::psi() const
{
	float result = 360;
	
	try
	{
		if (mIndex + 1 < mPolymer->size())
		{
			auto& next = mPolymer->operator[](mIndex + 1);
			if (mSeqID + 1 == next.mSeqID)
				result = DihedralAngle(N().location(), CAlpha().location(), C().location(), next.N().location()); 
		}
	}
	catch (const exception& ex)
	{
		if (cif::VERBOSE)
			cerr << ex.what() << endl;
	}

	return result;
}

float Monomer::alpha() const
{
	float result = 360;

	try
	{
		if (mIndex >= 1 and mIndex + 2 < mPolymer->size())
		{
			auto& prev = mPolymer->operator[](mIndex - 1);
			auto& next = mPolymer->operator[](mIndex + 1);
			auto& nextNext = mPolymer->operator[](mIndex + 2);
			
			result = DihedralAngle(prev.CAlpha().location(), CAlpha().location(), next.CAlpha().location(), nextNext.CAlpha().location());
		}
	}
	catch (const exception& ex)
	{
		if (cif::VERBOSE)
			cerr << ex.what() << endl;
	}

	return result;
}

float Monomer::kappa() const
{
	double result = 360;
	
	try
	{
		if (mIndex >= 2 and mIndex + 2 < mPolymer->size())
		{
			auto& prevPrev = mPolymer->operator[](mIndex - 2);
			auto& nextNext = mPolymer->operator[](mIndex + 2);
			
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
		if (cif::VERBOSE)
			cerr << "When trying to calculate kappa for " << asymID() << ':' << seqID() << ": "
				 << ex.what() << endl;
	}

	return result;
}

float Monomer::tco() const
{
	double result = 0.0;
	
	try
	{
		if (mIndex > 0)
		{
			auto& prev = mPolymer->operator[](mIndex - 1);
			if (prev.mSeqID + 1 == mSeqID)
				result = CosinusAngle(C().location(), O().location(), prev.C().location(), prev.O().location()); 
		}
	}
	catch (const exception& ex)
	{
		if (cif::VERBOSE)
			cerr << "When trying to calculate tco for " << asymID() << ':' << seqID() << ": "
				 << ex.what() << endl;
	}

	return result;
}

float Monomer::omega() const
{
	double result = 360;
	
	try
	{
		if (not is_last_in_chain())
			result = omega(*this, mPolymer->operator[](mIndex + 1));
	}
	catch (const exception& ex)
	{
		if (cif::VERBOSE)
			cerr << "When trying to calculate omega for " << asymID() << ':' << seqID() << ": "
				 << ex.what() << endl;
	}

	return result;
}

const map<string,vector<string>> kChiAtomsMap = {
   	{ "ASP", {"CG", "OD1"} },
	{ "ASN", {"CG", "OD1"} },
    { "ARG", {"CG", "CD", "NE", "CZ"} },
    { "HIS", {"CG", "ND1"} },
    { "GLN", {"CG", "CD", "OE1"} },
	{ "GLU", {"CG", "CD", "OE1"} },
    { "SER", {"OG"} },
    { "THR", {"OG1"} },
    { "LYS", {"CG", "CD", "CE", "NZ"} },
    { "TYR", {"CG", "CD1"} },
    { "PHE", {"CG", "CD1"} },
    { "LEU", {"CG", "CD1"} },
    { "TRP", {"CG", "CD1"} },
    { "CYS", {"SG"} },
    { "ILE", {"CG1", "CD1"} },
    { "MET", {"CG", "SD", "CE"} },
    { "MSE", {"CG", "SE", "CE"} },
    { "PRO", {"CG", "CD"} },
    { "VAL", {"CG1"} }
};

size_t Monomer::nrOfChis() const
{
	size_t result = 0;

	auto i = kChiAtomsMap.find(mCompoundID);
	if (i != kChiAtomsMap.end())
		result = i->second.size();
	
	return result;
}

float Monomer::chi(size_t nr) const
{
	float result = 0;

	try
	{
		auto i = kChiAtomsMap.find(mCompoundID);
		if (i != kChiAtomsMap.end() and nr < i->second.size())
		{
			vector<string> atoms{ "N", "CA", "CB" };

			atoms.insert(atoms.end(), i->second.begin(), i->second.end());

			// in case we have a positive chiral volume we need to swap atoms
			if (chiralVolume() > 0)
			{
				if (mCompoundID == "LEU") atoms.back() = "CD2";
				if (mCompoundID == "VAL") atoms.back() = "CG2";
			}

			result = DihedralAngle(
				atomByID(atoms[nr + 0]).location(),
				atomByID(atoms[nr + 1]).location(),
				atomByID(atoms[nr + 2]).location(),
				atomByID(atoms[nr + 3]).location()
			);
		}
	}
	catch(const std::exception& e)
	{
		if (cif::VERBOSE)
			std::cerr << e.what() << std::endl;
		result = 0;
	}

	return result;
}

bool Monomer::isCis() const
{
	bool result = false;
	
	if (mIndex + 1 < mPolymer->size())
	{
		auto& next = mPolymer->operator[](mIndex + 1);
		
		result = Monomer::isCis(*this, next);
	}
	
	return result;
}

bool Monomer::isComplete() const
{
	int seen = 0;
	for (auto& a: mAtoms)
	{
		     if (a.labelAtomID() == "CA")		seen |= 1;
		else if (a.labelAtomID() == "C")		seen |= 2;
		else if (a.labelAtomID() == "N")		seen |= 4;
		else if (a.labelAtomID() == "O")		seen |= 8;
		// else if (a.labelAtomID() == "OXT")		seen |= 16;
	}
	return seen == 15;
}

bool Monomer::hasAlternateBackboneAtoms() const
{
	bool result;
	for (auto& a: mAtoms)
	{
		if (not a.isAlternate())
			continue;

		auto atomID = a.labelAtomID();
		if (atomID == "CA" or atomID == "C" or atomID == "N" or atomID == "O")
		{
			result = true;
			break;
		}
	}
	return result;
}

float Monomer::chiralVolume() const
{
	float result = 0;

	if (mCompoundID == "LEU")
	{
		auto centre = atomByID("CG");
		auto atom1 = atomByID("CB");
		auto atom2 = atomByID("CD1");
		auto atom3 = atomByID("CD2");

		result = DotProduct(atom1.location() - centre.location(),
			CrossProduct(atom2.location() - centre.location(), atom3.location() - centre.location()));
	}
	else if (mCompoundID == "VAL")
	{
		auto centre = atomByID("CB");
		auto atom1 = atomByID("CA");
		auto atom2 = atomByID("CG1");
		auto atom3 = atomByID("CG2");

		result = DotProduct(atom1.location() - centre.location(),
			CrossProduct(atom2.location() - centre.location(), atom3.location() - centre.location()));
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

float Monomer::omega(const mmcif::Monomer& a, const mmcif::Monomer& b)
{
	float result = 360;
	
	try
	{
		result = DihedralAngle(
			a.atomByID("CA").location(),
			a.atomByID("C").location(),
			b.atomByID("N").location(),
			b.atomByID("CA").location());
	}
	catch (...) {}

	return result;
}

bool Monomer::isCis(const mmcif::Monomer& a, const mmcif::Monomer& b)
{
	return omega(a, b) < 30.0f;
}

// --------------------------------------------------------------------
// polymer
//
//Polymer::iterator::iterator(const Polymer& p, uint32_t index)
//	: mPolymer(&p), mIndex(index), mCurrent(p, index)
//{
//	auto& polySeq = mPolymer->mPolySeq;
//
//	if (index < polySeq.size())
//	{
//		int seqID;
//		string asymID, monID;
//		cif::tie(asymID, seqID, monID) =
//			polySeq[mIndex].get("asym_id", "seq_id", "mon_id");
//		
//		mCurrent = Monomer(*mPolymer, index, seqID, monID, "");
//	}
//}
//
//Monomer Polymer::operator[](size_t index) const
//{
//	if (index >= mPolySeq.size())
//		throw out_of_range("Invalid index for residue in polymer");
//
//	string compoundID;
//	int seqID;
//	
//	auto r = mPolySeq[index];
//	
//	cif::tie(seqID, compoundID) =
//		r.get("seq_id", "mon_id");
//	
//	return Monomer(const_cast<Polymer&>(*this), index, seqID, compoundID, "");
//}
//
//Polymer::iterator::iterator(const iterator& rhs)
//	: mPolymer(rhs.mPolymer), mIndex(rhs.mIndex), mCurrent(rhs.mCurrent)
//{
//}
//
//Polymer::iterator& Polymer::iterator::operator++()
//{
//	auto& polySeq = mPolymer->mPolySeq;
//
//	if (mIndex < polySeq.size())
//		++mIndex;
//
//	if (mIndex < polySeq.size())
//	{
//		int seqID;
//		string asymID, monID;
//		cif::tie(asymID, seqID, monID) =
//			polySeq[mIndex].get("asym_id", "seq_id", "mon_id");
//		
//		mCurrent = Monomer(*mPolymer, mIndex, seqID, monID, "");
//	}
//	
//	return *this;
//}

//Polymer::Polymer(const Structure& s, const string& asymID)
//	: mStructure(const_cast<Structure*>(&s)), mAsymID(asymID)
//	, mPolySeq(s.category("pdbx_poly_seq_scheme").find(cif::Key("asym_id") == mAsymID))
//{
//	mEntityID = mPolySeq.front()["entity_id"].as<string>();
//
//#if DEBUG
//	for (auto r: mPolySeq)
//		assert(r["entity_id"] == mEntityID);
//#endif
//	
//}

//Polymer::Polymer(Polymer&& rhs)
//	: vector<Monomer>(move(rhs))
//	, mStructure(rhs.mStructure)
//	, mEntityID(move(rhs.mEntityID)), mAsymID(move(rhs.mAsymID)), mPolySeq(move(rhs.mPolySeq))
//{
//	rhs.mStructure = nullptr;
//}
//
//Polymer& Polymer::operator=(Polymer&& rhs)
//{
//	vector<Monomer>::operator=(move(rhs));
//	mStructure = rhs.mStructure;			rhs.mStructure = nullptr;
//	mEntityID = move(rhs.mEntityID);
//	mAsymID = move(rhs.mAsymID);
//	mPolySeq = move(rhs.mPolySeq);
//	return *this;
//}

Polymer::Polymer(const Structure& s, const string& entityID, const string& asymID)
	: mStructure(const_cast<Structure*>(&s)), mEntityID(entityID), mAsymID(asymID)
	, mPolySeq(s.category("pdbx_poly_seq_scheme").find(cif::Key("asym_id") == mAsymID and cif::Key("entity_id") == mEntityID))
{
	map<uint32_t,uint32_t> ix;

	reserve(mPolySeq.size());

	for (auto r: mPolySeq)
	{
		int seqID;
		string compoundID;
		cif::tie(seqID, compoundID) = r.get("seq_id", "mon_id");

		uint32_t index = size();
		
		// store only the first
		if (not ix.count(seqID))
		{
			ix[seqID] = index;
			emplace_back(*this, index, seqID, compoundID);
		}
		else if (cif::VERBOSE)
		{
			Monomer m{*this, index, seqID, compoundID};
			std::cerr << "Dropping alternate residue " << m << std::endl;
		}
	}
}

string Polymer::chainID() const
{
	return mPolySeq.front()["pdb_strand_id"].as<string>();
}

Monomer& Polymer::getBySeqID(int seqID)
{
	for (auto& m: *this)
		if (m.seqID() == seqID)
			return m;
	
	throw runtime_error("Monomer with seqID " + to_string(seqID) + " not found in polymer " + mAsymID);
}

const Monomer& Polymer::getBySeqID(int seqID) const
{
	for (auto& m: *this)
		if (m.seqID() == seqID)
			return m;
	
	throw runtime_error("Monomer with seqID " + to_string(seqID) + " not found in polymer " + mAsymID);
}

int Polymer::Distance(const Monomer& a, const Monomer& b) const
{
	int result = numeric_limits<int>::max();
	
	if (a.asymID() == b.asymID())
	{
		int ixa = numeric_limits<int>::max(), ixb = numeric_limits<int>::max();
	
		int ix = 0, f = 0;
		for (auto& m: *this)
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

File::File(const std::string& File)
	: mImpl(new FileImpl)
{
	load(File);
}

File::~File()
{
	delete mImpl;
}

void File::load(const std::string& p)
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

void File::save(const std::string& file)
{
	mImpl->save(file);
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

Structure::Structure(File& f, uint32_t modelNr, StructureOpenOptions options)
	: mFile(f), mModelNr(modelNr)
{
	auto& db = *mFile.impl().mDb;
	auto& atomCat = db["atom_site"];
	
	for (auto& a: atomCat)
	{
		std::string id, typeSymbol;
		std::optional<int> modelNr;

		cif::tie(id, typeSymbol, modelNr) = a.get("id", "type_symbol", "pdbx_PDB_model_num");

		if (modelNr and *modelNr != mModelNr)
			continue;
		
		if ((options bitand StructureOpenOptions::SkipHydrogen) and typeSymbol == "H")
			continue;
		
		mAtoms.emplace_back(new AtomImpl(f, id, a));
	}
	
	loadData();
}

Structure::Structure(const Structure& s)
	: mFile(s.mFile), mModelNr(s.mModelNr)
{
	mAtoms.reserve(s.mAtoms.size());
	for (auto& atom: s.mAtoms)
		mAtoms.emplace_back(atom.clone());
	
	loadData();
}

Structure::~Structure()
{
}

void Structure::loadData()
{
	updateAtomIndex();

	auto& polySeqScheme = category("pdbx_poly_seq_scheme");
	
	for (auto& r: polySeqScheme)
	{
		string asymID, entityID, seqID, monID;
		cif::tie(asymID, entityID, seqID, monID) =
			r.get("asym_id", "entity_id", "seq_id", "mon_id");
		
		if (mPolymers.empty() or mPolymers.back().asymID() != asymID or mPolymers.back().entityID() != entityID)
			mPolymers.emplace_back(*this, entityID, asymID);
	}

	auto& nonPolyScheme = category("pdbx_nonpoly_scheme");
	
	for (auto& r: nonPolyScheme)
	{
		string asymID, monID, pdbSeqNum;
		cif::tie(asymID, monID, pdbSeqNum) =
			r.get("asym_id", "mon_id", "pdb_seq_num");
		
		if (monID == "HOH")
			mNonPolymers.emplace_back(*this, monID, asymID, pdbSeqNum);
		else if (mNonPolymers.empty() or mNonPolymers.back().asymID() != asymID)
			mNonPolymers.emplace_back(*this, monID, asymID);
	}
}

void Structure::updateAtomIndex()
{
	mAtomIndex = vector<size_t>(mAtoms.size());
	
	iota(mAtomIndex.begin(), mAtomIndex.end(), 0);
	
	sort(mAtomIndex.begin(), mAtomIndex.end(), [this](size_t a, size_t b) { return mAtoms[a].id() < mAtoms[b].id(); });
}

void Structure::sortAtoms()
{
	sort(mAtoms.begin(), mAtoms.end(), [](auto& a, auto& b) { return a.compare(b) < 0; });

	int id = 1;
	for (auto& atom: mAtoms)
	{
		atom.setID(id);
		++id;
	}

	updateAtomIndex();
}

AtomView Structure::waters() const
{
	AtomView result;
	
	auto& db = datablock();
	
	// Get the entity id for water
	auto& entityCat = db["entity"];
	string waterEntityID;
	for (auto& e: entityCat)
	{
		string id, type;
		cif::tie(id, type) = e.get("id", "type");
		if (ba::iequals(type, "water"))
		{
			waterEntityID = id;
			break;
		}
	}

	for (auto& a: mAtoms)
	{
		if (a.property<string>("label_entity_id") == waterEntityID)
			result.push_back(a);
	}
	
	return result;
}

Atom Structure::getAtomByID(string id) const
{
	auto i = lower_bound(mAtomIndex.begin(), mAtomIndex.end(),
		id, [this](auto& a, auto& b) { return mAtoms[a].id() < b; });

//	auto i = find_if(mAtoms.begin(), mAtoms.end(),
//		[&id](auto& a) { return a.id() == id; });

	if (i == mAtomIndex.end() or mAtoms[*i].id() != id)
		throw out_of_range("Could not find atom with id " + id);

	return mAtoms[*i];
}

Atom Structure::getAtomByLabel(const string& atomID, const string& asymID, const string& compID, int seqID, const string& altID)
{
	for (auto& a: mAtoms)
	{
		if (a.labelAtomID() == atomID and
			a.labelAsymID() == asymID and
			a.labelCompID() == compID and
			a.labelSeqID() == seqID and
			a.labelAltID() == altID)
		{
			return a;
		}
	}

	throw out_of_range("Could not find atom with specified label");
}

File& Structure::getFile() const
{
	return mFile;
}

cif::Category& Structure::category(const char* name) const
{
	auto& db = datablock();
	return db[name];
}

tuple<char,int,char> Structure::MapLabelToAuth(
	const string& asymID, int seqID) const
{
	auto& db = *getFile().impl().mDb;
	
	tuple<char,int,char> result;
	bool found = false;
	
	for (auto r: db["pdbx_poly_seq_scheme"].find(
						cif::Key("asym_id") == asymID and
						cif::Key("seq_id") == seqID))
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
							cif::Key("asym_id") == asymID and
							cif::Key("seq_id") == seqID))
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
	const string& asymID, int seqID, const string& monID,
	const string& authSeqID) const
{
	auto& db = datablock();
	
	tuple<string,int,string,string> result;
	
	if (monID == "HOH")
	{
		for (auto r: db["pdbx_nonpoly_scheme"].find(
							cif::Key("asym_id") == asymID and
							cif::Key("pdb_seq_num") == authSeqID and
							cif::Key("mon_id") == monID))
		{
			result = r.get("pdb_strand_id", "pdb_seq_num", "pdb_mon_id", "pdb_ins_code");
			break;
		}
	}
	else
	{
		for (auto r: db["pdbx_poly_seq_scheme"].find(
							cif::Key("asym_id") == asymID and
							cif::Key("seq_id") == seqID and
							cif::Key("mon_id") == monID))
		{
			result = r.get("pdb_strand_id", "pdb_seq_num", "pdb_mon_id", "pdb_ins_code");
			break;
		}

		for (auto r: db["pdbx_nonpoly_scheme"].find(
							cif::Key("asym_id") == asymID and
							cif::Key("mon_id") == monID))
		{
			result = r.get("pdb_strand_id", "pdb_seq_num", "pdb_mon_id", "pdb_ins_code");
			break;
		}
	}

	return result;
}

tuple<string,int,string> Structure::MapPDBToLabel(const string& asymID, int seqID,
	const string& compID, const string& iCode) const
{
	auto& db = datablock();
	
	tuple<string,int,string> result;

	if (iCode.empty())
	{
		for (auto r: db["pdbx_poly_seq_scheme"].find(
							cif::Key("pdb_strand_id") == asymID and
							cif::Key("pdb_seq_num") == seqID and
							cif::Key("pdb_mon_id") == compID and
							cif::Key("pdb_ins_code") == cif::Empty()))
		{
			result = r.get("asym_id", "seq_id", "mon_id");
			break;
		}
							
		for (auto r: db["pdbx_nonpoly_scheme"].find(
							cif::Key("pdb_strand_id") == asymID and
							cif::Key("pdb_seq_num") == seqID and
							cif::Key("pdb_mon_id") == compID and
							cif::Key("pdb_ins_code") == cif::Empty()))
		{
			result = r.get("asym_id", "ndb_seq_num", "mon_id");
			break;
		}
		
	}
	else
	{
		for (auto r: db["pdbx_poly_seq_scheme"].find(
							cif::Key("pdb_strand_id") == asymID and
							cif::Key("pdb_seq_num") == seqID and
							cif::Key("pdb_mon_id") == compID and
							cif::Key("pdb_ins_code") == iCode))
		{
			result = r.get("asym_id", "seq_id", "mon_id");
			break;
		}
							
		for (auto r: db["pdbx_nonpoly_scheme"].find(
							cif::Key("pdb_strand_id") == asymID and
							cif::Key("pdb_seq_num") == seqID and
							cif::Key("pdb_mon_id") == compID and
							cif::Key("pdb_ins_code") == iCode))
		{
			result = r.get("asym_id", "ndb_seq_num", "mon_id");
			break;
		}
	}

	return result;
}

cif::Datablock& Structure::datablock() const
{
	return *mFile.impl().mDb;
}

void Structure::insertCompound(const string& compoundID, bool isEntity)
{
	auto compound = Compound::create(compoundID);
	if (compound == nullptr)
		throw runtime_error("Trying to insert unknown compound " + compoundID + " (not found in CCP4 monomers lib)");

	cif::Datablock& db = *mFile.impl().mDb;
	
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

// --------------------------------------------------------------------

Structure::residue_iterator::residue_iterator(const Structure* s, poly_iterator polyIter, size_t polyResIndex, size_t nonPolyIndex)
	: mStructure(*s), mPolyIter(polyIter), mPolyResIndex(polyResIndex), mNonPolyIndex(nonPolyIndex)
{
	while (mPolyIter != mStructure.mPolymers.end() and mPolyResIndex == mPolyIter->size())
		++mPolyIter;
}

auto Structure::residue_iterator::operator*() -> reference
{
	if (mPolyIter != mStructure.mPolymers.end())
		return (*mPolyIter)[mPolyResIndex];
	else
		return mStructure.mNonPolymers[mNonPolyIndex];
}

auto Structure::residue_iterator::operator->() -> pointer
{
	if (mPolyIter != mStructure.mPolymers.end())
		return &(*mPolyIter)[mPolyResIndex];
	else
		return &mStructure.mNonPolymers[mNonPolyIndex];
}
		
Structure::residue_iterator& Structure::residue_iterator::operator++()
{
	if (mPolyIter != mStructure.mPolymers.end())
	{
		++mPolyResIndex;
		if (mPolyResIndex >= mPolyIter->size())
		{
			++mPolyIter;
			mPolyResIndex = 0;
		}
	}
	else
		++mNonPolyIndex;
	
	return *this;
}

Structure::residue_iterator Structure::residue_iterator::operator++(int)
{
	auto result = *this;
	operator++();
	return result;
}
		
bool Structure::residue_iterator::operator==(const Structure::residue_iterator& rhs) const
{
	return mPolyIter == rhs.mPolyIter and mPolyResIndex == rhs.mPolyResIndex and mNonPolyIndex == rhs.mNonPolyIndex;
}

bool Structure::residue_iterator::operator!=(const Structure::residue_iterator& rhs) const
{
	return mPolyIter != rhs.mPolyIter or mPolyResIndex != rhs.mPolyResIndex or mNonPolyIndex != rhs.mNonPolyIndex;
}

// --------------------------------------------------------------------

void Structure::removeAtom(Atom& a)
{
	cif::Datablock& db = *mFile.impl().mDb;
	
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
	
	updateAtomIndex();
}

void Structure::swapAtoms(Atom& a1, Atom& a2)
{
	cif::Datablock& db = *mFile.impl().mDb;
	auto& atomSites = db["atom_site"];

	auto r1 = atomSites.find(cif::Key("id") == a1.id());
	auto r2 = atomSites.find(cif::Key("id") == a2.id());
	
	if (r1.size() != 1)
		throw runtime_error("Cannot swap atoms since the number of atoms with id " + a1.id() + " is " + to_string(r1.size()));
	
	if (r2.size() != 1)
		throw runtime_error("Cannot swap atoms since the number of atoms with id " + a2.id() + " is " + to_string(r2.size()));
	
	auto l1 = r1.front()["label_atom_id"];
	auto l2 = r2.front()["label_atom_id"];
	l1.swap(l2);
	
	a1.impl()->swapAtomLabels(*a2.impl());

	auto l3 = r1.front()["auth_atom_id"];
	auto l4 = r2.front()["auth_atom_id"];
	l3.swap(l4);
}

void Structure::moveAtom(Atom& a, Point p)
{
	a.location(p);
}

void Structure::changeResidue(const Residue& res, const string& newCompound,
		const vector<tuple<string,string>>& remappedAtoms)
{
	cif::Datablock& db = *mFile.impl().mDb;

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
		
		auto i = find_if(atoms.begin(), atoms.end(), [&](const Atom& a) { return a.labelAtomID() == a1; });
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

}
