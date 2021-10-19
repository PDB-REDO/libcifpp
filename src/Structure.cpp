/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "cif++/Structure.hpp"

#include <filesystem>
#include <fstream>
#include <numeric>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "cif++/Cif2PDB.hpp"
#include "cif++/CifParser.hpp"
#include "cif++/PDB2Cif.hpp"
// #include "cif++/AtomShape.hpp"

namespace fs = std::filesystem;
namespace ba = boost::algorithm;
namespace io = boost::iostreams;

extern int cif::VERBOSE;

namespace mmcif
{

// --------------------------------------------------------------------
// FileImpl

struct FileImpl
{
	cif::File mData;
	cif::Datablock *mDb = nullptr;

	void load_data(const char *data, size_t length);

	void load(const std::filesystem::path &path);
	void save(const std::filesystem::path &path);
};

void FileImpl::load_data(const char *data, size_t length)
{
	bool gzipped = length > 2 and data[0] == static_cast<char>(0x1f) and data[1] == static_cast<char>(0x8b);

	try
	{
		// First try mmCIF
		struct membuf : public std::streambuf
		{
			membuf(char *data, size_t length) { this->setg(data, data, data + length); }
		} buffer(const_cast<char *>(data), length);

		std::istream is(&buffer);

		io::filtering_stream<io::input> in;
		if (gzipped)
			in.push(io::gzip_decompressor());
		in.push(is);

		mData.load(in);
	}
	catch (const cif::CifParserError &e)
	{
		// First try mmCIF
		struct membuf : public std::streambuf
		{
			membuf(char *data, size_t length) { this->setg(data, data, data + length); }
		} buffer(const_cast<char *>(data), length);

		std::istream is(&buffer);

		io::filtering_stream<io::input> in;
		if (gzipped)
			in.push(io::gzip_decompressor());
		in.push(is);

		ReadPDBFile(in, mData);
	}

	// Yes, we've parsed the data. Now locate the datablock.
	mDb = &mData.firstDatablock();

	// And validate, otherwise lots of functionality won't work
	//	if (mData.getValidator() == nullptr)
	mData.loadDictionary("mmcif_pdbx_v50");
	if (not mData.isValid())
		std::cerr << "Invalid mmCIF file" << (cif::VERBOSE ? "." : " use --verbose option to see errors") << std::endl;
}

void FileImpl::load(const std::filesystem::path &path)
{
	std::ifstream inFile(path, std::ios_base::in | std::ios_base::binary);
	if (not inFile.is_open())
		throw std::runtime_error("No such file: " + path.string());

	io::filtering_stream<io::input> in;
	std::string ext = path.extension().string();

	if (path.extension() == ".gz")
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
					std::cerr << "unrecognized file extension, trying cif" << std::endl;

				mData.load(in);
			}
			catch (const cif::CifParserError &e)
			{
				if (cif::VERBOSE)
					std::cerr << "Not cif, trying plain old PDB" << std::endl;

				// pffft...
				in.reset();

				if (inFile.is_open())
					inFile.seekg(0);
				else
					inFile.open(path, std::ios_base::in | std::ios::binary);

				if (path.extension() == ".gz")
					in.push(io::gzip_decompressor());

				in.push(inFile);

				ReadPDBFile(in, mData);
			}
		}
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Error trying to load file " << path << std::endl;
		throw;
	}

	// Yes, we've parsed the data. Now locate the datablock.
	mDb = &mData.firstDatablock();

	// And validate, otherwise lots of functionality won't work
	//	if (mData.getValidator() == nullptr)
	mData.loadDictionary("mmcif_pdbx_v50");
	if (not mData.isValid())
		std::cerr << "Invalid mmCIF file" << (cif::VERBOSE ? "." : " use --verbose option to see errors") << std::endl;
}

void FileImpl::save(const std::filesystem::path &path)
{
	std::ofstream outFile(path, std::ios_base::out | std::ios_base::binary);
	io::filtering_stream<io::output> out;

	if (path.extension() == ".gz")
		out.push(io::gzip_compressor());

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
	AtomImpl(const AtomImpl &i)
		: mDb(i.mDb)
		, mID(i.mID)
		, mType(i.mType)
		, mAtomID(i.mAtomID)
		, mCompID(i.mCompID)
		, mAsymID(i.mAsymID)
		, mSeqID(i.mSeqID)
		, mAltID(i.mAltID)
		, mLocation(i.mLocation)
		, mRefcount(1)
		, mRow(i.mRow)
		, mCachedRefs(i.mCachedRefs)
		, mCompound(i.mCompound)
		, mRadius(i.mRadius)
		, mSymmetryCopy(i.mSymmetryCopy)
		, mClone(true)
	// , mRTop(i.mRTop), mD(i.mD)
	{
	}

	AtomImpl(cif::Datablock &db, const std::string &id)
		: mDb(db)
		, mID(id)
		, mRefcount(1)
		, mCompound(nullptr)
	{
		auto &cat = db["atom_site"];

		mRow = cat[cif::Key("id") == mID];

		prefetch();
	}

	AtomImpl(cif::Datablock &db, cif::Row &row)
		: mDb(db)
		, mID(row["id"].as<std::string>())
		, mRefcount(1)
		, mRow(row)
		, mCompound(nullptr)
	{
		prefetch();
	}

	AtomImpl(cif::Datablock &db, const std::string &id, cif::Row row)
		: mDb(db)
		, mID(id)
		, mRefcount(1)
		, mRow(row)
		, mCompound(nullptr)
	{
		prefetch();
	}

	AtomImpl(const AtomImpl &impl, const Point &loc, const std::string &sym_op)
		: mDb(impl.mDb)
		, mID(impl.mID)
		, mType(impl.mType)
		, mAtomID(impl.mAtomID)
		, mCompID(impl.mCompID)
		, mAsymID(impl.mAsymID)
		, mSeqID(impl.mSeqID)
		, mAltID(impl.mAltID)
		, mLocation(loc)
		, mRefcount(1)
		, mRow(impl.mRow)
		, mCachedRefs(impl.mCachedRefs)
		, mCompound(impl.mCompound)
		, mRadius(impl.mRadius)
		, mSymmetryCopy(true)
		, mSymmetryOperator(sym_op)
	{
	}

	void prefetch()
	{
		// Prefetch some data
		std::string symbol;
		cif::tie(symbol, mAtomID, mCompID, mAsymID, mSeqID, mAltID) =
			mRow.get("type_symbol", "label_atom_id", "label_comp_id", "label_asym_id", "label_seq_id", "label_alt_id");

		if (symbol != "X")
			mType = AtomTypeTraits(symbol).type();

		float x, y, z;
		cif::tie(x, y, z) = mRow.get("Cartn_x", "Cartn_y", "Cartn_z");

		mLocation = Point(x, y, z);

		std::string compID;
		cif::tie(compID) = mRow.get("label_comp_id");

		// mCompound = CompoundFactory::instance().create(compID);
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
		bool result = false;

		auto cat = mDb.get("atom_site_anisotrop");
		if (cat)
		{
			auto r = cat->find1(cif::Key("id") == mID);

			if (not r.empty())
			{
				result = true;
				cif::tie(anisou[0], anisou[1], anisou[2], anisou[3], anisou[4], anisou[5]) =
					r.get("U[1][1]", "U[1][2]", "U[1][3]", "U[2][2]", "U[2][3]", "U[3][3]");
			}
		}

		return result;
	}

	void moveTo(const Point &p)
	{
		assert(not mSymmetryCopy);
		if (mSymmetryCopy)
			throw std::runtime_error("Moving symmetry copy");

		if (not mClone)
		{
			property("Cartn_x", std::to_string(p.getX()));
			property("Cartn_y", std::to_string(p.getY()));
			property("Cartn_z", std::to_string(p.getZ()));
		}

		//		boost::format kPosFmt("%.3f");
		//
		//		mRow["Cartn_x"] = (kPosFmt % p.getX()).str();
		//		mRow["Cartn_y"] = (kPosFmt % p.getY()).str();
		//		mRow["Cartn_z"] = (kPosFmt % p.getZ()).str();

		mLocation = p;
	}

	const Compound &comp() const
	{
		if (mCompound == nullptr)
		{
			std::string compID;
			cif::tie(compID) = mRow.get("label_comp_id");

			mCompound = CompoundFactory::instance().create(compID);

			if (cif::VERBOSE and mCompound == nullptr)
				std::cerr << "Compound not found: '" << compID << '\'' << std::endl;
		}

		if (mCompound == nullptr)
			throw std::runtime_error("no compound");

		return *mCompound;
	}

	bool isWater() const
	{
		// mCompound may still be null here, and besides, this check is not that exciting anyway
		return mCompID == "HOH" or mCompID == "H2O" or mCompID == "WAT";
	}

	float radius() const
	{
		return mRadius;
	}

	const std::string property(const std::string_view name) const
	{
		for (auto &&[tag, ref] : mCachedRefs)
		{
			if (tag == name)
				return ref.as<std::string>();
		}

		mCachedRefs.emplace_back(name, mRow[name]);
		return std::get<1>(mCachedRefs.back()).as<std::string>();
	}

	void property(const std::string_view name, const std::string &value)
	{
		for (auto &&[tag, ref] : mCachedRefs)
		{
			if (tag != name)
				continue;
			
			ref = value;
			return;
		}

		mCachedRefs.emplace_back(name, mRow[name]);
		std::get<1>(mCachedRefs.back()) = value;
	}

	int compare(const AtomImpl &b) const
	{
		int d = mAsymID.compare(b.mAsymID);
		if (d == 0)
			d = mSeqID - b.mSeqID;
		if (d == 0)
			d = mAtomID.compare(b.mAtomID);
		return d;
	}

	void swapAtomLabels(AtomImpl &b)
	{
		std::swap(mAtomID, b.mAtomID);
	}

	const cif::Datablock &mDb;
	std::string mID;
	AtomType mType;

	std::string mAtomID;
	std::string mCompID;
	std::string mAsymID;
	int mSeqID;
	std::string mAltID;

	Point mLocation;
	int mRefcount;
	cif::Row mRow;

	mutable std::vector<std::tuple<std::string,cif::detail::ItemReference>> mCachedRefs;

	mutable const Compound *mCompound = nullptr;
	float mRadius = std::nanf("4");

	bool mSymmetryCopy = false;
	bool mClone = false;

	std::string mSymmetryOperator = "1_555";
	// clipper::RTop_orth	mRTop;
	// Point				mD;
	// int32_t				mRTix;
};

//Atom::Atom(const File& f, const std::string& id)
//	: mImpl(new AtomImpl(f, id))
//{
//}
//

Atom::Atom()
	: mImpl_(nullptr)
{
}

Atom::Atom(AtomImpl *impl)
	: mImpl_(impl)
{
}

Atom::Atom(cif::Datablock &db, cif::Row &row)
	: mImpl_(new AtomImpl(db, row))
{
}

AtomImpl *Atom::impl()
{
	if (mImpl_ == nullptr)
		throw std::runtime_error("atom is not set");
	return mImpl_;
}

const AtomImpl *Atom::impl() const
{
	if (mImpl_ == nullptr)
		throw std::runtime_error("atom is not set");
	return mImpl_;
}

Atom Atom::clone() const
{
	return Atom(mImpl_ ? new AtomImpl(*mImpl_) : nullptr);
}

Atom::Atom(const Atom &rhs, const Point &loc, const std::string &sym_op)
	: mImpl_(new AtomImpl(*rhs.mImpl_, loc, sym_op))
{
}

Atom::Atom(const Atom &rhs)
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

Atom &Atom::operator=(const Atom &rhs)
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

const cif::Row Atom::getRow() const
{
	return mImpl_->mRow;
}

const cif::Row Atom::getRowAniso() const
{
	auto &db = mImpl_->mDb;
	auto cat = db.get("atom_site_anisotrop");
	if (not cat)
		return {};
	else
		return cat->find1(cif::Key("id") == mImpl_->mID);
}

template <>
std::string Atom::property<std::string>(const std::string_view name) const
{
	return impl()->property(name);
}

template <>
int Atom::property<int>(const std::string_view name) const
{
	auto v = impl()->property(name);
	return v.empty() ? 0 : stoi(v);
}

template <>
float Atom::property<float>(const std::string_view name) const
{
	return stof(impl()->property(name));
}

void Atom::property(const std::string_view name, const std::string &value)
{
	impl()->property(name, value);
}

const std::string &Atom::id() const
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

float Atom::uIso() const
{
	float result;

	if (not property<std::string>("U_iso_or_equiv").empty())
		result = property<float>("U_iso_or_equiv");
	else if (not property<std::string>("B_iso_or_equiv").empty())
		result = property<float>("B_iso_or_equiv") / static_cast<float>(8 * kPI * kPI);
	else
		throw std::runtime_error("Missing B_iso or U_iso");

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

std::string Atom::labelAtomID() const
{
	return impl()->mAtomID;
}

std::string Atom::labelCompID() const
{
	return impl()->mCompID;
}

std::string Atom::labelAsymID() const
{
	return impl()->mAsymID;
}

std::string Atom::labelEntityID() const
{
	return property<std::string>("label_entity_id");
}

std::string Atom::labelAltID() const
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

std::string Atom::authAsymID() const
{
	return property<std::string>("auth_asym_id");
}

std::string Atom::authAtomID() const
{
	return property<std::string>("auth_atom_id");
}

std::string Atom::pdbxAuthAltID() const
{
	return property<std::string>("pdbx_auth_alt_id");
}

std::string Atom::pdbxAuthInsCode() const
{
	return property<std::string>("pdbx_PDB_ins_code");
}

std::string Atom::authCompID() const
{
	return property<std::string>("auth_comp_id");
}

std::string Atom::authSeqID() const
{
	return property<std::string>("auth_seq_id");
}

std::string Atom::labelID() const
{
	return property<std::string>("label_comp_id") + '_' + impl()->mAsymID + '_' + std::to_string(impl()->mSeqID) + ':' + impl()->mAtomID;
}

std::string Atom::pdbID() const
{
	return property<std::string>("auth_comp_id") + '_' +
	       property<std::string>("auth_asym_id") + '_' +
	       property<std::string>("auth_seq_id") +
	       property<std::string>("pdbx_PDB_ins_code");
}

Point Atom::location() const
{
	return impl()->mLocation;
}

void Atom::location(Point p)
{
	impl()->moveTo(p);
}

void Atom::translate(Point t)
{
	auto loc = location();
	loc += t;
	location(loc);
}

void Atom::rotate(Quaternion q)
{
	auto loc = location();
	loc.rotate(q);
	location(loc);
}

// Atom Atom::symmetryCopy(const Point& d, const clipper::RTop_orth& rt)
// {
// 	return Atom(new AtomImpl(*impl(), d, rt));
// }

// bool Atom::isSymmetryCopy() const
// {
// 	return impl()->mSymmetryCopy;
// }

// std::string Atom::symmetry() const
// {
// 	return clipper::Symop(impl()->mRTop).format() + "\n" + impl()->mRTop.format();
// }

bool Atom::isSymmetryCopy() const
{
	return mImpl_->mSymmetryCopy;
}

std::string Atom::symmetry() const
{
	return mImpl_->mSymmetryOperator;
}

// const clipper::RTop_orth& Atom::symop() const
// {
// 	return impl()->mRTop;
// }

const Compound &Atom::comp() const
{
	return impl()->comp();
}

bool Atom::isWater() const
{
	return impl()->isWater();
}

bool Atom::operator==(const Atom &rhs) const
{
	return impl() == rhs.impl() or
	       (&impl()->mDb == &rhs.impl()->mDb and impl()->mID == rhs.impl()->mID);
}

// clipper::Atom Atom::toClipper() const
// {
// 	return impl()->toClipper();
// }

// void Atom::calculateRadius(float resHigh, float resLow, float perc)
// {
// 	AtomShape shape(*this, resHigh, resLow, false);
// 	impl()->mRadius = shape.radius();

// 	// verbose
// 	if (cif::VERBOSE > 1)
// 		cout << "Calculated radius for " << AtomTypeTraits(impl()->mType).name() << " with charge " << charge() << " is " << impl()->mRadius << std::endl;
// }

float Atom::radius() const
{
	return impl()->mRadius;
}

int Atom::compare(const Atom &b) const
{
	return impl() == b.impl() ? 0 : impl()->compare(*b.impl());
}

void Atom::setID(int id)
{
	impl()->mID = std::to_string(id);
}

std::ostream &operator<<(std::ostream &os, const Atom &atom)
{
	os << atom.labelCompID() << ' ' << atom.labelAsymID() << ':' << atom.labelSeqID() << ' ' << atom.labelAtomID();

	if (atom.isAlternate())
		os << '(' << atom.labelAltID() << ')';
	if (atom.authAsymID() != atom.labelAsymID() or atom.authSeqID() != std::to_string(atom.labelSeqID()) or atom.pdbxAuthInsCode().empty() == false)
		os << " [" << atom.authAsymID() << ':' << atom.authSeqID() << atom.pdbxAuthInsCode() << ']';

	return os;
}

// --------------------------------------------------------------------
// residue

// First constructor used to be for waters only, but now accepts sugars as well.

Residue::Residue(const Structure &structure, const std::string &compoundID,
	const std::string &asymID, const std::string &authSeqID)
	: mStructure(&structure)
	, mCompoundID(compoundID)
	, mAsymID(asymID)
	, mAuthSeqID(authSeqID)
{
	for (auto &a : mStructure->atoms())
	{
		if (a.labelAsymID() != mAsymID or
			a.labelCompID() != mCompoundID)
			continue;

		if (compoundID == "HOH")
		{
			if (not mAuthSeqID.empty() and a.authSeqID() != mAuthSeqID)
				continue;
		}
		else
		{
			if (mSeqID > 0 and a.labelSeqID() != mSeqID)
				continue;
		}

		mAtoms.push_back(a);
	}

	assert(not mAtoms.empty());
}

Residue::Residue(const Structure &structure, const std::string &compoundID, const std::string &asymID)
	: Residue(structure, compoundID, asymID, 0, {})
{
}

Residue::Residue(const Structure &structure, const std::string &compoundID,
	const std::string &asymID, int seqID, const std::string &authSeqID)
	: mStructure(&structure)
	, mCompoundID(compoundID)
	, mAsymID(asymID)
	, mSeqID(seqID)
	, mAuthSeqID(authSeqID)
{
	assert(mCompoundID != "HOH");

	for (auto &a : mStructure->atoms())
	{
		if (mSeqID > 0 and a.labelSeqID() != mSeqID)
			continue;

		if (a.labelAsymID() != mAsymID or
			a.labelCompID() != mCompoundID)
			continue;

		mAtoms.push_back(a);
	}
}

Residue::Residue(Residue &&rhs)
	: mStructure(rhs.mStructure)
	, mCompoundID(std::move(rhs.mCompoundID))
	, mAsymID(std::move(rhs.mAsymID))
	, mSeqID(rhs.mSeqID)
	, mAuthSeqID(rhs.mAuthSeqID)
	, mAtoms(std::move(rhs.mAtoms))
{
	//std::cerr << "move constructor residue" << std::endl;
	rhs.mStructure = nullptr;
}

Residue &Residue::operator=(Residue &&rhs)
{
	//std::cerr << "move assignment residue" << std::endl;
	mStructure = rhs.mStructure;
	rhs.mStructure = nullptr;
	mCompoundID = std::move(rhs.mCompoundID);
	mAsymID = std::move(rhs.mAsymID);
	mSeqID = rhs.mSeqID;
	mAuthSeqID = rhs.mAuthSeqID;
	mAtoms = std::move(rhs.mAtoms);

	return *this;
}

Residue::~Residue()
{
	//std::cerr << "~Residue" << std::endl;
}

std::string Residue::entityID() const
{
	return mAtoms.empty() ? "" : mAtoms.front().labelEntityID();
}

std::string Residue::authInsCode() const
{
	assert(mStructure);

	std::string result;

	try
	{
		char iCode;
		tie(std::ignore, std::ignore, iCode) = mStructure->MapLabelToAuth(mAsymID, mSeqID);

		result = std::string{iCode};
		ba::trim(result);
	}
	catch (...)
	{
	}

	return result;
}

std::string Residue::authAsymID() const
{
	assert(mStructure);

	std::string result;

	try
	{
		tie(result, std::ignore, std::ignore) = mStructure->MapLabelToAuth(mAsymID, mSeqID);
	}
	catch (...)
	{
		result = mAsymID;
	}

	return result;
}

std::string Residue::authSeqID() const
{
	assert(mStructure);

	std::string result;

	try
	{
		int seqID;
		tie(std::ignore, seqID, std::ignore) = mStructure->MapLabelToAuth(mAsymID, mSeqID);
		result = std::to_string(seqID);
	}
	catch (...)
	{
	}

	return result;
}

const Compound &Residue::compound() const
{
	auto result = CompoundFactory::instance().create(mCompoundID);
	if (result == nullptr)
		throw std::runtime_error("Failed to create compound " + mCompoundID);
	return *result;
}

const AtomView &Residue::atoms() const
{
	if (mStructure == nullptr)
		throw std::runtime_error("Invalid Residue object");

	return mAtoms;
}

std::string Residue::unique_alt_id() const
{
	if (mStructure == nullptr)
		throw std::runtime_error("Invalid Residue object");

	auto firstAlt = std::find_if(mAtoms.begin(), mAtoms.end(), [](auto &a)
		{ return not a.labelAltID().empty(); });

	return firstAlt != mAtoms.end() ? firstAlt->labelAltID() : "";
}

AtomView Residue::unique_atoms() const
{
	if (mStructure == nullptr)
		throw std::runtime_error("Invalid Residue object");

	AtomView result;
	std::string firstAlt;

	for (auto &atom : mAtoms)
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

std::set<std::string> Residue::getAlternateIDs() const
{
	std::set<std::string> result;

	for (auto a : mAtoms)
	{
		auto alt = a.labelAltID();
		if (not alt.empty())
			result.insert(alt);
	}

	return result;
}

Atom Residue::atomByID(const std::string &atomID) const
{
	Atom result;

	for (auto &a : mAtoms)
	{
		if (a.labelAtomID() == atomID)
		{
			result = a;
			break;
		}
	}

	if (not result and cif::VERBOSE > 1)
		std::cerr << "Atom with atom_id " << atomID << " not found in residue " << mAsymID << ':' << mSeqID << std::endl;

	return result;
}

// Residue is a single entity if the atoms for the asym with mAsymID is equal
// to the number of atoms in this residue...  hope this is correct....
bool Residue::isEntity() const
{
	auto &db = mStructure->datablock();

	auto a1 = db["atom_site"].find(cif::Key("label_asym_id") == mAsymID);
	//	auto a2 = atoms();
	auto &a2 = mAtoms;

	return a1.size() == a2.size();
}

std::string Residue::authID() const
{
	std::string result;

	try
	{
		char chainID, iCode;
		int seqNum;

		std::tie(chainID, seqNum, iCode) = mStructure->MapLabelToAuth(mAsymID, mSeqID);

		result = chainID + std::to_string(seqNum);
		if (iCode != ' ' and iCode != 0)
			result += iCode;
	}
	catch (...)
	{
		result = mAsymID + std::to_string(mSeqID);
	}

	return result;
}

std::string Residue::labelID() const
{
	if (mCompoundID == "HOH")
		return mAsymID + mAuthSeqID;
	else
		return mAsymID + std::to_string(mSeqID);
}

std::tuple<Point, float> Residue::centerAndRadius() const
{
	std::vector<Point> pts;
	for (auto &a : mAtoms)
		pts.push_back(a.location());

	auto center = Centroid(pts);
	float radius = 0;

	for (auto &pt : pts)
	{
		float d = static_cast<float>(Distance(pt, center));
		if (radius < d)
			radius = d;
	}

	return std::make_tuple(center, radius);
}

bool Residue::hasAlternateAtoms() const
{
	return std::find_if(mAtoms.begin(), mAtoms.end(), [](const Atom &atom)
			   { return atom.isAlternate(); }) != mAtoms.end();
}

std::set<std::string> Residue::getAtomIDs() const
{
	std::set<std::string> ids;
	for (auto a : mAtoms)
		ids.insert(a.labelAtomID());

	return ids;
}

AtomView Residue::getAtomsByID(const std::string &atomID) const
{
	AtomView atoms;
	for (auto a : mAtoms)
	{
		if (a.labelAtomID() == atomID)
			atoms.push_back(a);
	}
	return atoms;
}

std::ostream &operator<<(std::ostream &os, const Residue &res)
{
	os << res.compoundID() << ' ' << res.asymID() << ':' << res.seqID();

	if (res.authAsymID() != res.asymID() or res.authSeqID() != std::to_string(res.seqID()))
		os << " [" << res.authAsymID() << ':' << res.authSeqID() << ']';

	return os;
}

// --------------------------------------------------------------------
// monomer

//Monomer::Monomer(Monomer&& rhs)
//	: Residue(std::move(rhs)), mPolymer(rhs.mPolymer), mIndex(rhs.mIndex)
//{
//}

Monomer::Monomer(const Polymer &polymer, size_t index, int seqID, const std::string &authSeqID, const std::string &compoundID)
	: Residue(*polymer.structure(), compoundID, polymer.asymID(), seqID, authSeqID)
	, mPolymer(&polymer)
	, mIndex(index)
{
}

Monomer::Monomer(Monomer &&rhs)
	: Residue(std::move(rhs))
	, mPolymer(rhs.mPolymer)
	, mIndex(rhs.mIndex)
{
	std::cerr << "move constructor monomer" << std::endl;

	//	mStructure = rhs.mStructure;			rhs.mStructure = nullptr;
	//	mCompoundID = std::move(rhs.mCompoundID);
	//	mAsymID = std::move(rhs.mAsymID);
	//	mSeqID = rhs.mSeqID;
	//	mAtoms = std::move(rhs.mAtoms);
	//
	//	mPolymer = rhs.mPolymer; rhs.mPolymer = nullptr;
	//	mIndex = rhs.mIndex;
	rhs.mPolymer = nullptr;
}

Monomer &Monomer::operator=(Monomer &&rhs)
{
	std::cerr << "move assignment monomer" << std::endl;

	Residue::operator=(std::move(rhs));
	mPolymer = rhs.mPolymer;
	rhs.mPolymer = nullptr;
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
	return mIndex >= 1 and mIndex + 2 < mPolymer->size();
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
			auto &prev = mPolymer->operator[](mIndex - 1);
			if (prev.mSeqID + 1 == mSeqID)
				result = static_cast<float>(DihedralAngle(prev.C().location(), N().location(), CAlpha().location(), C().location()));
		}
	}
	catch (const std::exception &ex)
	{
		if (cif::VERBOSE)
			std::cerr << ex.what() << std::endl;
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
			auto &next = mPolymer->operator[](mIndex + 1);
			if (mSeqID + 1 == next.mSeqID)
				result = static_cast<float>(DihedralAngle(N().location(), CAlpha().location(), C().location(), next.N().location()));
		}
	}
	catch (const std::exception &ex)
	{
		if (cif::VERBOSE)
			std::cerr << ex.what() << std::endl;
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
			auto &prev = mPolymer->operator[](mIndex - 1);
			auto &next = mPolymer->operator[](mIndex + 1);
			auto &nextNext = mPolymer->operator[](mIndex + 2);

			result = static_cast<float>(DihedralAngle(prev.CAlpha().location(), CAlpha().location(), next.CAlpha().location(), nextNext.CAlpha().location()));
		}
	}
	catch (const std::exception &ex)
	{
		if (cif::VERBOSE)
			std::cerr << ex.what() << std::endl;
	}

	return result;
}

float Monomer::kappa() const
{
	float result = 360;

	try
	{
		if (mIndex >= 2 and mIndex + 2 < mPolymer->size())
		{
			auto &prevPrev = mPolymer->operator[](mIndex - 2);
			auto &nextNext = mPolymer->operator[](mIndex + 2);

			if (prevPrev.mSeqID + 4 == nextNext.mSeqID)
			{
				double ckap = CosinusAngle(CAlpha().location(), prevPrev.CAlpha().location(), nextNext.CAlpha().location(), CAlpha().location());
				double skap = std::sqrt(1 - ckap * ckap);
				result = static_cast<float>(std::atan2(skap, ckap) * 180 / kPI);
			}
		}
	}
	catch (const std::exception &ex)
	{
		if (cif::VERBOSE)
			std::cerr << "When trying to calculate kappa for " << asymID() << ':' << seqID() << ": "
					  << ex.what() << std::endl;
	}

	return result;
}

float Monomer::tco() const
{
	float result = 0.0;

	try
	{
		if (mIndex > 0)
		{
			auto &prev = mPolymer->operator[](mIndex - 1);
			if (prev.mSeqID + 1 == mSeqID)
				result = static_cast<float>(CosinusAngle(C().location(), O().location(), prev.C().location(), prev.O().location()));
		}
	}
	catch (const std::exception &ex)
	{
		if (cif::VERBOSE)
			std::cerr << "When trying to calculate tco for " << asymID() << ':' << seqID() << ": "
					  << ex.what() << std::endl;
	}

	return result;
}

float Monomer::omega() const
{
	float result = 360;

	try
	{
		if (not is_last_in_chain())
			result = omega(*this, mPolymer->operator[](mIndex + 1));
	}
	catch (const std::exception &ex)
	{
		if (cif::VERBOSE)
			std::cerr << "When trying to calculate omega for " << asymID() << ':' << seqID() << ": "
					  << ex.what() << std::endl;
	}

	return result;
}

const std::map<std::string, std::vector<std::string>> kChiAtomsMap = {
	{"ASP", {"CG", "OD1"}},
	{"ASN", {"CG", "OD1"}},
	{"ARG", {"CG", "CD", "NE", "CZ"}},
	{"HIS", {"CG", "ND1"}},
	{"GLN", {"CG", "CD", "OE1"}},
	{"GLU", {"CG", "CD", "OE1"}},
	{"SER", {"OG"}},
	{"THR", {"OG1"}},
	{"LYS", {"CG", "CD", "CE", "NZ"}},
	{"TYR", {"CG", "CD1"}},
	{"PHE", {"CG", "CD1"}},
	{"LEU", {"CG", "CD1"}},
	{"TRP", {"CG", "CD1"}},
	{"CYS", {"SG"}},
	{"ILE", {"CG1", "CD1"}},
	{"MET", {"CG", "SD", "CE"}},
	{"MSE", {"CG", "SE", "CE"}},
	{"PRO", {"CG", "CD"}},
	{"VAL", {"CG1"}}};

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
			std::vector<std::string> atoms{"N", "CA", "CB"};

			atoms.insert(atoms.end(), i->second.begin(), i->second.end());

			// in case we have a positive chiral volume we need to swap atoms
			if (chiralVolume() > 0)
			{
				if (mCompoundID == "LEU")
					atoms.back() = "CD2";
				if (mCompoundID == "VAL")
					atoms.back() = "CG2";
			}

			result = static_cast<float>(DihedralAngle(
				atomByID(atoms[nr + 0]).location(),
				atomByID(atoms[nr + 1]).location(),
				atomByID(atoms[nr + 2]).location(),
				atomByID(atoms[nr + 3]).location()));
		}
	}
	catch (const std::exception &e)
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
		auto &next = mPolymer->operator[](mIndex + 1);

		result = Monomer::isCis(*this, next);
	}

	return result;
}

bool Monomer::isComplete() const
{
	int seen = 0;
	for (auto &a : mAtoms)
	{
		if (a.labelAtomID() == "CA")
			seen |= 1;
		else if (a.labelAtomID() == "C")
			seen |= 2;
		else if (a.labelAtomID() == "N")
			seen |= 4;
		else if (a.labelAtomID() == "O")
			seen |= 8;
		// else if (a.labelAtomID() == "OXT")		seen |= 16;
	}
	return seen == 15;
}

bool Monomer::hasAlternateBackboneAtoms() const
{
	bool result = false;

	for (auto &a : mAtoms)
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

bool Monomer::areBonded(const Monomer &a, const Monomer &b, float errorMargin)
{
	bool result = false;

	try
	{
		Point atoms[4] = {
			a.atomByID("CA").location(),
			a.atomByID("C").location(),
			b.atomByID("N").location(),
			b.atomByID("CA").location()};

		auto distanceCACA = Distance(atoms[0], atoms[3]);
		double omega = DihedralAngle(atoms[0], atoms[1], atoms[2], atoms[3]);

		bool cis = abs(omega) <= 30.0;
		float maxCACADistance = cis ? 3.0f : 3.8f;

		result = abs(distanceCACA - maxCACADistance) < errorMargin;
	}
	catch (...)
	{
	}

	return result;
}

float Monomer::omega(const mmcif::Monomer &a, const mmcif::Monomer &b)
{
	float result = 360;

	try
	{
		result = static_cast<float>(DihedralAngle(
			a.atomByID("CA").location(),
			a.atomByID("C").location(),
			b.atomByID("N").location(),
			b.atomByID("CA").location()));
	}
	catch (...)
	{
	}

	return result;
}

bool Monomer::isCis(const mmcif::Monomer &a, const mmcif::Monomer &b)
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
//		std::string asymID, monID;
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
//	std::string compoundID;
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
//		std::string asymID, monID;
//		cif::tie(asymID, seqID, monID) =
//			polySeq[mIndex].get("asym_id", "seq_id", "mon_id");
//
//		mCurrent = Monomer(*mPolymer, mIndex, seqID, monID, "");
//	}
//
//	return *this;
//}

//Polymer::Polymer(const Structure& s, const std::string& asymID)
//	: mStructure(const_cast<Structure*>(&s)), mAsymID(asymID)
//	, mPolySeq(s.category("pdbx_poly_seq_scheme").find(cif::Key("asym_id") == mAsymID))
//{
//	mEntityID = mPolySeq.front()["entity_id"].as<std::string>();
//
//#if DEBUG
//	for (auto r: mPolySeq)
//		assert(r["entity_id"] == mEntityID);
//#endif
//
//}

//Polymer::Polymer(Polymer&& rhs)
//	: std::vector<Monomer>(std::move(rhs))
//	, mStructure(rhs.mStructure)
//	, mEntityID(std::move(rhs.mEntityID)), mAsymID(std::move(rhs.mAsymID)), mPolySeq(std::move(rhs.mPolySeq))
//{
//	rhs.mStructure = nullptr;
//}
//
//Polymer& Polymer::operator=(Polymer&& rhs)
//{
//	std::vector<Monomer>::operator=(std::move(rhs));
//	mStructure = rhs.mStructure;			rhs.mStructure = nullptr;
//	mEntityID = std::move(rhs.mEntityID);
//	mAsymID = std::move(rhs.mAsymID);
//	mPolySeq = std::move(rhs.mPolySeq);
//	return *this;
//}

Polymer::Polymer(const Structure &s, const std::string &entityID, const std::string &asymID)
	: mStructure(const_cast<Structure *>(&s))
	, mEntityID(entityID)
	, mAsymID(asymID)
	, mPolySeq(s.category("pdbx_poly_seq_scheme"), cif::Key("asym_id") == mAsymID and cif::Key("entity_id") == mEntityID)
{
	std::map<size_t, size_t> ix;

	reserve(mPolySeq.size());

	for (auto r : mPolySeq)
	{
		int seqID;
		std::string compoundID, authSeqID;
		cif::tie(seqID, authSeqID, compoundID) = r.get("seq_id", "auth_seq_num", "mon_id");

		size_t index = size();

		// store only the first
		if (not ix.count(seqID))
		{
			ix[seqID] = index;
			emplace_back(*this, index, seqID, authSeqID, compoundID);
		}
		else if (cif::VERBOSE)
		{
			Monomer m{*this, index, seqID, authSeqID, compoundID};
			std::cerr << "Dropping alternate residue " << m << std::endl;
		}
	}
}

std::string Polymer::chainID() const
{
	return mPolySeq.front()["pdb_strand_id"].as<std::string>();
}

Monomer &Polymer::getBySeqID(int seqID)
{
	for (auto &m : *this)
		if (m.seqID() == seqID)
			return m;

	throw std::runtime_error("Monomer with seqID " + std::to_string(seqID) + " not found in polymer " + mAsymID);
}

const Monomer &Polymer::getBySeqID(int seqID) const
{
	for (auto &m : *this)
		if (m.seqID() == seqID)
			return m;

	throw std::runtime_error("Monomer with seqID " + std::to_string(seqID) + " not found in polymer " + mAsymID);
}

int Polymer::Distance(const Monomer &a, const Monomer &b) const
{
	int result = std::numeric_limits<int>::max();

	if (a.asymID() == b.asymID())
	{
		int ixa = std::numeric_limits<int>::max(), ixb = std::numeric_limits<int>::max();

		int ix = 0, f = 0;
		for (auto &m : *this)
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

File::File(const char *data, size_t length)
	: mImpl(new FileImpl)
{
	mImpl->load_data(data, length);
}

File::File(const std::filesystem::path &File)
	: mImpl(new FileImpl)
{
	load(File);
}

File::~File()
{
	delete mImpl;
}

cif::Datablock& File::createDatablock(const std::string_view name)
{
	auto db = new cif::Datablock(name);

	mImpl->mData.append(db);
	mImpl->mDb = db;

	return *mImpl->mDb;
}

void File::load(const std::filesystem::path &p)
{
	mImpl->load(p);
}

void File::save(const std::filesystem::path &file)
{
	mImpl->save(file);
}

cif::Datablock &File::data()
{
	assert(mImpl);
	assert(mImpl->mDb);

	if (mImpl == nullptr or mImpl->mDb == nullptr)
		throw std::runtime_error("No data loaded");

	return *mImpl->mDb;
}

cif::File &File::file()
{
	assert(mImpl);

	if (mImpl == nullptr)
		throw std::runtime_error("No data loaded");

	return mImpl->mData;
}

// --------------------------------------------------------------------
//	Structure

Structure::Structure(File &f, size_t modelNr, StructureOpenOptions options)
	: mFile(f)
	, mModelNr(modelNr)
{
	auto db = mFile.impl().mDb;
	if (db == nullptr)
		throw std::logic_error("Empty file!");

	auto &atomCat = (*db)["atom_site"];

	loadAtomsForModel(options);

	// Check to see if we should actually load another model?
	if (mAtoms.empty() and mModelNr == 1)
	{
		std::optional<size_t> model_nr;
		cif::tie(model_nr) = atomCat.front().get("pdbx_PDB_model_num");
		if (model_nr and *model_nr != mModelNr)
		{
			if (cif::VERBOSE)
				std::cerr << "No atoms loaded for model 1, trying model " << *model_nr << std::endl;
			mModelNr = *model_nr;
			loadAtomsForModel(options);
		}
	}

	if (mAtoms.empty())
		std::cerr << "Warning: no atoms loaded" << std::endl;
	else
		loadData();
}

void Structure::loadAtomsForModel(StructureOpenOptions options)
{
	auto db = mFile.impl().mDb;
	auto &atomCat = (*db)["atom_site"];

	for (auto &a : atomCat)
	{
		std::string id, type_symbol;
		std::optional<size_t> model_nr;

		cif::tie(id, type_symbol, model_nr) = a.get("id", "type_symbol", "pdbx_PDB_model_num");

		if (model_nr and *model_nr != mModelNr)
			continue;

		if ((options bitand StructureOpenOptions::SkipHydrogen) and type_symbol == "H")
			continue;

		mAtoms.emplace_back(new AtomImpl(*db, id, a));
	}
}


Structure::Structure(const Structure &s)
	: mFile(s.mFile)
	, mModelNr(s.mModelNr)
{
	mAtoms.reserve(s.mAtoms.size());
	for (auto &atom : s.mAtoms)
		mAtoms.emplace_back(atom.clone());

	loadData();
}

Structure::~Structure()
{
}

void Structure::loadData()
{
	updateAtomIndex();

	auto &polySeqScheme = category("pdbx_poly_seq_scheme");

	for (auto &r : polySeqScheme)
	{
		std::string asymID, entityID, seqID, monID;
		cif::tie(asymID, entityID, seqID, monID) =
			r.get("asym_id", "entity_id", "seq_id", "mon_id");

		if (mPolymers.empty() or mPolymers.back().asymID() != asymID or mPolymers.back().entityID() != entityID)
			mPolymers.emplace_back(*this, entityID, asymID);
	}

	auto &nonPolyScheme = category("pdbx_nonpoly_scheme");

	for (auto &r : nonPolyScheme)
	{
		std::string asymID, monID, pdbSeqNum;
		cif::tie(asymID, monID, pdbSeqNum) =
			r.get("asym_id", "mon_id", "pdb_seq_num");

		if (monID == "HOH")
			mNonPolymers.emplace_back(*this, monID, asymID, pdbSeqNum);
		else if (mNonPolymers.empty() or mNonPolymers.back().asymID() != asymID)
			mNonPolymers.emplace_back(*this, monID, asymID);
	}

	auto &branchScheme = category("pdbx_branch_scheme");

	for (auto &r : branchScheme)
	{
		std::string asymID, monID, num;
		cif::tie(asymID, monID, num) =
			r.get("asym_id", "mon_id", "num");

		mBranchResidues.emplace_back(*this, monID, asymID, num);
	}
}

void Structure::updateAtomIndex()
{
	mAtomIndex = std::vector<size_t>(mAtoms.size());

	iota(mAtomIndex.begin(), mAtomIndex.end(), 0);

	sort(mAtomIndex.begin(), mAtomIndex.end(), [this](size_t a, size_t b)
		{ return mAtoms[a].id() < mAtoms[b].id(); });
}

void Structure::sortAtoms()
{
	sort(mAtoms.begin(), mAtoms.end(), [](auto &a, auto &b)
		{ return a.compare(b) < 0; });

	int id = 1;
	for (auto &atom : mAtoms)
	{
		atom.setID(id);
		++id;
	}

	updateAtomIndex();
}

AtomView Structure::waters() const
{
	AtomView result;

	auto &db = datablock();

	// Get the entity id for water
	auto &entityCat = db["entity"];
	std::string waterEntityID;
	for (auto &e : entityCat)
	{
		std::string id, type;
		cif::tie(id, type) = e.get("id", "type");
		if (ba::iequals(type, "water"))
		{
			waterEntityID = id;
			break;
		}
	}

	for (auto &a : mAtoms)
	{
		if (a.property<std::string>("label_entity_id") == waterEntityID)
			result.push_back(a);
	}

	return result;
}

Atom Structure::getAtomByID(std::string id) const
{
	auto i = std::lower_bound(mAtomIndex.begin(), mAtomIndex.end(),
		id, [this](auto &a, auto &b)
		{ return mAtoms[a].id() < b; });

	//	auto i = find_if(mAtoms.begin(), mAtoms.end(),
	//		[&id](auto& a) { return a.id() == id; });

	if (i == mAtomIndex.end() or mAtoms[*i].id() != id)
		throw std::out_of_range("Could not find atom with id " + id);

	return mAtoms[*i];
}

Atom Structure::getAtomByLabel(const std::string &atomID, const std::string &asymID, const std::string &compID, int seqID, const std::string &altID)
{
	for (auto &a : mAtoms)
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

	throw std::out_of_range("Could not find atom with specified label");
}

Atom Structure::getAtomByPosition(Point p) const
{
	double distance = std::numeric_limits<double>::max();
	size_t index = std::numeric_limits<size_t>::max();

	for (size_t i = 0; i < mAtoms.size(); ++i)
	{
		auto &a = mAtoms.at(i);

		auto d = Distance(a.location(), p);
		if (d < distance)
		{
			distance = d;
			index = i;
		}
	}

	if (index < mAtoms.size())
		return mAtoms.at(index);

	return {};
}

Atom Structure::getAtomByPositionAndType(Point p, std::string_view type, std::string_view res_type) const
{
	double distance = std::numeric_limits<double>::max();
	size_t index = std::numeric_limits<size_t>::max();

	for (size_t i = 0; i < mAtoms.size(); ++i)
	{
		auto &a = mAtoms.at(i);

		if (a.labelCompID() != res_type)
			continue;
		
		if (a.labelAtomID() != type)
			continue;

		auto d = Distance(a.location(), p);
		if (d < distance)
		{
			distance = d;
			index = i;
		}
	}

	if (index < mAtoms.size())
		return mAtoms.at(index);

	return {};
}

const Residue &Structure::getResidue(const std::string &asymID, const std::string &compID, int seqID) const
{
	for (auto &poly : mPolymers)
	{
		if (poly.asymID() != asymID)
			continue;

		for (auto &res : poly)
		{
			if (res.seqID() == seqID and res.compoundID() == compID)
				return res;
		}
	}

	if (seqID == 0)
	{
		for (auto &res : mNonPolymers)
		{
			if (res.asymID() != asymID or res.compoundID() != compID)
				continue;

			return res;
		}
	}

	for (auto &res : mBranchResidues)
	{
		if (res.asymID() != asymID or res.compoundID() != compID or res.seqID() != seqID)
			continue;

		return res;
	}

	throw std::out_of_range("Could not find residue " + asymID + '/' + std::to_string(seqID));
}

const Residue &Structure::getResidue(const std::string &asymID) const
{
	for (auto &res : mNonPolymers)
	{
		if (res.asymID() != asymID)
			continue;

		return res;
	}

	throw std::out_of_range("Could not find residue " + asymID);
}

File &Structure::getFile() const
{
	return mFile;
}

cif::Category &Structure::category(std::string_view name) const
{
	auto &db = datablock();
	return db[name];
}

std::tuple<char, int, char> Structure::MapLabelToAuth(
	const std::string &asymID, int seqID) const
{
	auto &db = *getFile().impl().mDb;

	std::tuple<char, int, char> result;
	bool found = false;

	for (auto r : db["pdbx_poly_seq_scheme"].find(
			 cif::Key("asym_id") == asymID and
			 cif::Key("seq_id") == seqID))
	{
		std::string auth_asym_id, pdb_ins_code, pdb_seq_num, auth_seq_num;

		cif::tie(auth_asym_id, auth_seq_num, pdb_seq_num, pdb_ins_code) =
			r.get("pdb_strand_id", "auth_seq_num", "pdb_seq_num", "pdb_ins_code");

		if (auth_seq_num.empty())
			auth_seq_num = pdb_seq_num;

		try
		{
			result = std::make_tuple(auth_asym_id.front(), std::stoi(auth_seq_num),
				pdb_ins_code.empty() ? ' ' : pdb_ins_code.front());
		}
		catch (const std::exception &ex)
		{
			result = std::make_tuple(auth_asym_id.front(), 0,
				pdb_ins_code.empty() ? ' ' : pdb_ins_code.front());
		}

		found = true;
		break;
	}

	if (not found)
	{
		auto r = db["pdbx_nonpoly_scheme"].find(cif::Key("asym_id") == asymID);
		if (r.size() == 1)
		{
			std::string pdb_strand_id, pdb_ins_code;
			int pdb_seq_num;

			cif::tie(pdb_strand_id, pdb_seq_num, pdb_ins_code) =
				r.front().get("pdb_strand_id", "pdb_seq_num", "pdb_ins_code");

			result = std::make_tuple(pdb_strand_id.front(), pdb_seq_num,
				pdb_ins_code.empty() ? ' ' : pdb_ins_code.front());

			found = true;
		}
	}

	return result;
}

std::tuple<std::string, int, std::string, std::string> Structure::MapLabelToPDB(
	const std::string &asymID, int seqID, const std::string &monID,
	const std::string &authSeqID) const
{
	auto &db = datablock();

	std::tuple<std::string, int, std::string, std::string> result;

	if (monID == "HOH")
	{
		for (auto r : db["pdbx_nonpoly_scheme"].find(
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
		for (auto r : db["pdbx_poly_seq_scheme"].find(
				 cif::Key("asym_id") == asymID and
				 cif::Key("seq_id") == seqID and
				 cif::Key("mon_id") == monID))
		{
			result = r.get("pdb_strand_id", "pdb_seq_num", "pdb_mon_id", "pdb_ins_code");
			break;
		}

		for (auto r : db["pdbx_nonpoly_scheme"].find(
				 cif::Key("asym_id") == asymID and
				 cif::Key("mon_id") == monID))
		{
			result = r.get("pdb_strand_id", "pdb_seq_num", "pdb_mon_id", "pdb_ins_code");
			break;
		}
	}

	return result;
}

std::tuple<std::string, int, std::string> Structure::MapPDBToLabel(const std::string &asymID, int seqID,
	const std::string &compID, const std::string &iCode) const
{
	auto &db = datablock();

	std::tuple<std::string, int, std::string> result;

	if (iCode.empty())
	{
		for (auto r : db["pdbx_poly_seq_scheme"].find(
				 cif::Key("pdb_strand_id") == asymID and
				 cif::Key("pdb_seq_num") == seqID and
				 cif::Key("pdb_mon_id") == compID and
				 cif::Key("pdb_ins_code") == cif::Empty()))
		{
			result = r.get("asym_id", "seq_id", "mon_id");
			break;
		}

		for (auto r : db["pdbx_nonpoly_scheme"].find(
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
		for (auto r : db["pdbx_poly_seq_scheme"].find(
				 cif::Key("pdb_strand_id") == asymID and
				 cif::Key("pdb_seq_num") == seqID and
				 cif::Key("pdb_mon_id") == compID and
				 cif::Key("pdb_ins_code") == iCode))
		{
			result = r.get("asym_id", "seq_id", "mon_id");
			break;
		}

		for (auto r : db["pdbx_nonpoly_scheme"].find(
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

cif::Datablock &Structure::datablock() const
{
	return *mFile.impl().mDb;
}

std::string Structure::insertCompound(const std::string &compoundID, bool isEntity)
{
	using namespace cif::literals;

	auto compound = CompoundFactory::instance().create(compoundID);
	if (compound == nullptr)
		throw std::runtime_error("Trying to insert unknown compound " + compoundID + " (not found in CCD)");

	cif::Datablock &db = *mFile.impl().mDb;

	auto &chemComp = db["chem_comp"];
	auto r = chemComp.find(cif::Key("id") == compoundID);
	if (r.empty())
	{
		chemComp.emplace({
			{"id", compoundID},
			{"name", compound->name()},
			{"formula", compound->formula()},
			{"formula_weight", compound->formulaWeight()},
			{"type", compound->type()}});
	}

	std::string entity_id;

	if (isEntity)
	{
		auto &pdbxEntityNonpoly = db["pdbx_entity_nonpoly"];
		try
		{
			entity_id = pdbxEntityNonpoly.find1<std::string>("comp_id"_key == compoundID, "entity_id");
		}
		catch(const std::exception& ex)
		{
			auto &entity = db["entity"];
			entity_id = entity.getUniqueID("");

			entity.emplace({
				{"id", entity_id},
				{"type", "non-polymer"},
				{"pdbx_description", compound->name()},
				{"formula_weight", compound->formulaWeight()}});

			pdbxEntityNonpoly.emplace({
				{"entity_id", entity_id},
				{"name", compound->name()},
				{"comp_id", compoundID}});
		}
	}

	return entity_id;
}

// // --------------------------------------------------------------------

// Structure::residue_iterator::residue_iterator(const Structure* s, poly_iterator polyIter, size_t polyResIndex, size_t nonPolyIndex)
// 	: mStructure(*s), mPolyIter(polyIter), mPolyResIndex(polyResIndex), mNonPolyIndex(nonPolyIndex)
// {
// 	while (mPolyIter != mStructure.mPolymers.end() and mPolyResIndex == mPolyIter->size())
// 		++mPolyIter;
// }

// auto Structure::residue_iterator::operator*() -> reference
// {
// 	if (mPolyIter != mStructure.mPolymers.end())
// 		return (*mPolyIter)[mPolyResIndex];
// 	else
// 		return mStructure.mNonPolymers[mNonPolyIndex];
// }

// auto Structure::residue_iterator::operator->() -> pointer
// {
// 	if (mPolyIter != mStructure.mPolymers.end())
// 		return &(*mPolyIter)[mPolyResIndex];
// 	else
// 		return &mStructure.mNonPolymers[mNonPolyIndex];
// }

// Structure::residue_iterator& Structure::residue_iterator::operator++()
// {
// 	if (mPolyIter != mStructure.mPolymers.end())
// 	{
// 		++mPolyResIndex;
// 		if (mPolyResIndex >= mPolyIter->size())
// 		{
// 			++mPolyIter;
// 			mPolyResIndex = 0;
// 		}
// 	}
// 	else
// 		++mNonPolyIndex;

// 	return *this;
// }

// Structure::residue_iterator Structure::residue_iterator::operator++(int)
// {
// 	auto result = *this;
// 	operator++();
// 	return result;
// }

// bool Structure::residue_iterator::operator==(const Structure::residue_iterator& rhs) const
// {
// 	return mPolyIter == rhs.mPolyIter and mPolyResIndex == rhs.mPolyResIndex and mNonPolyIndex == rhs.mNonPolyIndex;
// }

// bool Structure::residue_iterator::operator!=(const Structure::residue_iterator& rhs) const
// {
// 	return mPolyIter != rhs.mPolyIter or mPolyResIndex != rhs.mPolyResIndex or mNonPolyIndex != rhs.mNonPolyIndex;
// }

// --------------------------------------------------------------------

void Structure::removeAtom(Atom &a)
{
	cif::Datablock &db = *mFile.impl().mDb;

	auto &atomSites = db["atom_site"];

	for (auto i = atomSites.begin(); i != atomSites.end(); ++i)
	{
		std::string id;
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

void Structure::swapAtoms(Atom &a1, Atom &a2)
{
	cif::Datablock &db = *mFile.impl().mDb;
	auto &atomSites = db["atom_site"];

	auto rs1 = atomSites.find(cif::Key("id") == a1.id());
	auto rs2 = atomSites.find(cif::Key("id") == a2.id());

	if (rs1.size() != 1)
		throw std::runtime_error("Cannot swap atoms since the number of atoms with id " + a1.id() + " is " + std::to_string(rs1.size()));

	if (rs2.size() != 1)
		throw std::runtime_error("Cannot swap atoms since the number of atoms with id " + a2.id() + " is " + std::to_string(rs2.size()));

	auto r1 = rs1.front();
	auto r2 = rs2.front();

	auto l1 = r1["label_atom_id"];
	auto l2 = r2["label_atom_id"];
	l1.swap(l2);

	a1.impl()->swapAtomLabels(*a2.impl());

	auto l3 = r1["auth_atom_id"];
	auto l4 = r2["auth_atom_id"];
	l3.swap(l4);
}

void Structure::moveAtom(Atom &a, Point p)
{
	a.location(p);
}

void Structure::changeResidue(const Residue &res, const std::string &newCompound,
	const std::vector<std::tuple<std::string, std::string>> &remappedAtoms)
{
	using namespace cif::literals;

	cif::Datablock &db = *mFile.impl().mDb;
	std::string asymID = res.asymID();

	const auto compound = CompoundFactory::instance().create(newCompound);
	if (not compound)
		throw std::runtime_error("Unknown compound " + newCompound);

	// First make sure the compound is already known or insert it.
	// And if the residue is an entity, we must make sure it exists

	std::string entityID;
	if (res.isEntity())
	{
		// create a copy of the entity first
		auto &entity = db["entity"];

		try
		{
			entityID = entity.find1<std::string>("type"_key == "non-polymer" and "pdbx_description"_key == compound->name(), "id");
		}
		catch (const std::exception &ex)
		{
			entityID = entity.getUniqueID("");
			entity.emplace({{"id", entityID},
				{"type", "non-polymer"},
				{"pdbx_description", compound->name()},
				{"formula_weight", compound->formulaWeight()}});
		}

		auto &pdbxEntityNonpoly = db["pdbx_entity_nonpoly"];
		pdbxEntityNonpoly.emplace({{"entity_id", entityID},
			{"name", compound->name()},
			{"comp_id", newCompound}});

		auto &pdbxNonPolyScheme = db["pdbx_nonpoly_scheme"];
		for (auto &nps : pdbxNonPolyScheme.find("asym_id"_key == asymID))
		{
			nps.assign("mon_id", newCompound, true);
			nps.assign("auth_mon_id", newCompound, true);
			nps.assign("entity_id", entityID, true);
		}

		// create rest
		auto &chemComp = db["chem_comp"];
		if (not chemComp.exists(cif::Key("id") == newCompound))
		{
			chemComp.emplace({{"id", newCompound},
				{"name", compound->name()},
				{"formula", compound->formula()},
				{"formula_weight", compound->formulaWeight()},
				{"type", compound->type()}});
		}

		// update the struct_asym for the new entity
		db["struct_asym"].update_value("id"_key == asymID, "entity_id", entityID);
	}
	else
		insertCompound(newCompound, false);

	auto &atomSites = db["atom_site"];
	auto atoms = res.atoms();

	for (auto &a : remappedAtoms)
	{
		std::string a1, a2;
		tie(a1, a2) = a;

		auto i = find_if(atoms.begin(), atoms.end(), [&](const Atom &a)
			{ return a.labelAtomID() == a1; });
		if (i == atoms.end())
		{
			std::cerr << "Missing atom for atom ID " << a1 << std::endl;
			continue;
		}

		auto r = atomSites.find(cif::Key("id") == i->id());

		if (r.size() != 1)
			continue;

		if (a2.empty() or a2 == ".")
			atomSites.erase(cif::Key("id") == i->id());
		else if (a1 != a2)
		{
			auto ra = r.front();
			ra["label_atom_id"] = a2;
			ra["auth_atom_id"] = a2;
			ra["type_symbol"] = AtomTypeTraits(compound->getAtomByID(a2).typeSymbol).symbol();
		}
	}

	for (auto a : atoms)
	{
		atomSites.update_value(cif::Key("id") == a.id(), "label_comp_id", newCompound);
		atomSites.update_value(cif::Key("id") == a.id(), "auth_comp_id", newCompound);
	}
}

std::string Structure::createNonPolyEntity(const std::string &comp_id)
{
	return insertCompound(comp_id, true);
}

std::string Structure::createNonpoly(const std::string &entity_id, const std::vector<mmcif::Atom> &atoms)
{
	using namespace cif::literals;

	cif::Datablock &db = *mFile.impl().mDb;
	auto &struct_asym = db["struct_asym"];
	std::string asym_id = struct_asym.getUniqueID();

	struct_asym.emplace({
		{ "id", asym_id },
		{ "pdbx_blank_PDB_chainid_flag", "N" },
		{ "pdbx_modified", "N" },
		{ "entity_id", entity_id },
		{ "details", "?" }
	});

	std::string comp_id = db["pdbx_entity_nonpoly"].find1<std::string>("entity_id"_key == entity_id, "comp_id");

	auto &atom_site = db["atom_site"];

	for (auto& atom: atoms)
	{
		auto atom_id = atom_site.getUniqueID("");

		auto &&[row, inserted ] = atom_site.emplace({
			{ "group_PDB",			atom.property<std::string>("group_PDB") },
			{ "id",					atom_id },
			{ "type_symbol",		atom.property<std::string>("type_symbol") },
			{ "label_atom_id",		atom.property<std::string>("label_atom_id") },
			{ "label_alt_id",		atom.property<std::string>("label_alt_id") },
			{ "label_comp_id",		comp_id },
			{ "label_asym_id",		asym_id },
			{ "label_entity_id",	entity_id },
			{ "label_seq_id",		"." },
			{ "pdbx_PDB_ins_code",	"" },
			{ "Cartn_x",			atom.property<std::string>("Cartn_x") },
			{ "Cartn_y",			atom.property<std::string>("Cartn_y") },
			{ "Cartn_z",			atom.property<std::string>("Cartn_z") },
			{ "occupancy",			atom.property<std::string>("occupancy") },
			{ "B_iso_or_equiv",		atom.property<std::string>("B_iso_or_equiv") },
			{ "pdbx_formal_charge",	atom.property<std::string>("pdbx_formal_charge") },
			{ "auth_seq_id",		"" },
			{ "auth_comp_id",		comp_id },
			{ "auth_asym_id",		asym_id },
			{ "auth_atom_id",		atom.property<std::string>("label_atom_id") },
			{ "pdbx_PDB_model_num",	1 }
		});

		mAtoms.emplace_back(new AtomImpl(db, atom_id, row));
	}

	mNonPolymers.emplace_back(*this, comp_id, asym_id);

	return asym_id;
}

void Structure::cleanupEmptyCategories()
{
	using namespace cif::literals;

	cif::Datablock &db = *mFile.impl().mDb;

	auto &atomSite = db["atom_site"];

	// Remove chem_comp's for which there are no atoms at all
	auto &chem_comp = db["chem_comp"];
	cif::RowSet obsoleteChemComps(chem_comp);

	for (auto chemComp : chem_comp)
	{
		std::string compID = chemComp["id"].as<std::string>();
		if (atomSite.exists("label_comp_id"_key == compID or "auth_comp_id"_key == compID))
			continue;

		obsoleteChemComps.push_back(chemComp);
	}

	for (auto chemComp : obsoleteChemComps)
		chem_comp.erase(chemComp);

	// similarly, remove entities not referenced by any atom

	auto &entities = db["entity"];
	cif::RowSet obsoleteEntities(entities);

	for (auto entity : entities)
	{
		std::string entityID = entity["id"].as<std::string>();
		if (atomSite.exists("label_entity_id"_key == entityID))
			continue;

		obsoleteEntities.push_back(entity);
	}

	for (auto entity : obsoleteEntities)
		entities.erase(entity);

	// the rest?

	for (const char *cat : {"pdbx_entity_nonpoly"})
	{
		auto &category = db[cat];

		cif::RowSet empty(category);
		for (auto row : category)
		{
			if (not category.hasChildren(row) and not category.hasParents(row))
				empty.push_back(row);
		}

		for (auto row : empty)
			category.erase(row);
	}

	// count molecules
	for (auto entity : entities)
	{
		std::string type, id;
		cif::tie(type, id) = entity.get("type", "id");

		std::optional<size_t> count;
		if (type == "polymer")
			count = db["entity_poly"].find("entity_id"_key == id).size();
		else if (type == "non-polymer" or type == "water")
			count = db["pdbx_nonpoly_scheme"].find("entity_id"_key == id).size();
		else if (type == "branched")
		{
			// is this correct?
			std::set<std::string> asym_ids;
			for (const auto &[ asym_id ] : db["pdbx_branch_scheme"].find<std::string>("entity_id"_key == id, "asym_id"))
				asym_ids.insert(asym_id);
			count = asym_ids.size();
		}

		entity["pdbx_number_of_molecules"] = count;
	}
}

void Structure::translate(Point t)
{
	for (auto& a: mAtoms)
		a.translate(t);
}

void Structure::rotate(Quaternion q)
{
	for (auto& a: mAtoms)
		a.rotate(q);
}

} // namespace mmcif
