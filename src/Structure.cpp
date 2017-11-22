// Lib for working with structures as contained in file and PDB files

#include "cif++/Structure.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "cif++/PDB2Cif.h"
#include "cif++/CifParser.h"
#include "cif++/Cif2PDB.h"

using namespace std;

namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;

extern int VERBOSE;

namespace libcif
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
	string ext;
	
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
	
	// Yes, we've parsed the data. Now locate the datablock.
	mDb = &mData.firstDatablock();
	
	// And validate, otherwise lots of functionality won't work
//	if (mData.getValidator() == nullptr)
		mData.loadDictionary("mmcif_pdbx");
	mData.validate();
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
	
	void prefetch()
	{
		// Prefetch some data
		string symbol;
		cif::tie(symbol) = mRow.get("type_symbol");
		
		mType = AtomTypeTraits(symbol).type();

		float x, y, z;
		cif::tie(x, y, z) = mRow.get("Cartn_x", "Cartn_y", "Cartn_z");
		
		mLocation = Point(x, y, z);
		
		try
		{
			comp();
		}
		catch (...) {}
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
	
	const Compound& comp()
	{
		if (mCompound == nullptr)
		{
			string compId;
			cif::tie(compId) = mRow.get("label_comp_id");
			
			mCompound = Compound::create(compId);
		}
		
		if (mCompound == nullptr)
			throw runtime_error("no compound");
	
		return *mCompound;
	}
	
	bool isWater() const
	{
		return mCompound != nullptr and mCompound->isWater();
	}

	const File&			mFile;
	string				mId;
	int					mRefcount;
	cif::Row			mRow;
	const Compound*		mCompound;
	Point				mLocation;
	AtomType			mType;
	
//	const entity&		mEntity;
//	std::string			mAsymId;
//	std::string			mAtomId;
//	Point				mLoc;
//	propertyList		mProperties;
};

Atom::Atom(const File& f, const string& id)
	: mImpl(new AtomImpl(f, id))
{
}

Atom::Atom(AtomImpl* impl)
	: mImpl(impl)
{
}

Atom::Atom(const Atom& rhs)
	: mImpl(rhs.mImpl)
{
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
		mImpl->release();
		mImpl = rhs.mImpl;
		mImpl->reference();
	}

	return *this;
}

string Atom::id() const
{
	return mImpl->mId;
}

AtomType Atom::type() const
{
	return mImpl->mType;
}

int Atom::charge() const
{
	int charge;
	cif::tie(charge) = mImpl->mRow.get("pdbx_formal_charge");
	
	return charge;
}

string Atom::labelAtomId() const
{
	string atomId;
	cif::tie(atomId) = mImpl->mRow.get("label_atom_id");
	
	return atomId;
}

string Atom::labelCompId() const
{
	string compId;
	cif::tie(compId) = mImpl->mRow.get("label_comp_id");
	
	return compId;
}

string Atom::labelAsymId() const
{
	string asymId;
	cif::tie(asymId) = mImpl->mRow.get("label_asym_id");
	
	return asymId;
}

int Atom::labelSeqId() const
{
	int seqId;
	cif::tie(seqId) = mImpl->mRow.get("label_seq_id");
	
	return seqId;
}

string Atom::authAsymId() const
{
	string asymId;
	cif::tie(asymId) = mImpl->mRow.get("auth_asym_id");
	
	return asymId;
}

int Atom::authSeqId() const
{
	int seqId;
	cif::tie(seqId) = mImpl->mRow.get("auth_seq_id");
	
	return seqId;
}

Point Atom::location() const
{
	return mImpl->mLocation;
}

const Compound& Atom::comp() const
{
	return mImpl->comp();
}

bool Atom::isWater() const
{
	return mImpl->isWater();
}

boost::any Atom::property(const std::string& name) const
{
	string s = mImpl->mRow[name].as<string>();
	
	return boost::any(s);
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

// --------------------------------------------------------------------
// residue

//AtomView residue::Atoms()
//{
//	assert(false);
//}

// --------------------------------------------------------------------
// monomer

// --------------------------------------------------------------------
// polymer

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
	
	File*			mFile;
	uint32			mModelNr;
	AtomView		mAtoms;
};

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
	
	auto& db = *getFile().impl().mDb;
	
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
		if (boost::any_cast<string>(a.property("label_entity_id")) == waterEntityId)
			result.push_back(a);
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

//tuple<string,string> Structure::MapLabelToAuth(
//	const string& asymId, int seqId)
//{
//	auto& db = *getFile().impl().mDb;
//	
//	tuple<string,int,string,string> result;
//	bool found = false;
//	
//	for (auto r: db["pdbx_poly_seq_scheme"].find(
//						cif::Key("asym_id") == asym_id and
//						cif::Key("seq_id") == seq_id))
//	{
//		string auth_asym_id, pdb_mon_id, pdb_ins_code;
//		int pdb_seq_num;
//		
//		cif::tie(pdb_strand_id, pdb_seq_num, pdb_mon_id, pdb_ins_code) =
//			r.get("pdb_strand_id", "pdb_seq_num", "pdb_mon_id", "pdb_ins_code");
//
//		result = make_tuple(pdb_strand_id, pdb_seq_num, pdb_mon_id, pdb_ins_code);
//
//		found = true;
//		break;
//	}
//						
//	for (auto r: db["pdbx_nonpoly_scheme"].find(
//						cif::Key("asym_id") == asym_id and
//						cif::Key("seq_id") == seq_id and
//						cif::Key("mon_id") == mon_id))
//	{
//		string pdb_strand_id, pdb_mon_id, pdb_ins_code;
//		int pdb_seq_num;
//		
//		cif::tie(pdb_strand_id, pdb_seq_num, pdb_mon_id, pdb_ins_code) =
//			r.get("pdb_strand_id", "pdb_seq_num", "pdb_mon_id", "pdb_ins_code");
//
//		result = make_tuple(pdb_strand_id, pdb_seq_num, pdb_mon_id, pdb_ins_code);
//
//		found = true;
//		break;
//	}
//
//	return result;
//}

tuple<string,int,string,string> Structure::MapLabelToPDB(
	const string& asymId, int seqId, const string& monId)
{
	auto& db = *getFile().impl().mDb;
	
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

// --------------------------------------------------------------------
// actions

void Structure::removeAtom(Atom& a)
{
	mImpl->removeAtom(a);
}

}
