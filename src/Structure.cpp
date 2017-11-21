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
	cif::datablock*		mDb = nullptr;
	
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
		catch (const cif::cifParserError& e)
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
// atom

struct atomImpl
{
	atomImpl(const File& f, const string& id)
		: mFile(f), mId(id), mRefcount(1), mCompound(nullptr)
	{
		auto& db = *mFile.impl().mDb;
		auto& cat = db["atom_site"];
		
		mRow = cat[cif::key("id") == mId];

		prefetch();
	}
	
	atomImpl(const file& f, const string& id, cif::row row)
		: mFile(f), mId(id), mRefcount(1), mRow(row), mCompound(nullptr)
	{
		prefetch();
	}
	
	void prefetch()
	{
		// Prefetch some data
		string symbol;
		cif::tie(symbol) = mRow.get("type_symbol");
		
		mType = atomTypeTraits(symbol).type();

		float x, y, z;
		cif::tie(x, y, z) = mRow.get("Cartn_x", "Cartn_y", "Cartn_z");
		
		mLocation = point(x, y, z);
		
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
	
	const compound& comp()
	{
		if (mCompound == nullptr)
		{
			string compId;
			cif::tie(compId) = mRow.get("label_comp_id");
			
			mCompound = compound::create(compId);
		}
		
		if (mCompound == nullptr)
			throw runtime_error("no compound");
	
		return *mCompound;
	}
	
	bool isWater() const
	{
		return mCompound != nullptr and mCompound->isWater();
	}

	const file&			mFile;
	string				mId;
	int					mRefcount;
	cif::row			mRow;
	const compound*		mCompound;
	point				mLocation;
	atomType			mType;
	
//	const entity&		mEntity;
//	std::string			mAsymId;
//	std::string			mAtomId;
//	point				mLoc;
//	propertyList		mProperties;
};

atom::atom(const file& f, const string& id)
	: mImpl(new atomImpl(f, id))
{
}

atom::atom(atomImpl* impl)
	: mImpl(impl)
{
}

atom::atom(const atom& rhs)
	: mImpl(rhs.mImpl)
{
	mImpl->reference();
}

atom::~atom()
{
	if (mImpl)
		mImpl->release();
}

atom& atom::operator=(const atom& rhs)
{
	if (this != &rhs)
	{
		mImpl->release();
		mImpl = rhs.mImpl;
		mImpl->reference();
	}

	return *this;
}

string atom::id() const
{
	return mImpl->mId;
}

atomType atom::type() const
{
	return mImpl->mType;
}

int atom::charge() const
{
	int charge;
	cif::tie(charge) = mImpl->mRow.get("pdbx_formal_charge");
	
	return charge;
}

string atom::labelAtomId() const
{
	string atomId;
	cif::tie(atomId) = mImpl->mRow.get("label_atom_id");
	
	return atomId;
}

string atom::labelCompId() const
{
	string compId;
	cif::tie(compId) = mImpl->mRow.get("label_comp_id");
	
	return compId;
}

string atom::labelAsymId() const
{
	string asymId;
	cif::tie(asymId) = mImpl->mRow.get("label_asym_id");
	
	return asymId;
}

int atom::labelSeqId() const
{
	int seqId;
	cif::tie(seqId) = mImpl->mRow.get("label_seq_id");
	
	return seqId;
}

string atom::authAsymId() const
{
	string asymId;
	cif::tie(asymId) = mImpl->mRow.get("auth_asym_id");
	
	return asymId;
}

int atom::authSeqId() const
{
	int seqId;
	cif::tie(seqId) = mImpl->mRow.get("auth_seq_id");
	
	return seqId;
}

point atom::location() const
{
	return mImpl->mLocation;
}

const compound& atom::comp() const
{
	return mImpl->comp();
}

bool atom::isWater() const
{
	return mImpl->isWater();
}

boost::any atom::property(const std::string& name) const
{
	string s = mImpl->mRow[name].as<string>();
	
	return boost::any(s);
}

bool atom::operator==(const atom& rhs) const
{
	return mImpl == rhs.mImpl or
		(&mImpl->mFile == &rhs.mImpl->mFile and mImpl->mId == rhs.mImpl->mId); 	
}

const file& atom::getFile() const
{
	assert(mImpl);
	return mImpl->mFile;
}

// --------------------------------------------------------------------
// residue

//atomView residue::atoms()
//{
//	assert(false);
//}

// --------------------------------------------------------------------
// monomer

// --------------------------------------------------------------------
// polymer

// --------------------------------------------------------------------
// file

file::file()
	: mImpl(new FileImpl)
{
}

File::file(fs::path file)
	: mImpl(new FileImpl)
{
	load(file);
}

file::~file()
{
	delete mImpl;
}

void file::load(fs::path p)
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
//		atomPtr ap(new atom(this, atom_site));
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

void file::save(boost::filesystem::path file)
{
	mImpl->save(file);
}

cif::datablock& file::data()
{
	assert(mImpl);
	assert(mImpl->mDb);
	
	if (mImpl == nullptr or mImpl->mDb == nullptr)
		throw runtime_error("No data loaded");
	
	return *mImpl->mDb;
}

// --------------------------------------------------------------------
//	structure

struct structureImpl
{
	structureImpl(structure& s, file& f, uint32 modelNr)
		: mFile(&f), mModelNr(modelNr)
	{
		auto& db = *mFile->impl().mDb;
		auto& atomCat = db["atom_site"];
		
		for (auto& a: atomCat)
		{
			auto modelNr = a["pdbx_PDB_model_num"];
			
			if (modelNr.empty() or modelNr.as<uint32>() == mModelNr)
				mAtoms.emplace_back(new atomImpl(f, a["id"].as<string>(), a));
		}
	}
	
	void removeAtom(atom& a);
	
	file*			mFile;
	uint32			mModelNr;
	atomView		mAtoms;
};

void structureImpl::removeAtom(atom& a)
{
	cif::datablock& db = *mFile->impl().mDb;
	
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

structure::structure(file& f, uint32 modelNr)
	: mImpl(new structureImpl(*this, f, modelNr))
{
}

structure::~structure()
{
	delete mImpl;
}

atomView structure::atoms() const
{
	return mImpl->mAtoms;
}

atomView structure::waters() const
{
	atomView result;
	
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

atom structure::getAtomById(string id) const
{
	for (auto& a: mImpl->mAtoms)
	{
		if (a.id() == id)
			return a;
	}
	
	throw out_of_range("Could not find atom with id " + id);
}

file& structure::getFile() const
{
	return *mImpl->mFile;
}

//tuple<string,string> structure::MapLabelToAuth(
//	const string& asymId, int seqId)
//{
//	auto& db = *getFile().impl().mDb;
//	
//	tuple<string,int,string,string> result;
//	bool found = false;
//	
//	for (auto r: db["pdbx_poly_seq_scheme"].find(
//						cif::key("asym_id") == asym_id and
//						cif::key("seq_id") == seq_id))
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
//						cif::key("asym_id") == asym_id and
//						cif::key("seq_id") == seq_id and
//						cif::key("mon_id") == mon_id))
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

tuple<string,int,string,string> structure::MapLabelToPDB(
	const string& asymId, int seqId, const string& monId)
{
	auto& db = *getFile().impl().mDb;
	
	tuple<string,int,string,string> result;
	
	for (auto r: db["pdbx_poly_seq_scheme"].find(
						cif::key("asym_id") == asymId and
						cif::key("seq_id") == seqId and
						cif::key("mon_id") == monId))
	{
		string pdbStrandId, pdbMonId, pdbInsCode;
		int pdbSeqNum;
		
		cif::tie(pdbStrandId, pdbSeqNum, pdbMonId, pdbInsCode) =
			r.get("pdb_strand_id", "pdb_seq_num", "pdb_mon_id", "pdb_ins_code");

		result = make_tuple(pdbStrandId, pdbSeqNum, pdbMonId, pdbInsCode);

		break;
	}
						
	for (auto r: db["pdbx_nonpoly_scheme"].find(
						cif::key("asym_id") == asymId and
						cif::key("seq_id") == seqId and
						cif::key("mon_id") == monId))
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

void structure::removeAtom(atom& a)
{
	mImpl->removeAtom(a);
}

}
