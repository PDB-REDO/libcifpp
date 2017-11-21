// Lib for working with structures as contained in file and PDB files

#include "libcif/structure.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "pdb2cif.h"
#include "libcif/cif-parser.h"
#include "cif2pdb.h"

using namespace std;

namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;

extern int VERBOSE;

namespace libcif
{
	
// --------------------------------------------------------------------
// file_impl
	
struct file_impl
{
	cif::file			m_data;
	cif::datablock*		m_db = nullptr;
	
	void load(fs::path p);
	void save(fs::path p);
};

void file_impl::load(fs::path p)
{
	fs::ifstream infile(p, ios_base::in | ios_base::binary);
	if (not infile.is_open())
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
	
	in.push(infile);

	// OK, we've got the file, now create a protein
	if (ext == ".cif")
		m_data.load(in);
	else if (ext == ".pdb" or ext == ".ent")
		ReadPDBFile(in, m_data);
	else
	{
		try
		{
			if (VERBOSE)
				cerr << "unrecognized file extension, trying cif" << endl;

			m_data.load(in);
		}
		catch (const cif::cif_parser_error& e)
		{
			if (VERBOSE)
				cerr << "Not cif, trying plain old PDB" << endl;

			// pffft...
			in.reset();

			if (infile.is_open())
				infile.seekg(0);
			else
				infile.open(p, ios_base::in | ios::binary);

			if (p.extension() == ".bz2")
				in.push(io::bzip2_decompressor());
			else if (p.extension() == ".gz")
				in.push(io::gzip_decompressor());
			
			in.push(infile);

			ReadPDBFile(in, m_data);
		}
	}
	
	// Yes, we've parsed the data. Now locate the datablock.
	m_db = &m_data.first_datablock();
	
	// And validate, otherwise lots of functionality won't work
//	if (m_data.get_validator() == nullptr)
		m_data.load_dictionary("mmcif_pdbx");
	m_data.validate();
}

void file_impl::save(fs::path p)
{
	fs::ofstream outfile(p, ios_base::out | ios_base::binary);
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
	
	out.push(outfile);
	
	if (p.extension() == ".pdb")
		WritePDBFile(out, m_data);
	else
		m_data.save(out);
}


// --------------------------------------------------------------------
// atom

struct atom_impl
{
	atom_impl(const file& f, const string& id)
		: m_file(f), m_id(id), m_refcount(1), m_compound(nullptr)
	{
		auto& db = *m_file.impl().m_db;
		auto& cat = db["atom_site"];
		
		m_row = cat[cif::key("id") == m_id];

		prefetch();
	}
	
	atom_impl(const file& f, const string& id, cif::row row)
		: m_file(f), m_id(id), m_refcount(1), m_row(row), m_compound(nullptr)
	{
		prefetch();
	}
	
	void prefetch()
	{
		// Prefetch some data
		string symbol;
		cif::tie(symbol) = m_row.get("type_symbol");
		
		m_type = atom_type_traits(symbol).type();

		float x, y, z;
		cif::tie(x, y, z) = m_row.get("Cartn_x", "Cartn_y", "Cartn_z");
		
		m_location = point(x, y, z);
		
		try
		{
			comp();
		}
		catch (...) {}
	}

	void reference()
	{
		++m_refcount;
	}
	
	void release()
	{
		if (--m_refcount < 0)
			delete this;
	}
	
	const compound& comp()
	{
		if (m_compound == nullptr)
		{
			string comp_id;
			cif::tie(comp_id) = m_row.get("label_comp_id");
			
			m_compound = compound::create(comp_id);
		}
		
		if (m_compound == nullptr)
			throw runtime_error("no compound");
	
		return *m_compound;
	}
	
	bool is_water() const
	{
		return m_compound != nullptr and m_compound->is_water();
	}

	const file&			m_file;
	string				m_id;
	int					m_refcount;
	cif::row			m_row;
	const compound*		m_compound;
	point				m_location;
	atom_type			m_type;
	
//	const entity&		m_entity;
//	std::string			m_asym_id;
//	std::string			m_atom_id;
//	point				m_loc;
//	property_list		m_properties;
};

atom::atom(const file& f, const string& id)
	: m_impl(new atom_impl(f, id))
{
}

atom::atom(atom_impl* impl)
	: m_impl(impl)
{
}

atom::atom(const atom& rhs)
	: m_impl(rhs.m_impl)
{
	m_impl->reference();
}

atom::~atom()
{
	if (m_impl)
		m_impl->release();
}

atom& atom::operator=(const atom& rhs)
{
	if (this != &rhs)
	{
		m_impl->release();
		m_impl = rhs.m_impl;
		m_impl->reference();
	}

	return *this;
}

string atom::id() const
{
	return m_impl->m_id;
}

atom_type atom::type() const
{
	return m_impl->m_type;
}

int atom::charge() const
{
	int charge;
	cif::tie(charge) = m_impl->m_row.get("pdbx_formal_charge");
	
	return charge;
}

string atom::label_atom_id() const
{
	string atom_id;
	cif::tie(atom_id) = m_impl->m_row.get("label_atom_id");
	
	return atom_id;
}

string atom::label_comp_id() const
{
	string comp_id;
	cif::tie(comp_id) = m_impl->m_row.get("label_comp_id");
	
	return comp_id;
}

string atom::label_asym_id() const
{
	string asym_id;
	cif::tie(asym_id) = m_impl->m_row.get("label_asym_id");
	
	return asym_id;
}

int atom::label_seq_id() const
{
	int seq_id;
	cif::tie(seq_id) = m_impl->m_row.get("label_seq_id");
	
	return seq_id;
}

string atom::auth_asym_id() const
{
	string asym_id;
	cif::tie(asym_id) = m_impl->m_row.get("auth_asym_id");
	
	return asym_id;
}

int atom::auth_seq_id() const
{
	int seq_id;
	cif::tie(seq_id) = m_impl->m_row.get("auth_seq_id");
	
	return seq_id;
}

point atom::location() const
{
	return m_impl->m_location;
}

const compound& atom::comp() const
{
	return m_impl->comp();
}

bool atom::is_water() const
{
	return m_impl->is_water();
}

boost::any atom::property(const std::string& name) const
{
	string s = m_impl->m_row[name].as<string>();
	
	return boost::any(s);
}

bool atom::operator==(const atom& rhs) const
{
	return m_impl == rhs.m_impl or
		(&m_impl->m_file == &rhs.m_impl->m_file and m_impl->m_id == rhs.m_impl->m_id); 	
}

const file& atom::get_file() const
{
	assert(m_impl);
	return m_impl->m_file;
}

// --------------------------------------------------------------------
// residue

//atom_view residue::atoms()
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
	: m_impl(new file_impl)
{
}

file::file(fs::path file)
	: m_impl(new file_impl)
{
	load(file);
}

file::~file()
{
	delete m_impl;
}

void file::load(fs::path p)
{
	m_impl->load(p);
	
//	// all data is now in m_file, construct atoms and others
//	
//	auto& db = m_file.first_datablock();
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
//	auto& atom_sites = db["atom_site"];
//	for (auto& atom_site: atom_sites)
//	{
//		atom_ptr ap(new atom(this, atom_site));
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
	m_impl->save(file);
}

cif::datablock& file::data()
{
	assert(m_impl);
	assert(m_impl->m_db);
	
	if (m_impl == nullptr or m_impl->m_db == nullptr)
		throw runtime_error("No data loaded");
	
	return *m_impl->m_db;
}

// --------------------------------------------------------------------
//	structure

struct structure_impl
{
	structure_impl(structure& s, file& f, uint32 model_nr)
		: m_file(&f), m_model_nr(model_nr)
	{
		auto& db = *m_file->impl().m_db;
		auto& atom_cat = db["atom_site"];
		
		for (auto& a: atom_cat)
		{
			auto model_nr = a["pdbx_PDB_model_num"];
			
			if (model_nr.empty() or model_nr.as<uint32>() == m_model_nr)
				m_atoms.emplace_back(new atom_impl(f, a["id"].as<string>(), a));
		}
	}
	
	void remove_atom(atom& a);
	
	file*			m_file;
	uint32			m_model_nr;
	atom_view		m_atoms;
};

void structure_impl::remove_atom(atom& a)
{
	cif::datablock& db = *m_file->impl().m_db;
	
	auto& atom_sites = db["atom_site"];
	
	for (auto i = atom_sites.begin(); i != atom_sites.end(); ++i)
	{
		string id;
		cif::tie(id) = i->get("id");
		
		if (id == a.id())
		{
			atom_sites.erase(i);
			break;
		}
	}
	
	m_atoms.erase(remove(m_atoms.begin(), m_atoms.end(), a), m_atoms.end());
}

structure::structure(file& f, uint32 model_nr)
	: m_impl(new structure_impl(*this, f, model_nr))
{
}

structure::~structure()
{
	delete m_impl;
}

atom_view structure::atoms() const
{
	return m_impl->m_atoms;
}

atom_view structure::waters() const
{
	atom_view result;
	
	auto& db = *get_file().impl().m_db;
	
	// Get the entity id for water
	auto& entity_cat = db["entity"];
	string water_entity_id;
	for (auto& e: entity_cat)
	{
		string id, type;
		cif::tie(id, type) = e.get("id", "type");
		if (ba::iequals(type, "water"))
		{
			water_entity_id = id;
			break;
		}
	}

	for (auto& a: m_impl->m_atoms)
	{
		if (boost::any_cast<string>(a.property("label_entity_id")) == water_entity_id)
			result.push_back(a);
	}
	
	return result;
}

atom structure::get_atom_by_id(string id) const
{
	for (auto& a: m_impl->m_atoms)
	{
		if (a.id() == id)
			return a;
	}
	
	throw out_of_range("Could not find atom with id " + id);
}

file& structure::get_file() const
{
	return *m_impl->m_file;
}

//tuple<string,string> structure::MapLabelToAuth(
//	const string& asym_id, int seq_id)
//{
//	auto& db = *get_file().impl().m_db;
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
	const string& asym_id, int seq_id, const string& mon_id)
{
	auto& db = *get_file().impl().m_db;
	
	tuple<string,int,string,string> result;
	
	for (auto r: db["pdbx_poly_seq_scheme"].find(
						cif::key("asym_id") == asym_id and
						cif::key("seq_id") == seq_id and
						cif::key("mon_id") == mon_id))
	{
		string pdb_strand_id, pdb_mon_id, pdb_ins_code;
		int pdb_seq_num;
		
		cif::tie(pdb_strand_id, pdb_seq_num, pdb_mon_id, pdb_ins_code) =
			r.get("pdb_strand_id", "pdb_seq_num", "pdb_mon_id", "pdb_ins_code");

		result = make_tuple(pdb_strand_id, pdb_seq_num, pdb_mon_id, pdb_ins_code);

		break;
	}
						
	for (auto r: db["pdbx_nonpoly_scheme"].find(
						cif::key("asym_id") == asym_id and
						cif::key("seq_id") == seq_id and
						cif::key("mon_id") == mon_id))
	{
		string pdb_strand_id, pdb_mon_id, pdb_ins_code;
		int pdb_seq_num;
		
		cif::tie(pdb_strand_id, pdb_seq_num, pdb_mon_id, pdb_ins_code) =
			r.get("pdb_strand_id", "pdb_seq_num", "pdb_mon_id", "pdb_ins_code");

		result = make_tuple(pdb_strand_id, pdb_seq_num, pdb_mon_id, pdb_ins_code);

		break;
	}

	return result;
}

// --------------------------------------------------------------------
// actions

void structure::remove_atom(atom& a)
{
	m_impl->remove_atom(a);
}

}
