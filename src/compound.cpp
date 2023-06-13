/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2020-2022 NKI/AVL, Netherlands Cancer Institute
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

#include "cif++.hpp"

#include <filesystem>
#include <fstream>
#include <map>
#include <mutex>
#include <numeric>
#include <shared_mutex>

namespace fs = std::filesystem;

namespace cif
{

// --------------------------------------------------------------------

std::string to_string(bond_type bondType)
{
	switch (bondType)
	{
		case bond_type::sing: return "sing";
		case bond_type::doub: return "doub";
		case bond_type::trip: return "trip";
		case bond_type::quad: return "quad";
		case bond_type::arom: return "arom";
		case bond_type::poly: return "poly";
		case bond_type::delo: return "delo";
		case bond_type::pi: return "pi";
	}
	throw std::invalid_argument("Invalid bondType");
}

bond_type from_string(const std::string &bondType)
{
	if (cif::iequals(bondType, "sing"))
		return bond_type::sing;
	if (cif::iequals(bondType, "doub"))
		return bond_type::doub;
	if (cif::iequals(bondType, "trip"))
		return bond_type::trip;
	if (cif::iequals(bondType, "quad"))
		return bond_type::quad;
	if (cif::iequals(bondType, "arom"))
		return bond_type::arom;
	if (cif::iequals(bondType, "poly"))
		return bond_type::poly;
	if (cif::iequals(bondType, "delo"))
		return bond_type::delo;
	if (cif::iequals(bondType, "pi"))
		return bond_type::pi;
	throw std::invalid_argument("Invalid bondType: " + bondType);
}

// --------------------------------------------------------------------
// compound helper classes

struct compound_atom_less
{
	bool operator()(const compound_atom &a, const compound_atom &b) const
	{
		int d = a.id.compare(b.id);
		if (d == 0)
			d = a.type_symbol - b.type_symbol;
		return d < 0;
	}
};

struct compound_bond_less
{
	bool operator()(const compound_bond &a, const compound_bond &b) const
	{
		int d = a.atom_id[0].compare(b.atom_id[0]);
		if (d == 0)
			d = a.atom_id[1].compare(b.atom_id[1]);
		if (d == 0)
			d = static_cast<int>(a.type) - static_cast<int>(b.type);
		return d < 0;
	}
};

// --------------------------------------------------------------------
// compound

compound::compound(cif::datablock &db)
{
	auto &chemComp = db["chem_comp"];

	if (chemComp.size() != 1)
		throw std::runtime_error("Invalid compound file, chem_comp should contain a single row");

	cif::tie(m_id, m_name, m_type, m_formula, m_formula_weight, m_formal_charge) =
		chemComp.front().get("id", "name", "type", "formula", "formula_weight", "pdbx_formal_charge");

	// The name should not contain newline characters since that triggers validation errors later on
	cif::replace_all(m_name, "\n", "");

	m_group = "non-polymer";

	auto &chemCompAtom = db["chem_comp_atom"];
	for (auto row : chemCompAtom)
	{
		compound_atom atom;
		std::string type_symbol;
		cif::tie(atom.id, type_symbol, atom.charge, atom.aromatic, atom.leaving_atom, atom.stereo_config, atom.x, atom.y, atom.z) =
			row.get("atom_id", "type_symbol", "charge", "pdbx_aromatic_flag", "pdbx_leaving_atom_flag", "pdbx_stereo_config",
				"model_Cartn_x", "model_Cartn_y", "model_Cartn_z");
		atom.type_symbol = atom_type_traits(type_symbol).type();
		m_atoms.push_back(std::move(atom));
	}

	auto &chemCompBond = db["chem_comp_bond"];
	for (auto row : chemCompBond)
	{
		compound_bond bond;
		std::string valueOrder;
		cif::tie(bond.atom_id[0], bond.atom_id[1], valueOrder, bond.aromatic, bond.stereo_config) = row.get("atom_id_1", "atom_id_2", "value_order", "pdbx_aromatic_flag", "pdbx_stereo_config");
		bond.type = from_string(valueOrder);
		m_bonds.push_back(std::move(bond));
	}
}

compound::compound(cif::datablock &db, const std::string &id, const std::string &name, const std::string &type, const std::string &group)
	: m_id(id)
	, m_name(name)
	, m_type(type)
	, m_group(group)
{
	auto &chemCompAtom = db["chem_comp_atom"];
	for (auto row : chemCompAtom)
	{
		compound_atom atom;
		std::string type_symbol;
		cif::tie(atom.id, type_symbol, atom.charge, atom.x, atom.y, atom.z) =
			row.get("atom_id", "type_symbol", "charge", "x", "y", "z");
		atom.type_symbol = atom_type_traits(type_symbol).type();

		m_formal_charge += atom.charge;
		m_formula_weight += atom_type_traits(atom.type_symbol).weight();

		m_atoms.push_back(std::move(atom));
	}

	auto &chemCompBond = db["chem_comp_bond"];
	for (auto row : chemCompBond)
	{
		compound_bond bond;
		std::string btype;
		cif::tie(bond.atom_id[0], bond.atom_id[1], btype, bond.aromatic) = row.get("atom_id_1", "atom_id_2", "type", "aromatic");

		using cif::iequals;

		if (iequals(btype, "single"))
			bond.type = bond_type::sing;
		else if (iequals(btype, "double"))
			bond.type = bond_type::doub;
		else if (iequals(btype, "triple"))
			bond.type = bond_type::trip;
		else if (iequals(btype, "deloc") or iequals(btype, "aromat") or iequals(btype, "aromatic"))
			bond.type = bond_type::delo;
		else
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Unimplemented chem_comp_bond.type " << btype << " in " << id << std::endl;
			bond.type = bond_type::sing;
		}
		m_bonds.push_back(std::move(bond));
	}
}

compound_atom compound::get_atom_by_atom_id(const std::string &atom_id) const
{
	compound_atom result = {};
	for (auto &a : m_atoms)
	{
		if (a.id == atom_id)
		{
			result = a;
			break;
		}
	}

	if (result.id != atom_id)
		throw std::out_of_range("No atom " + atom_id + " in compound " + m_id);

	return result;
}

bool compound::atoms_bonded(const std::string &atomId_1, const std::string &atomId_2) const
{
	auto i = find_if(m_bonds.begin(), m_bonds.end(),
		[&](const compound_bond &b)
		{
			return (b.atom_id[0] == atomId_1 and b.atom_id[1] == atomId_2) or (b.atom_id[0] == atomId_2 and b.atom_id[1] == atomId_1);
		});

	return i != m_bonds.end();
}

float compound::bond_length(const std::string &atomId_1, const std::string &atomId_2) const
{
	auto i = find_if(m_bonds.begin(), m_bonds.end(),
		[&](const compound_bond &b)
		{
			return (b.atom_id[0] == atomId_1 and b.atom_id[1] == atomId_2) or (b.atom_id[0] == atomId_2 and b.atom_id[1] == atomId_1);
		});

	float result = std::numeric_limits<float>::max();

	if (i != m_bonds.end())
	{
		auto a = get_atom_by_atom_id(atomId_1);
		auto b = get_atom_by_atom_id(atomId_2);

		result = distance(point{a.x, a.y, a.z}, point{b.x, b.y, b.z});
	}

	return result;
}


// --------------------------------------------------------------------
// known amino acids and bases

const std::map<std::string, char> compound_factory::kAAMap{
	{ "ALA", 'A' },
	{ "ARG", 'R' },
	{ "ASN", 'N' },
	{ "ASP", 'D' },
	{ "CYS", 'C' },
	{ "GLN", 'Q' },
	{ "GLU", 'E' },
	{ "GLY", 'G' },
	{ "HIS", 'H' },
	{ "ILE", 'I' },
	{ "LEU", 'L' },
	{ "LYS", 'K' },
	{ "MET", 'M' },
	{ "PHE", 'F' },
	{ "PRO", 'P' },
	{ "SER", 'S' },
	{ "THR", 'T' },
	{ "TRP", 'W' },
	{ "TYR", 'Y' },
	{ "VAL", 'V' },
	{ "GLX", 'Z' },
	{ "ASX", 'B' }
};

const std::map<std::string, char> compound_factory::kBaseMap{
	{ "A", 'A' },
	{ "C", 'C' },
	{ "G", 'G' },
	{ "T", 'T' },
	{ "U", 'U' },
	{ "DA", 'A' },
	{ "DC", 'C' },
	{ "DG", 'G' },
	{ "DT", 'T' }
};

// --------------------------------------------------------------------
// a factory class to generate compounds

class compound_factory_impl : public std::enable_shared_from_this<compound_factory_impl>
{
  public:
	compound_factory_impl(std::shared_ptr<compound_factory_impl> next);

	compound_factory_impl(const fs::path &file, std::shared_ptr<compound_factory_impl> next);

	virtual ~compound_factory_impl()
	{
		for (auto c : m_compounds)
			delete c;
	}

	compound *get(std::string id)
	{
		cif::to_upper(id);

		std::shared_lock lock(mMutex);

		compound *result = nullptr;

		// walk the list, see if any of us has the compound already
		for (auto impl = shared_from_this(); impl; impl = impl->m_next)
		{
			for (auto cmp : impl->m_compounds)
			{
				if (iequals(cmp->id(), id))
				{
					result = cmp;
					break;
				}
			}

			if (result)
				break;
		}

		if (result == nullptr and m_missing.count(id) == 0)
		{
			for (auto impl = shared_from_this(); impl; impl = impl->m_next)
			{
				result = impl->create(id);
				if (result != nullptr)
					break;
			}

			if (result == nullptr)
				m_missing.insert(id);
		}

		return result;
	}

	std::shared_ptr<compound_factory_impl> next() const
	{
		return m_next;
	}

	bool is_known_peptide(const std::string &resName)
	{
		return m_known_peptides.count(resName) or
		       (m_next and m_next->is_known_peptide(resName));
	}

	bool is_known_base(const std::string &resName)
	{
		return m_known_bases.count(resName) or
		       (m_next and m_next->is_known_base(resName));
	}

  protected:
	virtual compound *create(const std::string &id)
	{
		// For the base class we assume every compound is preloaded
		return nullptr;
	}

	std::shared_timed_mutex mMutex;

	std::vector<compound *> m_compounds;
	std::set<std::string> m_known_peptides;
	std::set<std::string> m_known_bases;
	std::set<std::string> m_missing;
	std::shared_ptr<compound_factory_impl> m_next;
};

// --------------------------------------------------------------------

compound_factory_impl::compound_factory_impl(std::shared_ptr<compound_factory_impl> next)
	: m_next(next)
{
	for (const auto &[key, value] : compound_factory::kAAMap)
		m_known_peptides.insert(key);

	for (const auto &[key, value] : compound_factory::kBaseMap)
		m_known_bases.insert(key);
}

compound_factory_impl::compound_factory_impl(const fs::path &file, std::shared_ptr<compound_factory_impl> next)
	: m_next(next)
{
	cif::file cifFile(file);

	if (cifFile.contains("comp_list")) // So this is a CCP4 restraints file, special handling
	{
		auto &compList = cifFile["comp_list"];
		auto &chemComp = compList["chem_comp"];

		for (const auto &[id, name, group] : chemComp.rows<std::string, std::string, std::string>("id", "name", "group"))
		{
			std::string type;

			// known groups are (counted from ccp4 monomer dictionary)

			//	D-pyranose
			//	DNA
			//	L-PEPTIDE LINKING
			//	L-SACCHARIDE
			//	L-peptide
			//	L-pyranose
			//	M-peptide
			//	NON-POLYMER
			//	P-peptide
			//	RNA
			//	furanose
			//	non-polymer
			//	non_polymer
			//	peptide
			//	pyranose
			//	saccharide

			if (cif::iequals(id, "gly"))
				type = "peptide linking";
			else if (cif::iequals(group, "l-peptide") or cif::iequals(group, "L-peptide linking") or cif::iequals(group, "peptide") or cif::iequals(group, "p-peptide"))
				type = "L-peptide linking";
			else if (cif::iequals(group, "DNA"))
				type = "DNA linking";
			else if (cif::iequals(group, "RNA"))
				type = "RNA linking";
			else
				type = "non-polymer";

			auto &db = cifFile["comp_" + id];

			m_compounds.push_back(new compound(db, id, name, type, group));
		}
	}
	else
	{
		// A CCD components file, validate it first
		try
		{
			cifFile.load_dictionary("mmcif_pdbx.dic");

			if (not cifFile.is_valid())
			{
				std::cerr << "The components file " << file << " is not valid" << std::endl;
				if (cif::VERBOSE < 1)
					std::cerr << "(use --verbose to see why)" << std::endl;
			}
		}
		catch (const std::exception &e)
		{
			std::cerr << "When trying to load the components file " << file << " there was an exception:" << std::endl
					  << e.what() << std::endl;
		}

		for (auto &db : cifFile)
			m_compounds.push_back(new compound(db));
	}
}

// --------------------------------------------------------------------
// Version for the default compounds, based on the cached components.cif file from CCD

class CCD_compound_factory_impl : public compound_factory_impl
{
  public:
	CCD_compound_factory_impl(std::shared_ptr<compound_factory_impl> next, const fs::path &file)
		: compound_factory_impl(next)
		, mCompoundsFile(file)
	{
	}

	CCD_compound_factory_impl(std::shared_ptr<compound_factory_impl> next)
		: compound_factory_impl(next)
	{
	}

	compound *create(const std::string &id) override;

	cif::parser::datablock_index mIndex;
	fs::path mCompoundsFile;
};

compound *CCD_compound_factory_impl::create(const std::string &id)
{
	compound *result = nullptr;

	std::unique_ptr<std::istream> ccd;

	if (mCompoundsFile.empty())
	{
		ccd = cif::load_resource("components.cif");
		if (not ccd)
		{
			std::cerr << "Could not locate the CCD components.cif file, please make sure the software is installed properly and/or use the update-libcifpp-data to fetch the data." << std::endl;
			return nullptr;
		}
	}
	else
		ccd.reset(new std::ifstream(mCompoundsFile));

	cif::file file;

	if (mIndex.empty())
	{
		if (cif::VERBOSE > 1)
		{
			std::cout << "Creating component index "
					  << "...";
			std::cout.flush();
		}

		cif::parser parser(*ccd, file);
		mIndex = parser.index_datablocks();

		if (cif::VERBOSE > 1)
			std::cout << " done" << std::endl;

		// reload the resource, perhaps this should be improved...
		if (mCompoundsFile.empty())
		{
			ccd = cif::load_resource("components.cif");
			if (not ccd)
				throw std::runtime_error("Could not locate the CCD components.cif file, please make sure the software is installed properly and/or use the update-libcifpp-data to fetch the data.");
		}
		else
			ccd.reset(new std::ifstream(mCompoundsFile));
	}

	if (cif::VERBOSE > 1)
	{
		std::cout << "Loading component " << id << "...";
		std::cout.flush();
	}

	cif::parser parser(*ccd, file);
	parser.parse_single_datablock(id, mIndex);

	if (cif::VERBOSE > 1)
		std::cout << " done" << std::endl;

	if (not file.empty())
	{
		auto &db = file.front();
		if (db.name() == id)
		{
			result = new compound(db);

			std::shared_lock lock(mMutex);
			m_compounds.push_back(result);
		}
	}

	if (result == nullptr and cif::VERBOSE > 0)
		std::cerr << "Could not locate compound " << id << " in the CCD components file" << std::endl;

	return result;
}

// --------------------------------------------------------------------
// Version for the default compounds, based on the data found in CCP4's monomers lib

class CCP4_compound_factory_impl : public compound_factory_impl
{
  public:
	CCP4_compound_factory_impl(const fs::path &clibd_mon, std::shared_ptr<compound_factory_impl> next = nullptr);

	compound *create(const std::string &id) override;

  private:
	cif::file m_file;
	fs::path m_CLIBD_MON;
};

CCP4_compound_factory_impl::CCP4_compound_factory_impl(const fs::path &clibd_mon, std::shared_ptr<compound_factory_impl> next)
	: compound_factory_impl(next)
	, m_file((clibd_mon / "list" / "mon_lib_list.cif").string())
	, m_CLIBD_MON(clibd_mon)
{
	const std::regex peptideRx("(?:[lmp]-)?peptide", std::regex::icase);

	auto &chemComps = m_file["comp_list"]["chem_comp"];

	for (const auto &[group, comp_id] : chemComps.rows<std::string, std::string>("group", "id"))
	{
		if (std::regex_match(group, peptideRx))
			m_known_peptides.insert(comp_id);
		else if (cif::iequals(group, "DNA") or cif::iequals(group, "RNA"))
			m_known_bases.insert(comp_id);
	}
}

compound *CCP4_compound_factory_impl::create(const std::string &id)
{
	compound *result = nullptr;

	auto &cat = m_file["comp_list"]["chem_comp"];

	auto rs = cat.find(cif::key("id") == id);

	if (rs.size() == 1)
	{
		auto row = rs.front();

		std::string name, group;
		uint32_t numberAtomsAll, numberAtomsNh;
		cif::tie(name, group, numberAtomsAll, numberAtomsNh) =
			row.get("name", "group", "number_atoms_all", "number_atoms_nh");

		fs::path resFile = m_CLIBD_MON / cif::to_lower_copy(id.substr(0, 1)) / (id + ".cif");

		if (not fs::exists(resFile) and (id == "COM" or id == "CON" or "PRN")) // seriously...
			resFile = m_CLIBD_MON / cif::to_lower_copy(id.substr(0, 1)) / (id + '_' + id + ".cif");

		if (fs::exists(resFile))
		{
			cif::file cf(resFile.string());

			// locate the datablock
			auto &db = cf["comp_" + id];

			std::string type;

			// known groups are (counted from ccp4 monomer dictionary)

			//	D-pyranose
			//	DNA
			//	L-PEPTIDE LINKING
			//	L-SACCHARIDE
			//	L-peptide
			//	L-pyranose
			//	M-peptide
			//	NON-POLYMER
			//	P-peptide
			//	RNA
			//	furanose
			//	non-polymer
			//	non_polymer
			//	peptide
			//	pyranose
			//	saccharide

			if (cif::iequals(id, "gly"))
				type = "peptide linking";
			else if (cif::iequals(group, "l-peptide") or cif::iequals(group, "L-peptide linking") or cif::iequals(group, "peptide") or cif::iequals(group, "p-peptide"))
				type = "L-peptide linking";
			else if (cif::iequals(group, "DNA"))
				type = "DNA linking";
			else if (cif::iequals(group, "RNA"))
				type = "RNA linking";
			else
				type = "non-polymer";

			m_compounds.push_back(new compound(db, id, name, type, group));
			result = m_compounds.back();
		}
	}

	return result;
}

// --------------------------------------------------------------------

std::unique_ptr<compound_factory> compound_factory::s_instance;
thread_local std::unique_ptr<compound_factory> compound_factory::tl_instance;
bool compound_factory::s_use_thread_local_instance;

void compound_factory::init(bool useThreadLocalInstanceOnly)
{
	s_use_thread_local_instance = useThreadLocalInstanceOnly;
}

compound_factory::compound_factory()
	: m_impl(nullptr)
{
	auto ccd = cif::load_resource("components.cif");
	if (ccd)
		m_impl = std::make_shared<CCD_compound_factory_impl>(m_impl);
	else if (cif::VERBOSE > 0)
		std::cerr << "CCD components.cif file was not found" << std::endl;

	const char *clibd_mon = getenv("CLIBD_MON");
	if (clibd_mon != nullptr and fs::is_directory(clibd_mon))
		m_impl = std::make_shared<CCP4_compound_factory_impl>(clibd_mon, m_impl);
	else if (cif::VERBOSE > 0)
		std::cerr << "CCP4 monomers library not found, CLIBD_MON is not defined" << std::endl;
}

compound_factory::~compound_factory()
{
}

compound_factory &compound_factory::instance()
{
	if (s_use_thread_local_instance)
	{
		if (not tl_instance)
			tl_instance.reset(new compound_factory());
		return *tl_instance;
	}
	else
	{
		if (not s_instance)
			s_instance.reset(new compound_factory());
		return *s_instance;
	}
}

void compound_factory::clear()
{
	if (s_use_thread_local_instance)
		tl_instance.reset(nullptr);
	else
		s_instance.reset();
}

void compound_factory::set_default_dictionary(const fs::path &inDictFile)
{
	if (not fs::exists(inDictFile))
		throw std::runtime_error("file not found: " + inDictFile.string());

	try
	{
		m_impl.reset(new CCD_compound_factory_impl(m_impl, inDictFile));
	}
	catch (const std::exception &)
	{
		std::throw_with_nested(std::runtime_error("Error loading dictionary " + inDictFile.string()));
	}
}

void compound_factory::push_dictionary(const fs::path &inDictFile)
{
	if (not fs::exists(inDictFile))
		throw std::runtime_error("file not found: " + inDictFile.string());

	try
	{
		m_impl.reset(new compound_factory_impl(inDictFile, m_impl));
	}
	catch (const std::exception &)
	{
		std::throw_with_nested(std::runtime_error("Error loading dictionary " + inDictFile.string()));
	}
}

void compound_factory::pop_dictionary()
{
	if (m_impl)
		m_impl = m_impl->next();
}

const compound *compound_factory::create(std::string id)
{
	return m_impl ? m_impl->get(id) : nullptr;
}

bool compound_factory::is_known_peptide(const std::string &resName) const
{
	return m_impl ? m_impl->is_known_peptide(resName) : kAAMap.count(resName) > 0;
}

bool compound_factory::is_known_base(const std::string &resName) const
{
	return m_impl ? m_impl->is_known_base(resName) : kBaseMap.count(resName) > 0;
}

} // namespace cif
