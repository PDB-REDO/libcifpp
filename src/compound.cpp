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

bond_type parse_bond_type_from_string(const std::string &bondType)
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

std::string to_string(stereo_config_type stereoConfig)
{
	switch (stereoConfig)
	{
		case stereo_config_type::N: return "N";
		case stereo_config_type::R: return "R";
		case stereo_config_type::S: return "S";
	}
	throw std::invalid_argument("Invalid stereoConfig");
}

stereo_config_type parse_stereo_config_from_string(const std::string &stereoConfig)
{
	if (cif::iequals(stereoConfig, "N"))
		return stereo_config_type::N;
	if (cif::iequals(stereoConfig, "R"))
		return stereo_config_type::R;
	if (cif::iequals(stereoConfig, "S"))
		return stereo_config_type::S;
	throw std::invalid_argument("Invalid stereoConfig: " + stereoConfig);
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

	std::string one_letter_code;

	cif::tie(m_id, m_name, m_type, m_formula, m_formula_weight, m_formal_charge, one_letter_code, m_parent_id) =
		chemComp.front().get("id", "name", "type", "formula", "formula_weight", "pdbx_formal_charge", "one_letter_code", "mon_nstd_parent_comp_id");
	
	if (one_letter_code.length() == 1)
		m_one_letter_code = one_letter_code.front();

	// The name should not contain newline characters since that triggers validation errors later on
	cif::replace_all(m_name, "\n", "");

	auto &chemCompAtom = db["chem_comp_atom"];
	for (auto row : chemCompAtom)
	{
		compound_atom atom;
		std::string type_symbol, stereo_config;
		cif::tie(atom.id, type_symbol, atom.charge, atom.aromatic, atom.leaving_atom, stereo_config, atom.x, atom.y, atom.z) =
			row.get("atom_id", "type_symbol", "charge", "pdbx_aromatic_flag", "pdbx_leaving_atom_flag", "pdbx_stereo_config",
				"model_Cartn_x", "model_Cartn_y", "model_Cartn_z");
		atom.type_symbol = atom_type_traits(type_symbol).type();
		if (stereo_config.empty())
			atom.stereo_config = stereo_config_type::N;
		else
		atom.stereo_config = parse_stereo_config_from_string(stereo_config);
		m_atoms.push_back(std::move(atom));
	}

	auto &chemCompBond = db["chem_comp_bond"];
	for (auto row : chemCompBond)
	{
		compound_bond bond;
		std::string valueOrder;
		cif::tie(bond.atom_id[0], bond.atom_id[1], valueOrder, bond.aromatic, bond.stereo_config) = row.get("atom_id_1", "atom_id_2", "value_order", "pdbx_aromatic_flag", "pdbx_stereo_config");
		if (valueOrder.empty())
			bond.type = bond_type::sing;
		else
		bond.type = parse_bond_type_from_string(valueOrder);
		m_bonds.push_back(std::move(bond));
	}
}

compound::compound(cif::datablock &db, int)
{
	auto &chemComp = db["chem_comp"];

	if (chemComp.size() != 1)
		throw std::runtime_error("Invalid compound file, chem_comp should contain a single row");

	cif::tie(m_id, m_name) =
		chemComp.front().get("id", "name");
	
	cif::trim(m_name);

	m_type = "NON-POLYMER";

	auto &chemCompAtom = db["chem_comp_atom"];
	for (auto row : chemCompAtom)
	{
		compound_atom atom;
		std::string type_symbol;
		cif::tie(atom.id, type_symbol, atom.charge, atom.x, atom.y, atom.z) =
			row.get("atom_id", "type_symbol", "charge", "x", "y", "z");
		atom.type_symbol = atom_type_traits(type_symbol).type();

		m_formal_charge += atom.charge;

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
				std::cerr << "Unimplemented chem_comp_bond.type " << btype << " in " << db.name() << '\n';
			bond.type = bond_type::sing;
		}
		m_bonds.push_back(std::move(bond));
	}

	// reconstruct a formula and weight

	m_formula_weight = 0;

	std::map<atom_type, int> f;
	for (auto &atom : m_atoms)
		f[atom.type_symbol] += 1;
	
	if (f.count(atom_type::C))
	{
		atom_type_traits att(atom_type::C);
		m_formula += att.symbol() + std::to_string(f[atom_type::C]) + ' ';
		m_formula_weight += att.weight() * f[atom_type::C];
	}

	for (const auto &[type, count] : f)
	{
		if (type == atom_type::C)
			continue;

		atom_type_traits att(type);
		m_formula += att.symbol() + std::to_string(count) + ' ';
		m_formula_weight += att.weight() * count;
	}
	
	if (not m_formula.empty())
		m_formula.pop_back();
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

		result = distance(point{ a.x, a.y, a.z }, point{ b.x, b.y, b.z });
	}

	return result;
}

// --------------------------------------------------------------------

bool compound::is_peptide() const
{
	return iequals(m_type, "l-peptide linking")	or iequals(m_type, "peptide linking");
}

bool compound::is_base() const
{
	return iequals(m_type, "dna linking")	or iequals(m_type, "rna linking");
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
	compound_factory_impl();
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

		// walk the list, see if any of the implementations has the compound already
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

	std::shared_ptr<compound_factory_impl> next()
	{
		return m_next;
	}

	void describe(std::ostream &os)
	{
		if (m_file.empty())
			os << "CCD components.cif resource\n";
		else
			os << "CCD components file: " << std::quoted(m_file.string()) << '\n';

		if (m_next)
			m_next->describe(os);
	}

  protected:
	compound_factory_impl(std::shared_ptr<compound_factory_impl> next);

	virtual compound *create(const std::string &id);

	std::shared_timed_mutex mMutex;

	fs::path m_file;
	cif::parser::datablock_index m_index;

	std::vector<compound *> m_compounds;
	std::set<std::string> m_missing;
	std::shared_ptr<compound_factory_impl> m_next;
};

compound_factory_impl::compound_factory_impl()
{
}

compound_factory_impl::compound_factory_impl(std::shared_ptr<compound_factory_impl> next)
	: m_next(next)
{
}

compound_factory_impl::compound_factory_impl(const fs::path &file, std::shared_ptr<compound_factory_impl> next)
	: compound_factory_impl(next)
{
	m_file = file;
}

compound *compound_factory_impl::create(const std::string &id)
{
	compound *result = nullptr;

	std::unique_ptr<std::istream> ccd;

	if (m_file.empty())
	{
		ccd = cif::load_resource("components.cif");
		if (not ccd)
		{
			std::cerr << "Could not locate the CCD components.cif file, please make sure the software is installed properly and/or use the update-libcifpp-data to fetch the data.\n";
			return nullptr;
		}
	}
	else
		ccd.reset(new std::ifstream(m_file));

	cif::file file;

	if (m_index.empty())
	{
		if (cif::VERBOSE > 1)
		{
			std::cout << "Creating component index "
					  << "...";
			std::cout.flush();
		}

		cif::parser parser(*ccd, file);
		m_index = parser.index_datablocks();

		if (cif::VERBOSE > 1)
			std::cout << " done" << std::endl;

		// reload the resource, perhaps this should be improved...
		if (m_file.empty())
		{
			ccd = cif::load_resource("components.cif");
			if (not ccd)
				throw std::runtime_error("Could not locate the CCD components.cif file, please make sure the software is installed properly and/or use the update-libcifpp-data to fetch the data.");
		}
		else
			ccd.reset(new std::ifstream(m_file));
	}

	if (cif::VERBOSE > 1)
	{
		std::cout << "Loading component " << id << "...";
		std::cout.flush();
	}

	cif::parser parser(*ccd, file);
	parser.parse_single_datablock(id, m_index);

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

	return result;
}

// --------------------------------------------------------------------

class local_compound_factory_impl : public compound_factory_impl
{
  public:
	local_compound_factory_impl(const cif::file &file, std::shared_ptr<compound_factory_impl> next)
		: compound_factory_impl(next)
		, m_local_file(file)
	{
		const std::regex peptideRx("(?:[lmp]-)?peptide", std::regex::icase);

		for (const auto &[id, name, threeLetterCode, group] :
			file["comp_list"]["chem_comp"].rows<std::string, std::string, std::string, std::string>("id", "name", "three_letter_code", "group"))
		{
			auto &rdb = m_local_file["comp_" + id];
			if (rdb.empty())
			{
				std::cerr << "Missing data in restraint file for id " + id + '\n';
				continue;
			}

			construct_compound(rdb, id, name, threeLetterCode, group);
		}
	}

	compound *create(const std::string &id) override;

  private:

	compound *construct_compound(const datablock &db, const std::string &id, const std::string &name, const std::string &three_letter_code, const std::string &group);

	cif::file m_local_file;
};

compound *local_compound_factory_impl::create(const std::string &id)
{
	compound *result = nullptr;

	for (auto &db : m_local_file)
	{
		if (db.name() == "comp_" + id)
		{
			auto chem_comp = db.get("chem_comp");
			if (not chem_comp)
				break;

			try
			{
				const auto &[id, name, threeLetterCode, group] =
					chem_comp->front().get<std::string, std::string, std::string, std::string>("id", "name", "three_letter_code", "group");

				result = construct_compound(db, id, name, threeLetterCode, group);
			}
			catch (const std::exception &ex)
			{
				std::throw_with_nested(std::runtime_error("Error loading compound " + id));
			}

			break;
		}
	}

	return result;
}

compound *local_compound_factory_impl::construct_compound(const datablock &rdb, const std::string &id,
	const std::string &name, const std::string &three_letter_code, const std::string &group)
{
	cif::datablock db(id);

	float formula_weight = 0;
	int formal_charge = 0;
	std::map<std::string,size_t> formula_data;

	for (size_t ord = 1; const auto &[atom_id, type_symbol, type, charge, x, y, z] :
		rdb["chem_comp_atom"].rows<std::string, std::string, std::string, int, float, float, float>(
			"atom_id", "type_symbol", "type", "charge", "x", "y", "z"))
	{
		auto atom = cif::atom_type_traits(type_symbol);
		formula_weight += atom.weight();

		formula_data[type_symbol] += 1;

		db["chem_comp_atom"].emplace({
			{ "comp_id", id },
			{ "atom_id", atom_id },
			{ "type_symbol", type_symbol },
			{ "charge", charge },
			{ "model_Cartn_x",  x, 3 },
			{ "model_Cartn_y",  y, 3 },
			{ "model_Cartn_z",  z, 3 },
			{ "pdbx_ordinal", ord++ }
		});

		formal_charge += charge;
	}

	for (size_t ord = 1; const auto &[atom_id_1, atom_id_2, type, aromatic] :
		rdb["chem_comp_bond"].rows<std::string, std::string, std::string, bool>("atom_id_1", "atom_id_2", "type", "aromatic"))
	{
		std::string value_order("SING");

		if (cif::iequals(type, "single") or cif::iequals(type, "sing"))
			value_order = "SING";
		else if (cif::iequals(type, "double") or cif::iequals(type, "doub"))
			value_order = "DOUB";
		else if (cif::iequals(type, "triple") or cif::iequals(type, "trip"))
			value_order = "TRIP";

		db["chem_comp_bond"].emplace({
			{ "comp_id", id },
			{ "atom_id_1", atom_id_1 },
			{ "atom_id_2", atom_id_2 },
			{ "value_order", value_order },
			{ "pdbx_aromatic_flag", aromatic },
			// TODO: fetch stereo_config info from chem_comp_chir
			{ "pdbx_ordinal", ord++ }
		});
	}

	db.emplace_back(rdb["pdbx_chem_comp_descriptor"]);

	std::string formula;
	for (bool first = true; const auto &[symbol, count]: formula_data)
	{
		if (std::exchange(first, false))
			formula += ' ';
		formula += symbol;
		if (count > 1)
			formula += std::to_string(count);
	}

	std::string type;
	if (cif::iequals(group, "peptide") or cif::iequals(group, "l-peptide") or cif::iequals(group, "l-peptide linking"))
		type = "L-PEPTIDE LINKING";
	else if (cif::iequals(group, "dna"))
		type = "DNA LINKING";
	else if (cif::iequals(group, "rna"))
		type = "RNA LINKING";
	else
		type = "NON-POLYMER";

	db["chem_comp"].emplace({
		{ "id", id },
		{ "name", name },
		{ "type", type },
		{ "formula", formula },
		{ "pdbx_formal_charge", formal_charge },
		{ "formula_weight", formula_weight },
		{ "three_letter_code", three_letter_code }
	});

	std::shared_lock lock(mMutex);

	auto result = new compound(db);
	m_compounds.push_back(result);
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
		m_impl = std::make_shared<compound_factory_impl>();
	else if (cif::VERBOSE > 0)
		std::cerr << "CCD components.cif resource was not found\n";
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
		m_impl.reset(new compound_factory_impl(inDictFile, m_impl));
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

void compound_factory::push_dictionary(const cif::file &inDictFile)
{
	try
	{
		m_impl.reset(new local_compound_factory_impl(inDictFile, m_impl));
	}
	catch (const std::exception &)
	{
		std::throw_with_nested(std::runtime_error("Error loading dictionary from local mmCIF file"));
	}
}

void compound_factory::pop_dictionary()
{
	if (m_impl)
		m_impl = m_impl->next();
}

const compound *compound_factory::create(std::string_view id)
{
	auto result = m_impl ? m_impl->get(std::string{ id }) : nullptr;
	if (not result)
		report_missing_compound(id);
	return result;
}

bool compound_factory::is_known_peptide(const std::string &res_name) const
{
	return kAAMap.count(res_name) > 0;
}

bool compound_factory::is_known_base(const std::string &res_name) const
{
	return kBaseMap.count(res_name) > 0;
}

/// Return whether @a res_name is a peptide
bool compound_factory::is_peptide(std::string_view res_name) const
{
	bool result = is_std_peptide(res_name);
	if (not result and m_impl)
	{
		auto compound = const_cast<compound_factory&>(*this).create(res_name);
		result = compound != nullptr and compound->is_peptide();
	}
	return result;
}

/// Return whether @a res_name is a base
bool compound_factory::is_base(std::string_view res_name) const
{
	bool result = is_std_base(res_name);
	if (not result and m_impl)
	{
		auto compound = const_cast<compound_factory&>(*this).create(res_name);
		result = compound != nullptr and compound->is_base();
	}
	return result;
}

/// Return whether @a res_name is one of the standard peptides
bool compound_factory::is_std_peptide(std::string_view res_name) const
{
	return kAAMap.count(std::string{ res_name }) > 0;
}

/// Return whether @a res_name is one of the standard bases
bool compound_factory::is_std_base(std::string_view res_name) const
{
	return kBaseMap.count(std::string{ res_name }) > 0;
}

/// Return whether @a res_name is a monomer (either base or peptide)
bool compound_factory::is_monomer(std::string_view res_name) const
{
	return is_peptide(res_name) or is_base(res_name);
}

void compound_factory::report_missing_compound(std::string_view compound_id)
{
	static bool s_reported = false;
	if (std::exchange(s_reported, true) == false)
	{
		using namespace cif::colour;

		std::clog << "\n"
				  << cif::coloured("Configuration error:", white, red) << "\n\n"
				  << "The attempt to retrieve compound information for " << std::quoted(compound_id) << " failed.\n\n"
				  << "This information is searched for in a CCD file called components.cif or\n"
				  << "components.cif.gz which should be located in one of the following directories:\n\n";

		cif::list_data_directories(std::clog);

		std::clog << "\n(Note that you can add a directory to the search paths by setting the \n"
				  << "LIBCIFPP_DATA_DIR environmental variable)\n\n";

#if defined(CACHE_DIR)
		std::clog << "On Linux an optional cron script might have been installed that automatically updates\n"
				  << "components.cif and mmCIF dictionary files. This script only works when the file\n"
				  << "libcifpp.conf contains an uncommented line with the text:\n\n"
				  << "update=true\n\n"
				  << "If you do not have a working cron script, you can manually update the files\n"
				  << "in /var/cache/libcifpp using the following commands:\n\n"
				  << "curl -o " << CACHE_DIR << "/components.cif https://files.wwpdb.org/pub/pdb/data/monomers/components.cif\n"
				  << "curl -o " << CACHE_DIR << "/mmcif_pdbx.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic\n"
				  << "curl -o " << CACHE_DIR << "/mmcif_ma.dic https://github.com/ihmwg/ModelCIF/raw/master/dist/mmcif_ma.dic\n\n";
#endif

		if (m_impl)
		{
			std::clog << "The current order of compound factory objects is:\n\n";
			m_impl->describe(std::clog);
		}
		else
			std::clog << "No compound factory objects are created since none of the data sources is found.\n";

		cif::list_file_resources(std::clog);

		std::clog.flush();
	}
}

} // namespace cif
