/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2024 NKI/AVL, Netherlands Cancer Institute
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

// --------------------------------------------------------------------

namespace cif::pdb
{


void checkAtomRecords(datablock &db)
{
	using namespace literals;

	auto &cf = compound_factory::instance();

	auto &atom_site = db["atom_site"];
	auto &atom_type = db["atom_type"];
	auto &chem_comp = db["chem_comp"];

	for (auto row : atom_site)
	{
		const auto &[symbol, label_asym_id, auth_asym_id, label_comp_id, auth_comp_id, label_seq_id, auth_seq_id, label_atom_id, auth_atom_id] =
			row.get<std::string, std::optional<std::string>, std::optional<std::string>, std::optional<std::string>, std::optional<std::string>,
				std::optional<int>, std::optional<std::string>, std::optional<std::string>, std::optional<std::string>>(
				"type_symbol", "label_asym_id", "auth_asym_id", "label_comp_id", "auth_comp_id", "label_seq_id", "auth_seq_id", "label_atom_id", "auth_atom_id");

		if (symbol.empty())
			throw std::runtime_error("Missing type symbol in atom_site record");

		if (atom_type.count("symbol"_key == symbol) == 0)
			atom_type.emplace({ { "symbol", symbol } });

		if (not(label_asym_id.has_value() or auth_asym_id.has_value()))
			throw std::runtime_error("atom_site records does not have a label_asym_id nor an auth_asym_id, cannot continue");

		if (not(label_comp_id.has_value() or auth_comp_id.has_value()))
			throw std::runtime_error("atom_site records does not have a label_comp_id nor an auth_comp_id, cannot continue");

		if (not(label_atom_id.has_value() or auth_atom_id.has_value()))
			throw std::runtime_error("atom_site records does not have a label_atom_id nor an auth_atom_id, cannot continue");

		std::string asym_id = label_asym_id.value_or(*auth_asym_id);
		std::string comp_id = label_comp_id.value_or(*auth_comp_id);

		bool is_peptide = cf.is_known_peptide(comp_id);
		auto compound = cf.create(comp_id);

		if (not compound)
			throw std::runtime_error("Missing compound information for " + comp_id);

		std::string mon_nstd_flag(".");
		if (is_peptide)
		{
			if (compound_factory::kAAMap.find(comp_id) != compound_factory::kAAMap.end())
				mon_nstd_flag = "y";
			else
				mon_nstd_flag = "n";
		}

		auto chem_comp_entry = chem_comp.find_first("id"_key == comp_id);

		if (not chem_comp_entry)
		{
			chem_comp.emplace({ //
				{ "id", comp_id },
				{ "type", compound->type() },
				{ "mon_nstd_flag", mon_nstd_flag },
				{ "name", compound->name() },
				{ "formula", compound->formula() },
				{ "formula_weight", compound->formula_weight() } });
		}
		else
		{
			std::vector<item> items;

			if (not chem_comp_entry["type"])
				items.emplace_back(item{ "type", compound->type() });
			if (not chem_comp_entry["mon_nstd_flag"])
				items.emplace_back(item{ "mon_nstd_flag", mon_nstd_flag });
			if (not chem_comp_entry["name"])
				items.emplace_back(item{ "name", compound->name() });
			if (not chem_comp_entry["formula"])
				items.emplace_back(item{ "formula", compound->formula() });
			if (not chem_comp_entry["formula_weight"])
				items.emplace_back(item{ "formula_weight", compound->formula_weight() });

			if (not items.empty())
				chem_comp_entry.assign(std::move(items));
		}

		if (is_peptide and not(label_seq_id.has_value() or auth_seq_id.has_value()))
			throw std::runtime_error("atom_site record has peptide comp_id but no sequence number, cannot continue");

		std::string seq_id;
		if (label_seq_id.has_value())
			seq_id = std::to_string(*label_seq_id);
		else if (auth_seq_id.has_value())
			seq_id = *auth_seq_id;

		row.assign({ //
			{ "auth_asym_id", auth_asym_id.value_or(*label_asym_id) },
			{ "auth_seq_id", auth_seq_id.value_or(std::to_string(*label_seq_id)) },
			{ "auth_comp_id", auth_comp_id.value_or(*label_comp_id) },
			{ "auth_atom_id", auth_atom_id.value_or(*label_atom_id) } });
	}
}

void createStructAsym(datablock &db)
{
	auto &atom_site = db["atom_site"];
	auto &struct_asym = db["struct_asym"];

	for (auto label_asym_id : atom_site.rows<std::string>("label_asym_id"))
	{
		if (label_asym_id.empty())
			throw std::runtime_error("File contains atom_site records without a label_asym_id");
		if (struct_asym.count(key("id") == label_asym_id) == 0)
		{
			struct_asym.emplace({ //
				{ "id", label_asym_id } });
		}
	}
}

void createEntity(datablock &db)
{
	using namespace literals;

	auto &cf = compound_factory::instance();

	auto &atom_site = db["atom_site"];
	atom_site.add_column("label_entity_id");

	auto &struct_asym = db["struct_asym"];
	struct_asym.add_column("entity_id");

	std::map<std::string,std::vector<std::tuple<std::string,int>>> asyms;

	for (auto asym_id : db["struct_asym"].rows<std::string>("id"))
	{
		int last_seq_id = -1;

		for (const auto &[comp_id, seq_id] : atom_site.find<std::string,int>("label_asym_id"_key == asym_id, "label_comp_id", "label_seq_id"))
		{
			if (seq_id == last_seq_id)
				continue;
			
			last_seq_id = seq_id;

			asyms[asym_id].emplace_back(comp_id, last_seq_id);
		}
	}

	auto less = [](const std::vector<std::tuple<std::string,int>> &a, const std::vector<std::tuple<std::string,int>> &b)
	{
		int d = static_cast<int>(a.size()) - static_cast<int>(b.size());
		return d == 0 ? a > b : d > 0;
	};

	std::set<std::vector<std::tuple<std::string,int>>,decltype(less)> entities(less);

	for (const auto &[asym_id, content] : asyms)
		entities.emplace(content);
	
	auto water_weight = cf.create("HOH")->formula_weight();

	int poly_count = 0;

	auto &entity = db["entity"];
	for (auto &content : entities)
	{
		auto entity_id = entity.get_unique_id("");

		std::string type, desc;
		float weight = 0;
		int count = 0;

		auto first_comp_id = std::get<0>(content.front());

		if (first_comp_id == "HOH")
		{
			type = "water";
			desc = "water";
			weight = water_weight;
		}
		else if (content.size() == 1)
		{
			auto c = cf.create(first_comp_id);

			type = "non-polymer";
			desc = c->name();
			weight = c->formula_weight();
		}
		else
		{
			type = "polymer";
			desc = "polymer-" + std::to_string(++poly_count);

			weight = water_weight;
			for (const auto &[comp_id, seq_id] : content)
				weight += cf.create(comp_id)->formula_weight() - water_weight;
		}

		for (const auto &[asym_id, ac] : asyms)
		{
			if (ac != content)
				continue;
			
			atom_site.update_value("label_asym_id"_key == asym_id, "label_entity_id", entity_id);
			struct_asym.update_value("id"_key == asym_id, "entity_id", entity_id);

			if (type != "water")
				++count;
			else
				count = atom_site.count("label_asym_id"_key == asym_id and "label_atom_id"_key == "O");
		}

		entity.emplace({ // 
			{ "id", entity_id },
			{ "type", type },
			{ "pdbx_description", desc },
			{ "formula_weight", weight },
			{ "pdbx_number_of_molecules", count }
		});
	}
}

void createEntityPoly(datablock &db)
{
	using namespace literals;

	auto &cf = compound_factory::instance();

	auto &atom_site = db["atom_site"];
	auto &entity_poly = db["entity_poly"];

	for (auto entity_id : db["entity"].find<std::string>("type"_key == "polymer", "id"))
	{
		std::string type;
		int last_seq_id = -1;
		std::string seq, seq_can;
		bool non_std_monomer = false;
		bool non_std_linkage = false;
		std::string pdb_strand_id;

		for (const auto &[comp_id, seq_id, auth_asym_id] : atom_site.find<std::string,int,std::string>("label_entity_id"_key == entity_id, "label_comp_id", "label_seq_id", "auth_asym_id"))
		{
			if (seq_id == last_seq_id)
				continue;
			
			last_seq_id = seq_id;

			auto c = cf.create(comp_id);

			std::string letter;
			char letter_can;

			// TODO: Perhaps we should improve this... 
			if (type != "other")
			{
				std::string c_type;
				if (cf.is_known_base(comp_id))
				{
					c_type = "polydeoxyribonucleotide";
					letter = letter_can = compound_factory::kBaseMap.at(comp_id);
				}
				else if (cf.is_known_peptide(comp_id))
				{
					c_type = "polypeptide(L)";
					letter = letter_can = compound_factory::kAAMap.at(comp_id);
				}
				else if (iequals(c->type(), "D-PEPTIDE LINKING"))
				{
					c_type = "polypeptide(D)";

					letter_can = c->one_letter_code();
					if (letter_can == 0)
						letter_can = 'X';
					
					letter = '(' + comp_id + ')';

					non_std_linkage = true;
					non_std_monomer = true;
				}
				else if (iequals(c->type(), "L-PEPTIDE LINKING") or iequals(c->type(), "PEPTIDE LINKING"))
				{
					c_type = "polypeptide(L)";

					letter_can = c->one_letter_code();
					if (letter_can == 0)
						letter_can = 'X';

					letter = '(' + comp_id + ')';

					non_std_monomer = true;
				}

				if (type.empty())
					type = c_type;
				else if (type != c_type)
					type = "other";
			}

			seq += letter;
			seq_can += letter_can;

			pdb_strand_id = auth_asym_id;
		}

		for (auto i = seq.begin() + 80; i < seq.end(); i += 80)
			i = seq.insert(i, '\n') + 1;
		
		for (auto i = seq_can.begin() + 76; i < seq_can.end(); i += 76)
		{
			auto j = i;
			while (j < i + 4 and j < seq_can.end())
			{
				if (*j == '(')
					break;
				++j;
			}

			if (j < seq_can.end())
				i = seq_can.insert(j, '\n') + 1;
			else
				i = j;
		}

		entity_poly.emplace({ // 
			{ "entity_id", entity_id },
			{ "type", type },
			{ "nstd_linkage", non_std_linkage },
			{ "nstd_monomer", non_std_monomer },
			{ "pdbx_seq_one_letter_code", seq },
			{ "pdbx_seq_one_letter_code_can", seq_can },
			{ "pdbx_strand_id", pdb_strand_id }
		});
	}
}

void createEntityPolySeq(datablock &db)
{
	if (db.get("entity_poly") == nullptr)
		createEntityPoly(db);

	using namespace literals;

	auto &atom_site = db["atom_site"];
	auto &entity_poly = db["entity_poly"];
	auto &entity_poly_seq = db["entity_poly_seq"];
	auto &struct_asym = db["struct_asym"];

	for (auto entity_id : entity_poly.rows<std::string>("entity_id"))
	{
		int last_seq_id = -1;
		std::string last_comp_id;
		std::string asym_id = struct_asym.find_first<std::string>("entity_id"_key == entity_id, "id");

		for (const auto &[comp_id, seq_id] : atom_site.find<std::string,int>("label_entity_id"_key == entity_id and "label_asym_id"_key == asym_id, "label_comp_id", "label_seq_id"))
		{
			bool hetero = false;

			if (seq_id == last_seq_id)
			{
				if (last_comp_id != comp_id)
					hetero = true;
				else
					continue;
			}

			if (hetero)
			{
				entity_poly_seq.back().assign({
					{ "hetero", true }
				});
			}

			entity_poly_seq.emplace({ // 
				{ "entity_id", entity_id },
				{ "num", seq_id },
				{ "mon_id", comp_id },
				{ "hetero", hetero }
			});
			
			last_seq_id = seq_id;
			last_comp_id = comp_id;
		}

		// you cannot assume this is correct...
		entity_poly_seq.sort([](row_handle a, row_handle b)
		{
			return a.get<int>("num") < b.get<int>("num");
		});
	}
}

void createPdbxPolySeqScheme(datablock &db)
{
	if (db.get("entity_poly_seq") == nullptr)
		createEntityPolySeq(db);

	using namespace literals;

	auto &atom_site = db["atom_site"];
	auto &entity_poly = db["entity_poly"];
	auto &entity_poly_seq = db["entity_poly_seq"];
	auto &struct_asym = db["struct_asym"];
	auto &pdbx_poly_seq_scheme = db["pdbx_poly_seq_scheme"];

	for (const auto &[entity_id, pdb_strand_id] : entity_poly.rows<std::string, std::string>("entity_id", "pdbx_strand_id"))
	{
		for (auto asym_id : struct_asym.find<std::string>("entity_id"_key == entity_id, "id"))
		{
			for (const auto &[comp_id, num, hetero] : entity_poly_seq.find<std::string,int,bool>("entity_id"_key == entity_id, "mon_id", "num", "hetero"))
			{
				const auto &[auth_seq_num, auth_mon_id, ins_code] =
					atom_site.find_first<std::string,std::string,std::optional<std::string>>(
						"label_asym_id"_key == asym_id and "label_seq_id"_key == num,
						"auth_seq_id", "auth_comp_id", "pdbx_PDB_ins_code"
					);
				
				pdbx_poly_seq_scheme.emplace({ //
					{ "asym_id", asym_id },
					{ "entity_id", entity_id  },
					{ "seq_id", num },
					{ "mon_id", comp_id },
					{ "ndb_seq_num", num },
					{ "pdb_seq_num", auth_seq_num },
					{ "auth_seq_num", auth_seq_num },
					{ "pdb_mon_id", auth_mon_id },
					{ "auth_mon_id", auth_mon_id },
					{ "pdb_strand_id", pdb_strand_id },
					{ "pdb_ins_code", ins_code },
					{ "hetero", hetero }
				});
			}
		}
	}
}

void reconstruct_pdbx(file &file, std::string_view dictionary)
{
	if (file.empty())
		throw std::runtime_error("Cannot reconstruct PDBx, file seems to be empty");

	// assuming the first datablock contains the entry ...
	auto &db = file.front();

	// ... and any additional datablock will contain compound information
	cif::compound_source cs(file);

	if (db.get("atom_site") == nullptr)
		throw std::runtime_error("Cannot reconstruct PDBx file, atom data missing");

	auto &validator = validator_factory::instance()[dictionary];

	std::string entry_id;

	// Phenix files do not have an entry record
	if (db.get("entry") == nullptr)
	{
		entry_id = db.name();
		category entry("entry");
		entry.emplace({ { "id", entry_id } });
		db.emplace_back(std::move(entry));
	}
	else
	{
		auto &entry = db["entry"];
		if (entry.size() != 1)
			throw std::runtime_error("Unexpected size of entry category");

		entry_id = entry.front().get<std::string>("id");
	}

	for (auto &cat : db)
	{
		auto cv = validator.get_validator_for_category(cat.name());
		if (not cv)
			continue;

		for (auto link : validator.get_links_for_child(cat.name()))
		{
			if (link->m_parent_category != "entry")
				continue;

			// So, this cat should have a link to the entry

			auto pk = find(link->m_parent_keys.begin(), link->m_parent_keys.end(), "id");
			if (pk == link->m_parent_keys.end())
				continue;

			auto ix = pk - link->m_parent_keys.begin();
			auto key = link->m_child_keys[ix];

			for (auto row : cat)
			{
				row.assign({ { key, entry_id } });
			}
		}

		// See if all categories that need a key do have a value
		if (cv->m_keys.size() == 1)
		{
			auto key = cv->m_keys.front();
			for (auto row : cat)
			{
				auto ord = row.get<std::string>(key.c_str());
				if (ord.empty())
					row.assign({ //
						{ key, cat.get_unique_id([](int nr)
								   { return std::to_string(nr); }) } });
			}
		}
	}

	file.load_dictionary(dictionary);

	// Now create any missing categories

	// First, see if atom records make sense at all
	// Will take care of atom_type and chem_comp as well.
	checkAtomRecords(db);

	// Next make sure we have struct_asym records
	if (db.get("struct_asym") == nullptr)
		createStructAsym(db);
	
	if (db.get("entity") == nullptr)
		createEntity(db);

	if (db.get("pdbx_poly_seq_scheme") == nullptr)
		createPdbxPolySeqScheme(db);
}

} // namespace cif::pdb
