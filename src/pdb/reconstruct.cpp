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

using residue_key_type = std::tuple<
	std::optional<std::string>,
	std::optional<int>,
	std::optional<std::string>,
	std::optional<std::string>,
	std::optional<int>,
	std::optional<std::string>>;

template <typename T>
auto get_either_or(std::optional<T> &a, std::optional<T> &b)
{
	if (a.has_value())
		return a.value();
	else if (b.has_value())
		return b.value();
	else
		return T{};
}

inline std::string get_asym_id(residue_key_type &k)
{
	return get_either_or(std::get<0>(k), std::get<3>(k));
}

inline int get_seq_id(residue_key_type &k)
{
	return get_either_or(std::get<1>(k), std::get<4>(k));
}

inline std::string get_comp_id(residue_key_type &k)
{
	return get_either_or(std::get<2>(k), std::get<5>(k));
}

inline bool has_asym_id(residue_key_type &k)
{
	return std::get<0>(k).has_value() or std::get<3>(k).has_value();
}

inline bool has_seq_id(residue_key_type &k)
{
	return std::get<1>(k).has_value() or std::get<4>(k).has_value();
}

inline bool has_comp_id(residue_key_type &k)
{
	return std::get<2>(k).has_value() or std::get<5>(k).has_value();
}

condition get_condition(residue_key_type &k)
{
	return key("auth_asym_id") == std::get<0>(k) and
	       key("auth_seq_id") == std::get<1>(k) and
	       key("auth_comp_id") == std::get<2>(k) and
	       key("label_asym_id") == std::get<3>(k) and
	       key("label_seq_id") == std::get<4>(k) and
	       key("label_comp_id") == std::get<5>(k);
}

// --------------------------------------------------------------------

void checkEntities(datablock &db)
{
	using namespace cif::literals;

	auto &cf = cif::compound_factory::instance();

	for (auto entity : db["entity"].find("formula_weight"_key == null or "formula_weight"_key == 0))
	{
		const auto &[entity_id, type] = entity.get<std::string, std::string>("id", "type");

		float formula_weight = 0;

		if (type == "polymer")
		{
			int n = 0;

			for (std::string comp_id : db["pdbx_poly_seq_scheme"].find<std::string>("entity_id"_key == entity_id, "mon_id"))
			{
				auto compound = cf.create(comp_id);
				assert(compound);
				if (not compound)
					throw std::runtime_error("missing information for compound " + comp_id);
				formula_weight += compound->formula_weight();
				++n;
			}

			formula_weight -= (n - 1) * 18.015;
		}
		else if (type == "water")
			formula_weight = 18.015;
		else if (type == "branched")
		{
			int n = 0;

			for (std::string comp_id : db["pdbx_entity_branch_list"].find<std::string>("entity_id"_key == entity_id, "comp_id"))
			{
				auto compound = cf.create(comp_id);
				assert(compound);
				if (not compound)
					throw std::runtime_error("missing information for compound " + comp_id);
				formula_weight += compound->formula_weight();
				++n;
			}

			formula_weight -= (n - 1) * 18.015;
		}
		else if (type == "non-polymer")
		{
			auto comp_id = db["pdbx_nonpoly_scheme"].find_first<std::optional<std::string>>("entity_id"_key == entity_id, "mon_id");
			if (comp_id.has_value())
			{
				auto compound = cf.create(*comp_id);
				if (not compound)
				{
					std::cerr << "missing information for compound " << *comp_id << "\n";
					continue;
				}
				formula_weight = compound->formula_weight();
			}
		}

		if (formula_weight > 0)
			entity.assign({ { "formula_weight", formula_weight, 3 } });
	}
}

void createEntityIDs(datablock &db)
{
	// Suppose the file does not have entity ID's. We have to make up some

	// walk the atoms. For each auth_asym_id we have a new struct_asym.
	// Within the same auth_asym_id's check for a break between polymer and
	// non-polymer atoms. If found, create new entity
	// Each residue with separate seq_id in asym with same auth_asym_id and
	// of type non-polymer is a separate struct asym
	//
	// that should cover it

	auto &atom_site = db["atom_site"];
	auto &cf = compound_factory::instance();

	std::vector<std::vector<residue_key_type>> entities;

	std::string lastAsymID;
	int lastSeqID = -1;
	std::vector<residue_key_type> waters;

	for (residue_key_type k : atom_site.rows<std::optional<std::string>,
							  std::optional<int>,
							  std::optional<std::string>,
							  std::optional<std::string>,
							  std::optional<int>,
							  std::optional<std::string>>(
			 "auth_asym_id", "auth_seq_id", "auth_comp_id",
			 "label_asym_id", "label_seq_id", "label_comp_id"))
	{
		std::string comp_id = get_comp_id(k);

		if (cf.is_water(comp_id))
		{
			waters.emplace_back(k);
			continue;
		}

		std::string asym_id = get_asym_id(k);
		int seq_id = get_seq_id(k);

		bool is_monomer = cf.is_monomer(comp_id);

		if (lastAsymID == asym_id and lastSeqID == seq_id and not is_monomer)
			continue;

		if (asym_id != lastAsymID or (not is_monomer and lastSeqID != seq_id))
			entities.push_back({});

		entities.back().emplace_back(k);

		lastAsymID = asym_id;
		lastSeqID = seq_id;
	}

	std::map<std::size_t, std::string> entity_ids;

	atom_site.add_item("label_entity_id");

	for (std::size_t i = 0; i < entities.size(); ++i)
	{
		if (entity_ids.contains(i))
			continue;

		auto entity_id = std::to_string(i + 1);
		entity_ids[i] = entity_id;

		for (std::size_t j = i + 1; j < entities.size(); ++j)
		{
			if (entities[i] == entities[j])
				entity_ids[j] = entity_id;
		}
	}

	for (std::size_t ix = 0; auto &e : entities)
	{
		auto k = e.front();
		const auto &entity_id = entity_ids[ix++];

		std::string comp_id = get_comp_id(k);

		for (auto &k : e)
			atom_site.update_value(get_condition(k), "label_entity_id", entity_id);
	}

	if (not waters.empty())
	{
		std::string waterEntityID = std::to_string(entities.size() + 1);
		for (auto &k : waters)
			atom_site.update_value(get_condition(k), "label_entity_id", waterEntityID);
	}
}

void fillLabelAsymID(category &atom_site)
{
	std::map<std::tuple<std::string, std::string>, std::string> mapAuthAsymIDAndEntityToLabelAsymID;

	// pray that label_entity_id is filled in and use that to discriminate between asyms

	if (atom_site.has_item("label_asym_id"))
	{
		for (const auto &[label_entity_id, auth_asym_id, label_asym_id] :
			atom_site.find<std::optional<std::string>, std::string, std::string>(
				key("label_asym_id") != cif::null, "label_entity_id", "auth_asym_id", "label_asym_id"))
		{
			if (not label_entity_id.has_value())
				continue;

			auto key = make_tuple(auth_asym_id, *label_entity_id);
			auto i = mapAuthAsymIDAndEntityToLabelAsymID.find(key);

			if (i == mapAuthAsymIDAndEntityToLabelAsymID.end())
				mapAuthAsymIDAndEntityToLabelAsymID.emplace(make_pair(key, label_asym_id));
			else if (i->second != label_asym_id)
			{
				if (cif::VERBOSE > 0)
					std::clog << "Inconsistent assignment of label_asym_id for the tuple entity_id: " << *label_entity_id << " and auth_asym_id: " << auth_asym_id << '\n';

				mapAuthAsymIDAndEntityToLabelAsymID.clear();
				break;
			}
		}
	}
	else
	{
		// horror scenario..
		// We filled in entity_ids, right? use those along with auth_asym_id
		// to come up with new label_asym_ids

		atom_site.add_item("label_asym_id");

		for (auto key : atom_site.rows<std::string, std::string>(
				 "auth_asym_id", "label_entity_id"))
		{
			if (not mapAuthAsymIDAndEntityToLabelAsymID.contains(key))
			{
				std::string asym_id = cif_id_for_number(mapAuthAsymIDAndEntityToLabelAsymID.size());
				mapAuthAsymIDAndEntityToLabelAsymID[key] = asym_id;
			}
		}
	}

	for (const auto &[key, value] : mapAuthAsymIDAndEntityToLabelAsymID)
	{
		const auto &[auth_asym_id, label_entity_id] = key;

		atom_site.update_value(
			cif::key("label_asym_id") == null and
				cif::key("auth_asym_id") == auth_asym_id and
				cif::key("label_entity_id") == label_entity_id,
			"label_asym_id", value);
	}

	// Check to see if we're done
	if (atom_site.contains(key("label_asym_id") == cif::null))
	{
		// nope, not yet.
		throw std::runtime_error("atom_site category still contains records with empty label_asym_id, don't know how to continue");
	}
}

void fixNegativeSeqID(category &atom_site)
{
	std::set<std::string> asymsWithNegativeSeqID;
	for (auto asym_id : atom_site.find<std::string>(key("label_seq_id") < 0, "label_asym_id"))
		asymsWithNegativeSeqID.emplace(asym_id);

	for (auto asym_id : asymsWithNegativeSeqID)
	{
		// create a pseudo entity_poly_seq first

		std::vector<std::tuple<std::string, int>> poly_seq;
		for (auto key : atom_site.find<std::string, int>(key("label_asym_id") == asym_id, "auth_seq_id", "label_seq_id"))
		{
			if (poly_seq.empty() or poly_seq.back() != key)
				poly_seq.emplace_back(key);
		}

		// simply renumber all items, but only if it is really a poly (i.e. size > 1)
		if (poly_seq.size() > 1)
		{
			int seq_id = 1;
			for (const auto &[auth_seq_id, label_seq_id] : poly_seq)
			{
				for (auto row : atom_site.find(key("label_asym_id") == asym_id and
											   key("auth_seq_id") == auth_seq_id and
											   key("label_seq_id") == label_seq_id))
				{
					row.assign("label_seq_id", std::to_string(seq_id), false, false);
				}

				++seq_id;
			}
		}
		else if (poly_seq.size() == 1) // a monomer?
		{
			const auto &[auth_seq_id, label_seq_id] = poly_seq.front();

			for (auto row : atom_site.find(key("label_asym_id") == asym_id and
										   key("auth_seq_id") == auth_seq_id and
										   key("label_seq_id") == label_seq_id))
			{
				row.assign("label_seq_id", ".", false, false);
			}
		}
	}
}

void checkChemCompRecords(datablock &db)
{
	auto &cf = compound_factory::instance();
	auto &chem_comp = db["chem_comp"];

	for (auto chem_comp_entry : chem_comp)
	{
		auto compound = cf.create(chem_comp_entry["id"].text());

		if (not compound)
			std::cerr << "Unknown compound: " << chem_comp_entry["id"].text() << '\n';
		else
		{
			std::vector<item> items;

			if (not chem_comp_entry["type"])
				items.emplace_back(item{ "type", compound->type() });
			if (not chem_comp_entry["name"])
				items.emplace_back(item{ "name", compound->name() });
			if (not chem_comp_entry["formula"])
				items.emplace_back(item{ "formula", compound->formula() });
			if (not chem_comp_entry["formula_weight"])
				items.emplace_back(item{ "formula_weight", compound->formula_weight() });

			if (not items.empty())
				chem_comp_entry.assign(items);
		}
	}
}

void checkAtomRecords(datablock &db)
{
	using namespace literals;

	auto &cf = compound_factory::instance();

	auto &atom_site = db["atom_site"];
	auto &atom_type = db["atom_type"];
	auto &chem_comp = db["chem_comp"];

	// Some common errors: missing label_asym_id for some of the atom records
	if (atom_site.contains(key("label_asym_id") == cif::null))
		fillLabelAsymID(atom_site);

	// And negative seq_id values
	if (atom_site.contains(key("label_seq_id") < 0))
		fixNegativeSeqID(atom_site);

	std::set<int> polymer_entities;
	for (int id : db["entity"].find<int>("type"_key == "polymer", "id"))
		polymer_entities.insert(id);

	std::set<std::string> missingCompounds;

	for (auto row : atom_site)
	{
		residue_key_type k = row.get<std::optional<std::string>,
			std::optional<int>,
			std::optional<std::string>,
			std::optional<std::string>,
			std::optional<int>,
			std::optional<std::string>>(
			"auth_asym_id", "auth_seq_id", "auth_comp_id",
			"label_asym_id", "label_seq_id", "label_comp_id");

		if (row["type_symbol"].empty())
			throw std::runtime_error("Missing type symbol in atom_site record");

		std::string symbol{ row["type_symbol"].text() };
		if (atom_type.count("symbol"_key == symbol) == 0)
			atom_type.emplace({ { "symbol", symbol } });

		if (not has_asym_id(k))
			throw std::runtime_error("atom_site records does not have a label_asym_id nor an auth_asym_id, cannot continue");

		if (not has_comp_id(k))
			throw std::runtime_error("atom_site records does not have a label_comp_id nor an auth_comp_id, cannot continue");

		if (not has_seq_id(k))
			throw std::runtime_error("atom_site records does not have a label_atom_id nor an auth_atom_id, cannot continue");

		std::string asym_id = get_asym_id(k);
		std::string comp_id = get_comp_id(k);

		if (missingCompounds.contains(comp_id))
			continue;

		bool is_polymer = polymer_entities.contains(row["label_entity_id"].as<int>());
		auto compound = cf.create(comp_id);

		if (not compound)
		{
			missingCompounds.insert(comp_id);
			std::cerr << "Missing compound information for " << comp_id << "\n";
			continue;
		}

		auto chem_comp_entry = chem_comp.find_first("id"_key == comp_id);

		std::optional<bool> non_std;
		if  (cf.is_monomer(comp_id))
			non_std = cf.is_std_monomer(comp_id);

		if (not chem_comp_entry)
		{
			chem_comp.emplace({ //
				{ "id", comp_id },
				{ "type", compound->type() },
				{ "mon_nstd_flag", non_std },
				{ "name", compound->name() },
				{ "formula", compound->formula() },
				{ "formula_weight", compound->formula_weight() } });
		}
		else
		{
			std::vector<item> items;

			if (not chem_comp_entry["type"])
				items.emplace_back(item{ "type", compound->type() });
			if (not chem_comp_entry["mon_nstd_flag"] and non_std.has_value())
				items.emplace_back(item{ "mon_nstd_flag", non_std });
			if (not chem_comp_entry["name"])
				items.emplace_back(item{ "name", compound->name() });
			if (not chem_comp_entry["formula"])
				items.emplace_back(item{ "formula", compound->formula() });
			if (not chem_comp_entry["formula_weight"])
				items.emplace_back(item{ "formula_weight", compound->formula_weight() });

			if (not items.empty())
				chem_comp_entry.assign(items);
		}

		if (is_polymer and not has_seq_id(k))
			throw std::runtime_error("atom_site record has peptide comp_id but no sequence number, cannot continue");

		int seq_id = get_seq_id(k);

		if (is_polymer and row["label_seq_id"].empty() and cf.is_monomer(comp_id))
			row["label_seq_id"] = std::to_string(seq_id);

		if (row["label_atom_id"].empty())
			row["label_atom_id"] = row["auth_atom_id"].text();
		if (row["label_asym_id"].empty())
			row["label_asym_id"] = row["auth_asym_id"].text();
		if (row["label_comp_id"].empty())
			row["label_comp_id"] = row["auth_comp_id"].text();
		if (row["label_atom_id"].empty())
			row["label_atom_id"] = row["auth_atom_id"].text();

		// Rewrite the coordinates and other items that look better in a fixed format
		// Be careful not to nuke invalidly formatted data here
		for (auto [item_name, prec] : std::vector<std::tuple<std::string_view, std::string::size_type>>{
				 { "cartn_x", 3 },
				 { "cartn_y", 3 },
				 { "cartn_z", 3 },
				 { "occupancy", 2 },
				 { "b_iso_or_equiv", 2 } })
		{
			if (row[item_name].empty())
				continue;

			float v;
			auto s = row.get<std::string>(item_name);
			if (auto [ptr, ec] = cif::from_chars(s.data(), s.data() + s.length(), v); (bool)ec)
				continue;

			if (s.length() < prec + 1 or s[s.length() - prec - 1] != '.')
			{
				char b[12];

				if (auto [ptr, ec] = cif::to_chars(b, b + sizeof(b), v, cif::chars_format::fixed, prec); (bool)ec)
					row.assign(item_name, { b, static_cast<std::string::size_type>(ptr - b) }, false, false);
			}
		}
	}

	// auto *cv = atom_site.get_cat_validator();
	// if (cv)
	// {
	// 	// See if there are items that are no longer known
	// 	for (auto item_name : atom_site.get_items())
	// 	{
	// 		if (cv->get_validator_for_item(item_name) != nullptr)
	// 			continue;

	// 		auto r = atom_site.find_first(key(item_name) != null);
	// 		if (not r)
	// 		{
	// 			if (cif::VERBOSE > 0)
	// 				std::clog << "Dropping unknown item " << item_name << '\n';

	// 			atom_site.remove_item(item_name);
	// 		}
	// 		else if (cif::VERBOSE > 0)
	// 			std::clog << "Keeping unknown item " << std::quoted(item_name) << " in atom_site since it is not empty\n";
	// 	}
	// }
}

void checkAtomAnisotropRecords(datablock &db)
{
	using namespace literals;

	auto &atom_site = db["atom_site"];
	auto &atom_site_anisotrop = db["atom_site_anisotrop"];

	// auto m_validator = db.get_validator();
	// if (not m_validator)
	// 	return;

	std::vector<row_handle> to_be_deleted;

	bool warnReplaceTypeSymbol = true;
	for (auto row : atom_site_anisotrop)
	{
		auto parents = atom_site_anisotrop.get_parents(row, atom_site);
		if (parents.size() != 1)
		{
			to_be_deleted.emplace_back(row);
			continue;
		}

		// this happens sometimes (Phenix):

		auto parent = parents.front();

		if (row["type_symbol"].empty())
			row["type_symbol"] = parent["type_symbol"].text();
		else if (row["type_symbol"].text() != parent["type_symbol"].text())
		{
			if (cif::VERBOSE and std::exchange(warnReplaceTypeSymbol, false))
				std::clog << "Replacing type_symbol in atom_site_anisotrop record(s)\n";
			row["type_symbol"] = parent["type_symbol"].text();
		}

		if (row["pdbx_auth_alt_id"].empty() and not parent["pdbx_auth_alt_id"].empty())
			row["pdbx_auth_alt_id"] = parent["pdbx_auth_alt_id"].text();
		if (row["pdbx_label_seq_id"].empty() and not parent["pdbx_label_seq_id"].empty())
			row["pdbx_label_seq_id"] = parent["label_seq_id"].text();
		if (row["pdbx_label_asym_id"].empty() and not parent["pdbx_label_asym_id"].empty())
			row["pdbx_label_asym_id"] = parent["label_asym_id"].text();
		if (row["pdbx_label_atom_id"].empty() and not parent["pdbx_label_atom_id"].empty())
			row["pdbx_label_atom_id"] = parent["label_atom_id"].text();
		if (row["pdbx_label_comp_id"].empty() and not parent["pdbx_label_comp_id"].empty())
			row["pdbx_label_comp_id"] = parent["label_comp_id"].text();
		// if (row["pdbx_PDB_model_num"].empty() and not parent["pdbx_PDB_model_num"].empty())
		// 	row["pdbx_PDB_model_num"] = parent["pdbx_PDB_model_num"].text();
	}

	if (not to_be_deleted.empty())
	{
		if (cif::VERBOSE > 0)
			std::clog << "Dropped " << to_be_deleted.size() << " anisotrop records since they did not have exactly one parent\n";

		for (auto row : to_be_deleted)
			atom_site_anisotrop.erase(row);
	}
}

void createStructAsym(datablock &db)
{
	auto &atom_site = db["atom_site"];
	auto &struct_asym = db["struct_asym"];

	for (const auto &[label_asym_id, entity_id] : atom_site.rows<std::string, std::string>("label_asym_id", "label_entity_id"))
	{
		if (label_asym_id.empty())
			throw std::runtime_error("File contains atom_site records without a label_asym_id");
		if (struct_asym.count(key("id") == label_asym_id) == 0)
		{
			struct_asym.emplace({
				// clang-format off
				{ "id", label_asym_id },
				{ "entity_id", entity_id }
				//clang-format on
			});
		}
	}
}

void createEntity(datablock &db)
{
	using namespace literals;

	auto &cf = compound_factory::instance();

	auto &atom_site = db["atom_site"];
	atom_site.add_item("label_entity_id");

	auto &struct_asym = db["struct_asym"];
	struct_asym.add_item("entity_id");

	std::map<std::string, std::vector<std::tuple<std::string, int>>> asyms;

	for (auto asym_id : db["struct_asym"].rows<std::string>("id"))
	{
		int last_seq_id = -1;

		for (const auto &[comp_id, seq_id] : atom_site.find<std::string, int>("label_asym_id"_key == asym_id, "label_comp_id", "label_seq_id"))
		{
			if (seq_id == last_seq_id)
				continue;

			last_seq_id = seq_id;

			asyms[asym_id].emplace_back(comp_id, last_seq_id);
		}
	}

	auto less = [](const std::vector<std::tuple<std::string, int>> &a, const std::vector<std::tuple<std::string, int>> &b)
	{
		int d = static_cast<int>(a.size()) - static_cast<int>(b.size());
		return d == 0 ? a > b : d > 0;
	};

	std::set<std::vector<std::tuple<std::string, int>>, decltype(less)> entities(less);

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
			{ "pdbx_number_of_molecules", count } });
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
		std::map<std::string, std::string> seq, seq_can;
		bool non_std_monomer = false;
		bool non_std_linkage = false;
		std::vector<std::string> pdb_strand_ids;

		for (const auto &[comp_id, seq_id, auth_asym_id] : atom_site.find<std::string, int, std::string>(
				 "label_entity_id"_key == entity_id, "label_comp_id", "label_seq_id", "auth_asym_id"))
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
				if (cf.is_base(comp_id))
				{
					c_type = "polydeoxyribonucleotide";
					letter_can = compound_factory::kBaseMap.at(comp_id);
					if (comp_id.length() == 1)
						letter = letter_can;
					else
						letter = '(' + letter_can + ')';
				}
				else if (cf.is_peptide(comp_id))
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
				else
				{
					// c_type = "other";

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

			seq[auth_asym_id] += letter;
			seq_can[auth_asym_id] += letter_can;

			if (find(pdb_strand_ids.begin(), pdb_strand_ids.end(), auth_asym_id) == pdb_strand_ids.end())
				pdb_strand_ids.emplace_back(auth_asym_id);
		}

		// sanity check, each seq should be the same

		std::string entity_seq;
		std::string entity_seq_can;

		for (const auto &[auth_asym_id_1, seq_1] : seq)
		{
			if (entity_seq.empty())
			{
				entity_seq = seq_1;
				entity_seq_can = seq_can[auth_asym_id_1];
			}

			for (const auto &[auth_asym_id_2, seq_2] : seq)
			{
				if (auth_asym_id_1 != auth_asym_id_2 and seq_1 != seq_2)
					throw std::runtime_error("Inconsistent sequences for auth_asym_id " + auth_asym_id_1 + " and " + auth_asym_id_2);
			}
		}

		for (auto i = entity_seq.begin() + 80; i < entity_seq.end(); i += 80)
			i = entity_seq.insert(i, '\n') + 1;

		for (auto i = entity_seq_can.begin() + 76; i < entity_seq_can.end(); i += 76)
		{
			auto j = i;
			while (j < i + 4 and j < entity_seq_can.end())
			{
				if (*j == '(')
					break;
				++j;
			}

			if (j < entity_seq_can.end())
				i = entity_seq_can.insert(j, '\n') + 1;
			else
				i = j;
		}

		entity_poly.emplace({ //
			{ "entity_id", entity_id },
			{ "type", type },
			{ "nstd_linkage", non_std_linkage },
			{ "nstd_monomer", non_std_monomer },
			{ "pdbx_seq_one_letter_code", entity_seq },
			{ "pdbx_seq_one_letter_code_can", entity_seq_can },
			{ "pdbx_strand_id", join(pdb_strand_ids, ",") } });
	}
}

void createEntityPolySeq(datablock &db)
{
	if (auto cat = db.get("entity_poly"); cat == nullptr or cat->empty())
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

		for (const auto &[comp_id, seq_id] : atom_site.find<std::string, int>("label_entity_id"_key == entity_id and "label_asym_id"_key == asym_id, "label_comp_id", "label_seq_id"))
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
				entity_poly_seq.back().assign({ { "hetero", true } });
			}

			entity_poly_seq.emplace({ //
				{ "entity_id", entity_id },
				{ "num", seq_id },
				{ "mon_id", comp_id },
				{ "hetero", hetero } });

			last_seq_id = seq_id;
			last_comp_id = comp_id;
		}

		// you cannot assume this is correct...
		entity_poly_seq.sort([](row_handle a, row_handle b)
			{ return a.get<int>("num") < b.get<int>("num"); });
	}
}

void createPdbxPolySeqScheme(datablock &db)
{
	if (auto cat = db.get("entity_poly"); cat == nullptr or cat->empty())
		createEntityPoly(db);

	if (auto cat = db.get("entity_poly_seq"); cat == nullptr or cat->empty())
		createEntityPolySeq(db);

	using namespace literals;

	auto &atom_site = db["atom_site"];
	auto &entity_poly = db["entity_poly"];
	auto &entity_poly_seq = db["entity_poly_seq"];
	auto &struct_asym = db["struct_asym"];
	auto &pdbx_poly_seq_scheme = db["pdbx_poly_seq_scheme"];

	// Find the mapping between asym_id and pdb_strand_id first
	std::map<std::string, std::string> asym_id_to_pdb_strand_map;

	for (const auto &[entity_id, pdb_strand_ids] : entity_poly.rows<std::string, std::string>("entity_id", "pdbx_strand_id"))
	{
		for (auto pdb_strand_id : split(pdb_strand_ids, ","))
		{
			auto asym_id = atom_site.find_first<std::string>(key("label_entity_id") == entity_id and key("auth_asym_id") == pdb_strand_id, "label_asym_id");
			asym_id_to_pdb_strand_map[asym_id] = pdb_strand_id;
		}
	}

	for (auto &entity_id : entity_poly.rows<std::string>("entity_id"))
	{
		for (auto asym_id : struct_asym.find<std::string>("entity_id"_key == entity_id, "id"))
		{
			for (const auto &[comp_id, num, hetero] : entity_poly_seq.find<std::string, int, bool>("entity_id"_key == entity_id, "mon_id", "num", "hetero"))
			{
				const auto &[auth_seq_num, auth_mon_id, ins_code] =
					atom_site.find_first<std::string, std::string, std::optional<std::string>>(
						"label_asym_id"_key == asym_id and "label_seq_id"_key == num,
						"auth_seq_id", "auth_comp_id", "pdbx_PDB_ins_code");

				pdbx_poly_seq_scheme.emplace({ //
					{ "asym_id", asym_id },
					{ "entity_id", entity_id },
					{ "seq_id", num },
					{ "mon_id", comp_id },
					{ "ndb_seq_num", num },
					{ "pdb_seq_num", auth_seq_num },
					{ "auth_seq_num", auth_seq_num },
					{ "pdb_mon_id", auth_mon_id },
					{ "auth_mon_id", auth_mon_id },
					{ "pdb_strand_id", asym_id_to_pdb_strand_map[asym_id] },
					{ "pdb_ins_code", ins_code },
					{ "hetero", hetero } });
			}
		}
	}
}

// Some programs write out a ndb_poly_seq_scheme, which has been replaced by pdbx_poly_seq_scheme
void comparePolySeqSchemes(datablock &db)
{
	auto &ndb_poly_seq_scheme = db["ndb_poly_seq_scheme"];
	auto &pdbx_poly_seq_scheme = db["pdbx_poly_seq_scheme"];

	// Since often ndb_poly_seq_scheme only contains an id and mon_id item
	// we assume that it should match the accompanying pdbx_poly_seq

	std::vector<std::string> asym_ids_ndb, asym_ids_pdbx;

	for (auto asym_id : ndb_poly_seq_scheme.rows<std::string>("id"))
	{
		auto i = std::lower_bound(asym_ids_ndb.begin(), asym_ids_ndb.end(), asym_id);
		if (i == asym_ids_ndb.end() or *i != asym_id)
			asym_ids_ndb.insert(i, asym_id);
	}

	for (auto asym_id : pdbx_poly_seq_scheme.rows<std::string>("asym_id"))
	{
		auto i = std::lower_bound(asym_ids_pdbx.begin(), asym_ids_pdbx.end(), asym_id);
		if (i == asym_ids_pdbx.end() or *i != asym_id)
			asym_ids_pdbx.insert(i, asym_id);
	}

	// If we have different Asym ID's assume the ndb is invalid.
	if (asym_ids_ndb != asym_ids_pdbx)
	{
		if (cif::VERBOSE > 0)
			std::clog << "The asym ID's of ndb_poly_seq_scheme and pdbx_poly_seq_scheme are not equal, dropping ndb_poly_seq_scheme\n";
		ndb_poly_seq_scheme.clear();
	}
	else
	{
		for (const auto &asym_id : asym_ids_ndb)
		{
			bool valid = true;

			auto ndb_range = ndb_poly_seq_scheme.find(key("id") == asym_id);
			auto pdbx_range = pdbx_poly_seq_scheme.find(key("asym_id") == asym_id);

			for (auto ndb_i = ndb_range.begin(), pdbx_i = pdbx_range.begin();
				 ndb_i != ndb_range.end() or pdbx_i != pdbx_range.end(); ++ndb_i, ++pdbx_i)
			{
				if (ndb_i == ndb_range.end() or pdbx_i == pdbx_range.end())
				{
					if (cif::VERBOSE > 0)
						std::clog << "The sequences in ndb_poly_seq_scheme and pdbx_poly_seq_scheme are unequal in size for asym ID " << asym_id << '\n';
					valid = false;
					break;
				}

				auto ndb_mon_id = ndb_i->get<std::string>("mon_id");
				auto pdbx_mon_id = pdbx_i->get<std::string>("mon_id");

				if (ndb_mon_id != pdbx_mon_id)
				{
					if (cif::VERBOSE > 0)
						std::clog << "The sequences in ndb_poly_seq_scheme and pdbx_poly_seq_scheme contain different mon ID's for asym ID " << asym_id << '\n';
					valid = false;
					break;
				}
			}

			if (not valid)
			{
				if (cif::VERBOSE > 0)
					std::clog << "Dropping asym ID " << asym_id << " from ndb_poly_seq_scheme\n";
				ndb_poly_seq_scheme.erase(key("id") == asym_id);
			}
		}
	}
}

bool reconstruct_pdbx(file &file, std::string_view dictionary)
{
	if (file.empty())
		throw std::runtime_error("Cannot reconstruct PDBx, file seems to be empty");

	// assuming the first datablock contains the entry ...
	auto &db = file.front();

	// ... and any additional datablock will contain compound information
	cif::compound_source cs(file);

	if (auto cat = db.get("atom_site"); cat == nullptr or cat->empty())
		throw std::runtime_error("Cannot reconstruct PDBx file, atom data missing");

	auto &validator = validator_factory::instance()[dictionary];

	std::string entry_id;

	// Phenix files do not have an entry record
	if (auto cat = db.get("entry"); cat == nullptr or cat->empty())
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

	// Start with chem_comp, it is often missing many fields
	// that can easily be filled in.
	checkChemCompRecords(db);

	// If the data is really horrible, it might not contain entities
	if (not db["atom_site"].find_first(key("label_entity_id") != null))
		createEntityIDs(db);

	// Now see if atom records make sense at all
	checkAtomRecords(db);

	std::vector<std::string> invalidCategories;

	// clean up each category
	for (auto &cat : db)
	{
		try
		{
			auto cv = validator.get_validator_for_category(cat.name());
			if (not cv)
				continue;

			// Start by renaming items that may have old names based on alias info

			for (auto item_name : cat.get_items())
			{
				auto iv = cv->get_validator_for_item(item_name);
				if (iv) // know, must be OK then`
					continue;

				iv = cv->get_validator_for_aliased_item(item_name);
				if (not iv)
					continue;

				if (cif::VERBOSE > 0)
					std::clog << "Renaming " << item_name << " to " << iv->m_item_name << " in category " << cat.name() << '\n';
				cat.rename_item(item_name, iv->m_item_name);
			}

			// In case a single ID field is missing, add it
			if (cv->m_keys.size() == 1 and not cat.has_item(cv->m_keys.front()))
			{
				std::string key = cv->m_keys.front();

				auto iv = cv->get_validator_for_item(key);
				bool number = iv != nullptr and
				              iv->m_type != nullptr and
				              iv->m_type->m_primitive_type == cif::DDL_PrimitiveType::Numb;

				for (std::size_t ix = 0; auto row : cat)
				{
					if (number)
						row.assign(key, std::to_string(++ix), false, false);
					else
						row.assign(key, cif::cif_id_for_number(ix++), false, false);
				}
			}

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

			// Fill in all mandatory items
			for (auto item : cv->m_mandatory_items)
			{
				if (not cat.has_item(item))
				{
					if (cif::VERBOSE > 0)
						std::clog << "Adding mandatory item " << item << " to category " << cat.name() << '\n';

					cat.add_item(item);

					cat.update_value(all(), item, "?");
				}
			}

			// validate all values, and if they do not validate replace the content with an unknown flag
			for (auto item_name : cat.get_items())
			{
				auto iv = cv->get_validator_for_item(item_name);
				if (not iv)
				{
					// Drop this item
					cat.remove_item(item_name);
					continue;
				}

				auto ix = cat.get_item_ix(item_name);

				for (auto row : cat)
				{
					std::error_code ec;
					std::string_view value = row[ix].text();

					if (not iv->validate_value(value, ec))
					{
						if (cif::VERBOSE > 0)
							std::clog << "Replacing value (" << std::quoted(value) << ") for item " << item_name << " in category " << cat.name() << " since it does not validate\n";

						row[ix] = "?";
					}
				}
			}

			enum class State
			{
				Start,
				MissingKeys,
				DuplicateKeys
			} state = State::Start;

			for (;;)
			{
				// See if we can build an index
				try
				{
					cat.set_validator(&validator, db);
				}
				catch (const missing_key_error &ex)
				{
					if (state == State::MissingKeys)
					{
						if (cif::VERBOSE > 0)
							std::clog << "Repairing failed for category " << cat.name() << ", missing keys remain: " << ex.what() << '\n';

						throw;
					}

					state = State::MissingKeys;

					auto key = ex.get_key();

					if (cif::VERBOSE > 0)
						std::clog << "Need to add key " << key << " to category " << cat.name() << '\n';

					for (auto row : cat)
					{
						auto ord = row.get<std::string>(key.c_str());
						if (ord.empty())
							row.assign({ //
								{ key, cat.get_unique_value(key) } });
					}

					continue;
				}
				catch (const duplicate_key_error &ex)
				{
					if (state == State::DuplicateKeys)
					{
						if (cif::VERBOSE > 0)
							std::clog << "Repairing failed for category " << cat.name() << ", duplicate keys remain: " << ex.what() << '\n';

						throw;
					}

					state = State::DuplicateKeys;

					if (cif::VERBOSE > 0)
						std::clog << "Attempt to fix " << cat.name() << " failed: " << ex.what() << '\n';

					// replace items that do not define a relation to a parent

					std::set<std::string> replaceableKeys;
					for (auto key : cv->m_keys)
					{
						bool replaceable = true;
						for (auto lv : validator.get_links_for_child(cat.name()))
						{
							if (find(lv->m_child_keys.begin(), lv->m_child_keys.end(), key) != lv->m_child_keys.end())
							{
								replaceable = false;
								break;
							}
						}

						if (replaceable)
							replaceableKeys.insert(key);
					}

					if (replaceableKeys.empty())
						throw std::runtime_error("Cannot repair category " + cat.name() + " since it contains duplicate keys that cannot be replaced");

					for (auto key : replaceableKeys)
					{
						for (auto row : cat)
							row.assign(key, cat.get_unique_value(key), false, false);
					}

					continue;
				}

				break;
			}
		}
		catch (const std::exception &ex)
		{
			if (cif::VERBOSE > 0)
				std::clog << ex.what() << '\n';

			std::clog << "Will drop category " << cat.name() << " since it cannot be repaired\n";

			invalidCategories.emplace_back(cat.name());
		}
	}

	for (auto cat_name : invalidCategories)
	{
		auto i = find_if(db.begin(), db.end(), [cat_name](const category &cat)
			{ return cat.name() == cat_name; });
		if (i != db.end())
			db.erase(i);
	}

	db["chem_comp"].reorder_by_index();

	file.load_dictionary(dictionary);

	if (db.get("atom_site_anisotrop"))
		checkAtomAnisotropRecords(db);

	// Now create any missing categories
	// Next make sure we have struct_asym records
	if (auto cat = db.get("struct_asym"); cat == nullptr or cat->empty())
		createStructAsym(db);

	if (auto cat = db.get("entity"); cat == nullptr or cat->empty())
		createEntity(db);

	// fill in missing formula_weight, e.g.
	checkEntities(db);

	if (auto cat = db.get("pdbx_poly_seq_scheme"); cat == nullptr or cat->empty())
		createPdbxPolySeqScheme(db);

	if (auto cat = db.get("ndb_poly_seq_scheme"); cat == nullptr or cat->empty())
		comparePolySeqSchemes(db);

	// skip unknown categories for now
	bool valid = true;
	for (auto &cat : db)
		valid = valid and (cat.get_cat_validator() == nullptr or cat.is_valid());

	return valid and is_valid_pdbx_file(file, dictionary);
}

} // namespace cif::pdb
