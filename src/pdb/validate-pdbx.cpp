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

namespace cif::pdb
{

condition get_parents_condition(const validator &validator, row_handle rh, const category &parentCat)
{
	condition result;

	auto &childCat = rh.get_category();
	auto childName = childCat.name();
	auto parentName = parentCat.name();

	auto links = validator.get_links_for_child(childName);
	links.erase(remove_if(links.begin(), links.end(), [n = parentName](auto &l)
					{ return l->m_parent_category != n; }),
		links.end());

	if (not links.empty())
	{
		for (auto &link : links)
		{
			condition cond;

			for (std::size_t ix = 0; ix < link->m_child_keys.size(); ++ix)
			{
				auto childValue = rh[link->m_child_keys[ix]];

				if (childValue.empty())
					continue;

				cond = std::move(cond) and key(link->m_parent_keys[ix]) == childValue.text();
			}

			result = std::move(result) or std::move(cond);
		}
	}
	else if (cif::VERBOSE > 0)
		std::cerr << "warning: no child to parent links were found for child " << childName << " and parent " << parentName << '\n';

	return result;
}

bool is_valid_pdbx_file(const file &file, std::string_view dictionary)
{
	std::error_code ec;
	bool result = is_valid_pdbx_file(file, dictionary, ec);
	return result and not (bool)ec;
}

bool is_valid_pdbx_file(const file &file, std::error_code &ec)
{
	bool result = false;

	if (file.empty())
		ec = make_error_code(validation_error::empty_file);
	else
	{
		std::string dictionary = "mmcif_pdbx";

		for (auto &db : file)
		{
			auto audit_conform = db.get("audit_conform");
			if (audit_conform == nullptr)
				continue;
			
			if (not audit_conform->empty())
			{
				auto specified_dict = audit_conform->front()["dict_name"];
				if (not specified_dict.empty())
					dictionary = specified_dict.as<std::string>();
			}

			break;
		}

		result = is_valid_pdbx_file(file, dictionary, ec);
	}
	
	return result;
}

bool is_valid_pdbx_file(const file &file, std::string_view dictionary, std::error_code &ec)
{
	using namespace cif::literals;

	bool result = true;

	try
	{
		auto &cf = cif::compound_factory::instance();
		auto &validator = cif::validator_factory::instance().operator[](dictionary);

		if (file.empty())
			throw std::runtime_error("Empty file");

		auto &db = file.front();

		if (db.empty())
			throw std::runtime_error("Empty datablock");

		auto &atom_site = db["atom_site"];
		if (atom_site.empty())
			throw std::runtime_error("Empty or missing atom_site category");

		auto &pdbx_poly_seq_scheme = db["pdbx_poly_seq_scheme"];

		std::string last_asym_id;
		int last_seq_id = -1;
		for (auto r : atom_site)
		{
			auto seq_id = r.get<std::optional<int>>("label_seq_id");
			if (not seq_id.has_value()) // not a residue in a polymer
				continue;

			if (*seq_id == last_seq_id)
				continue;

			last_seq_id = *seq_id;

			auto comp_id = r.get<std::string>("label_comp_id");
			if (not cf.is_monomer(comp_id))
				continue;

			auto p = pdbx_poly_seq_scheme.find(get_parents_condition(validator, r, pdbx_poly_seq_scheme));
			if (p.size() != 1)
			{
				if (cif::VERBOSE > 0)
					std::clog << "In atom_site record: " << r["id"].text() << '\n';
				throw std::runtime_error("For each monomer in atom_site there should be exactly one pdbx_poly_seq_scheme record");
			}
		}

		auto &entity = db["entity"];
		if (entity.empty())
			throw std::runtime_error("Entity category is missing or empty");

		auto &entity_poly = db["entity_poly"];
		if (entity_poly.empty())
			throw std::runtime_error("Entity_poly category is missing or empty");

		auto &entity_poly_seq = db["entity_poly_seq"];
		if (entity_poly_seq.empty())
			throw std::runtime_error("Entity_poly_seq category is missing or empty");

		auto &struct_asym = db["struct_asym"];
		if (struct_asym.empty())
			throw std::runtime_error("struct_asym category is missing or empty");

		for (auto entity_id : entity.find<std::string>("type"_key == "polymer", "id"))
		{
			if (entity_poly.count("entity_id"_key == entity_id) != 1)
				throw std::runtime_error("There should be exactly one entity_poly record per polymer entity");

			const auto entity_poly_type = entity_poly.find1<std::string>("entity_id"_key == entity_id, "type");

			std::map<int,std::set<std::string>> mon_per_seq_id;

			for (const auto &[num, mon_id, hetero] : entity_poly_seq.find<int, std::string, bool>("entity_id"_key == entity_id, "num", "mon_id", "hetero"))
			{
				mon_per_seq_id[num].emplace(mon_id);

				for (auto asym_id : struct_asym.find<std::string>("entity_id"_key == entity_id, "id"))
				{
					if (pdbx_poly_seq_scheme.count(
							"entity_id"_key == entity_id and
							"asym_id"_key == asym_id and
							"mon_id"_key == mon_id and
							"seq_id"_key == num and
							"hetero"_key == hetero) != 1)
					{
						throw std::runtime_error("For each entity_poly_seq record there should be exactly one pdbx_poly_seq record");
					}
				}
			}

			for (const auto &[seq_id, mon_id, hetero] : pdbx_poly_seq_scheme.find<int, std::string, bool>("entity_id"_key == entity_id, "seq_id", "mon_id", "hetero"))
			{
				if (entity_poly_seq.count(
						"entity_id"_key == entity_id and
						"mon_id"_key == mon_id and
						"num"_key == seq_id and
						"hetero"_key == hetero) != 1)
				{
					throw std::runtime_error("For each pdbx_poly_seq/struct_asym record there should be exactly one entity_poly_seq record");
				}

				if ((mon_per_seq_id[seq_id].size() > 1) != hetero)
					throw std::runtime_error("Mismatch between the hetero flag in the poly seq schemes and the number residues per seq_id");
			}

			for (const auto &[seq_id, mon_ids] : mon_per_seq_id)
			{
				for (auto asym_id : struct_asym.find<std::string>("entity_id"_key == entity_id, "id"))
				{
					condition cond;
					
					for (auto mon_id : mon_ids)
						cond = std::move(cond) or "label_comp_id"_key == mon_id;

					cond = "label_entity_id"_key == entity_id and
						"label_asym_id"_key == asym_id and
						"label_seq_id"_key == seq_id and not std::move(cond);
					
					if (atom_site.contains(std::move(cond)))
						throw std::runtime_error("An atom_site record exists that has no parent in the poly seq scheme categories");
				}
			}

			auto &&[seq, seq_can] = entity_poly.find1<std::optional<std::string>, std::optional<std::string>>("entity_id"_key == entity_id,
				"pdbx_seq_one_letter_code", "pdbx_seq_one_letter_code_can");
			
			std::string::const_iterator si, sci, se, sce;

			auto seq_match = [&](bool can, std::string::const_iterator si, std::string::const_iterator se)
			{
				for (const auto &[seq_id, comp_ids] : mon_per_seq_id)
				{
					if (si == se)
						return false;

					bool match = false;

					for (auto comp_id : comp_ids)
					{
						std::string letter;

						if (can)
						{
							if (compound_factory::kBaseMap.contains(comp_id))
								letter = compound_factory::kBaseMap.at(comp_id);
							else
							{
								auto c = cf.create(comp_id);
								if (c and c->one_letter_code())
									letter = c->one_letter_code();
								else
									letter = "X";
							}
						}
						else
						{
							if (compound_factory::kAAMap.contains(comp_id))
								letter = compound_factory::kAAMap.at(comp_id);
							else if (comp_id.length() == 1 and compound_factory::kBaseMap.contains(comp_id))
								letter = compound_factory::kBaseMap.at(comp_id);
							else
								letter = '(' + comp_id + ')';
						}
						
						if (iequals(std::string{si, si + letter.length()}, letter))
						{
							match = true;
							si += letter.length();
							break;
						}
						else
							return false;
					}

					if (not match)
						break;
				}

				return si == se;
			};

			if (not seq.has_value())
			{
				if (cif::VERBOSE > 0)
					std::clog << "Warning: entity_poly has no sequence for entity_id " << entity_id << '\n';
			}
			else
			{
				seq->erase(std::remove_if(seq->begin(), seq->end(), [](char ch) { return std::isspace(ch); }), seq->end());

				if (not seq_match(false, seq->begin(), seq->end()))
					throw std::runtime_error("Sequences do not match for entity " + entity_id);
			}

			if (not seq_can.has_value())
			{
				if (cif::VERBOSE > 0)
					std::clog << "Warning: entity_poly has no sequence for entity_id " << entity_id << '\n';
			}
			else
			{
				seq_can->erase(std::remove_if(seq_can->begin(), seq_can->end(), [](char ch) { return std::isspace(ch); }), seq_can->end());

				if (not seq_match(true, seq_can->begin(), seq_can->end()))
					throw std::runtime_error("Canonical sequences do not match for entity " + entity_id);
			}
		}

		result = true;
	}
	catch (const std::exception &ex)
	{
		result = false;
		if (cif::VERBOSE > 0)
			std::clog << ex.what() << '\n';
		ec = make_error_code(validation_error::not_valid_pdbx);
	}

	if (not result and (bool)ec)
		ec = make_error_code(validation_error::not_valid_pdbx);

	return result;
}

} // namespace cif::pdb
  