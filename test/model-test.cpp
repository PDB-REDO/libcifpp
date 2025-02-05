/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute
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

#include "test-main.hpp"

#include <stdexcept>

#include <cif++.hpp>

// --------------------------------------------------------------------

cif::file operator""_cf(const char *text, std::size_t length)
{
	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(text), length);

	std::istream is(&buffer);
	return cif::file(is);
}

// --------------------------------------------------------------------

TEST_CASE("create_nonpoly_1")
{
	cif::VERBOSE = 1;

	cif::file file;
	file.load_dictionary("mmcif_pdbx.dic");
	file.emplace("TEST"); // create a datablock

	cif::mm::structure structure(file);

	std::string entity_id = structure.create_non_poly_entity("HEM");

	auto atoms = R"(
data_HEM
loop_
_atom_site.id
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
1 HETATM C  CHA . ? -5.248  39.769 -0.250  1.00 7.67  ?
2 HETATM C  CHB . ? -3.774  36.790 3.280   1.00 7.05  ?
3 HETATM C  CHC . ? -2.879  33.328 0.013   1.00 7.69  ?
4 HETATM C  CHD . ? -4.342  36.262 -3.536  1.00 8.00  ?
# that's enough to test with
)"_cf;

	atoms.load_dictionary("mmcif_pdbx.dic");

	auto &hem_data = atoms["HEM"];
	auto &atom_site = hem_data["atom_site"];

	auto hem_atoms = atom_site.rows();
	std::vector<cif::mm::atom> atom_data;
	for (auto hem_atom : hem_atoms)
		atom_data.emplace_back(hem_data, hem_atom);

	structure.create_non_poly(entity_id, atom_data);

	auto expected = R"(
data_TEST
# 
_pdbx_nonpoly_scheme.asym_id         A 
_pdbx_nonpoly_scheme.ndb_seq_num     1 
_pdbx_nonpoly_scheme.entity_id       1 
_pdbx_nonpoly_scheme.mon_id          HEM 
_pdbx_nonpoly_scheme.pdb_seq_num     1 
_pdbx_nonpoly_scheme.auth_seq_num    1 
_pdbx_nonpoly_scheme.pdb_mon_id      HEM 
_pdbx_nonpoly_scheme.auth_mon_id     HEM 
_pdbx_nonpoly_scheme.pdb_strand_id   A 
_pdbx_nonpoly_scheme.pdb_ins_code    . 
#
loop_
_atom_site.id
_atom_site.auth_asym_id
_atom_site.label_alt_id
_atom_site.label_asym_id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.type_symbol
_atom_site.group_PDB
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
1 A ? A CHA HEM 1 . C HETATM ? -5.248 39.769 -0.250 1.00 7.67 ? 1 HEM CHA 1
2 A ? A CHB HEM 1 . C HETATM ? -3.774 36.790 3.280  1.00 7.05 ? 1 HEM CHB 1
3 A ? A CHC HEM 1 . C HETATM ? -2.879 33.328 0.013  1.00 7.69 ? 1 HEM CHC 1
4 A ? A CHD HEM 1 . C HETATM ? -4.342 36.262 -3.536 1.00 8.00 ? 1 HEM CHD 1
#
_chem_comp.id               HEM
_chem_comp.type             NON-POLYMER
_chem_comp.name             'PROTOPORPHYRIN IX CONTAINING FE'
_chem_comp.formula          'C34 H32 Fe N4 O4'
_chem_comp.formula_weight   616.487000
#
_pdbx_entity_nonpoly.entity_id   1
_pdbx_entity_nonpoly.name        'PROTOPORPHYRIN IX CONTAINING FE'
_pdbx_entity_nonpoly.comp_id     HEM
#
_entity.id                 1
_entity.type               non-polymer
_entity.pdbx_description   'PROTOPORPHYRIN IX CONTAINING FE'
_entity.formula_weight     616.487000
#
_struct_asym.id                            A
_struct_asym.entity_id                     1
_struct_asym.pdbx_blank_PDB_chainid_flag   N
_struct_asym.pdbx_modified                 N
_struct_asym.details                       ?
#
_atom_type.symbol   C
)"_cf;

	expected.load_dictionary("mmcif_pdbx.dic");

	if (not(expected.front() == structure.get_datablock()))
	{
		std::cerr << expected.front() << '\n'
				  << '\n'
				  << structure.get_datablock() << '\n';
		REQUIRE(false);
	}
}

// --------------------------------------------------------------------

TEST_CASE("create_nonpoly_2")
{
	cif::VERBOSE = 1;

	cif::file file;
	file.load_dictionary("mmcif_pdbx.dic");
	file.emplace("TEST"); // create a datablock

	cif::mm::structure structure(file);

	cif::file lig(gTestDir / "HEM.cif");
	auto &chem_comp_atom = lig["HEM"]["chem_comp_atom"];

	std::vector<cif::row_initializer> atoms;

	for (const auto &[type_symbol, label_atom_id, Cartn_x, Cartn_y, Cartn_z] :
		chem_comp_atom.rows<std::string, std::string, float, float, float>(
			"type_symbol", "atom_id", "model_Cartn_x", "model_Cartn_y", "model_Cartn_z"))
	{
		atoms.emplace_back(cif::row_initializer{
			{ "type_symbol", type_symbol },
			{ "label_atom_id", label_atom_id },
			{ "auth_atom_id", label_atom_id },
			{ "Cartn_x", Cartn_x },
			{ "Cartn_y", Cartn_y },
			{ "Cartn_z", Cartn_z } });

		if (atoms.size() == 4)
			break;
	}

	std::string entity_id = structure.create_non_poly_entity("HEM");
	structure.create_non_poly(entity_id, atoms);

	auto expected = R"(
data_TEST
# 
_pdbx_nonpoly_scheme.asym_id         A 
_pdbx_nonpoly_scheme.ndb_seq_num     1 
_pdbx_nonpoly_scheme.entity_id       1 
_pdbx_nonpoly_scheme.mon_id          HEM 
_pdbx_nonpoly_scheme.pdb_seq_num     1 
_pdbx_nonpoly_scheme.auth_seq_num    1 
_pdbx_nonpoly_scheme.pdb_mon_id      HEM 
_pdbx_nonpoly_scheme.auth_mon_id     HEM 
_pdbx_nonpoly_scheme.pdb_strand_id   A 
_pdbx_nonpoly_scheme.pdb_ins_code    . 
#
loop_
_atom_site.id
_atom_site.auth_asym_id
_atom_site.label_alt_id
_atom_site.label_asym_id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.type_symbol
_atom_site.group_PDB
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
1 A ? A CHA HEM 1 . C HETATM ? 2.748 -19.531 39.896 1.00 ? 1 HEM CHA 1
2 A ? A CHB HEM 1 . C HETATM ? 3.258 -17.744 35.477 1.00 ? 1 HEM CHB 1
3 A ? A CHC HEM 1 . C HETATM ? 1.703 -21.9   33.637 1.00 ? 1 HEM CHC 1
4 A ? A CHD HEM 1 . C HETATM ? 1.149 -23.677 38.059 1.00 ? 1 HEM CHD 1
#
_chem_comp.id               HEM
_chem_comp.type             NON-POLYMER
_chem_comp.name             'PROTOPORPHYRIN IX CONTAINING FE'
_chem_comp.formula          'C34 H32 Fe N4 O4'
_chem_comp.formula_weight   616.487000
#
_pdbx_entity_nonpoly.entity_id   1
_pdbx_entity_nonpoly.name        'PROTOPORPHYRIN IX CONTAINING FE'
_pdbx_entity_nonpoly.comp_id     HEM
#
_entity.id                 1
_entity.type               non-polymer
_entity.pdbx_description   'PROTOPORPHYRIN IX CONTAINING FE'
_entity.formula_weight     616.487000
#
_struct_asym.id                            A
_struct_asym.entity_id                     1
_struct_asym.pdbx_blank_PDB_chainid_flag   N
_struct_asym.pdbx_modified                 N
_struct_asym.details                       ?
#
_atom_type.symbol   C
)"_cf;

	expected.load_dictionary("mmcif_pdbx.dic");

	REQUIRE(expected.front() == structure.get_datablock());

	if (not(expected.front() == structure.get_datablock()))
	{
		// REQUIRE(false);
		std::cout << expected.front() << '\n'
				  << '\n'
				  << structure.get_datablock() << '\n';

		expected.save("/tmp/a");
		file.save("/tmp/b");
	}
}

// --------------------------------------------------------------------

TEST_CASE("test_atom_id")
{
	auto data = R"(
data_TEST
# 
_pdbx_nonpoly_scheme.asym_id         A 
_pdbx_nonpoly_scheme.ndb_seq_num     1 
_pdbx_nonpoly_scheme.entity_id       1 
_pdbx_nonpoly_scheme.mon_id          HEM 
_pdbx_nonpoly_scheme.pdb_seq_num     1 
_pdbx_nonpoly_scheme.auth_seq_num    1 
_pdbx_nonpoly_scheme.pdb_mon_id      HEM 
_pdbx_nonpoly_scheme.auth_mon_id     HEM 
_pdbx_nonpoly_scheme.pdb_strand_id   A 
_pdbx_nonpoly_scheme.pdb_ins_code    . 
#
loop_
_atom_site.id
_atom_site.auth_asym_id
_atom_site.label_alt_id
_atom_site.label_asym_id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.type_symbol
_atom_site.group_PDB
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
1 A ? A CHA HEM 1 . C HETATM ? -5.248 39.769 -0.250 1.00 7.67 ? 1 HEM CHA 1
3 A ? A CHB HEM 1 . C HETATM ? -3.774 36.790 3.280  1.00 7.05 ? 1 HEM CHB 1
2 A ? A CHC HEM 1 . C HETATM ? -2.879 33.328 0.013  1.00 7.69 ? 1 HEM CHC 1
4 A ? A CHD HEM 1 . C HETATM ? -4.342 36.262 -3.536 1.00 8.00 ? 1 HEM CHD 1
#
_chem_comp.id               HEM
_chem_comp.type             NON-POLYMER
_chem_comp.name             'PROTOPORPHYRIN IX CONTAINING FE'
_chem_comp.formula          'C34 H32 Fe N4 O4'
_chem_comp.formula_weight   616.487000
#
_pdbx_entity_nonpoly.entity_id   1
_pdbx_entity_nonpoly.name        'PROTOPORPHYRIN IX CONTAINING FE'
_pdbx_entity_nonpoly.comp_id     HEM
#
_entity.id                 1
_entity.type               non-polymer
_entity.pdbx_description   'PROTOPORPHYRIN IX CONTAINING FE'
_entity.formula_weight     616.487000
#
_struct_asym.id                            A
_struct_asym.entity_id                     1
_struct_asym.pdbx_blank_PDB_chainid_flag   N
_struct_asym.pdbx_modified                 N
_struct_asym.details                       ?
#
)"_cf;

	data.load_dictionary("mmcif_pdbx.dic");

	cif::mm::structure s(data);

	REQUIRE(s.get_atom_by_id("1").get_label_atom_id() == "CHA");
	REQUIRE(s.get_atom_by_id("2").get_label_atom_id() == "CHC");
	REQUIRE(s.get_atom_by_id("3").get_label_atom_id() == "CHB");
	REQUIRE(s.get_atom_by_id("4").get_label_atom_id() == "CHD");
}

// --------------------------------------------------------------------

TEST_CASE("atom_numbers_1")
{
	const std::filesystem::path test1(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file(test1.string());
	cif::mm::structure structure(file);

	auto &db = file.front();

	auto &atoms = structure.atoms();
	auto ai = atoms.begin();

	for (const auto &[id, label_asym_id, label_seq_id, label_atom_id, auth_seq_id, label_comp_id] :
		db["atom_site"].rows<std::string, std::string, int, std::string, std::string, std::string>("id", "label_asym_id", "label_seq_id", "label_atom_id", "auth_seq_id", "label_comp_id"))
	{
		auto atom = structure.get_atom_by_id(id);

		REQUIRE(atom.get_label_asym_id() == label_asym_id);
		REQUIRE(atom.get_label_seq_id() == label_seq_id);
		REQUIRE(atom.get_label_atom_id() == label_atom_id);
		REQUIRE(atom.get_auth_seq_id() == auth_seq_id);
		REQUIRE(atom.get_label_comp_id() == label_comp_id);

		REQUIRE(ai != atoms.end());

		REQUIRE(ai->id() == id);
		++ai;
	}

	REQUIRE(ai == atoms.end());
}
// --------------------------------------------------------------------

TEST_CASE("test_load_2")
{
	using namespace cif::literals;

	const std::filesystem::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file(example.string());

	auto &db = file.front();

	cif::mm::structure s(file);

	REQUIRE(s.polymers().size() == 1UL);

	auto &pdbx_poly_seq_scheme = db["pdbx_poly_seq_scheme"];

	for (auto &poly : s.polymers())
	{
		REQUIRE(poly.size() == pdbx_poly_seq_scheme.find("asym_id"_key == poly.get_asym_id()).size());
	}
}

TEST_CASE("remove_residue_1")
{
	using namespace cif::literals;

	const std::filesystem::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file(example.string());

	cif::mm::structure s(file);
	s.remove_residue(s.get_residue("B"));

	REQUIRE_NOTHROW(s.validate_atoms());
}
