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

#include "test-main.hpp"

#include <cif++.hpp>

#include <stdexcept>

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

TEST_CASE("test-1")
{
	auto f = R"(data_1CBS
# 
_entry.id   1CBS 
# 
_entity.id                     1
_entity.type                   polymer
# 
_entity_poly.entity_id                      1 
_entity_poly.type                           'polypeptide(L)' 
_entity_poly.nstd_linkage                   no 
_entity_poly.nstd_monomer                   no 
_entity_poly.pdbx_seq_one_letter_code       
;PNFSG
;
_entity_poly.pdbx_seq_one_letter_code_can   
;PNFSG
;
_entity_poly.pdbx_strand_id                 A 
_entity_poly.pdbx_target_identifier         ? 
# 
loop_
_entity_poly_seq.entity_id 
_entity_poly_seq.num 
_entity_poly_seq.mon_id 
_entity_poly_seq.hetero 
1 1   PRO n 
1 2   ASN n 
1 3   PHE n 
1 4   SER n 
1 5   GLY n 
#
loop_
_struct_asym.id 
_struct_asym.pdbx_blank_PDB_chainid_flag 
_struct_asym.pdbx_modified 
_struct_asym.entity_id 
_struct_asym.details 
A N N 1 ? 
# 
loop_
_atom_type.symbol 
C 
N 
O 
S 
# 
loop_
_atom_site.group_PDB 
_atom_site.id 
_atom_site.type_symbol 
_atom_site.label_atom_id 
_atom_site.label_alt_id 
_atom_site.label_comp_id 
_atom_site.label_asym_id 
_atom_site.label_entity_id 
_atom_site.label_seq_id 
_atom_site.pdbx_PDB_ins_code 
_atom_site.Cartn_x 
_atom_site.Cartn_y 
_atom_site.Cartn_z 
_atom_site.occupancy 
_atom_site.B_iso_or_equiv 
_atom_site.pdbx_formal_charge 
_atom_site.auth_seq_id 
_atom_site.auth_comp_id 
_atom_site.auth_asym_id 
_atom_site.auth_atom_id 
_atom_site.pdbx_PDB_model_num 
ATOM   2    C CA  . PRO A 1 1   ? 18.150 13.525 43.680 1.00 28.82 ? 1   PRO A CA  1 
ATOM   9    C CA  . ASN A 1 2   ? 20.576 16.457 43.578 1.00 20.79 ? 2   ASN A CA  1 
ATOM   17   C CA  . PHE A 1 3   ? 21.144 17.838 40.087 1.00 12.62 ? 3   PHE A CA  1 
ATOM   28   C CA  . SER A 1 4   ? 23.170 20.780 41.464 1.00 11.30 ? 4   SER A CA  1 
ATOM   34   C CA  . GLY A 1 5   ? 26.628 21.486 40.103 1.00 10.86 ? 5   GLY A CA  1 
# 
loop_
_pdbx_poly_seq_scheme.asym_id 
_pdbx_poly_seq_scheme.entity_id 
_pdbx_poly_seq_scheme.seq_id 
_pdbx_poly_seq_scheme.mon_id 
_pdbx_poly_seq_scheme.ndb_seq_num 
_pdbx_poly_seq_scheme.pdb_seq_num 
_pdbx_poly_seq_scheme.auth_seq_num 
_pdbx_poly_seq_scheme.pdb_mon_id 
_pdbx_poly_seq_scheme.auth_mon_id 
_pdbx_poly_seq_scheme.pdb_strand_id 
_pdbx_poly_seq_scheme.pdb_ins_code 
_pdbx_poly_seq_scheme.hetero 
A 1 1   PRO 1   1   1   PRO PRO A . n 
A 1 2   ASN 2   2   2   ASN ASN A . n 
A 1 3   PHE 3   3   3   PHE PHE A . n 
A 1 4   SER 4   4   4   SER SER A . n 
A 1 5   GLY 5   5   5   GLY GLY A . n 
# 
)"_cf;

	SECTION("Plain file")
	{
		REQUIRE(cif::pdb::is_valid_pdbx_file(f));
	}

	SECTION("Delete one atom_site")
	{
		auto &db = f.front();
		auto n = db["atom_site"].erase(cif::key("id") == 2);

		REQUIRE(n == 1);

		REQUIRE(cif::pdb::is_valid_pdbx_file(f));
	}

	SECTION("Delete a pdbx_poly_seq_scheme record")
	{
		auto &db = f.front();
		auto n = db["pdbx_poly_seq_scheme"].erase(cif::key("seq_id") == 2);

		REQUIRE(n == 1);

		REQUIRE_FALSE(cif::pdb::is_valid_pdbx_file(f));
	}

	SECTION("Delete an entity_poly_seq record")
	{
		auto &db = f.front();
		auto n = db["entity_poly_seq"].erase(cif::key("num") == 2);

		REQUIRE(n == 1);

		REQUIRE_FALSE(cif::pdb::is_valid_pdbx_file(f));
	}

	SECTION("Delete an entity_poly record")
	{
		auto &db = f.front();
		auto n = db["entity_poly"].erase(cif::key("entity_id") == 1);

		REQUIRE(n == 1);

		REQUIRE_FALSE(cif::pdb::is_valid_pdbx_file(f));
	}

	SECTION("Mutate an atom_site record")
	{
		auto &db = f.front();
		auto r = db["atom_site"].find1(cif::key("id") == 9);
		r.assign({
			{ "label_comp_id", "ALA" },
			{ "auth_comp_id", "ALA" }
		});

		REQUIRE_FALSE(cif::pdb::is_valid_pdbx_file(f));
	}

	SECTION("Hetero consistency")
	{
		auto &db = f.front();
		db["entity_poly_seq"].emplace({ //
			{ "entity_id", 1 },
			{ "num", 1 },
			{ "mon_id", "ALA" },
			{ "hetero", "n" }
		});

		db["pdbx_poly_seq_scheme"].emplace({ //
			{ "asym_id", "A" },
			{ "entity_id", "1" },
			{ "seq_id", "1" },
			{ "mon_id", "ALA" },
			{ "ndb_seq_num", "1" },
			{ "pdb_seq_num", "1" },
			{ "auth_seq_num", "1" },
			{ "pdb_mon_id", "ALA" },
			{ "auth_mon_id", "ALA" },
			{ "pdb_strand_id", "A" },
			{ "pdb_ins_code", "." },
			{ "hetero", "n" }
		});

		REQUIRE_FALSE(cif::pdb::is_valid_pdbx_file(f));
	}

	SECTION("Missing hetero for record in atom_site")
	{
		auto &db = f.front();
		
		auto r1 = db["atom_site"].front();
		cif::row_initializer cr(r1);
		cr.set_value("id", "3");
		cr.set_value("label_comp_id", "ALA");

		db["atom_site"].emplace(std::move(cr));

		REQUIRE_FALSE(cif::pdb::is_valid_pdbx_file(f));
	}

	SECTION("Missing letter in entity_poly.pdbx_seq_one_letter_code")
	{
		auto &db = f.front();
		auto &entity_poly = db["entity_poly"];

		entity_poly.front().assign({
			{ "pdbx_seq_one_letter_code", "PNSG" }
		});

		REQUIRE_FALSE(cif::pdb::is_valid_pdbx_file(f));
	}

	SECTION("Too many letters in entity_poly.pdbx_seq_one_letter_code")
	{
		auto &db = f.front();
		auto &entity_poly = db["entity_poly"];

		entity_poly.front().assign({
			{ "pdbx_seq_one_letter_code", "PNFSGX" }
		});

		REQUIRE_FALSE(cif::pdb::is_valid_pdbx_file(f));
	}

	SECTION("Mismatch in entity_poly.pdbx_seq_one_letter_code")
	{
		auto &db = f.front();
		auto &entity_poly = db["entity_poly"];

		entity_poly.front().assign({
			{ "pdbx_seq_one_letter_code", "PNASG" }
		});

		REQUIRE_FALSE(cif::pdb::is_valid_pdbx_file(f));
	}


}
