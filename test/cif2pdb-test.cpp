/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2026 NKI/AVL, Netherlands Cancer Institute
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

#include <catch2/catch_test_macros.hpp>
#include <cif++/pdb.hpp>

#include "../src/pdb/pdb_record.hpp"

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

// --------------------------------------------------------------------

cif::file operator""_cf(const char *text, std::size_t length)
{
	return cif::file(text, length);
}

// cif2pdbDate is an internal helper in src/pdb/cif2pdb.cpp. It is not part of
// the public API, but since it has external linkage we can test it directly.
namespace cif::pdb
{
	std::string cif2pdbDate(const std::string &d);

	std::size_t WriteCitation(std::ostream &pdbFile, const datablock &db, const_row_handle r, int reference);
	void WriteTitle(std::ostream &pdbFile, const datablock &db);
	int WriteHeterogen(std::ostream &pdbFile, const datablock &db);
	void WriteConnectivity(std::ostream &pdbFile, const datablock &db);
	std::tuple<int, int> WriteCoordinate(std::ostream &pdbFile, const datablock &db);
}

TEST_CASE("cif2pdb-date")
{
	CHECK(cif::pdb::cif2pdbDate("2024-05-17") == "17-MAY-24");
	CHECK(cif::pdb::cif2pdbDate("2024-05") == "MAY-24");
	CHECK(cif::pdb::cif2pdbDate("2024-1") == "JAN-24");
	CHECK(cif::pdb::cif2pdbDate("2024-12") == "DEC-24");
	CHECK(cif::pdb::cif2pdbDate("1999-12-31") == "31-DEC-99");
	CHECK(cif::pdb::cif2pdbDate("2000-01-01") == "01-JAN-00");
	CHECK(cif::pdb::cif2pdbDate("2024-02-29") == "29-FEB-24");
	CHECK(cif::pdb::cif2pdbDate("2000-02-29") == "29-FEB-00");

	// Day numbers must be valid for the given month and year.
	CHECK(cif::pdb::cif2pdbDate("2023-02-29").empty());
	CHECK(cif::pdb::cif2pdbDate("1900-02-29").empty());
	CHECK(cif::pdb::cif2pdbDate("2024-02-30").empty());
	CHECK(cif::pdb::cif2pdbDate("2024-04-31").empty());
	CHECK(cif::pdb::cif2pdbDate("2024-05-00").empty());
	CHECK(cif::pdb::cif2pdbDate("2024-05-32").empty());

	// Month numbers outside 1..12 used to index kMonths out of bounds.
	CHECK(cif::pdb::cif2pdbDate("2024-00").empty());
	CHECK(cif::pdb::cif2pdbDate("2024-13").empty());
	CHECK(cif::pdb::cif2pdbDate("2024-99").empty());

	// Anything that is not a valid date results in an empty string.
	CHECK(cif::pdb::cif2pdbDate("").empty());
	CHECK(cif::pdb::cif2pdbDate("not a date").empty());
	CHECK(cif::pdb::cif2pdbDate("2024-05-1").empty());
}

namespace
{

std::unique_ptr<cif::pdb::PDBRecord> make_record(const std::string &name, const std::string &value)
{
	return std::unique_ptr<cif::pdb::PDBRecord>(new (value.length()) cif::pdb::PDBRecord(1, name, value));
}

} // namespace

TEST_CASE("pdb-record-bounds")
{
	// A record name longer than the 10-character field must not overflow it.
	CHECK_THROWS_AS(make_record("TOO-LONG-NAME", "value"), std::invalid_argument);

	auto rec = make_record("REMARK 200", "1234567890");

	// Valid columns map to mValue starting at column 7.
	CHECK(rec->vS(7, 10) == "1234");
	CHECK(rec->vC(7) == '1');
	CHECK(rec->vI(7, 10) == 1234);

	// Columns before the value area and reversed ranges are rejected.
	CHECK_THROWS_AS(rec->vS(1, 6), std::out_of_range);
	CHECK_THROWS_AS(rec->vC(6), std::out_of_range);
	CHECK_THROWS_AS(rec->vI(1, 6), std::out_of_range);
	CHECK_THROWS_AS(rec->vS(10, 7), std::out_of_range);
	CHECK_THROWS_AS(rec->vI(10, 7), std::out_of_range);

	// Reading past the end is clamped, not out of bounds.
	CHECK(rec->vS(7, 100) == "1234567890");
	CHECK(rec->vS(100, 200).empty());
	CHECK(rec->vC(100) == ' ');

	// The error path in vI must also stay within the record when the
	// requested range extends past a truncated line.
	auto truncated = make_record("ATOM  ", "xx");
	CHECK_THROWS_AS(truncated->vI(7, 100), std::runtime_error);
}

TEST_CASE("cif2pdb-citation-reference-newline")
{
	auto f = R"(
data_test
loop_
_citation.id
_citation.title
_citation.journal_abbrev
_citation.journal_volume
_citation.journal_id_ASTM
_citation.country
_citation.journal_id_ISSN
_citation.journal_id_CSD
_citation.book_publisher
_citation.pdbx_database_id_PubMed
_citation.pdbx_database_id_DOI
_citation.page_first
_citation.page_last
_citation.year
1 'A title' . . . . . . . . . . . 2024
)"_cf;

	std::ostringstream os;
	cif::pdb::WriteCitation(os, f.front(), f.front()["citation"].front(), 1);

	CHECK(os.str().starts_with("REMARK   1 REFERENCE 1\n"));
}

TEST_CASE("cif2pdb-revdat-single-line")
{
	auto f = R"(
data_TEST
loop_
_database_PDB_rev.num
_database_PDB_rev.mod_type
_database_PDB_rev.date
_database_PDB_rev.replaces
1 0 2024-05-17 ?
#
loop_
_database_PDB_rev_record.rev_num
_database_PDB_rev_record.type
1 HEADER
1 COMPND
1 SOURCE
1 KEYWDS
)"_cf;

	std::ostringstream os;
	cif::pdb::WriteTitle(os, f.front());

	auto text = os.str();
	auto pos = text.find("REVDAT");
	REQUIRE(pos != std::string::npos);
	auto eol = text.find('\n', pos);
	REQUIRE(eol != std::string::npos);
	auto line = text.substr(pos, eol - pos);

	CHECK(line.find("HEADER") != std::string::npos);
	CHECK(line.find("COMPND") != std::string::npos);
	CHECK(line.find("SOURCE") != std::string::npos);
	CHECK(line.find("KEYWDS") != std::string::npos);
}

TEST_CASE("cif2pdb-formul-long-formula")
{
	cif::datablock db("test");
	db["chem_comp"].emplace(cif::row_initializer{
		{ "id", "LIG" },
		{ "name", "LIGAND" },
		{ "mon_nstd_flag", "N" },
		{ "formula", std::string(51, 'C') + "Z" },
		{ "pdbx_synonyms", nullptr } });

	std::ostringstream os;
	cif::pdb::WriteHeterogen(os, db);

	auto text = os.str();
	CHECK(text.find("FORMUL") != std::string::npos);
	CHECK(text.find('Z') != std::string::npos);
}

TEST_CASE("cif2pdb-het-model-count")
{
	cif::datablock db("test");
	db["chem_comp"].emplace(cif::row_initializer{
		{ "id", "LIG" },
		{ "name", "LIGAND" },
		{ "mon_nstd_flag", "N" },
		{ "formula", "C" },
		{ "pdbx_synonyms", nullptr } });

	auto &atom_site = db["atom_site"];
	for (int model : { 1, 2 })
	{
		atom_site.emplace(cif::row_initializer{
			{ "group_PDB", "HETATM" },
			{ "id", std::to_string(model) },
			{ "label_entity_id", "2" },
			{ "label_comp_id", "LIG" },
			{ "auth_comp_id", "LIG" },
			{ "auth_asym_id", "A" },
			{ "auth_seq_id", 1 },
			{ "pdbx_PDB_ins_code", nullptr },
			{ "label_atom_id", "C1" },
			{ "type_symbol", "C" },
			{ "pdbx_PDB_model_num", model } });
	}

	std::ostringstream os;
	cif::pdb::WriteHeterogen(os, db);

	auto text = os.str();
	auto pos = text.find("HET    ");
	REQUIRE(pos != std::string::npos);

	auto eol = text.find('\n', pos);
	REQUIRE(eol != std::string::npos);

	auto line = text.substr(pos, eol - pos);
	CHECK(line.ends_with("    1"));
	CHECK(text.find("HET    ", pos + 1) == std::string::npos);
}

TEST_CASE("cif2pdb-het-insertion-code")
{
	cif::datablock db("test");
	db["chem_comp"].emplace(cif::row_initializer{
		{ "id", "LIG" },
		{ "name", "LIGAND" },
		{ "mon_nstd_flag", "N" },
		{ "formula", "C" },
		{ "pdbx_synonyms", nullptr } });

	auto &atom_site = db["atom_site"];
	for (const auto &icode : { std::string{}, std::string{ "A" } })
	{
		atom_site.emplace(cif::row_initializer{
			{ "group_PDB", "HETATM" },
			{ "id", icode.empty() ? "1" : "2" },
			{ "label_entity_id", "2" },
			{ "label_comp_id", "LIG" },
			{ "auth_comp_id", "LIG" },
			{ "auth_asym_id", "A" },
			{ "auth_seq_id", 1 },
			{ "pdbx_PDB_ins_code", icode.empty() ? cif::item_value_type::MISSING : cif::item_value(icode) },
			{ "label_atom_id", "C1" },
			{ "type_symbol", "C" },
			{ "pdbx_PDB_model_num", 1 } });
	}

	std::ostringstream os;
	cif::pdb::WriteHeterogen(os, db);

	auto text = os.str();
	std::size_t count = 0;
	for (auto pos = text.find("HET    "); pos != std::string::npos; pos = text.find("HET    ", pos + 1))
		++count;

	CHECK(count == 2);
}

TEST_CASE("cif2pdb-link-distance")
{
	cif::datablock db("test");

	db["struct_conn"].emplace(cif::row_initializer{
		{ "conn_type_id", "covale" },
		{ "ptnr1_label_atom_id", "CA" },
		{ "pdbx_ptnr1_label_alt_id", nullptr },
		{ "ptnr1_label_comp_id", "ALA" },
		{ "ptnr1_auth_asym_id", "A" },
		{ "ptnr1_auth_seq_id", 1 },
		{ "pdbx_ptnr1_PDB_ins_code", nullptr },
		{ "ptnr2_label_atom_id", "CA" },
		{ "pdbx_ptnr2_label_alt_id", nullptr },
		{ "ptnr2_label_comp_id", "ALA" },
		{ "ptnr2_auth_asym_id", "A" },
		{ "ptnr2_auth_seq_id", 2 },
		{ "pdbx_ptnr2_PDB_ins_code", nullptr },
		{ "ptnr1_symmetry", "1_555" },
		{ "ptnr2_symmetry", "1_555" },
		{ "pdbx_dist_value", "1.23" } });

	std::ostringstream os;
	cif::pdb::WriteConnectivity(os, db);

	auto text = os.str();
	auto pos = text.find("LINK");
	REQUIRE(pos != std::string::npos);

	auto eol = text.find('\n', pos);
	REQUIRE(eol != std::string::npos);

	CHECK(text.substr(pos, eol - pos).find("1.23") != std::string::npos);
}

TEST_CASE("cif2pdb-model-two-coordinate-count")
{
	cif::datablock db("test");

	db["atom_site"].emplace(cif::row_initializer{
		{ "group_PDB", "ATOM" },
		{ "id", "1" },
		{ "label_atom_id", "CA" },
		{ "label_alt_id", nullptr },
		{ "auth_comp_id", "ALA" },
		{ "auth_asym_id", "A" },
		{ "auth_seq_id", 1 },
		{ "pdbx_PDB_ins_code", nullptr },
		{ "Cartn_x", 1.0 },
		{ "Cartn_y", 2.0 },
		{ "Cartn_z", 3.0 },
		{ "occupancy", 1.0 },
		{ "B_iso_or_equiv", 20.0 },
		{ "type_symbol", "C" },
		{ "pdbx_formal_charge", 0 },
		{ "label_comp_id", "ALA" },
		{ "label_entity_id", "1" },
		{ "pdbx_PDB_model_num", 2 } });

	std::ostringstream os;
	auto [numCoord, numTer] = cif::pdb::WriteCoordinate(os, db);

	CHECK(numCoord == 1);
	CHECK(os.str().find("ATOM") != std::string::npos);
}

TEST_CASE("pdb-remark3-truncated-at-eof")
{
	std::istringstream is(
		"HEADER    TEST\n"
		"ATOM      1  N   PRO A   1      16.979  13.301  44.555  1.00 30.05           N\n"
		"REMARK   3   REFINEMENT.\n"
		"REMARK   3   PROGRAM     : X-PLOR\n"
		"REMARK   3   DATA USED IN REFINEMENT.\n"
		"REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 1.80\n");

	CHECK_NOTHROW(cif::pdb::read(is));
}

TEST_CASE("pdb-remark3-wrapped-value")
{
	std::string pdb =
		"HEADER    TEST\n"
		"ATOM      1  N   PRO A   1      16.979  13.301  44.555  1.00 30.05           N\n"
		"REMARK   3   REFINEMENT.\n"
		"REMARK   3   PROGRAM     : X-PLOR\n"
		"REMARK   3   DATA USED IN REFINEMENT.\n"
		"REMARK   3   FREE R VALUE TEST SET SELECTION : RANDOM\n";

	pdb += "REMARK   3" + std::string(35, ' ') + "SELECTION\n";
	pdb += "REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 8.00\n";
	pdb += "ATOM      2  CA  PRO A   1      18.150  13.525  43.680  1.00 28.82           C\n";

	std::istringstream is(pdb);
	auto f = cif::pdb::read(is);
	auto &refine = f.front()["refine"];

	REQUIRE(refine.size() == 1);
	CHECK(refine.front()["pdbx_R_Free_selection_details"].str() == "RANDOM");
}
