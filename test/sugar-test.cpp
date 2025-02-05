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

cif::file operator""_cf(const char* text, std::size_t length)
{
    struct membuf : public std::streambuf
    {
        membuf(char* text, std::size_t length)
        {
            this->setg(text, text, text + length);
        }
    } buffer(const_cast<char*>(text), length);

    std::istream is(&buffer);
    return cif::file(is);
}

// --------------------------------------------------------------------

TEST_CASE("sugar_name_1")
{
	using namespace cif::literals;

	const std::filesystem::path example(gTestDir / "1juh.cif.gz");
	cif::file file(example.string());
	cif::mm::structure s(file);

	auto &db = s.get_datablock();
	auto &entity = db["entity"];

	auto &branches = s.branches();

	REQUIRE(branches.size() == 4UL);

	for (auto &branch : branches)
	{
		auto entityID = branch.front().get_entity_id();

		auto name = entity.find1<std::string>("id"_key == entityID, "pdbx_description");
		REQUIRE(branch.name() == name);
	}
}

// // --------------------------------------------------------------------

// TEST_CASE("create_sugar_1")
// {
// 	using namespace cif::literals;

// 	const std::filesystem::path example(gTestDir / "1juh.cif.gz");
// 	cif::file file(example.string());
// 	cif::mm::structure s(file);

// 	// collect atoms from asym L first
// 	auto &NAG = s.get_residue("L");
// 	auto nagAtoms = NAG.atoms();

// 	std::vector<cif::row_initializer> ai;

// 	auto &db = s.get_datablock();
// 	auto &as = db["atom_site"];

// 	// NOTE, row_initializer does not actually hold the data, so copy it first
// 	// before it gets destroyed by remove_residue

// 	for (auto r : as.find("label_asym_id"_key == "L"))
// 		/*auto &ri = */ai.emplace_back(r);

// 	s.remove_residue(NAG);

// 	auto &branch = s.create_branch(ai);

// 	REQUIRE(branch.name() == "2-acetamido-2-deoxy-beta-D-glucopyranose");
// 	REQUIRE(branch.size() == 1);

// 	REQUIRE(branch[0].atoms().size() == nagAtoms.size());

// 	REQUIRE(file.is_valid());

// 	file.save(gTestDir / "test-create_sugar_1.cif");
// }

// // --------------------------------------------------------------------

// TEST_CASE("create_sugar_2")
// {
// 	using namespace cif::literals;

// 	const std::filesystem::path example(gTestDir / "1juh.cif.gz");
// 	cif::file file(example.string());
// 	cif::mm::structure s(file);

// 	// Get branch for H
// 	auto &bH = s.get_branch_by_asym_id("H");

// 	REQUIRE(bH.size() == 2);

// 	std::vector<cif::row_initializer> ai[2];

// 	auto &db = s.get_datablock();
// 	auto &as = db["atom_site"];

// 	for (std::size_t i = 0; i < 2; ++i)
// 	{
// 		for (auto r : as.find("label_asym_id"_key == "H" and "auth_seq_id"_key == i + 1))
// 			/*auto &ri = */ai[i].emplace_back(r);
// 	}

// 	s.remove_branch(bH);

// 	REQUIRE(file.is_valid());

// 	auto &bN = s.create_branch(ai[0]);
// 	s.extend_branch(bN.get_asym_id(), ai[1], 1, "O4");

// 	REQUIRE(bN.name() == "2-acetamido-2-deoxy-beta-D-glucopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose");
// 	REQUIRE(bN.size() == 2);

// 	REQUIRE(file.is_valid());

// 	file.save(gTestDir / "test-create_sugar_2.cif");

// 	REQUIRE_NO_THROW(cif::mm::structure s2(file));
// }

// --------------------------------------------------------------------

TEST_CASE("delete_sugar_1")
{
	using namespace cif::literals;

	const std::filesystem::path example(gTestDir / "1juh.cif.gz");
	cif::file file(example.string());
	cif::mm::structure s(file);

	// Get branch for H
	auto &bG = s.get_branch_by_asym_id("G");

	REQUIRE(bG.size() == 4UL);

	s.remove_residue(bG[1]);

	REQUIRE(bG.size() == 1UL);

	auto &bN = s.get_branch_by_asym_id("G");

	REQUIRE(bN.name() == "2-acetamido-2-deoxy-beta-D-glucopyranose");
	REQUIRE(bN.size() == 1UL);

	REQUIRE(file.is_valid());

	// file.save(gTestDir / "test-create_sugar_3.cif");

	cif::mm::structure s2(file);
}
