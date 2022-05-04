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

#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include <stdexcept>

#include "cif++/Cif++.hpp"
#include "cif++/Structure.hpp"
#include "cif++/CifValidator.hpp"

// --------------------------------------------------------------------

cif::File operator""_cf(const char* text, size_t length)
{
    struct membuf : public std::streambuf
    {
        membuf(char* text, size_t length)
        {
            this->setg(text, text, text + length);
        }
    } buffer(const_cast<char*>(text), length);

    std::istream is(&buffer);
    return cif::File(is);
}

// --------------------------------------------------------------------

std::filesystem::path gTestDir = std::filesystem::current_path();

bool init_unit_test()
{
    cif::VERBOSE = 1;

	// not a test, just initialize test dir
	if (boost::unit_test::framework::master_test_suite().argc == 2)
		gTestDir = boost::unit_test::framework::master_test_suite().argv[1];

	// do this now, avoids the need for installing
	cif::addFileResource("mmcif_pdbx_v50.dic", gTestDir / ".." / "rsrc" / "mmcif_pdbx_v50.dic");

	// initialize CCD location
	cif::addFileResource("components.cif", gTestDir / ".." / "data" / "ccd-subset.cif");

	mmcif::CompoundFactory::instance().pushDictionary(gTestDir / "HEM.cif");

	return true;
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(sugar_name_1)
{
	using namespace cif::literals;

	const std::filesystem::path example(gTestDir / "1juh.cif.gz");
	mmcif::File file(example.string());
	mmcif::Structure s(file);

	auto &db = s.datablock();
	auto &entity = db["entity"];

	auto &branches = s.branches();

	BOOST_CHECK_EQUAL(branches.size(), 4);

	for (auto &branch : branches)
	{
		auto entityID = branch.front().entityID();

		auto name = entity.find1<std::string>("id"_key == entityID, "pdbx_description");
		BOOST_CHECK_EQUAL(branch.name(), name);
	}
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(create_sugar_1)
{
	using namespace cif::literals;

	const std::filesystem::path example(gTestDir / "1juh.cif.gz");
	mmcif::File file(example.string());
	mmcif::Structure s(file);

	// collect atoms from asym L first
	auto &NAG = s.getResidue("L");
	auto nagAtoms = NAG.atoms();

	std::vector<std::vector<cif::Item>> ai;

	auto &db = s.datablock();
	auto &as = db["atom_site"];

	for (auto r : as.find("label_asym_id"_key == "L"))
		ai.emplace_back(r.begin(), r.end());

	s.removeResidue(NAG);

	auto &branch = s.createBranch(ai);

	BOOST_CHECK_EQUAL(branch.name(), "2-acetamido-2-deoxy-beta-D-glucopyranose");
	BOOST_CHECK_EQUAL(branch.size(), 1);

	BOOST_CHECK_EQUAL(branch[0].atoms().size(), nagAtoms.size());
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(create_sugar_2)
{
	using namespace cif::literals;

	const std::filesystem::path example(gTestDir / "1juh.cif.gz");
	mmcif::File file(example.string());
	mmcif::Structure s(file);

	// Get branch for H
	auto &bH = s.getBranchByAsymID("H");

	BOOST_CHECK_EQUAL(bH.size(), 2);

	std::vector<std::vector<cif::Item>> ai[2];

	auto &db = s.datablock();
	auto &as = db["atom_site"];

	for (size_t i = 0; i < 2; ++i)
	{
		for (auto r : as.find("label_asym_id"_key == "H" and "auth_seq_id"_key == i + 1))
			ai[i].emplace_back(r.begin(), r.end());
	}

	s.removeBranch(bH);

	auto &bN = s.createBranch(ai[0]);
	s.extendBranch(bN.asymID(), ai[1], 1, "O4");

	BOOST_CHECK_EQUAL(bN.name(), "2-acetamido-2-deoxy-beta-D-glucopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose");
	BOOST_CHECK_EQUAL(bN.size(), 2);

	file.save(gTestDir / "test-create_sugar_2.cif");

	BOOST_CHECK_NO_THROW(mmcif::Structure s2(file));
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(delete_sugar_1)
{
	using namespace cif::literals;

	const std::filesystem::path example(gTestDir / "1juh.cif.gz");
	mmcif::File file(example.string());
	mmcif::Structure s(file);

	// Get branch for H
	auto &bG = s.getBranchByAsymID("G");

	BOOST_CHECK_EQUAL(bG.size(), 4);

	s.removeResidue(bG[1]);

	BOOST_CHECK_EQUAL(bG.size(), 1);

	auto &bN = s.getBranchByAsymID("G");

	BOOST_CHECK_EQUAL(bN.name(), "2-acetamido-2-deoxy-beta-D-glucopyranose");
	BOOST_CHECK_EQUAL(bN.size(), 1);

	file.save(gTestDir / "test-create_sugar_3.cif");

	BOOST_CHECK_NO_THROW(mmcif::Structure s2(file));
}
