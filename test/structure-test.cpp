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

#define BOOST_TEST_MODULE Structure_Test
#include <boost/test/included/unit_test.hpp>

#include <stdexcept>

#include "cif++/Cif++.hpp"
#include "cif++/Structure.hpp"

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

BOOST_AUTO_TEST_CASE(create_nonpoly_1)
{
    cif::VERBOSE = 1;

	// do this now, avoids the need for installing
	cif::addFileResource("mmcif_pdbx_v50.dic", "../rsrc/mmcif_pdbx_v50.dic");

	auto expected = R"(
data_TEST
loop_
_entity.id
_entity.type
_entity.src_method
_entity.pdbx_description
_entity.formula_weight
1 non-polymer syn 'PROTOPORPHYRIN IX CONTAINING FE' 616.487
	)"_cf;

	expected.loadDictionary("mmcif_pdbx_v50.dic");

	mmcif::File file;
	file.file().loadDictionary("mmcif_pdbx_v50.dic");
	file.createDatablock("TEST");	// create a datablock
	
	mmcif::Structure structure(file);

	structure.createEntityNonPoly({
		{ "src_method", "syn" },
		{ "pdbx_description", "PROTOPORPHYRIN IX CONTAINING FE" },
		{ "formula_weight", 616.487 }
	}, "HEM" );

	BOOST_TEST(expected.firstDatablock() == structure.getFile().data());

	std::cout << expected.firstDatablock() << std::endl
			  << std::endl
			  << structure.getFile().data() << std::endl;

//     // using namespace mmcif;

//     auto f = R"(data_TEST
// #
// loop_
// _test.id
// _test.name
// 1 aap
// 2 noot
// 3 mies
//     )"_cf;

//     auto& db = f.firstDatablock();

//     BOOST_CHECK(db.getName() == "TEST");
    
//     auto& test = db["test"];
//     BOOST_CHECK(test.size() == 3);

//     // wrong! the next lines will crash. And that's OK, don't do that
//     // for (auto r: test)
//     // 	test.erase(r);
    
//     // BOOST_CHECK(test.empty());

//     // test.purge();

//     auto n = test.erase(cif::Key("id") == 1, [](const cif::Row& r) {
//         BOOST_CHECK_EQUAL(r["id"].as<int>(), 1);
//         BOOST_CHECK_EQUAL(r["name"].as<std::string>(), "aap");
//       });

//     BOOST_CHECK_EQUAL(n, 1);
}
