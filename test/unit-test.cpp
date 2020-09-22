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

#define BOOST_TEST_MODULE LibCifPP_Test
#include <boost/test/included/unit_test.hpp>

#include "cif++/DistanceMap.hpp"

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

BOOST_AUTO_TEST_CASE(ut1)
{
	using namespace mmcif;

	auto f = R"(data_TEST
#
loop_
_test.id
_test.name
1 aap
2 noot
3 mies
	)"_cf;

	auto& db = f.firstDatablock();

	BOOST_CHECK(db.getName() == "TEST");
	
	auto& test = db["test"];
	BOOST_CHECK(test.size() == 3);

	// wrong! the next lines will crash. And that's OK, don't do that
	// for (auto r: test)
	// 	test.erase(r);
	
	// BOOST_CHECK(test.empty());

	// test.purge();

	auto n = test.erase(cif::Key("id") == 1, [](const cif::Row& r) {
		BOOST_CHECK_EQUAL(r["id"].as<int>(), 1);
		BOOST_CHECK_EQUAL(r["name"].as<std::string>(), "aap");
  	});

	BOOST_CHECK_EQUAL(n, 1);
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(ut2)
{
	using namespace mmcif;

	auto f = R"(data_TEST
#
loop_
_test.id
_test.name
1 aap
2 noot
3 mies
	)"_cf;

	auto& db = f.firstDatablock();

	BOOST_CHECK(db.getName() == "TEST");
	
	auto& test = db["test"];
	BOOST_CHECK(test.size() == 3);

	int n = 0;
	for (auto r: test.find(cif::Key("name") == "aap"))
	{
		BOOST_CHECK(++n == 1);
		BOOST_CHECK(r["id"].as<int>() == 1);
		BOOST_CHECK(r["name"].as<std::string>() == "aap");
	}

	auto t = test.find(cif::Key("id") == 1);
	BOOST_CHECK(not t.empty());
	BOOST_CHECK(t.front()["name"].as<std::string>() == "aap");
}