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

#include <string>

// cif2pdbDate is an internal helper in src/pdb/cif2pdb.cpp. It is not part of
// the public API, but since it has external linkage we can test it directly.
namespace cif::pdb
{
	std::string cif2pdbDate(const std::string &d);
}

TEST_CASE("cif2pdb-date")
{
	CHECK(cif::pdb::cif2pdbDate("2024-05-17") == "17-MAY-24");
	CHECK(cif::pdb::cif2pdbDate("2024-05") == "MAY-24");
	CHECK(cif::pdb::cif2pdbDate("2024-1") == "JAN-24");
	CHECK(cif::pdb::cif2pdbDate("2024-12") == "DEC-24");
	CHECK(cif::pdb::cif2pdbDate("1999-12-31") == "31-DEC-99");
	CHECK(cif::pdb::cif2pdbDate("2000-01-01") == "01-JAN-00");

	// Month numbers outside 1..12 used to index kMonths out of bounds.
	CHECK(cif::pdb::cif2pdbDate("2024-00").empty());
	CHECK(cif::pdb::cif2pdbDate("2024-13").empty());
	CHECK(cif::pdb::cif2pdbDate("2024-99").empty());

	// Anything that is not a valid date results in an empty string.
	CHECK(cif::pdb::cif2pdbDate("").empty());
	CHECK(cif::pdb::cif2pdbDate("not a date").empty());
	CHECK(cif::pdb::cif2pdbDate("2024-05-1").empty());
}
