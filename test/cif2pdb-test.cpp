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
#include <stdexcept>
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
