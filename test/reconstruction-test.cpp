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

#include "cif++/utilities.hpp"
#include "test-main.hpp"

#include <cif++/cif++.hpp>

#include <filesystem>
#include <iostream>
#include <fstream>

TEST_CASE("reconstruct")
{
	cif::VERBOSE = 1;

	cif::compound_factory::instance().push_dictionary(gTestDir / "REA.cif");

	for (std::filesystem::directory_iterator i(gTestDir / "reconstruct"); i != std::filesystem::directory_iterator{}; ++i)
	{
		std::cout << i->path() << '\n';

		if (i->path().extension() == ".pdb")
		{
			cif::file f = cif::pdb::read(i->path());

			std::error_code ec;

			if (not cif::pdb::is_valid_pdbx_file(f, ec))
				CHECK(cif::pdb::reconstruct_pdbx(f));
		}
		else
		{
			cif::file f(i->path());

			std::error_code ec;
			CHECK_FALSE(cif::pdb::is_valid_pdbx_file(f, ec));
			CHECK(ec != std::errc{});

			auto valid = cif::pdb::reconstruct_pdbx(f);

			CHECK(valid);

			if (not valid)
			{
				std::ofstream of(std::filesystem::temp_directory_path() / i->path().filename());
				of << f;
				of.close();
			}
		}
	}
}