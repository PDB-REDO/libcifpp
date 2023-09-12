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

#include <cif++.hpp>

#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{
	cif::VERBOSE = 3;

	try
	{
		std::filesystem::path testdir = std::filesystem::current_path();

		if (argc == 3)
			testdir = argv[2];
		else
		{
			while (not testdir.empty() and not std::filesystem::exists(testdir / "test"))
				testdir = testdir.parent_path();
			testdir /= "test";
		}

		if (std::filesystem::exists(testdir / ".." / "data" / "ccd-subset.cif"))
			cif::add_file_resource("components.cif", testdir / ".." / "data" / "ccd-subset.cif");

		if (std::filesystem::exists(testdir / ".." / "rsrc" / "mmcif_pdbx.dic"))
			cif::add_file_resource("mmcif_pdbx.dic", testdir / ".." / "rsrc" / "mmcif_pdbx.dic");

		cif::compound_factory::instance().push_dictionary(testdir / "REA.cif");
		cif::compound_factory::instance().push_dictionary(testdir / "RXA.cif");

		cif::file f(testdir / ".."/"examples"/"1cbs.cif.gz");
		cif::mm::structure structure(f);

		auto &res = structure.get_residue("B");
		structure.change_residue(res, "RXA", {});

		structure.cleanup_empty_categories();

		f.save(std::cout);

		if (not f.is_valid())
			throw std::runtime_error("Invalid");

		f.save(std::cout);
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		exit(1);
	}
	
	return 0;	
}
