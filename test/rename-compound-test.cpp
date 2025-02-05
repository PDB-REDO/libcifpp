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

#include <cif++.hpp>

#include <iostream>
#include <fstream>

TEST_CASE("rename")
{
	cif::VERBOSE = 3;

	if (std::filesystem::exists(gTestDir / ".." / "rsrc" / "ccd-subset.cif"))
		cif::add_file_resource("components.cif", gTestDir / ".." / "rsrc" / "ccd-subset.cif");

	if (std::filesystem::exists(gTestDir / ".." / "rsrc" / "mmcif_pdbx.dic"))
		cif::add_file_resource("mmcif_pdbx.dic", gTestDir / ".." / "rsrc" / "mmcif_pdbx.dic");

	cif::compound_factory::instance().push_dictionary(gTestDir / "REA.cif");
	cif::compound_factory::instance().push_dictionary(gTestDir / "RXA.cif");

	cif::file f(gTestDir / ".."/"examples"/"1cbs.cif.gz");
	cif::mm::structure structure(f);

	auto &res = structure.get_residue("B");
	structure.change_residue(res, "RXA", {});

	structure.cleanup_empty_categories();

	f.save(std::cout);

	if (not f.is_valid())
		throw std::runtime_error("Invalid");

	f.save(std::cout);
}
