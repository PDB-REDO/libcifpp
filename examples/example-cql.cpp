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

#include <cif++/cif++.hpp>
#include <filesystem>
#include <iostream>

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		std::cerr << "Usage: example <inputfile>\n";
		exit(1);
	}

	try
	{
		cif::file file(argv[1]);

		if (file.empty())
		{
			std::cerr << "Empty file\n";
			exit(1);
		}

		auto &db = file.front();
		cif::cql::connection c(db);
		cif::cql::transaction tx(c);

		auto N = tx.exec("SELECT COUNT(*) FROM atom_site").one_field().get<std::size_t>();
		auto M = tx.exec("SELECT COUNT(*) FROM atom_site WHERE label_atom_id = 'OXT'").one_field().get<std::size_t>();

		std::cout << "File contains " << N << " atoms of which " << M << (M == 1 ? " is" : " are") << " OXT\n"
				  << "residues with an OXT are:\n";

		for (const auto &[asym, comp, seqnr] : tx.stream<std::string, std::string, int>(
				 "SELECT label_asym_id, label_comp_id, label_seq_id FROM atom_site WHERE label_atom_id = 'OXT'"))
		{
			std::cout << asym << ' ' << comp << ' ' << seqnr << '\n';
		}
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Exception: {}" << ex.what() << '\n';
		return 1;
	}

	return 0;
}
