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


#include "test-main.hpp"

#include <stdexcept>

#include <cif++.hpp>

// --------------------------------------------------------------------

TEST_CASE("fmt_1")
{
	std::ostringstream os;

	std::string world("world");
	os << cif::format("Hello, %-10.10s, the magic number is %d and pi is %g", world, 42, cif::kPI);
	REQUIRE(os.str() == "Hello, world     , the magic number is 42 and pi is 3.14159");

	REQUIRE(cif::format("Hello, %-10.10s, the magic number is %d and pi is %g", world, 42, cif::kPI).str() ==
		"Hello, world     , the magic number is 42 and pi is 3.14159");
}

// --------------------------------------------------------------------

TEST_CASE("clr_1")
{
	using namespace cif::colour;

	std::cout << "Hello, " << cif::coloured("world!", white, red, cif::colour::regular) << '\n'
			  << "Hello, " << cif::coloured("world!", white, red, bold) << '\n'
			  << "Hello, " << cif::coloured("world!", black, red) << '\n'
			  << "Hello, " << cif::coloured("world!", white, green) << '\n'
			  << "Hello, " << cif::coloured("world!", white, blue) << '\n'
			  << "Hello, " << cif::coloured("world!", blue, white) << '\n'
			  << "Hello, " << cif::coloured("world!", red, white, bold) << '\n';
}