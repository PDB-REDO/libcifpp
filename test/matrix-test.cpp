/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2025 NKI/AVL, Netherlands Cancer Institute
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

#include "cif++/matrix.hpp"
#include "test-main.hpp"

#include <catch2/catch_test_macros.hpp>
#include <cif++/cif++.hpp>

TEST_CASE("m1")
{
	cif::matrix3x3<int> m = cif::identity_matrix<int>(3);

	CHECK(cif::determinant(m) == 1);
}

TEST_CASE("m2")
{
	cif::matrix4x4<int> m = cif::identity_matrix<int>(4);

	cif::sub_matrix<cif::matrix4x4<int>> ms(m, 1, 1);
	CHECK(ms == cif::identity_matrix<int>(3));
}

TEST_CASE("m3")
{
	cif::matrix4x4<int> m{
		{ 1, 2, 3, 4,      //
			5, 6, 7, 8,    //
			9, 10, 11, 12, //
			13, 14, 15, 16 }
	};
	cif::sub_matrix<cif::matrix4x4<int>> ms(m, 1, 1);

	cif::matrix3x3<int> t{
		{ 1, 3, 4, 9, 11, 12, 13, 15, 16 }
	};

	CHECK(ms == t);
}

TEST_CASE("m4")
{
	cif::matrix4x4<int> m{
		{
			-2,
			3,
			1,
			0,
			4,
			1,
			-3,
			2,
			0,
			-1,
			2,
			5,
			3,
			2,
			0,
			-4,
		}
	};

	std::cout << m << "\n\n";

	// std::cout << cif::matrix3x3<int>(cif::sub_matrix<decltype(m)>(m, 0, 0)) << "\n\n";
	// std::cout << cif::matrix3x3<int>(cif::sub_matrix<decltype(m)>(m, 0, 1)) << "\n\n";
	// std::cout << cif::matrix3x3<int>(cif::sub_matrix<decltype(m)>(m, 0, 2)) << "\n\n";
	// std::cout << cif::matrix3x3<int>(cif::sub_matrix<decltype(m)>(m, 0, 3)) << "\n\n";

	// std::cout << cif::determinant(cif::matrix3x3<int>(cif::sub_matrix<decltype(m)>(m, 0, 0))) << "\n\n";
	// std::cout << cif::determinant(cif::matrix3x3<int>(cif::sub_matrix<decltype(m)>(m, 0, 1))) << "\n\n";
	// std::cout << cif::determinant(cif::matrix3x3<int>(cif::sub_matrix<decltype(m)>(m, 0, 2))) << "\n\n";
	// std::cout << cif::determinant(cif::matrix3x3<int>(cif::sub_matrix<decltype(m)>(m, 0, 3))) << "\n\n";

	CHECK(cif::determinant(m) == 332);
}

// --------------------------------------------------------------------

TEST_CASE("m5")
{
    cif::matrix4x4<float> m = cif::identity_matrix<float>(4);
    cif::matrix_fixed<float, 1, 4> v({ 0, 0.5f, 0, 1.0f });

    auto mv = v * m;

    CHECK(mv == v);

}

TEST_CASE("m6")
{
	cif::matrix_fixed<float, 1, 2> a({ 1, 2 });
	cif::matrix_fixed<float, 2, 1> b({ 1, 2 });

	auto c = a * b;

	CHECK(c.dim_m() == 1);
	CHECK(c.dim_n() == 1);
	CHECK(c(0, 0) == 5);

	auto d = b * a;

	CHECK(d.dim_m() == 2);
	CHECK(d.dim_n() == 2);

	CHECK(d(0, 0) == 1);
	CHECK(d(0, 1) == 2);
	CHECK(d(1, 0) == 2);
	CHECK(d(1, 1) == 4);
}