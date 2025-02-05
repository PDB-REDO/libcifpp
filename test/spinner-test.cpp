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

#include "test-main.hpp"

#include "cif++/utilities.hpp"

#include <random>
#include <thread>

TEST_CASE("test_one")
{
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(100, 1000);

	cif::progress_bar pb(10, "test");

	for (int i = 0; i < 10; ++i)
	{
		std::this_thread::sleep_for(std::chrono::milliseconds(distrib(gen)));

		pb.message("step " + std::to_string(i));
		pb.consumed(1);
	}
}

TEST_CASE("test_two")
{
	cif::progress_bar pb(10, "test");


	for (int i = 0; i < 5; ++i)
		pb.consumed(1);
}

TEST_CASE("test_three")
{
	using namespace std::literals;

	cif::progress_bar pb(10, "test");
	pb.consumed(10);

	std::this_thread::sleep_for(100ms);
}
