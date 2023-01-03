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

#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include <stdexcept>

#include <cif++.hpp>

#include <cif++/dictionary_parser.hpp>
#include <cif++/parser.hpp>

namespace tt = boost::test_tools;

std::filesystem::path gTestDir = std::filesystem::current_path(); // filled in first test

// --------------------------------------------------------------------

cif::file operator""_cf(const char *text, size_t length)
{
	struct membuf : public std::streambuf
	{
		membuf(char *text, size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(text), length);

	std::istream is(&buffer);
	return cif::file(is);
}

// --------------------------------------------------------------------

bool init_unit_test()
{
	cif::VERBOSE = 1;

	// not a test, just initialize test dir
	if (boost::unit_test::framework::master_test_suite().argc == 2)
		gTestDir = boost::unit_test::framework::master_test_suite().argv[1];

	// do this now, avoids the need for installing
	cif::add_file_resource("mmcif_pdbx.dic", gTestDir / ".." / "rsrc" / "mmcif_pdbx.dic");

	// initialize CCD location
	cif::add_file_resource("components.cif", gTestDir / ".." / "data" / "ccd-subset.cif");

	return true;
}

// --------------------------------------------------------------------

// 3d tests

BOOST_AUTO_TEST_CASE(t1)
{
	// std::random_device rnd;
	// std::mt19937 gen(rnd());
	// std::uniform_real_distribution<float> dis(0, 1);

	// Quaternion q{ dis(gen), dis(gen), dis(gen), dis(gen) };
	// q = Normalize(q);

	// Quaternion q{ 0.1, 0.2, 0.3, 0.4 };
	cif::quaternion q{ 0.5, 0.5, 0.5, 0.5 };
	q = normalize(q);

	const auto &&[angle0, axis0] = cif::quaternion_to_angle_axis(q);

	std::vector<cif::point> p1{
		{ 16.979, 13.301, 44.555 },
		{ 18.150, 13.525, 43.680 },
		{ 18.656, 14.966, 43.784 },
		{ 17.890, 15.889, 44.078 },
		{ 17.678, 13.270, 42.255 },
		{ 16.248, 13.734, 42.347 },
		{ 15.762, 13.216, 43.724 }
	};

	auto p2 = p1;

	cif::center_points(p1);

	for (auto &p : p2)
		p.rotate(q);

	cif::center_points(p2);

	auto q2 = cif::align_points(p1, p2);

	const auto &&[angle, axis] = cif::quaternion_to_angle_axis(q2);

	BOOST_TEST(std::fmod(360 + angle, 360) == std::fmod(360 - angle0, 360), tt::tolerance(0.01));

	for (auto &p : p1)
		p.rotate(q2);

	float rmsd = cif::RMSd(p1, p2);

	BOOST_TEST(rmsd < 1e-5);

	// std::cout << "rmsd: " << RMSd(p1, p2) << std::endl;
}

BOOST_AUTO_TEST_CASE(t2)
{
	cif::point p[] = {
		{ 1, 1, 0 },
		{ 2, 1, 0 },
		{ 1, 2, 0 }
	};

	cif::point xp = cif::cross_product(p[1] - p[0], p[2] - p[0]);

	auto q = cif::construct_from_angle_axis(45, xp); // mmcif::Normalize(Quaternion{45 * mmcif::kPI / 180, xp.mX, xp.mY, xp.mZ});

	auto &&[angle, axis] = cif::quaternion_to_angle_axis(q);

	BOOST_TEST(angle == 45, tt::tolerance(0.01));
}

BOOST_AUTO_TEST_CASE(t3)
{
	cif::point p[] = {
		{ 1, 1, 0 },
		{ 2, 1, 0 },
		{ 1, 2, 0 }
	};

	cif::point xp = cif::cross_product(p[1] - p[0], p[2] - p[0]);

	auto q = cif::construct_from_angle_axis(45, xp); // mmcif::Normalize(Quaternion{45 * mmcif::kPI / 180, xp.mX, xp.mY, xp.mZ});

	auto v = p[1];
	v -= p[0];
	v.rotate(q);
	v += p[0];

	std::cout << v << std::endl;

	double a = cif::angle(v, p[0], p[1]);

	BOOST_TEST(a == 45, tt::tolerance(0.01));
}

BOOST_AUTO_TEST_CASE(dh_q_1)
{
	struct
	{
		float angle;
		cif::point pts[4];
	} tests[] = {
		{ -97.5,
			{ { 68.8649979, -7.34800005, 54.3769989 },
				{ 68.1350021, -8.18700027, 53.6489983 },
				{ 68.7760239, -9.07335377, 52.7140236 },
				{ 68.9000015, -10.3944235, 53.2217026 } } },
		{ 80.3,
			{ { 0.304512024, 0.531184196, 2.25860214 },
				{ 0.956512451, 0.0321846008, 1.07460022 },
				{ 0, 0, 0 },
				{ 0.21336633, -1.09552193, -0.878999829 } } },
		{ -97.5,
			{ { 0.088973999, 1.72535372, 1.66297531 },
				{ -0.641021729, 0.886353493, 0.93497467 },
				{ 0, 0, 0 },
				{ 1.29433727, -0.395142615, 0.432300746 } } },
		{ -97.5,
			{
				{ 0.088973999, 1.72535372, 1.66297531 },
				{ -0.641021729, 0.886353493, 0.93497467 },
				{ 0, 0, 0 },
				{ 1.33983064, 0.384027064, -0.275154471 },

			} }
	};

	for (auto &&[angle, pts] : tests)
	{
		auto q = cif::construct_for_dihedral_angle(pts[0], pts[1], pts[2], pts[3], angle, 1);

		pts[3] -= pts[2];
		pts[3].rotate(q);
		pts[3] += pts[2];

		auto dh = cif::dihedral_angle(pts[0], pts[1], pts[2], pts[3]);
		BOOST_TEST(dh == angle, tt::tolerance(0.1f));
	}
}