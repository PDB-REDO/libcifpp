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

namespace tt = boost::test_tools;
namespace utf = boost::unit_test;

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
		{ 16.979f, 13.301f, 44.555f },
		{ 18.150f, 13.525f, 43.680f },
		{ 18.656f, 14.966f, 43.784f },
		{ 17.890f, 15.889f, 44.078f },
		{ 17.678f, 13.270f, 42.255f },
		{ 16.248f, 13.734f, 42.347f },
		{ 15.762f, 13.216f, 43.724f }
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

	auto rmsd = cif::RMSd(p1, p2);

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

BOOST_AUTO_TEST_CASE(dh_q_0)
{
	cif::point axis(1, 0, 0);

	cif::point p(1, 1, 0);
	
	cif::point t[3] =
	{
		{ 0, 1, 0 },
		{ 0, 0, 0 },
		{ 1, 0, 0 }
	};

	auto a = cif::dihedral_angle(t[0], t[1], t[2], p);
	BOOST_TEST(a == 0, tt::tolerance(0.01f));

	auto q = cif::construct_from_angle_axis(90, axis);

	p.rotate(q);

	BOOST_TEST(p.m_x == 1, tt::tolerance(0.01f));
	BOOST_TEST(p.m_y == 0, tt::tolerance(0.01f));
	BOOST_TEST(p.m_z == 1, tt::tolerance(0.01f));

	a = cif::dihedral_angle(t[0], t[1], t[2], p);
	BOOST_TEST(a == 90, tt::tolerance(0.01f));

	q = cif::construct_from_angle_axis(-90, axis);

	p.rotate(q);

	BOOST_TEST(p.m_x == 1, tt::tolerance(0.01f));
	BOOST_TEST(p.m_y == 1, tt::tolerance(0.01f));
	BOOST_TEST(p.m_z == 0, tt::tolerance(0.01f));

	a = cif::dihedral_angle(t[0], t[1], t[2], p);
	BOOST_TEST(a == 0, tt::tolerance(0.01f));

}

BOOST_AUTO_TEST_CASE(dh_q_1)
{
	struct
	{
		float angle;
		cif::point pts[4];
	} tests[] = {
		{ -97.5f,
			{ { 68.8649979f, -7.34800005f, 54.3769989f },
				{ 68.1350021f, -8.18700027f, 53.6489983f },
				{ 68.7760239f, -9.07335377f, 52.7140236f },
				{ 68.9000015f, -10.3944235f, 53.2217026f } } },
		{ 80.3f,
			{ { 0.304512024f, 0.531184196f, 2.25860214f },
				{ 0.956512451f, 0.0321846008f, 1.07460022f },
				{ 0, 0, 0 },
				{ 0.21336633f, -1.09552193f, -0.878999829f } } },
		{ -97.5f,
			{ { 0.088973999f, 1.72535372f, 1.66297531f },
				{ -0.641021729f, 0.886353493f, 0.93497467f },
				{ 0, 0, 0 },
				{ 1.29433727f, -0.395142615f, 0.432300746f } } },
		{ -97.5f,
			{
				{ 0.088973999f, 1.72535372f, 1.66297531f },
				{ -0.641021729f, 0.886353493f, 0.93497467f },
				{ 0, 0, 0 },
				{ 1.33983064f, 0.384027064f, -0.275154471f },

			} }
	};

	for (auto &&[angle, pts] : tests)
	{
		auto q = cif::construct_for_dihedral_angle(pts[0], pts[1], pts[2], pts[3], angle, 1);

		pts[3].rotate(q, pts[2]);

		auto dh = cif::dihedral_angle(pts[0], pts[1], pts[2], pts[3]);
		BOOST_TEST(dh == angle, tt::tolerance(0.1f));
	}
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(symm_1)
{
	cif::cell c(10, 10, 10);

	cif::point p{ 1, 1, 1 };

	cif::point f = fractional(p, c);

	BOOST_TEST(f.m_x == 0.1f, tt::tolerance(0.01));
	BOOST_TEST(f.m_y == 0.1f, tt::tolerance(0.01));
	BOOST_TEST(f.m_z == 0.1f, tt::tolerance(0.01));

	cif::point o = orthogonal(f, c);

	BOOST_TEST(o.m_x == 1.f, tt::tolerance(0.01));
	BOOST_TEST(o.m_y == 1.f, tt::tolerance(0.01));
	BOOST_TEST(o.m_z == 1.f, tt::tolerance(0.01));
}

BOOST_AUTO_TEST_CASE(symm_2)
{
	using namespace cif::literals;

	auto symop = "1_555"_symop;

	BOOST_TEST(symop.is_identity() == true);
}

BOOST_AUTO_TEST_CASE(symm_3)
{
	using namespace cif::literals;

	cif::spacegroup sg(18);

	BOOST_TEST(sg.size() == 4);
	BOOST_TEST(sg.get_name() == "P 21 21 2");
}

BOOST_AUTO_TEST_CASE(symm_4, *utf::tolerance(0.1f))
{
	using namespace cif::literals;

	// based on 2b8h
	auto sg = cif::spacegroup(154); // p 32 2 1
	auto c = cif::cell(107.516, 107.516, 338.487, 90.00, 90.00, 120.00);
	
	cif::point a{   -8.688,  79.351, 10.439 }; // O6 NAG A 500
	cif::point b{  -35.356,  33.693, -3.236 }; // CG2 THR D 400
	cif::point sb(  -6.916,   79.34,   3.236); // 4_565 copy of b

	BOOST_TEST(distance(a, sg(a, c, "1_455"_symop)) == static_cast<float>(c.get_a()));
	BOOST_TEST(distance(a, sg(a, c, "1_545"_symop)) == static_cast<float>(c.get_b()));
	BOOST_TEST(distance(a, sg(a, c, "1_554"_symop)) == static_cast<float>(c.get_c()));

	auto sb2 = sg(b, c, "4_565"_symop);
	BOOST_TEST(sb.m_x == sb2.m_x);
	BOOST_TEST(sb.m_y == sb2.m_y);
	BOOST_TEST(sb.m_z == sb2.m_z);

	BOOST_TEST(distance(a, sb2) == 7.42f);	
}

BOOST_AUTO_TEST_CASE(symm_2bi3_1, *utf::tolerance(0.1f))
{
	cif::file f(gTestDir / "2bi3.cif.gz");

	auto &db = f.front();
	cif::mm::structure s(db);

	cif::spacegroup sg(db);
	cif::cell c(db);

	auto struct_conn = db["struct_conn"];
	for (const auto &[
			asym1, seqid1, authseqid1, atomid1, symm1,
			asym2, seqid2, authseqid2, atomid2, symm2,
			dist] : struct_conn.find<
				std::string,int,std::string,std::string,std::string,
				std::string,int,std::string,std::string,std::string,
				float>(
			cif::key("ptnr1_symmetry") != "1_555" or cif::key("ptnr2_symmetry") != "1_555",
			"ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_auth_seq_id", "ptnr1_label_atom_id", "ptnr1_symmetry", 
			"ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr2_auth_seq_id", "ptnr2_label_atom_id", "ptnr2_symmetry", 
			"pdbx_dist_value"
		))
	{
		auto &r1 = s.get_residue(asym1, seqid1, authseqid1);
		auto &r2 = s.get_residue(asym2, seqid2, authseqid2);

		auto a1 = r1.get_atom_by_atom_id(atomid1);
		auto a2 = r2.get_atom_by_atom_id(atomid2);

		auto sa1 = symmetry_copy(a1.get_location(), sg, c, cif::sym_op(symm1));
		auto sa2 = symmetry_copy(a2.get_location(), sg, c, cif::sym_op(symm2));

		BOOST_TEST(cif::distance(sa1, sa2) == dist);

		auto pa1 = a1.get_location();

		const auto &[d, p, so] = cif::closest_symmetry_copy(sg, c, pa1, a2.get_location());

		BOOST_TEST(p.m_x == sa2.m_x);
		BOOST_TEST(p.m_y == sa2.m_y);
		BOOST_TEST(p.m_z == sa2.m_z);

		BOOST_TEST(d == dist);
		BOOST_TEST(so.string() == symm2);
	}
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(eigen_1, *utf::tolerance(0.1f))
{
	cif::symmetric_matrix4x4<float> m;

	m(0, 0) = 4;
	m(0, 1) = -30;
	m(0, 2) = 60;
	m(0, 3) = -35;
	m(1, 1) = 300;
	m(1, 2) = -675;
	m(1, 3) = 420;
	m(2, 2) = 1620;
	m(2, 3) = -1050;
	m(3, 3) = 700;

	cif::matrix4x4<float> m2;
	m2 = m;

	const auto &[ev, em] = cif::eigen(m2, true);

	BOOST_TEST(ev[0] == 0.1666428611718905f);
	BOOST_TEST(ev[1] == 1.4780548447781369f);
	BOOST_TEST(ev[2] == 37.1014913651276582f);
	BOOST_TEST(ev[3] == 2585.25381092892231f);






// 	=== Example ===

// Let 
// <math>
// 	S = \begin{pmatrix} 4 & -30 & 60 & -35 \\ -30 & 300 & -675 & 420 \\ 60 & -675 & 1620 & -1050 \\ -35 & 420 & -1050 & 700 \end{pmatrix}
// </math>

// Then ''jacobi'' produces the following eigenvalues and eigenvectors after 3 sweeps (19 iterations) :

// <math>
// 	e_1 = 2585.25381092892231
// </math>

// <math>
// 	E_1 = \begin{pmatrix}0.0291933231647860588\\ -0.328712055763188997\\ 0.791411145833126331\\ -0.514552749997152907\end{pmatrix}
// </math>

// <math>
// 	e_2 = 37.1014913651276582
// </math>

// <math>
// 	E_2 = \begin{pmatrix}-0.179186290535454826\\ 0.741917790628453435\\ -0.100228136947192199\\ -0.638282528193614892\end{pmatrix}
// </math>

// <math>
// 	e_3 = 1.4780548447781369
// </math>

// <math>
// 	E_3 = \begin{pmatrix}-0.582075699497237650\\ 0.370502185067093058\\ 0.509578634501799626\\ 0.514048272222164294\end{pmatrix}
// </math>

// <math>
// 	e_4 = 0.1666428611718905
// </math>

// <math>
// 	E_4 = \begin{pmatrix}0.792608291163763585\\ 0.451923120901599794\\ 0.322416398581824992\\ 0.252161169688241933\end{pmatrix}
// </math>

}