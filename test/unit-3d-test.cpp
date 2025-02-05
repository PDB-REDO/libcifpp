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

#include <Eigen/Eigenvalues>

// --------------------------------------------------------------------

cif::file operator""_cf(const char *text, std::size_t length)
{
	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(text), length);

	std::istream is(&buffer);
	return cif::file(is);
}

// --------------------------------------------------------------------

// 3d tests

TEST_CASE("t1")
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

	REQUIRE_THAT(std::fmod(360 + angle, 360), Catch::Matchers::WithinRel(std::fmod(360 - angle0, 360), 0.01));

	for (auto &p : p1)
		p.rotate(q2);

	auto rmsd = cif::RMSd(p1, p2);

	REQUIRE(rmsd < 1e-5);

	// std::cout << "rmsd: " << RMSd(p1, p2) << '\n';
}

TEST_CASE("t2")
{
	cif::point p[] = {
		{ 1, 1, 0 },
		{ 2, 1, 0 },
		{ 1, 2, 0 }
	};

	cif::point xp = cif::cross_product(p[1] - p[0], p[2] - p[0]);

	auto q = cif::construct_from_angle_axis(45, xp); // mmcif::Normalize(Quaternion{45 * mmcif::kPI / 180, xp.mX, xp.mY, xp.mZ});

	auto &&[angle, axis] = cif::quaternion_to_angle_axis(q);

	REQUIRE_THAT(angle, Catch::Matchers::WithinRel(45.f, 0.01f));
}

TEST_CASE("t3")
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

	std::cout << v << '\n';

	double a = cif::angle(v, p[0], p[1]);

	REQUIRE_THAT(a, Catch::Matchers::WithinRel(45.f, 0.01f));
}

TEST_CASE("dh_q_0")
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
	REQUIRE_THAT(a, Catch::Matchers::WithinRel(0.f, 0.01f));

	auto q = cif::construct_from_angle_axis(90, axis);

	p.rotate(q);

	REQUIRE(std::abs(p.m_x - 1.f) < 0.01f);
	REQUIRE(std::abs(p.m_y - 0.f) < 0.01f);
	REQUIRE(std::abs(p.m_z - 1.f) < 0.01f);

	a = cif::dihedral_angle(t[0], t[1], t[2], p);
	REQUIRE(std::abs(a - 90.f) < 0.01f);

	q = cif::construct_from_angle_axis(-90, axis);

	p.rotate(q);

	REQUIRE(std::abs(p.m_x - 1.f) < 0.01f);
	REQUIRE(std::abs(p.m_y - 1.f) < 0.01f);
	REQUIRE(std::abs(p.m_z - 0.f) < 0.01f);

	a = cif::dihedral_angle(t[0], t[1], t[2], p);
	REQUIRE(std::abs(a - 0.f) < 0.01f);

}

TEST_CASE("dh_q_1")
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
		REQUIRE_THAT(dh, Catch::Matchers::WithinRel(angle, 0.1f));
	}
}

// --------------------------------------------------------------------

TEST_CASE("m2q_0, *utf::tolerance(0.001f)")
{
	for (std::size_t i = 0; i < cif::kSymopNrTableSize; ++i)
	{
		auto d = cif::kSymopNrTable[i].symop().data();

		cif::matrix3x3<float> rot;
		float Qxx = rot(0, 0) = d[0];
		float Qxy = rot(0, 1) = d[1];
		float Qxz = rot(0, 2) = d[2];
		float Qyx = rot(1, 0) = d[3];
		float Qyy = rot(1, 1) = d[4];
		float Qyz = rot(1, 2) = d[5];
		float Qzx = rot(2, 0) = d[6];
		float Qzy = rot(2, 1) = d[7];
		float Qzz = rot(2, 2) = d[8];

		Eigen::Matrix4f em;

		em << Qxx - Qyy - Qzz, Qyx + Qxy, Qzx + Qxz, Qzy - Qyz,
		      Qyx + Qxy, Qyy - Qxx - Qzz, Qzy + Qyz, Qxz - Qzx,
			  Qzx + Qxz, Qzy + Qyz, Qzz - Qxx - Qyy, Qyx - Qxy,
			  Qzy - Qyz, Qxz - Qzx, Qyx - Qxy, Qxx + Qyy + Qzz;

		Eigen::EigenSolver<Eigen::Matrix4f> es(em / 3);

		auto ev = es.eigenvalues();

		std::size_t bestJ = 0;
		float bestEV = -1;

		for (std::size_t j = 0; j < 4; ++j)
		{
			if (bestEV < ev[j].real())
			{
				bestEV = ev[j].real();
				bestJ = j;
			}
		}

		if (std::abs(bestEV - 1) > 0.01)
			continue; // not a rotation matrix

		auto col = es.eigenvectors().col(bestJ);

		auto q = normalize(cif::quaternion{
			static_cast<float>(col(3).real()),
			static_cast<float>(col(0).real()),
			static_cast<float>(col(1).real()),
			static_cast<float>(col(2).real()) });
		
		cif::point p1{ 1, 1, 1 };
		cif::point p2 = p1;
		p2.rotate(q);

		cif::point p3 = rot * p1;

		REQUIRE_THAT(p2.m_x, Catch::Matchers::WithinRel(p3.m_x, 0.01f));
		REQUIRE_THAT(p2.m_y, Catch::Matchers::WithinRel(p3.m_y, 0.01f));
		REQUIRE_THAT(p2.m_z, Catch::Matchers::WithinRel(p3.m_z, 0.01f));
	}
}

// "TEST_CASE(m2q_1, *utf::tolerance(0.001f)")
// {
// 	for (std::size_t i = 0; i < cif::kSymopNrTableSize; ++i)
// 	{
// 		auto d = cif::kSymopNrTable[i].symop().data();

// 		cif::matrix3x3<float> rot;
// 		float Qxx = rot(0, 0) = d[0];
// 		float Qxy = rot(0, 1) = d[1];
// 		float Qxz = rot(0, 2) = d[2];
// 		float Qyx = rot(1, 0) = d[3];
// 		float Qyy = rot(1, 1) = d[4];
// 		float Qyz = rot(1, 2) = d[5];
// 		float Qzx = rot(2, 0) = d[6];
// 		float Qzy = rot(2, 1) = d[7];
// 		float Qzz = rot(2, 2) = d[8];

// 		cif::matrix4x4<float> m({
// 			Qxx - Qyy - Qzz, Qyx + Qxy, Qzx + Qxz, Qzy - Qyz,
// 			Qyx + Qxy, Qyy - Qxx - Qzz, Qzy + Qyz, Qxz - Qzx,
// 			Qzx + Qxz, Qzy + Qyz, Qzz - Qxx - Qyy, Qyx - Qxy,
// 			Qzy - Qyz, Qxz - Qzx, Qyx - Qxy, Qxx + Qyy + Qzz
// 		});

// 		auto &&[ev, em] = cif::eigen(m * (1/3.0f), false);

// 		std::size_t bestJ = 0;
// 		float bestEV = -1;

// 		for (std::size_t j = 0; j < 4; ++j)
// 		{
// 			if (bestEV < ev[j])
// 			{
// 				bestEV = ev[j];
// 				bestJ = j;
// 			}
// 		}

// 		if (std::abs(bestEV - 1) > 0.01)
// 			continue; // not a rotation matrix

// 		auto q = normalize(cif::quaternion{
// 			static_cast<float>(em(bestJ, 3)),
// 			static_cast<float>(em(bestJ, 0)),
// 			static_cast<float>(em(bestJ, 1)),
// 			static_cast<float>(em(bestJ, 2)) });
		
// 		cif::point p1{ 1, 1, 1 };
// 		cif::point p2 = p1;
// 		p2.rotate(q);

// 		cif::point p3 = rot * p1;

// 		REQUIRE(p2.m_x == p3.m_x);
// 		REQUIRE(p2.m_y == p3.m_y);
// 		REQUIRE(p2.m_z == p3.m_z);
// 	}
// }

// --------------------------------------------------------------------

TEST_CASE("symm_1")
{
	cif::cell c(10, 10, 10);

	cif::point p{ 1, 1, 1 };

	cif::point f = fractional(p, c);

	REQUIRE_THAT(f.m_x, Catch::Matchers::WithinRel(0.1f, 0.01f));
	REQUIRE_THAT(f.m_y, Catch::Matchers::WithinRel(0.1f, 0.01f));
	REQUIRE_THAT(f.m_z, Catch::Matchers::WithinRel(0.1f, 0.01f));

	cif::point o = orthogonal(f, c);

	REQUIRE_THAT(o.m_x, Catch::Matchers::WithinRel(1.f, 0.01f));
	REQUIRE_THAT(o.m_y, Catch::Matchers::WithinRel(1.f, 0.01f));
	REQUIRE_THAT(o.m_z, Catch::Matchers::WithinRel(1.f, 0.01f));
}

TEST_CASE("symm_2")
{
	using namespace cif::literals;

	auto symop = "1_555"_symop;

	REQUIRE(symop.is_identity() == true);
}

TEST_CASE("symm_3")
{
	using namespace cif::literals;

	cif::spacegroup sg(18);

	REQUIRE(sg.size() == 4UL);
	REQUIRE(sg.get_name() == "P 21 21 2");
}

TEST_CASE("symm_4, *utf::tolerance(0.1f)")
{
	using namespace cif::literals;

	// based on 2b8h
	auto sg = cif::spacegroup(154); // p 32 2 1
	auto c = cif::cell(107.516, 107.516, 338.487, 90.00, 90.00, 120.00);
	
	cif::point a{   -8.688,  79.351, 10.439 }; // O6 NAG A 500
	cif::point b{  -35.356,  33.693, -3.236 }; // CG2 THR D 400
	cif::point sb(  -6.916,   79.34,   3.236); // 4_565 copy of b

	REQUIRE_THAT(distance(a, sg(a, c, "1_455"_symop)), Catch::Matchers::WithinRel(static_cast<float>(c.get_a()), 0.01f));
	REQUIRE_THAT(distance(a, sg(a, c, "1_545"_symop)), Catch::Matchers::WithinRel(static_cast<float>(c.get_b()), 0.01f));
	REQUIRE_THAT(distance(a, sg(a, c, "1_554"_symop)), Catch::Matchers::WithinRel(static_cast<float>(c.get_c()), 0.01f));

	auto sb2 = sg(b, c, "4_565"_symop);
	REQUIRE_THAT(sb.m_x, Catch::Matchers::WithinRel(sb2.m_x, 0.01f));
	REQUIRE_THAT(sb.m_y, Catch::Matchers::WithinRel(sb2.m_y, 0.01f));
	REQUIRE_THAT(sb.m_z, Catch::Matchers::WithinRel(sb2.m_z, 0.01f));

	REQUIRE_THAT(distance(a, sb2), Catch::Matchers::WithinRel(7.42f, 0.01f));	
}

// --------------------------------------------------------------------

TEST_CASE("symm_4wvp_1, *utf::tolerance(0.1f)")
{
	using namespace cif::literals;

	cif::file f(gTestDir / "4wvp.cif.gz");

	auto &db = f.front();
	cif::mm::structure s(db);

	cif::crystal c(db);

	cif::point p{ -78.722, 98.528,  11.994 };
	auto a = s.get_residue("A", 10, "").get_atom_by_atom_id("O");

	auto sp1 = c.symmetry_copy(a.get_location(), "2_565"_symop);
	REQUIRE_THAT(sp1.m_x, Catch::Matchers::WithinAbs(p.m_x, 0.5f));
	REQUIRE_THAT(sp1.m_y, Catch::Matchers::WithinAbs(p.m_y, 0.5f));
	REQUIRE_THAT(sp1.m_z, Catch::Matchers::WithinAbs(p.m_z, 0.5f));

	const auto &[d, sp2, so] = c.closest_symmetry_copy(p, a.get_location());

	REQUIRE(d < 1);

	REQUIRE_THAT(sp2.m_x, Catch::Matchers::WithinAbs(p.m_x, 0.5f));
	REQUIRE_THAT(sp2.m_y, Catch::Matchers::WithinAbs(p.m_y, 0.5f));
	REQUIRE_THAT(sp2.m_z, Catch::Matchers::WithinAbs(p.m_z, 0.5f));

}

TEST_CASE("symm_2bi3_1, *utf::tolerance(0.1f)")
{
	cif::file f(gTestDir / "2bi3.cif.gz");

	auto &db = f.front();
	cif::mm::structure s(db);

	cif::crystal c(db);

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

		auto sa1 = c.symmetry_copy(a1.get_location(), cif::sym_op(symm1));
		auto sa2 = c.symmetry_copy(a2.get_location(), cif::sym_op(symm2));

		REQUIRE_THAT(cif::distance(sa1, sa2), Catch::Matchers::WithinAbs(dist, 0.5f));

		auto pa1 = a1.get_location();

		const auto &[d, p, so] = c.closest_symmetry_copy(pa1, a2.get_location());

		REQUIRE_THAT(p.m_x, Catch::Matchers::WithinAbs(sa2.m_x, 0.5f));
		REQUIRE_THAT(p.m_y, Catch::Matchers::WithinAbs(sa2.m_y, 0.5f));
		REQUIRE_THAT(p.m_z, Catch::Matchers::WithinAbs(sa2.m_z, 0.5f));

		REQUIRE_THAT(d, Catch::Matchers::WithinAbs(dist, 0.5f));
		REQUIRE(so.string() == symm2);
	}
}

TEST_CASE("symm_2bi3_1a, *utf::tolerance(0.1f)")
{
	using namespace cif::literals;

	cif::file f(gTestDir / "2bi3.cif.gz");

	auto &db = f.front();

	cif::crystal c(db);
	auto struct_conn = db["struct_conn"];
	auto atom_site = db["atom_site"];

	for (const auto &[
			asym1, seqid1, authseqid1, atomid1, symm1,
			asym2, seqid2, authseqid2, atomid2, symm2,
			dist] : struct_conn.find<
				std::string,std::optional<int>,std::string,std::string,std::string,
				std::string,std::optional<int>,std::string,std::string,std::string,
				float>(
			cif::key("ptnr1_symmetry") != "1_555" or cif::key("ptnr2_symmetry") != "1_555",
			"ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_auth_seq_id", "ptnr1_label_atom_id", "ptnr1_symmetry", 
			"ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr2_auth_seq_id", "ptnr2_label_atom_id", "ptnr2_symmetry", 
			"pdbx_dist_value"
		))
	{
		cif::point p1 = atom_site.find1<float,float,float>(
			"label_asym_id"_key == asym1 and "label_seq_id"_key == seqid1 and "auth_seq_id"_key == authseqid1 and "label_atom_id"_key == atomid1,
			"cartn_x", "cartn_y", "cartn_z");
		cif::point p2 = atom_site.find1<float,float,float>(
			"label_asym_id"_key == asym2 and "label_seq_id"_key == seqid2 and "auth_seq_id"_key == authseqid2 and "label_atom_id"_key == atomid2,
			"cartn_x", "cartn_y", "cartn_z");

		auto sa1 = c.symmetry_copy(p1, cif::sym_op(symm1));
		auto sa2 = c.symmetry_copy(p2, cif::sym_op(symm2));

		REQUIRE_THAT(cif::distance(sa1, sa2), Catch::Matchers::WithinAbs(dist, 0.5f));

		const auto &[d, p, so] = c.closest_symmetry_copy(p1, p2);

		REQUIRE_THAT(p.m_x, Catch::Matchers::WithinAbs(sa2.m_x, 0.5f));
		REQUIRE_THAT(p.m_y, Catch::Matchers::WithinAbs(sa2.m_y, 0.5f));
		REQUIRE_THAT(p.m_z, Catch::Matchers::WithinAbs(sa2.m_z, 0.5f));

		REQUIRE_THAT(d, Catch::Matchers::WithinAbs(dist, 0.5f));
		REQUIRE(so.string() == symm2);
	}
}

TEST_CASE("symm_3bwh_1, *utf::tolerance(0.1f)")
{
	cif::file f(gTestDir / "3bwh.cif.gz");

	auto &db = f.front();

	cif::crystal c(db);
	cif::mm::structure s(db);

	for (auto a1 : s.atoms())
	{
		for (auto a2 : s.atoms())
		{
			if (a1 == a2)
				continue;
			
			const auto&[ d, p, so ] = c.closest_symmetry_copy(a1.get_location(), a2.get_location());

			REQUIRE_THAT(d, Catch::Matchers::WithinAbs(distance(a1.get_location(), p), 0.5f));
		}
	}
}

TEST_CASE("volume_3bwh_1, *utf::tolerance(0.1f)")
{
	cif::file f(gTestDir / "1juh.cif.gz");

	auto &db = f.front();

	cif::crystal c(db);

	REQUIRE_THAT(c.get_cell().get_volume(), Catch::Matchers::WithinRel(741009.625f, 0.01f));
}

