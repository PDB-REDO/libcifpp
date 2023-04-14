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

#include "cif++/symmetry.hpp"
#include "cif++/datablock.hpp"
#include "cif++/point.hpp"

#include <stdexcept>

#include "symop_table_data.hpp"

namespace cif
{

// --------------------------------------------------------------------

cell::cell(float a, float b, float c, float alpha, float beta, float gamma)
	: m_a(a)
	, m_b(b)
	, m_c(c)
	, m_alpha(alpha)
	, m_beta(beta)
	, m_gamma(gamma)
{
	init();
}

cell::cell(const datablock &db)
{
	auto &_cell = db["cell"];

	tie(m_a, m_b, m_c, m_alpha, m_beta, m_gamma) =
		_cell.front().get("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma");

	init();
}

void cell::init()
{
	auto alpha = (m_alpha * kPI) / 180;
	auto beta = (m_beta * kPI) / 180;
	auto gamma = (m_gamma * kPI) / 180;

	auto alpha_star = std::acos((std::cos(gamma) * std::cos(beta) - std::cos(alpha)) / (std::sin(beta) * std::sin(gamma)));

	m_orthogonal = identity_matrix(3);

	m_orthogonal(0, 0) = m_a;
	m_orthogonal(0, 1) = m_b * std::cos(gamma);
	m_orthogonal(0, 2) = m_c * std::cos(beta);
	m_orthogonal(1, 1) = m_b * std::sin(gamma);
	m_orthogonal(1, 2) = -m_c * std::sin(beta) * std::cos(alpha_star);
	m_orthogonal(2, 2) = m_c * std::sin(beta) * std::sin(alpha_star);

	m_fractional = inverse(m_orthogonal);
}

// --------------------------------------------------------------------

sym_op::sym_op(std::string_view s)
{
	auto b = s.data();
	auto e = b + s.length();

	int rnri;

	auto r = std::from_chars(b, e, rnri);
	
	m_nr = rnri;
	m_ta = r.ptr[1] - '0';
	m_tb = r.ptr[2] - '0';
	m_tc = r.ptr[3] - '0';

	if (r.ec != std::errc() or rnri > 192 or r.ptr[0] != '_' or m_ta > 9 or m_tb > 9 or m_tc > 9)
		throw std::invalid_argument("Could not convert string into sym_op");
}

std::string sym_op::string() const
{
	char b[9];
	auto r = std::to_chars(b, b + sizeof(b), m_nr);
	if (r.ec != std::errc() or r.ptr > b + 4)
		throw std::runtime_error("Could not write out symmetry operation to string");
	
	*r.ptr++ = '_';
	*r.ptr++ = '0' + m_ta;
	*r.ptr++ = '0' + m_tb;
	*r.ptr++ = '0' + m_tc;
	*r.ptr = 0;

	return { b, static_cast<size_t>(r.ptr - b) };
}

// --------------------------------------------------------------------

transformation::transformation(const symop_data &data)
{
	const auto &d = data.data();

	m_rotation(0, 0) = d[0];
	m_rotation(0, 1) = d[1];
	m_rotation(0, 2) = d[2];
	m_rotation(1, 0) = d[3];
	m_rotation(1, 1) = d[4];
	m_rotation(1, 2) = d[5];
	m_rotation(2, 0) = d[6];
	m_rotation(2, 1) = d[7];
	m_rotation(2, 2) = d[8];

	m_translation.m_x = d[9] == 0 ? 0 : 1.0 * d[9] / d[10];
	m_translation.m_y = d[11] == 0 ? 0 : 1.0 * d[11] / d[12];
	m_translation.m_z = d[13] == 0 ? 0 : 1.0 * d[13] / d[14];
}

transformation operator*(const transformation &lhs, const transformation &rhs)
{
	auto r = lhs.m_rotation * rhs.m_rotation;
	auto t = lhs.m_rotation * rhs.m_translation;
	t = t + lhs.m_translation;

	return transformation(r, t);

	// return transformation(lhs.m_rotation * rhs.m_rotation, lhs.m_rotation * rhs.m_translation + lhs.m_translation);
}

transformation inverse(const transformation &t)
{
	auto inv_matrix = inverse(t.m_rotation);
	return { inv_matrix, -(inv_matrix * t.m_translation) };
}

// --------------------------------------------------------------------

spacegroup::spacegroup(int nr)
	: m_nr(nr)
{
	const size_t N = kSymopNrTableSize;
	int32_t L = 0, R = static_cast<int32_t>(N - 1);
	while (L <= R)
	{
		int32_t i = (L + R) / 2;
		if (kSymopNrTable[i].spacegroup() < m_nr)
			L = i + 1;
		else
			R = i - 1;
	}

	m_index = L;

	for (size_t i = L; i < N and kSymopNrTable[i].spacegroup() == m_nr; ++i)
		emplace_back(kSymopNrTable[i].symop().data());
}

std::string spacegroup::get_name() const
{
	for (auto &s : kSpaceGroups)
	{
		if (s.nr == m_nr)
			return s.name;
	}

	throw std::runtime_error("Spacegroup has an invalid number: " + std::to_string(m_nr));
}

point offsetToOrigin(const cell &c, const point &p)
{
	point d{};

	while (p.m_x + d.m_x < (c.get_a() / 2))
		d.m_x += c.get_a();
	while (p.m_x + d.m_x > (c.get_a() / 2))
		d.m_x -= c.get_a();

	while (p.m_y + d.m_y < (c.get_b() / 2))
		d.m_y += c.get_b();
	while (p.m_y + d.m_y > (c.get_b() / 2))
		d.m_y -= c.get_b();

	while (p.m_z + d.m_z < (c.get_c() / 2))
		d.m_z += c.get_c();
	while (p.m_z + d.m_z > (c.get_c() / 2))
		d.m_z -= c.get_c();

	return d;
};

std::tuple<int,int,int> offsetToOriginInt(const cell &c, const point &p)
{
	auto o = offsetToOrigin(c, p);
	return {
		std::rintf(o.m_x / c.get_a()),
		std::rintf(o.m_y / c.get_b()),
		std::rintf(o.m_z / c.get_c())
	};
}

point spacegroup::operator()(const point &pt, const cell &c, sym_op symop) const
{
	if (symop.m_nr < 1 or symop.m_nr > size())
		throw std::out_of_range("symmetry operator number out of range");
	
	transformation t = at(symop.m_nr - 1);

	t.m_translation.m_x += symop.m_ta - 5;
	t.m_translation.m_y += symop.m_tb - 5;
	t.m_translation.m_z += symop.m_tc - 5;

	auto t_orth = orthogonal(t, c);

	auto o = offsetToOrigin(c, pt);
	transformation tlo(identity_matrix<float>(3), o);
	auto itlo = cif::inverse(tlo);

	point result = pt + o;
	result = t_orth(result);
	result = itlo(result);

	return result;
}

point spacegroup::inverse(const point &pt, const cell &c, sym_op symop) const
{
	if (symop.m_nr < 1 or symop.m_nr > size())
		throw std::out_of_range("symmetry operator number out of range");
	
	transformation t = at(symop.m_nr - 1);

	t.m_translation.m_x += symop.m_ta - 5;
	t.m_translation.m_y += symop.m_tb - 5;
	t.m_translation.m_z += symop.m_tc - 5;

	auto t_orth = cif::inverse(orthogonal(t, c));

	auto o = offsetToOrigin(c, pt);
	transformation tlo(identity_matrix<float>(3), o);
	auto itlo = cif::inverse(tlo);

	point result = pt + o;
	result = t_orth(result);
	result = itlo(result);

	return result;
}

// --------------------------------------------------------------------
// Unfortunately, clipper has a different numbering scheme than PDB
// for rotation numbers. So we created a table to map those.
// Perhaps a bit over the top, but hey....

// --------------------------------------------------------------------

int get_space_group_number(std::string_view spacegroup)
{
	if (spacegroup == "P 21 21 2 A")
		spacegroup = "P 21 21 2 (a)";
	else if (spacegroup.empty())
		throw std::runtime_error("No spacegroup, cannot continue");

	int result = 0;

	const size_t N = kNrOfSpaceGroups;
	int32_t L = 0, R = static_cast<int32_t>(N - 1);
	while (L <= R)
	{
		int32_t i = (L + R) / 2;

		int d = spacegroup.compare(kSpaceGroups[i].name);

		if (d > 0)
			L = i + 1;
		else if (d < 0)
			R = i - 1;
		else
		{
			result = kSpaceGroups[i].nr;
			break;
		}
	}

	// not found, see if we can find a match based on xHM name
	if (result == 0)
	{
		for (size_t i = 0; i < kNrOfSpaceGroups; ++i)
		{
			auto &sp = kSpaceGroups[i];
			if (sp.xHM == spacegroup)
			{
				result = sp.nr;
				break;
			}
		}
	}

	if (result == 0)
		throw std::runtime_error("Spacegroup name " + std::string(spacegroup) + " was not found in table");

	return result;
}

// --------------------------------------------------------------------

int get_space_group_number(std::string_view spacegroup, space_group_name type)
{
	if (spacegroup == "P 21 21 2 A")
		spacegroup = "P 21 21 2 (a)";
	else if (spacegroup.empty())
		throw std::runtime_error("No spacegroup, cannot continue");

	int result = 0;

	if (type == space_group_name::full)
	{
		const size_t N = kNrOfSpaceGroups;
		int32_t L = 0, R = static_cast<int32_t>(N - 1);
		while (L <= R)
		{
			int32_t i = (L + R) / 2;

			int d = spacegroup.compare(kSpaceGroups[i].name);

			if (d > 0)
				L = i + 1;
			else if (d < 0)
				R = i - 1;
			else
			{
				result = kSpaceGroups[i].nr;
				break;
			}
		}
	}
	else if (type == space_group_name::xHM)
	{
		for (auto &sg : kSpaceGroups)
		{
			if (sg.xHM == spacegroup)
			{
				result = sg.nr;
				break;
			}
		}
	}
	else
	{
		for (auto &sg : kSpaceGroups)
		{
			if (sg.Hall == spacegroup)
			{
				result = sg.nr;
				break;
			}
		}
	}

	// not found, see if we can find a match based on xHM name
	if (result == 0)
		throw std::runtime_error("Spacegroup name " + std::string(spacegroup) + " was not found in table");

	return result;
}

int get_space_group_number(const datablock &db)
{
	auto &_symmetry = db["symmetry"];

	if (_symmetry.size() != 1)
		throw std::runtime_error("Could not find a unique symmetry in this mmCIF file");
	
	return _symmetry.front().get<int>("Int_Tables_number");
}

// --------------------------------------------------------------------

std::tuple<float,point,sym_op> closest_symmetry_copy(const spacegroup &sg, const cell &c, point a, point b)
{
	if (c.get_a() == 0 or c.get_b() == 0 or c.get_c() == 0)
		throw std::runtime_error("Invalid cell, contains a dimension that is zero");

	point result_p;
	float result_d = std::numeric_limits<float>::max();
	sym_op result_s;

	auto o = offsetToOrigin(c, a);

	a += o;
	b += o;

	auto fa = fractional(a, c);
	auto fb = fractional(b, c);

	for (size_t i = 0; i < sg.size(); ++i)
	{
		sym_op s(i + 1);
		auto &t = sg[i];

		auto fsb = t(fb);

		while (fsb.m_x - 0.5f > fa.m_x)
		{
			fsb.m_x -= 1;
			s.m_ta -= 1;
		}

		while (fsb.m_x + 0.5f < fa.m_x)
		{
			fsb.m_x += 1;
			s.m_ta += 1;			
		}

		while (fsb.m_y - 0.5f > fa.m_y)
		{
			fsb.m_y -= 1;
			s.m_tb -= 1;
		}

		while (fsb.m_y + 0.5f < fa.m_y)
		{
			fsb.m_y += 1;
			s.m_tb += 1;			
		}

		while (fsb.m_z - 0.5f > fa.m_z)
		{
			fsb.m_z -= 1;
			s.m_tc -= 1;
		}

		while (fsb.m_z + 0.5f < fa.m_z)
		{
			fsb.m_z += 1;
			s.m_tc += 1;			
		}

		auto p = orthogonal(fsb, c);
		auto dsq = distance_squared(a, p);
		if (result_d > dsq)
		{
			result_d = dsq;
			result_p = p;
			result_s = s;
		}
	}

	transformation tlo(identity_matrix<float>(3), o);
	auto itlo = cif::inverse(tlo);

	return { std::sqrt(result_d), itlo(result_p), result_s };
}

} // namespace cif
