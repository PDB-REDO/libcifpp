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

	cif::tie(m_a, m_b, m_c, m_alpha, m_beta, m_gamma) =
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
	m_translation.m_y = d[13] == 0 ? 0 : 1.0 * d[13] / d[14];
}

// --------------------------------------------------------------------

spacegroup::spacegroup(int nr)
	: m_nr(nr)
{
	const size_t N = cif::kSymopNrTableSize;
	int32_t L = 0, R = static_cast<int32_t>(N - 1);
	while (L <= R)
	{
		int32_t i = (L + R) / 2;
		if (cif::kSymopNrTable[i].spacegroup() < m_nr)
			L = i + 1;
		else
			R = i - 1;
	}

	m_index = L;

	for (size_t i = L; i < N and cif::kSymopNrTable[i].spacegroup() == m_nr; ++i)
		emplace_back(cif::kSymopNrTable[i].symop().data());
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

point spacegroup::operator()(const point &pt, const cell &c, sym_op symop)
{
	if (symop.m_nr < 1 or symop.m_nr > size())
		throw std::out_of_range("symmetry operator number out of range");
	
	transformation t = at(symop.m_nr - 1);

	return pt;
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

} // namespace cif
