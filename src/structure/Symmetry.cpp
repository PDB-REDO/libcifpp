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

#include <atomic>
#include <mutex>

#include <cif++/structure/Symmetry.hpp>
#include <cif++/utilities.hpp>

#include "./SymOpTable_data.hpp"

namespace mmcif
{

// --------------------------------------------------------------------
// Unfortunately, clipper has a different numbering scheme than PDB
// for rotation numbers. So we created a table to map those.
// Perhaps a bit over the top, but hey....

// --------------------------------------------------------------------

int GetSpacegroupNumber(std::string spacegroup)
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
			auto& sp = kSpaceGroups[i];
			if (sp.xHM == spacegroup)
			{
				result = sp.nr;
				break;
			}
		}
	}

	if (result == 0)
		throw std::runtime_error("Spacegroup name " + spacegroup + " was not found in table");
	
	return result;
}

// --------------------------------------------------------------------

int GetSpacegroupNumber(std::string spacegroup, SpacegroupName type)
{
	if (spacegroup == "P 21 21 2 A")
		spacegroup = "P 21 21 2 (a)";
	else if (spacegroup.empty())
		throw std::runtime_error("No spacegroup, cannot continue");

	int result = 0;

	if (type == SpacegroupName::full)
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
	else if (type == SpacegroupName::xHM)
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
		throw std::runtime_error("Spacegroup name " + spacegroup + " was not found in table");
	
	return result;
}

}
