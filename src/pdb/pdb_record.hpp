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

#pragma once

#include "cif++/file.hpp"

/// \file pdb_record.hpp

namespace cif::pdb
{

// --------------------------------------------------------------------

struct PDBRecord
{
	PDBRecord *mNext;
	uint32_t mLineNr;
	char mName[11];
	std::size_t mVlen;
	char mValue[1];

	PDBRecord(uint32_t lineNr, const std::string &name, const std::string &value);
	~PDBRecord();

	void *operator new(std::size_t);
	void *operator new(std::size_t size, std::size_t vLen);

	void operator delete(void *p);
	void operator delete(void *p, std::size_t vLen);

	bool is(const char *name) const;

	char vC(std::size_t column);
	std::string vS(std::size_t columnFirst, std::size_t columnLast = std::numeric_limits<std::size_t>::max());
	int vI(int columnFirst, int columnLast);
	std::string vF(std::size_t columnFirst, std::size_t columnLast);
};

} // namespace pdbx