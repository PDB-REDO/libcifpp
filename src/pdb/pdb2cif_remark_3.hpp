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

#include "pdb_record.hpp"

// --------------------------------------------------------------------

namespace cif::pdb
{

struct TemplateLine;

class Remark3Parser
{
  public:
	virtual ~Remark3Parser() {}

	static bool parse(const std::string &expMethod, PDBRecord *r, cif::datablock &db);

	virtual std::string program();
	virtual std::string version();

  protected:
	Remark3Parser(const std::string &name, const std::string &expMethod, PDBRecord *r, cif::datablock &db,
		const TemplateLine templatelines[], uint32_t templateLineCount, std::regex programVersion);

	virtual float parse();
	std::string nextLine();

	bool match(const char *expr, int nextState);
	void storeCapture(const char *category, std::initializer_list<const char *> items, bool createNew = false);
	void storeRefineLsRestr(const char *type, std::initializer_list<const char *> values);
	void updateRefineLsRestr(const char *type, std::initializer_list<const char *> values);

	virtual void fixup() {}

	std::string mName;
	std::string mExpMethod;
	PDBRecord *mRec;
	cif::datablock mDb;
	std::string mLine;
	std::smatch mM;
	uint32_t mState;

	const TemplateLine *mTemplate;
	uint32_t mTemplateCount;
	std::regex mProgramVersion;
};

} // namespace pdbx
