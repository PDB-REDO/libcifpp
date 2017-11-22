#pragma once

#include "cif++/Cif++.h"

// --------------------------------------------------------------------

struct PDBRecord
{
	PDBRecord*	mNext;
	uint32		mLineNr;
	char		mName[11];
	size_t		mVlen;
	char		mValue[0];

	PDBRecord(uint32 lineNr, const std::string& name, const std::string& value);
	~PDBRecord();
	
	void* operator new(size_t);
	void* operator new(size_t size, size_t vLen);
	
	void operator delete(void* p);

	bool is(const char* name) const;
	
	char vC(size_t column);
	std::string vS(size_t columnFirst, size_t columnLast = std::numeric_limits<size_t>::max());
	int vI(int columnFirst, int columnLast);
	std::string vF(size_t columnFirst, size_t columnLast);
};

// --------------------------------------------------------------------

void ReadPDBFile(std::istream& pdbFile, cif::File& cifFile);
