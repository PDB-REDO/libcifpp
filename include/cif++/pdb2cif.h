#pragma once

#include "cif++.h"

// --------------------------------------------------------------------

struct PDBRecord
{
	PDBRecord*	m_next;
	uint32		m_line_nr;
	char		m_name[11];
	size_t		m_vlen;
	char		m_value[0];

	PDBRecord(uint32 line_nr, const std::string& name, const std::string& value);
	~PDBRecord();
	
	void* operator new(size_t);
	void* operator new(size_t size, size_t v_len);
	
	void operator delete(void* p);

	bool is(const char* name) const;
	
	char v_c(size_t column);
	std::string v_s(size_t column_first, size_t column_last = std::numeric_limits<size_t>::max());
	int v_i(int column_first, int column_last);
	std::string v_f(size_t column_first, size_t column_last);
};

// --------------------------------------------------------------------

void ReadPDBFile(std::istream& pdbFile, cif::file& cifFile);
