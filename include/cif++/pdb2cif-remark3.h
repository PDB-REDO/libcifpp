#pragma once

#include "pdb2cif.h"

// --------------------------------------------------------------------

struct TemplateLine;

class Remark3Parser
{
  public:
	virtual ~Remark3Parser() {}

	static bool Parse(const std::string& expMethod, PDBRecord* r, cif::datablock& db);

	virtual std::string Program();
	virtual std::string Version();

  protected:

	Remark3Parser(const std::string& name, const std::string& expMethod, PDBRecord* r, cif::datablock& db,
			const TemplateLine templatelines[], uint32 templateLineCount, std::regex program_version);

	virtual float Parse();
	std::string NextLine();

	bool Match(const char* expr, int nextState);
	void StoreCapture(const char* category, std::initializer_list<const char*> items, bool createNew = false);
	void StoreRefineLsRestr(const char* type, std::initializer_list<const char*> values);
	void UpdateRefineLsRestr(const char* type, std::initializer_list<const char*> values);

	virtual void Fixup() {}

	std::string		m_name;
	std::string		m_expMethod;
	PDBRecord*		m_rec;
	cif::datablock	m_db;
	std::string		m_line;
	std::smatch		m_m;
	uint32			m_state;

	const TemplateLine*	m_template;
	uint32				m_templateCount;
	std::regex			m_program_version;
};


