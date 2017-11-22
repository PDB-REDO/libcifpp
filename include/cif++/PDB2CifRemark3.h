#pragma once

#include "cif++/PDB2Cif.h"

// --------------------------------------------------------------------

struct TemplateLine;

class Remark3Parser
{
  public:
	virtual ~Remark3Parser() {}

	static bool parse(const std::string& expMethod, PDBRecord* r, cif::Datablock& db);

	virtual std::string program();
	virtual std::string version();

  protected:

	Remark3Parser(const std::string& name, const std::string& expMethod, PDBRecord* r, cif::Datablock& db,
			const TemplateLine templatelines[], uint32 templateLineCount, std::regex programVersion);

	virtual float parse();
	std::string nextLine();

	bool match(const char* expr, int nextState);
	void storeCapture(const char* category, std::initializer_list<const char*> items, bool createNew = false);
	void storeRefineLsRestr(const char* type, std::initializer_list<const char*> values);
	void updateRefineLsRestr(const char* type, std::initializer_list<const char*> values);

	virtual void fixup() {}

	std::string		mName;
	std::string		mExpMethod;
	PDBRecord*		mRec;
	cif::Datablock	mDb;
	std::string		mLine;
	std::smatch		mM;
	uint32			mState;

	const TemplateLine*	mTemplate;
	uint32				mTemplateCount;
	std::regex			mProgramVersion;
};


