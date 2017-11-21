#pragma once

#include "pdb2cif.h"

// --------------------------------------------------------------------

struct TemplateLine;

class Remark3Parser
{
  public:
	virtual ~Remark3Parser() {}

	static bool parse(const std::string& expMethod, PDBRecord* r, cif::datablock& db);

	virtual std::string program();
	virtual std::string version();

  protected:

	Remark3Parser(const std::string& name, const std::string& expMethod, PDBRecord* r, cif::datablock& db,
			const TemplateLine templatelines[], uint32 templateLineCount, std::regex programVersion);

	virtual float Parse();
	std::string NextLine();

	bool Match(const char* expr, int nextState);
	void StoreCapture(const char* category, std::initializer_list<const char*> items, bool createNew = false);
	void StoreRefineLsRestr(const char* type, std::initializer_list<const char*> values);
	void UpdateRefineLsRestr(const char* type, std::initializer_list<const char*> values);

	virtual void Fixup() {}

	std::string		mName;
	std::string		mExpMethod;
	PDBRecord*		mRec;
	cif::datablock	mDb;
	std::string		mLine;
	std::smatch		mM;
	uint32			mState;

	const TemplateLine*	mTemplate;
	uint32				mTemplateCount;
	std::regex			mProgramVersion;
};


