#pragma once

#include <map>
#include <string>

#include <boost/filesystem/path.hpp>

extern const std::map<std::string,char> kAAMap, kBaseMap;

class PeptideDB
{
  public:
	static PeptideDB& Instance();
	
	void PushDictionary(boost::filesystem::path dict);
	void PopDictionary();
	
	bool IsKnownPeptide(const std::string& res_name) const;
	bool IsKnownBase(const std::string& res_name) const;

	std::string GetNameForResidue(const std::string& res_name) const;
	std::string GetFormulaForResidue(const std::string& res_name) const;
	std::string Unalias(const std::string& res_name) const;

  private:
	PeptideDB();
	~PeptideDB();
	
	PeptideDB(const PeptideDB&) = delete;
	PeptideDB& operator=(const PeptideDB&) = delete;

	struct PeptideDBImpl*	mImpl;
	static PeptideDB*		sInstance;
};
