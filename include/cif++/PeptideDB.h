#pragma once

#include <map>
#include <string>

#include <boost/filesystem/path.hpp>

extern const std::map<std::string,char> kAAMap, kBaseMap;

class PeptideDB
{
  public:
	static PeptideDB& Instance();
	
	void pushDictionary(boost::filesystem::path dict);
	void popDictionary();
	
	bool isKnownPeptide(const std::string& res_name) const;
	bool isKnownBase(const std::string& res_name) const;

	std::string nameForResidue(const std::string& res_name) const;
	std::string formulaForResidue(const std::string& res_name) const;
	std::string unalias(const std::string& res_name) const;

  private:
	PeptideDB();
	~PeptideDB();
	
	PeptideDB(const PeptideDB&) = delete;
	PeptideDB& operator=(const PeptideDB&) = delete;

	struct PeptideDBImpl*	mImpl;
	static PeptideDB*		sInstance;
};
