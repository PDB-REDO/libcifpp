// Generated revision file

#pragma once

#include <ostream>

const char kLibCIFPPProjectName[] = "cifpp";
const char kLibCIFPPVersionNumber[] = "5.0.1";
const char kLibCIFPPVersionGitTag[] = "64e40e7";
const char kLibCIFPPBuildInfo[] = "830*";
const char kLibCIFPPBuildDate[] = "2022-11-07T11:26:40Z";

inline void write_version_string(std::ostream &os, bool verbose)
{
	os << kLibCIFPPProjectName << " version " << kLibCIFPPVersionNumber << std::endl;
	if (verbose)
	{
		os << "build: " << kLibCIFPPBuildInfo << ' ' << kLibCIFPPBuildDate << std::endl;
		if (kLibCIFPPVersionGitTag[0] != 0)
			os << "git tag: " << kLibCIFPPVersionGitTag << std::endl;
	}
}
