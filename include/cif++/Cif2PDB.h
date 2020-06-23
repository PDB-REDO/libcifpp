#pragma once

#include "cif++/Cif++.h"

void WritePDBFile(std::ostream& pdbFile, cif::File& cifFile);

/// \brief Just the HEADER, COMPND, SOURCE and AUTHOR lines
void WritePDBHeaderLines(std::ostream& os, cif::File& cifFile);

std::string GetPDBHEADERLine(cif::File& cifFile, int truncate_at = 127);
std::string GetPDBCOMPNDLine(cif::File& cifFile, int truncate_at = 127);
std::string GetPDBSOURCELine(cif::File& cifFile, int truncate_at = 127);
std::string GetPDBAUTHORLine(cif::File& cifFile, int truncate_at = 127);
