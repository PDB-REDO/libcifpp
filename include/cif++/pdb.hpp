/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2023 NKI/AVL, Netherlands Cancer Institute
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

#include "cif++/file.hpp"

#include <system_error>

/**
 * @file pdb.hpp
 *
 * This file presents the API to read and write files in the
 * legacy and ancient PDB format.
 *
 * The code works on the basis of best effort since it is
 * impossible to have correct round trip fidelity.
 *
 */

namespace cif::pdb
{

/// --------------------------------------------------------------------
// PDB to mmCIF

/** @brief Read a file in either mmCIF or PDB format from file @a file,
 * compressed or not, depending on the content.
 */

file read(const std::filesystem::path &file);

/** @brief Read a file in either mmCIF or PDB format from std::istream @a is,
 * compressed or not, depending on the content.
 */

file read(std::istream &is);

/**
 * @brief Read a file in legacy PDB format from std::istream @a is and
 * put the data into @a cifFile
 */
file read_pdb_file(std::istream &pdbFile);

// mmCIF to PDB

/** @brief Write out the data in @a db in legacy PDB format
 * to std::ostream @a os
 */
void write(std::ostream &os, const datablock &db);

/** @brief Write out the data in @a f in legacy PDB format
 * to std::ostream @a os
 */
inline void write(std::ostream &os, const file &f)
{
	write(os, f.front());
}

/** @brief Write out the data in @a db to file @a file
 * in legacy PDB format or mmCIF format, depending on the
 * filename extension.
 *
 * If extension of @a file is *.gz* the resulting file will
 * be written in gzip compressed format.
 */
void write(const std::filesystem::path &file, const datablock &db);

/** @brief Write out the data in @a f to file @a file
 * in legacy PDB format or mmCIF format, depending on the
 * filename extension.
 *
 * If extension of @a file is *.gz* the resulting file will
 * be written in gzip compressed format.
 */
inline void write(const std::filesystem::path &p, const file &f)
{
	write(p, f.front());
}

// --------------------------------------------------------------------

/** \brief Reconstruct all missing categories for an assumed PDBx file.
 *
 * Some people believe that simply dumping some atom records is enough.
 *
 * \param file The cif::file that hopefully contains some valid data
 * \param dictionary The mmcif dictionary to use
 * \result Returns true if the resulting file is valid
 */

bool reconstruct_pdbx(file &pdbx_file, std::string_view dictionary = "mmcif_pdbx");

/** \brief This is an extension to cif::validator, use the logic in common
 * PDBx files to see if the file is internally consistent.
 *
 * This function for now checks if the following categories are consistent:
 *
 * atom_site -> pdbx_poly_seq_scheme -> entity_poly_seq -> entity_poly -> entity
 *
 * Use the common \ref cif::VERBOSE flag to turn on diagnostic messages.
 * 
 * This function throws a std::system_error in case of an error
 *
 * \param file The input file
 * \param dictionary The mmcif dictionary to use
 * \result Returns true if the file was valid and consistent
 */

bool is_valid_pdbx_file(const file &pdbx_file, std::string_view dictionary = "mmcif_pdbx");

/** \brief This is an extension to cif::validator, use the logic in common
 * PDBx files to see if the file is internally consistent.
 *
 * This function for now checks if the following categories are consistent:
 *
 * atom_site -> pdbx_poly_seq_scheme -> entity_poly_seq -> entity_poly -> entity
 *
 * Use the common \ref cif::VERBOSE flag to turn on diagnostic messages.
 * 
 * The dictionary is assumed to be specified in the file or to be the
 * default mmcif_pdbx.dic dictionary.
 *
 * \param file The input file
 * \param ec The error_code in case something was wrong
 * \result Returns true if the file was valid and consistent
 */

bool is_valid_pdbx_file(const file &pdbx_file, std::error_code &ec);

/** \brief This is an extension to cif::validator, use the logic in common
 * PDBx files to see if the file is internally consistent.
 *
 * This function for now checks if the following categories are consistent:
 *
 * atom_site -> pdbx_poly_seq_scheme -> entity_poly_seq -> entity_poly -> entity
 *
 * Use the common \ref cif::VERBOSE flag to turn on diagnostic messages.
 * 
 * \param file The input file
 * \param dictionary The dictionary to use
 * \param ec The error_code in case something was wrong
 * \result Returns true if the file was valid and consistent
 */

bool is_valid_pdbx_file(const file &pdbx_file, std::string_view dictionary,
	std::error_code &ec);

// --------------------------------------------------------------------
// Other I/O related routines

/** @brief Return the HEADER line for the data in @a data
 *
 * The line returned should be compatible with the legacy PDB
 * format and is e.g. used in the DSSP program.
 *
 * @param data The datablock to use as source for the requested data
 * @param truncate_at The maximum length of the line returned
 */

std::string get_HEADER_line(const datablock &data, std::string::size_type truncate_at = 127);
/** @brief Return the COMPND line for the data in @a data
 *
 * The line returned should be compatible with the legacy PDB
 * format and is e.g. used in the DSSP program.
 *
 * @param data The datablock to use as source for the requested data
 * @param truncate_at The maximum length of the line returned
 */

std::string get_COMPND_line(const datablock &data, std::string::size_type truncate_at = 127);
/** @brief Return the SOURCE line for the data in @a data
 *
 * The line returned should be compatible with the legacy PDB
 * format and is e.g. used in the DSSP program.
 *
 * @param data The datablock to use as source for the requested data
 * @param truncate_at The maximum length of the line returned
 */

std::string get_SOURCE_line(const datablock &data, std::string::size_type truncate_at = 127);
/** @brief Return the AUTHOR line for the data in @a data
 *
 * The line returned should be compatible with the legacy PDB
 * format and is e.g. used in the DSSP program.
 *
 * @param data The datablock to use as source for the requested data
 * @param truncate_at The maximum length of the line returned
 */

std::string get_AUTHOR_line(const datablock &data, std::string::size_type truncate_at = 127);

} // namespace cif::pdb
