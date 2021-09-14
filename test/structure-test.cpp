/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute
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

#define BOOST_TEST_MODULE Structure_Test
#include <boost/test/included/unit_test.hpp>

#include <stdexcept>

#include "cif++/Cif++.hpp"
#include "cif++/Structure.hpp"

// --------------------------------------------------------------------

cif::File operator""_cf(const char* text, size_t length)
{
    struct membuf : public std::streambuf
    {
        membuf(char* text, size_t length)
        {
            this->setg(text, text, text + length);
        }
    } buffer(const_cast<char*>(text), length);

    std::istream is(&buffer);
    return cif::File(is);
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(create_nonpoly_1)
{
    cif::VERBOSE = 1;

	// do this now, avoids the need for installing
	cif::addFileResource("mmcif_pdbx_v50.dic", "../rsrc/mmcif_pdbx_v50.dic");

	mmcif::File file;
	file.file().loadDictionary("mmcif_pdbx_v50.dic");
	file.createDatablock("TEST");	// create a datablock
	
	mmcif::Structure structure(file);

	std::string entity_id = structure.createNonPolyEntity("HEM");

	auto atoms = R"(
data_HEM
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
HETATM C  CHA . ? -5.248  39.769 -0.250  1.00 7.67  ?
HETATM C  CHB . ? -3.774  36.790 3.280   1.00 7.05  ?
HETATM C  CHC . ? -2.879  33.328 0.013   1.00 7.69  ?
HETATM C  CHD . ? -4.342  36.262 -3.536  1.00 8.00  ?
# that's enough to test with
)"_cf;

	auto &atom_site = atoms["HEM"]["atom_site"];

	auto hem_atoms = atom_site.rows();
	std::vector<cif::Row> atom_data(hem_atoms.begin(), hem_atoms.end());

	structure.createNonpoly(entity_id, atom_data);

	auto expected = R"(
data_TEST
#
loop_
_atom_site.id
_atom_site.auth_asym_id
_atom_site.label_alt_id
_atom_site.label_asym_id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.type_symbol
_atom_site.group_PDB
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
1 A ? A CHA HEM 1 . C HETATM ? -5.248 39.769 -0.250 1.00 7.67 ? ? HEM CHA 1
2 A ? A CHB HEM 1 . C HETATM ? -3.774 36.790 3.280  1.00 7.05 ? ? HEM CHB 1
3 A ? A CHC HEM 1 . C HETATM ? -2.879 33.328 0.013  1.00 7.69 ? ? HEM CHC 1
4 A ? A CHD HEM 1 . C HETATM ? -4.342 36.262 -3.536 1.00 8.00 ? ? HEM CHD 1
#
_chem_comp.id               HEM
_chem_comp.type             NON-POLYMER
_chem_comp.name             'PROTOPORPHYRIN IX CONTAINING FE'
_chem_comp.formula          'C34 H32 Fe N4 O4'
_chem_comp.formula_weight   616.487000
#
_pdbx_entity_nonpoly.entity_id   1
_pdbx_entity_nonpoly.name        'PROTOPORPHYRIN IX CONTAINING FE'
_pdbx_entity_nonpoly.comp_id     HEM
#
_entity.id                 1
_entity.type               non-polymer
_entity.pdbx_description   'PROTOPORPHYRIN IX CONTAINING FE'
_entity.formula_weight     616.487000
#
_struct_asym.id                            A
_struct_asym.entity_id                     1
_struct_asym.pdbx_blank_PDB_chainid_flag   N
_struct_asym.pdbx_modified                 N
_struct_asym.details                       ?
#
)"_cf;

	expected.loadDictionary("mmcif_pdbx_v50.dic");

	if (not (expected.firstDatablock() == structure.getFile().data()))
	{
		BOOST_TEST(false);
		std::cout << expected.firstDatablock() << std::endl
				<< std::endl
				<< structure.getFile().data() << std::endl;
	}
}
