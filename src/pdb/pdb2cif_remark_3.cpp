/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
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

#include "pdb2cif_remark_3.hpp"

#include <cif++.hpp>

#include <map>
#include <set>

namespace cif::pdb
{

// --------------------------------------------------------------------

struct TemplateLine
{
	const char *rx;
	int nextStateOffset;
	const char *category;
	std::initializer_list<const char *> items;
	const char *lsRestrType = nullptr;
	bool createNew;
};

// --------------------------------------------------------------------

const TemplateLine kBusterTNT_Template[] = {
	/* 0 */ { R"(DATA USED IN REFINEMENT\.)", 1 },
	/* 1 */ { R"(RESOLUTION RANGE HIGH \(ANGSTROMS\) :\s+(.+?))", 1, "refine", { "ls_d_res_high" } },
	/* 2 */ { R"(RESOLUTION RANGE LOW \(ANGSTROMS\) :\s+(.+?))", 1, "refine", { "ls_d_res_low" } },
	/* 3 */ { R"(DATA CUTOFF \(SIGMA\(F\)\) :\s+(.+?))", 1, "refine", { "pdbx_ls_sigma_F" } },
	/* 4 */ { R"(COMPLETENESS FOR RANGE \(%\) :\s+(.+?))", 1, "refine", { "ls_percent_reflns_obs" } },
	/* 5 */ { R"(NUMBER OF REFLECTIONS :\s+(.+?))", 1, "refine", { "ls_number_reflns_obs" } },
	/* 6 */ { R"(FIT TO DATA USED IN REFINEMENT\.)", 1 },
	/* 7 */ { R"(CROSS-VALIDATION METHOD :\s+(.+?))", 1, "refine", { "pdbx_ls_cross_valid_method" } },
	/* 8 */ { R"(FREE R VALUE TEST SET SELECTION :\s+(.+?))", 1, "refine", { "pdbx_R_Free_selection_details" } },
	/* 9 */ { R"(R VALUE \(WORKING ?\+ ?TEST SET\) :\s+(.+?))", 1, "refine", { "ls_R_factor_obs" } },
	/* 10 */ { R"(R VALUE \(WORKING SET\) :\s+(.+?))", 1, "refine", { "ls_R_factor_R_work" } },
	/* 11 */ { R"(FREE R VALUE :\s+(.+?))", 1, "refine", { "ls_R_factor_R_free" } },
	/* 12 */ { R"(FREE R VALUE TEST SET SIZE \(%\) :\s+(.+?))", 1, "refine", { "ls_percent_reflns_R_free" } },
	/* 13 */ { R"(FREE R VALUE TEST SET COUNT :\s+(.+?))", 1, "refine", { "ls_number_reflns_R_free" } },
	/* 14 */ { R"(ESTIMATED ERROR OF FREE R VALUE :\s+(.+?))", 1, "refine", { "ls_R_factor_R_free_error" } },
	/* 15 */ { R"(FIT IN THE HIGHEST RESOLUTION BIN\.)", 1 },
	/* 16 */ { R"(TOTAL NUMBER OF BINS USED :\s+(.+?))", 1, "refine_ls_shell", { "pdbx_total_number_of_bins_used" } },
	/* 17 */ { R"(BIN RESOLUTION RANGE HIGH \(A(?:NGSTROMS)?\) :\s+(.+?))", 1, "refine_ls_shell", { "d_res_high" } },
	/* 18 */ { R"(BIN RESOLUTION RANGE LOW \(A(?:NGSTROMS)?\) :\s+(.+?))", 1, "refine_ls_shell", { "d_res_low" } },
	/* 19 */ { R"(BIN COMPLETENESS \(WORKING\+TEST\) \(%\) :\s+(.+?))", 1, "refine_ls_shell", { "percent_reflns_obs" } },
	/* 20 */ { R"(REFLECTIONS IN BIN \(WORKING ?\+ ?TEST(?: SET)?\) :\s+(.+?))", 1, "refine_ls_shell", { "number_reflns_all" } },
	/* 21 */ { R"(BIN R VALUE \(WORKING ?\+ ?TEST(?: SET)?\) :\s+(.+?))", 1, "refine_ls_shell", { "R_factor_all" } },
	/* 22 */ { R"(REFLECTIONS IN BIN \(WORKING SET\) :\s+(.+?))", 1, "refine_ls_shell", { "number_reflns_R_work" } },
	/* 23 */ { R"(BIN R VALUE \(WORKING SET\) :\s+(.+?))", 1, "refine_ls_shell", { "R_factor_R_work" } },
	/* 24 */ { R"(BIN FREE R VALUE :\s+(.+?))", 1, "refine_ls_shell", { "R_factor_R_free" } },
	/* 25 */ { R"(BIN FREE R VALUE TEST SET SIZE \(%\) :\s+(.+?))", 1, "refine_ls_shell", { "percent_reflns_R_free" } },
	/* 26 */ { R"(BIN FREE R VALUE TEST SET COUNT :\s+(.+?))", 1, "refine_ls_shell", { "number_reflns_R_free" } },
	/* 27 */ { R"(ESTIMATED ERROR OF BIN FREE R VALUE :\s+(.+?))", 1, "refine_ls_shell", { "R_factor_R_free_error" } },
	/* 28 */ { R"(NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT\.)", 1 },
	/* 29 */ { R"(PROTEIN ATOMS :\s+(.+?))", 1, "refine_hist", { "pdbx_number_atoms_protein" } },
	/* 30 */ { R"(NUCLEIC ACID ATOMS :\s+(.+?))", 1, "refine_hist", { "pdbx_number_atoms_nucleic_acid" } },
	/* 31 */ { R"(HETEROGEN ATOMS :\s+(.+?))", 1, "refine_hist", { "pdbx_number_atoms_ligand" } },
	/* 32 */ { R"(SOLVENT ATOMS :\s+(.+?))", 1, "refine_hist", { "number_atoms_solvent" } },
	/* 33 */ { R"(B VALUES\.)", 1 },
	/* 34 */ { R"(B VALUE TYPE :\s+(.+?))", 1, "refine", { "pdbx_TLS_residual_ADP_flag" } },
	/* 35 */ { R"(FROM WILSON PLOT \(A\*\*2\) :\s+(.+?))", 1, "reflns", { "B_iso_Wilson_estimate" } },
	/* 36 */ { R"(MEAN B VALUE \(OVERALL, A\*\*2\) :\s+(.+?))", 1, "refine", { "B_iso_mean" } },
	/* 37 */ { R"(OVERALL ANISOTROPIC B VALUE\.)", 1 },
	/* 38 */ { R"(B11 \(A\*\*2\) :\s+(.+?))", 1, "refine", { "aniso_B[1][1]" } },
	/* 39 */ { R"(B22 \(A\*\*2\) :\s+(.+?))", 1, "refine", { "aniso_B[2][2]" } },
	/* 40 */ { R"(B33 \(A\*\*2\) :\s+(.+?))", 1, "refine", { "aniso_B[3][3]" } },
	/* 41 */ { R"(B12 \(A\*\*2\) :\s+(.+?))", 1, "refine", { "aniso_B[1][2]" } },
	/* 42 */ { R"(B13 \(A\*\*2\) :\s+(.+?))", 1, "refine", { "aniso_B[1][3]" } },
	/* 43 */ { R"(B23 \(A\*\*2\) :\s+(.+?))", 1, "refine", { "aniso_B[2][3]" } },
	/* 44 */ { R"(ESTIMATED COORDINATE ERROR\.)", 1 },
	/* 45 */ { R"(ESD FROM LUZZATI PLOT \(A\) :\s+(.+?))", 1, "refine_analyze", { "Luzzati_coordinate_error_obs" } },
	/* 46 */ { R"(DPI \(BLOW EQ-10\) BASED ON R VALUE \(A\) :\s+(.+?))", 1, "refine", { "pdbx_overall_SU_R_Blow_DPI" } },
	/* 47 */ { R"(DPI \(BLOW EQ-9\) BASED ON FREE R VALUE \(A\) :\s+(.+?))", 1, "refine", { "pdbx_overall_SU_R_free_Blow_DPI" } },
	/* 48 */ { R"(DPI \(CRUICKSHANK\) BASED ON R VALUE \(A\) :\s+(.+?))", 1, "refine", { "overall_SU_R_Cruickshank_DPI" } },
	/* 49 */ { R"(DPI \(CRUICKSHANK\) BASED ON FREE R VALUE \(A\) :\s+(.+?))", 1, "refine", { "pdbx_overall_SU_R_free_Cruickshank_DPI" } },
	/* 50 */ { R"(REFERENCES: BLOW.+)", 1 },
	/* 51 */ { R"(CORRELATION COEFFICIENTS\.)", 1 },
	/* 52 */ { R"(CORRELATION COEFFICIENT FO-FC :\s+(.+?))", 1, "refine", { "correlation_coeff_Fo_to_Fc" } },
	/* 53 */ { R"(CORRELATION COEFFICIENT FO-FC FREE :\s+(.+?))", 1, "refine", { "correlation_coeff_Fo_to_Fc_free" } },
	/* 54 */ { R"(NUMBER OF GEOMETRIC FUNCTION TERMS DEFINED : 15)", 1 },
	/* 55 */ { R"(TERM COUNT WEIGHT FUNCTION\.)", 1 },
	/* 56 */ { R"(BOND LENGTHS :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_bond_d", true },
	/* 57 */ { R"(BOND ANGLES :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_angle_deg", true },
	/* 58 */ { R"(TORSION ANGLES :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_dihedral_angle_d", true },
	/* 59 */ { R"(TRIGONAL CARBON PLANES :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_trig_c_planes", true },
	/* 60 */ { R"(GENERAL PLANES :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_gen_planes", true },
	/* 61 */ { R"(ISOTROPIC THERMAL FACTORS :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_it", true },
	/* 62 */ { R"(BAD NON-BONDED CONTACTS :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_nbd", true },
	/* 63 */ { R"(IMPROPER TORSIONS :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_improper_torsion", true },
	/* 64 */ { R"(PSEUDOROTATION ANGLES :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_pseud_angle", true },
	/* 65 */ { R"(CHIRAL IMPROPER TORSION :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_chiral_improper_torsion", true },
	/* 66 */ { R"(SUM OF OCCUPANCIES :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_sum_occupancies", true },
	/* 67 */ { R"(UTILITY DISTANCES :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_utility_distance", true },
	/* 68 */ { R"(UTILITY ANGLES :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_utility_angle", true },
	/* 69 */ { R"(UTILITY TORSION :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_utility_torsion", true },
	/* 70 */ { R"(IDEAL-DIST CONTACT TERM :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "number", "weight", "pdbx_restraint_function" }, "t_ideal_dist_contact", true },
	/* 71 */ { R"(RMS DEVIATIONS FROM IDEAL VALUES\.)", 1 },
	/* 72 */ { R"(BOND LENGTHS \(A\) :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_bond_d", false },
	/* 73 */ { R"(BOND ANGLES \(DEGREES\) :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_angle_deg", false },
	/* 74 */ { R"(TORSION ANGLES \(DEGREES\) :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_dihedral_angle_d", false },
	/* 75 */ { R"(PSEUDO ROTATION ANGLES \(DEGREES\) :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_pseud_angle", false },
	/* 76 */ { R"(TRIGONAL CARBON PLANES \(A\) :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_trig_c_planes", false },
	/* 77 */ { R"(GENERAL PLANES \(A\) :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_gen_planes", false },
	/* 78 */ { R"(ISOTROPIC THERMAL FACTORS \(A\*\*2\) :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_it", false },
	/* 79 */ { R"(NON-BONDED CONTACTS \(A\) :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_nbd", false },
	/* 80 */ { R"(PEPTIDE OMEGA TORSION ANGLES \(DEGREES\) :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_omega_torsion", false },
	/* 81 */ { R"(OTHER TORSION ANGLES \(DEGREES\) :\s+(.+?);\s+(.+?);\s+(.+?))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_other_torsion", false },
	/* 82 */ { R"(TLS DETAILS\.?)", 1 },
	/* 83 */ { R"(NUMBER OF TLS GROUPS :.+)", 1 },
	/* 84 */ { R"(TLS GROUP :\s*(\d+))", 1, "pdbx_refine_tls", { "id" }, nullptr, true },
	/* 85 */ { R"((?:SELECTION|SET) *:\s+(.+?))", 1, "pdbx_refine_tls_group", { "selection_details" }, nullptr, true },
	/* 86 */ { R"(ORIGIN FOR THE GROUP \(A\):\s+(.+?)\s+(.+?)\s+(.+?))", 1, "pdbx_refine_tls", { "origin_x", "origin_y", "origin_z" } },
	/* 87 */ { R"(T TENSOR)", 1 },
	/* 88 */ { R"(T11:\s+(.+?) T22:\s+(.+?))", 1, "pdbx_refine_tls", { "T[1][1]", "T[2][2]" } },
	/* 89 */ { R"(T33:\s+(.+?) T12:\s+(.+?))", 1, "pdbx_refine_tls", { "T[3][3]", "T[1][2]" } },
	/* 90 */ { R"(T13:\s+(.+?) T23:\s+(.+?))", 1, "pdbx_refine_tls", { "T[1][3]", "T[2][3]" } },
	/* 91 */ { R"(L TENSOR)", 1 },
	/* 92 */ { R"(L11:\s+(.+?) L22:\s+(.+?))", 1, "pdbx_refine_tls", { "L[1][1]", "L[2][2]" } },
	/* 93 */ { R"(L33:\s+(.+?) L12:\s+(.+?))", 1, "pdbx_refine_tls", { "L[3][3]", "L[1][2]" } },
	/* 94 */ { R"(L13:\s+(.+?) L23:\s+(.+?))", 1, "pdbx_refine_tls", { "L[1][3]", "L[2][3]" } },
	/* 95 */ { R"(S TENSOR)", 1 },
	/* 96 */ { R"(S11:\s+(.+?) S12:\s+(.+?) S13:\s+(.+?))", 1, "pdbx_refine_tls", { "S[1][1]", "S[1][2]", "S[1][3]" } },
	/* 97 */ { R"(S21:\s+(.+?) S22:\s+(.+?) S23:\s+(.+?))", 1, "pdbx_refine_tls", { "S[2][1]", "S[2][2]", "S[2][3]" } },
	/* 98 */ { R"(S31:\s+(.+?) S32:\s+(.+?) S33:\s+(.+?))", 84 - 98, "pdbx_refine_tls", { "S[3][1]", "S[3][2]", "S[3][3]" } },
};

class BUSTER_TNT_Remark3Parser : public Remark3Parser
{
  public:
	BUSTER_TNT_Remark3Parser(const std::string &name, const std::string &expMethod, PDBRecord *r, cif::datablock &db)
		: Remark3Parser(name, expMethod, r, db,
			  kBusterTNT_Template, sizeof(kBusterTNT_Template) / sizeof(TemplateLine),
			  std::regex(R"((BUSTER(?:-TNT)?)(?: (\d+(?:\..+)?))?)"))
	{
	}
};

const TemplateLine kCNS_Template[] = {
	/* 0 */ { R"(REFINEMENT TARGET\s*:\s*(.+))", 1, "refine", { "pdbx_stereochemistry_target_values" } },
	/* 1 */ { R"(DATA USED IN REFINEMENT\.)", 1 },
	/* 2 */ { R"(RESOLUTION RANGE HIGH \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_high" } },
	/* 3 */ { R"(RESOLUTION RANGE LOW \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_low" } },
	/* 4 */ { R"(DATA CUTOFF \(SIGMA\(F\)\)\s*:\s*(.+))", 1, "refine", { "pdbx_ls_sigma_F" } },
	/* 5 */ { R"(DATA CUTOFF HIGH \(ABS\(F\)\)\s*:\s*(.+))", 1, "refine", { "pdbx_data_cutoff_high_absF" } },
	/* 6 */ { R"(DATA CUTOFF LOW \(ABS\(F\)\)\s*:\s*(.+))", 1, "refine", { "pdbx_data_cutoff_low_absF" } },
	/* 7 */ { R"(COMPLETENESS \(WORKING\+TEST\) \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_obs" } },
	/* 8 */ { R"(NUMBER OF REFLECTIONS\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_obs" } },
	/* 9 */ { R"(FIT TO DATA USED IN REFINEMENT\.)", 1 },
	/* 10 */ { R"(CROSS-VALIDATION METHOD\s*:\s*(.+))", 1, "refine", { "pdbx_ls_cross_valid_method" } },
	/* 11 */ { R"(FREE R VALUE TEST SET SELECTION\s*:\s*(.+))", 1, "refine", { "pdbx_R_Free_selection_details" } },
	/* 12 */ { R"(R VALUE \(WORKING \+ TEST SET\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_obs" } },
	/* 13 */ { R"(R VALUE \(WORKING SET\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_work" } },
	/* 14 */ { R"(FREE R VALUE\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_free" } },
	/* 15 */ { R"(FREE R VALUE TEST SET SIZE \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_R_free" } },
	/* 16 */ { R"(FREE R VALUE TEST SET COUNT\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_R_free" } },
	/* 17 */ { R"(ESTIMATED ERROR OF FREE R VALUE\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_free_error" } },
	/* 18 */ { R"(FIT/AGREEMENT OF MODEL WITH ALL DATA\.)", 1 },
	/* 19 */ { R"(R VALUE \(WORKING \+ TEST SET, NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "R_factor_all_no_cutoff" } },
	/* 20 */ { R"(R VALUE \(WORKING SET, NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "R_factor_obs_no_cutoff" } },
	/* 21 */ { R"(FREE R VALUE \(NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_factor_no_cutoff" } },
	/* 22 */ { R"(FREE R VALUE TEST SET SIZE \(%, NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_val_test_set_size_perc_no_cutoff" } },
	/* 23 */ { R"(FREE R VALUE TEST SET COUNT \(NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_val_test_set_ct_no_cutoff" } },
	/* 24 */ { R"(ESTIMATED ERROR OF FREE R VALUE \(NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_error_no_cutoff" } },
	/* 25 */ { R"(TOTAL NUMBER OF REFLECTIONS \(NO CUTOFF\)\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_all" } },
	/* 26 */ { R"(FIT IN THE HIGHEST RESOLUTION BIN\.)", 1 },
	/* 27 */ { R"(TOTAL NUMBER OF BINS USED\s*:\s*(.+))", 1, "refine_ls_shell", { "pdbx_total_number_of_bins_used" } },
	/* 28 */ { R"(BIN RESOLUTION RANGE HIGH \(A\)\s*:\s*(.+))", 1, "refine_ls_shell", { "d_res_high" } },
	/* 29 */ { R"(BIN RESOLUTION RANGE LOW \(A\)\s*:\s*(.+))", 1, "refine_ls_shell", { "d_res_low" } },
	/* 30 */ { R"(BIN COMPLETENESS \(WORKING\+TEST\) \(%\)\s*:\s*(.+))", 1, "refine_ls_shell", { "percent_reflns_obs" } },
	/* 31 */ { R"(REFLECTIONS IN BIN \(WORKING SET\)\s*:\s*(.+))", 1, "refine_ls_shell", { "number_reflns_R_work" } },
	/* 32 */ { R"(BIN R VALUE \(WORKING SET\)\s*:\s*(.+))", 1, "refine_ls_shell", { "R_factor_R_work" } },
	/* 33 */ { R"(BIN FREE R VALUE\s*:\s*(.+))", 1, "refine_ls_shell", { "R_factor_R_free" } },
	/* 34 */ { R"(BIN FREE R VALUE TEST SET SIZE \(%\)\s*:\s*(.+))", 1, "refine_ls_shell", { "percent_reflns_R_free" } },
	/* 35 */ { R"(BIN FREE R VALUE TEST SET COUNT\s*:\s*(.+))", 1, "refine_ls_shell", { "number_reflns_R_free" } },
	/* 36 */ { R"(ESTIMATED ERROR OF BIN FREE R VALUE\s*:\s*(.+))", 1, "refine_ls_shell", { "R_factor_R_free_error" } },
	/* 37 */ { R"(NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT\.)", 1 },
	/* 38 */ { R"(PROTEIN ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_protein" } },
	/* 39 */ { R"(NUCLEIC ACID ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_nucleic_acid" } },
	/* 40 */ { R"(HETEROGEN ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_ligand" } },
	/* 41 */ { R"(SOLVENT ATOMS\s*:\s*(.+))", 1, "refine_hist", { "number_atoms_solvent" } },
	/* 42 */ { R"(B VALUES\.)", 1 },
	/* 43 */ { R"(B VALUE TYPE\s*:\s*(.+))", 1, "refine", { "pdbx_TLS_residual_ADP_flag" } },
	/* 44 */ { R"(FROM WILSON PLOT \(A\*\*2\)\s*:\s*(.+))", 1, "reflns", { "B_iso_Wilson_estimate" } },
	/* 45 */ { R"(MEAN B VALUE \(OVERALL, A\*\*2\)\s*:\s*(.+))", 1, "refine", { "B_iso_mean" } },
	/* 46 */ { R"(OVERALL ANISOTROPIC B VALUE\.)", 1 },
	/* 47 */ { R"(B11 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][1]" } },
	/* 48 */ { R"(B22 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[2][2]" } },
	/* 49 */ { R"(B33 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[3][3]" } },
	/* 50 */ { R"(B12 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][2]" } },
	/* 51 */ { R"(B13 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][3]" } },
	/* 52 */ { R"(B23 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[2][3]" } },
	/* 53 */ { R"(ESTIMATED COORDINATE ERROR\.)", 1 },
	/* 54 */ { R"(ESD FROM LUZZATI PLOT \(A\)\s*:\s*(.+))", 1, "refine_analyze", { "Luzzati_coordinate_error_obs" } },
	/* 55 */ { R"(ESD FROM SIGMAA \(A\)\s*:\s*(.+))", 1, "refine_analyze", { "Luzzati_sigma_a_obs" } },
	/* 56 */ { R"(LOW RESOLUTION CUTOFF \(A\)\s*:\s*(.+))", 1, "refine_analyze", { "Luzzati_d_res_low_obs" } },
	/* 57 */ { R"(CROSS-VALIDATED ESTIMATED COORDINATE ERROR\.)", 1 },
	/* 58 */ { R"(ESD FROM C-V LUZZATI PLOT \(A\)\s*:\s*(.+))", 1, "refine_analyze", { "Luzzati_coordinate_error_free" } },
	/* 59 */ { R"(ESD FROM C-V SIGMAA \(A\)\s*:\s*(.+))", 1, "refine_analyze", { "Luzzati_sigma_a_free" } },
	/* 60 */ { R"(RMS DEVIATIONS FROM IDEAL VALUES\.)", 1 },
	/* 61 */ { R"(BOND LENGTHS \(A\)\s*:\s*(.+))", 1, "refine_ls_restr", { "dev_ideal" }, "c_bond_d", false },
	/* 62 */ { R"(BOND ANGLES \(DEGREES\)\s*:\s*(.+))", 1, "refine_ls_restr", { "dev_ideal" }, "c_angle_deg", false },
	/* 63 */ { R"(DIHEDRAL ANGLES \(DEGREES\)\s*:\s*(.+))", 1, "refine_ls_restr", { "dev_ideal" }, "c_dihedral_angle_d", false },
	/* 64 */ { R"(IMPROPER ANGLES \(DEGREES\)\s*:\s*(.+))", 1, "refine_ls_restr", { "dev_ideal" }, "c_improper_angle_d", false },
	/* 65 */ { R"(ISOTROPIC THERMAL MODEL\s*:\s*(.+))", 1, "refine", { "pdbx_isotropic_thermal_model" } },
	/* 66 */ { R"(ISOTROPIC THERMAL FACTOR RESTRAINTS\. RMS SIGMA)", 1 },
	/* 67 */ { R"(MAIN-CHAIN BOND \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "c_mcbond_it", false },
	/* 68 */ { R"(MAIN-CHAIN ANGLE \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "c_mcangle_it", false },
	/* 69 */ { R"(SIDE-CHAIN BOND \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "c_scbond_it", false },
	/* 70 */ { R"(SIDE-CHAIN ANGLE \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "c_scangle_it", false },
	/* 71 */ { R"(BULK SOLVENT MODELING\.)", 1 },
	/* 72 */ { R"(METHOD USED\s*:\s*(.+))", 1, "refine", { "solvent_model_details" } },
	/* 73 */ { R"(KSOL\s*:\s*(.+))", 1, "refine", { "solvent_model_param_ksol" } },
	/* 74 */ { R"(BSOL\s*:\s*(.+))", 1, "refine", { "solvent_model_param_bsol" } },
	/* 75 */ { R"(NCS MODEL\s*:\s*(.+))", 1, /* "refine_ls_restr_ncs", { "ncs_model_details" } */ },
	/* 76 */ { R"(NCS RESTRAINTS\. RMS SIGMA/WEIGHT)", 1 },
	/* 77 */ { R"(GROUP (\d+) POSITIONAL \(A\)\s*:\s*(.+))", 1, /* "refine_ls_restr_ncs", { "dom_id", "rms_dev_position", "weight_position" } */ },
	/* 78 */ { R"(GROUP (\d+) B-FACTOR \(A\*\*2\)\s*:\s*(.+))", 1, /* "refine_ls_restr_ncs", { "dom_id", "rms_dev_B_iso", "weight_B_iso" } */ },
	/* 79 */ { R"(PARAMETER FILE (\d+) :\s+(.+))", 1, /* "pdbx_xplor_file", { "serial_no", "param_file" } */ },
	/* 80 */ { R"(TOPOLOGY FILE (\d+) :\s+(.+))", 1, /* "pdbx_xplor_file", { "serial_no", "topol_file" } */ },
};

class CNS_Remark3Parser : public Remark3Parser
{
  public:
	CNS_Remark3Parser(const std::string &name, const std::string &expMethod, PDBRecord *r, cif::datablock &db)
		: Remark3Parser(name, expMethod, r, db, kCNS_Template,
			  sizeof(kCNS_Template) / sizeof(TemplateLine), std::regex(R"((CN[SX])(?: (\d+(?:\.\d+)?))?)"))
	{
	}
};

const TemplateLine kPHENIX_Template[] = {
	/* 0 */ { R"(REFINEMENT TARGET\s*:\s*(.+))", 1, "refine", { "pdbx_stereochemistry_target_values" } },
	/* 1 */ { R"(DATA USED IN REFINEMENT\.)", 1 },
	/* 2 */ { R"(RESOLUTION RANGE HIGH \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_high" } },
	/* 3 */ { R"(RESOLUTION RANGE LOW \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_low" } },
	/* 4 */ { R"(MIN\(FOBS/SIGMA_FOBS\)\s*:\s*(.+))", 1, "refine", { "pdbx_ls_sigma_F" } },
	/* 5 */ { R"(COMPLETENESS FOR RANGE \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_obs" } },
	/* 6 */ { R"(NUMBER OF REFLECTIONS\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_obs" } },
	/* 7 */ { R"(FIT TO DATA USED IN REFINEMENT\.)", 1 },
	/* 8 */ { R"(R VALUE \(WORKING \+ TEST SET\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_obs" } },
	/* 9 */ { R"(R VALUE \(WORKING SET\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_work" } },
	/* 10 */ { R"(FREE R VALUE\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_free" } },
	/* 11 */ { R"(FREE R VALUE TEST SET SIZE \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_R_free" } },
	/* 12 */ { R"(FREE R VALUE TEST SET COUNT\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_R_free" } },
	/* 13 */ { R"(FIT TO DATA USED IN REFINEMENT \(IN BINS\)\.)", 1 },
	/* 14 */ { R"(BIN RESOLUTION RANGE COMPL\. NWORK NFREE RWORK RFREE)", 1 },
	/* 15 */ { R"(\d+ (\d+(?:\.\d+)?) - (\d+(?:\.\d+)?) (\d+(?:\.\d+)?) (\d+) (\d+) (\d+(?:\.\d+)?) (\d+(?:\.\d+)?))", 0, "refine_ls_shell", { "d_res_low", "d_res_high", "percent_reflns_obs", "number_reflns_R_work", "number_reflns_R_free", "R_factor_R_work", "R_factor_R_free" }, nullptr, true },
	/* 16 */ { R"(BULK SOLVENT MODELLING\.)", 1 },
	/* 17 */ { R"(METHOD USED\s*:\s*(.+))", 1, "refine", { "solvent_model_details" } },
	/* 18 */ { R"(SOLVENT RADIUS\s*:\s*(.+))", 1, "refine", { "pdbx_solvent_vdw_probe_radii" } },
	/* 19 */ { R"(SHRINKAGE RADIUS\s*:\s*(.+))", 1, "refine", { "pdbx_solvent_shrinkage_radii" } },
	/* 20 */ { R"(K_SOL\s*:\s*(.+))", 1, "refine", { "solvent_model_param_ksol" } },
	/* 21 */ { R"(B_SOL\s*:\s*(.+))", 1, "refine", { "solvent_model_param_bsol" } },
	/* 22 */ { R"(ERROR ESTIMATES\.)", 1 },
	/* 23 */ { R"(COORDINATE ERROR \(MAXIMUM-LIKELIHOOD BASED\)\s*:\s*(.+))", 1, "refine", { "overall_SU_ML" } },
	/* 24 */ { R"(PHASE ERROR \(DEGREES, MAXIMUM-LIKELIHOOD BASED\)\s*:\s*(.+))", 1, "refine", { "pdbx_overall_phase_error" } },
	/* 25 */ { R"(B VALUES\.)", 1 },
	/* 26 */ { R"(B VALUE TYPE\s*:\s*(.+))", 1, "refine", { "pdbx_TLS_residual_ADP_flag" } },
	/* 27 */ { R"(FROM WILSON PLOT \(A\*\*2\)\s*:\s*(.+))", 1, "reflns", { "B_iso_Wilson_estimate" } },
	/* 28 */ { R"(MEAN B VALUE \(OVERALL, A\*\*2\)\s*:\s*(.+))", 1, "refine", { "B_iso_mean" } },
	/* 29 */ { R"(OVERALL ANISOTROPIC B VALUE\.)", 1 },
	/* 30 */ { R"(B11 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][1]" } },
	/* 31 */ { R"(B22 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[2][2]" } },
	/* 32 */ { R"(B33 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[3][3]" } },
	/* 33 */ { R"(B12 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][2]" } },
	/* 34 */ { R"(B13 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][3]" } },
	/* 35 */ { R"(B23 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[2][3]" } },
	/* 36 */ { R"(TWINNING INFORMATION\.)", 1 },
	/* 37 */ { R"(FRACTION:\s*(.+))", 1, "pdbx_reflns_twin", { "fraction" } },
	/* 38 */ { R"(OPERATOR:\s*(.+))", 1, "pdbx_reflns_twin", { "operator" } },
	/* 39 */ { R"(DEVIATIONS FROM IDEAL VALUES\.)", 1 },
	/* 40 */ { R"(RMSD COUNT)", 1 },
	/* 41 */ { R"(BOND\s*:\s*(\d+(?:\.\d+))\s+(\d+))", 1, "refine_ls_restr", { "dev_ideal", "number" }, "f_bond_d", false },
	/* 42 */ { R"(ANGLE\s*:\s*(\d+(?:\.\d+))\s+(\d+))", 1, "refine_ls_restr", { "dev_ideal", "number" }, "f_angle_d", false },
	/* 43 */ { R"(CHIRALITY\s*:\s*(\d+(?:\.\d+))\s+(\d+))", 1, "refine_ls_restr", { "dev_ideal", "number" }, "f_chiral_restr", false },
	/* 44 */ { R"(PLANARITY\s*:\s*(\d+(?:\.\d+))\s+(\d+))", 1, "refine_ls_restr", { "dev_ideal", "number" }, "f_plane_restr", false },
	/* 45 */ { R"(DIHEDRAL\s*:\s*(\d+(?:\.\d+))\s+(\d+))", 1, "refine_ls_restr", { "dev_ideal", "number" }, "f_dihedral_angle_d", false },
	/* 46 */ { R"(TLS DETAILS)", 1 },
	/* 47 */ { R"(NUMBER OF TLS GROUPS\s*:\s*(.+))", 1 },
	/* 48 */ { R"(TLS GROUP\s*:\s*(.+))", 1, "pdbx_refine_tls", { "id" }, nullptr, true },
	/* 49 */ { R"(SELECTION:\s*(.+))", 1, "pdbx_refine_tls_group", { "selection_details" }, nullptr, true },
	/* 50 */ { R"(ORIGIN FOR THE GROUP(?:\s*\(A\))?\s*:\s*(\S+)\s+(\S+)\s+(\S+))", 1, "pdbx_refine_tls", { "origin_x", "origin_y", "origin_z" } },
	/* 51 */ { R"(T TENSOR)", 1 },
	/* 52 */ { R"(T11\s*:\s*(.+) T22\s*:\s*(.+))", 1, "pdbx_refine_tls", { "T[1][1]", "T[2][2]" } },
	/* 53 */ { R"(T33\s*:\s*(.+) T12\s*:\s*(.+))", 1, "pdbx_refine_tls", { "T[3][3]", "T[1][2]" } },
	/* 54 */ { R"(T13\s*:\s*(.+) T23\s*:\s*(.+))", 1, "pdbx_refine_tls", { "T[1][3]", "T[2][3]" } },
	/* 55 */ { R"(L TENSOR)", 1 },
	/* 56 */ { R"(L11\s*:\s*(.+) L22\s*:\s*(.+))", 1, "pdbx_refine_tls", { "L[1][1]", "L[2][2]" } },
	/* 57 */ { R"(L33\s*:\s*(.+) L12\s*:\s*(.+))", 1, "pdbx_refine_tls", { "L[3][3]", "L[1][2]" } },
	/* 58 */ { R"(L13\s*:\s*(.+) L23\s*:\s*(.+))", 1, "pdbx_refine_tls", { "L[1][3]", "L[2][3]" } },
	/* 59 */ { R"(S TENSOR)", 1 },
	/* 60 */ { R"(S11\s*:\s*(.+) S12\s*:\s*(.+) S13\s*:\s*(.+))", 1, "pdbx_refine_tls", { "S[1][1]", "S[1][2]", "S[1][3]" } },
	/* 61 */ { R"(S21\s*:\s*(.+) S22\s*:\s*(.+) S23\s*:\s*(.+))", 1, "pdbx_refine_tls", { "S[2][1]", "S[2][2]", "S[2][3]" } },
	/* 62 */ { R"(S31\s*:\s*(.+) S32\s*:\s*(.+) S33\s*:\s*(.+))", 48 - 62, "pdbx_refine_tls", { "S[3][1]", "S[3][2]", "S[3][3]" } },
	/* 63 */ { R"(ANOMALOUS SCATTERER GROUPS DETAILS\.)", 1 },
	/* 64 */ { R"(NUMBER OF ANOMALOUS SCATTERER GROUPS\s*:\s*\d+)", 1 },
	/* 65 */ { R"(ANOMALOUS SCATTERER GROUP\s*:\s*\d+)", 1 },
	/* 66 */ { R"(SELECTION: .+)", 1 },
	/* 67 */ { R"(fp\s*:\s*.+)", 1 },
	/* 68 */ { R"(fdp\s*:\s*.+)", 63 - 68 },
	/* 69 */ { R"(NCS DETAILS)", 1 },
	/* 70 */ { R"(NUMBER OF NCS GROUPS\s*:\s*(.+))", 1 },
};

class PHENIX_Remark3Parser : public Remark3Parser
{
  public:
	PHENIX_Remark3Parser(const std::string &name, const std::string &expMethod, PDBRecord *r, cif::datablock &db)
		: Remark3Parser(name, expMethod, r, db, kPHENIX_Template, sizeof(kPHENIX_Template) / sizeof(TemplateLine),
			  std::regex(R"((PHENIX)(?: \(PHENIX\.REFINE:) (\d+(?:\.[^)]+)?)\)?)"))
	{
	}

	virtual void fixup();
};

void PHENIX_Remark3Parser::fixup()
{
	for (auto r : mDb["refine_ls_shell"])
	{
		try
		{
			float val = r["percent_reflns_obs"].as<float>();
			int perc = static_cast<int>(val * 100);
			r["percent_reflns_obs"] = perc;
		}
		catch (...)
		{
		}
	}
}

const TemplateLine kNUCLSQ_Template[] = {
	/* 0 */ { R"(DATA USED IN REFINEMENT\.)", 1 },
	/* 1 */ { R"(RESOLUTION RANGE HIGH \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_high" } },
	/* 2 */ { R"(RESOLUTION RANGE LOW \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_low" } },
	/* 3 */ { R"(DATA CUTOFF \(SIGMA\(F\)\)\s*:\s*(.+))", 1, "refine", { "pdbx_ls_sigma_F" } },
	/* 4 */ { R"(COMPLETENESS FOR RANGE \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_obs" } },
	/* 5 */ { R"(NUMBER OF REFLECTIONS\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_obs" } },
	/* 6 */ { R"(FIT TO DATA USED IN REFINEMENT\.)", 1 },
	/* 7 */ { R"(CROSS-VALIDATION METHOD\s*:\s*(.+))", 1, "refine", { "pdbx_ls_cross_valid_method" } },
	/* 8 */ { R"(FREE R VALUE TEST SET SELECTION\s*:\s*(.+))", 1, "refine", { "pdbx_R_Free_selection_details" } },
	/* 9 */ { R"(R VALUE \(WORKING \+ TEST SET\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_obs" } },
	/* 10 */ { R"(R VALUE \(WORKING SET\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_work" } },
	/* 11 */ { R"(FREE R VALUE\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_free" } },
	/* 12 */ { R"(FREE R VALUE TEST SET SIZE \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_R_free" } },
	/* 13 */ { R"(FREE R VALUE TEST SET COUNT\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_R_free" } },
	/* 14 */ { R"(FIT/AGREEMENT OF MODEL WITH ALL DATA\.)", 1 },
	/* 15 */ { R"(R VALUE \(WORKING \+ TEST SET, NO CUTOFF\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_all" } },
	/* 16 */ { R"(R VALUE \(WORKING SET, NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "R_factor_obs_no_cutoff" } },
	/* 17 */ { R"(FREE R VALUE \(NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_factor_no_cutoff" } },
	/* 18 */ { R"(FREE R VALUE TEST SET SIZE \(%, NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_val_test_set_size_perc_no_cutoff" } },
	/* 19 */ { R"(FREE R VALUE TEST SET COUNT \(NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_val_test_set_ct_no_cutoff" } },
	/* 20 */ { R"(TOTAL NUMBER OF REFLECTIONS \(NO CUTOFF\)\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_all" } },
	/* 21 */ { R"(NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT\.)", 1 },
	/* 22 */ { R"(PROTEIN ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_protein" } },
	/* 23 */ { R"(NUCLEIC ACID ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_nucleic_acid" } },
	/* 24 */ { R"(HETEROGEN ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_ligand" } },
	/* 25 */ { R"(SOLVENT ATOMS\s*:\s*(.+))", 1, "refine_hist", { "number_atoms_solvent" } },
	/* 26 */ { R"(B VALUES\.)", 1 },
	/* 27 */ { R"(B VALUE TYPE\s*:\s*(.+))", 1, "refine", { "pdbx_TLS_residual_ADP_flag" } },
	/* 28 */ { R"(FROM WILSON PLOT \(A\*\*2\)\s*:\s*(.+))", 1, "reflns", { "B_iso_Wilson_estimate" } },
	/* 29 */ { R"(MEAN B VALUE \(OVERALL, A\*\*2\)\s*:\s*(.+))", 1, "refine", { "B_iso_mean" } },
	/* 30 */ { R"(OVERALL ANISOTROPIC B VALUE\.)", 1 },
	/* 31 */ { R"(B11 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][1]" } },
	/* 32 */ { R"(B22 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[2][2]" } },
	/* 33 */ { R"(B33 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[3][3]" } },
	/* 34 */ { R"(B12 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][2]" } },
	/* 35 */ { R"(B13 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][3]" } },
	/* 36 */ { R"(B23 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[2][3]" } },
	/* 37 */ { R"(ESTIMATED COORDINATE ERROR\.)", 1 },
	/* 38 */ { R"(ESD FROM LUZZATI PLOT \(A\)\s*:\s*(.+))", 1, "refine_analyze", { "Luzzati_coordinate_error_obs" } },
	/* 39 */ { R"(ESD FROM SIGMAA \(A\)\s*:\s*(.+))", 1, "refine_analyze", { "Luzzati_sigma_a_obs" } },
	/* 40 */ { R"(LOW RESOLUTION CUTOFF \(A\)\s*:\s*(.+))", 1, "refine_analyze", { "Luzzati_d_res_low_obs" } },
	/* 41 */ { R"(RMS DEVIATIONS FROM IDEAL VALUES\.)", 1 },
	/* 42 */ { R"(DISTANCE RESTRAINTS\. RMS SIGMA)", 1 },
	/* 43 */ { R"(SUGAR-BASE BOND DISTANCE \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "n_sugar_bond_d", false },
	/* 44 */ { R"(SUGAR-BASE BOND ANGLE DISTANCE \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "n_sugar_bond_angle_d", false },
	/* 45 */ { R"(PHOSPHATE BONDS DISTANCE \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "n_phos_bond_d", false },
	/* 46 */ { R"(PHOSPHATE BOND ANGLE, H-BOND \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "n_phos_bond_angle_d", false },
	/* 47 */ { R"(PLANE RESTRAINT \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "n_plane_restr", false },
	/* 48 */ { R"(CHIRAL-CENTER RESTRAINT \(A\*\*3\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "n_chiral_restr", false },
	/* 49 */ { R"(NON-BONDED CONTACT RESTRAINTS\.)", 1 },
	/* 50 */ { R"(SINGLE TORSION \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "n_singtor_nbd", false },
	/* 51 */ { R"(MULTIPLE TORSION \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "n_multtor_nbd", false },
	/* 59 */ { R"(ISOTROPIC THERMAL FACTOR RESTRAINTS\. RMS SIGMA)", 1 },
	/* 60 */ { R"(SUGAR-BASE BONDS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "n_sugar_bond_it", false },
	/* 61 */ { R"(SUGAR-BASE ANGLES \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "n_sugar_angle_it", false },
	/* 62 */ { R"(PHOSPHATE BONDS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "n_phos_bond_it", false },
	/* 63 */ { R"(PHOSPHATE BOND ANGLE, H-BOND \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "n_phos_angle_it", false },
};

class NUCLSQ_Remark3Parser : public Remark3Parser
{
  public:
	NUCLSQ_Remark3Parser(const std::string &name, const std::string &expMethod, PDBRecord *r, cif::datablock &db)
		: Remark3Parser(name, expMethod, r, db, kNUCLSQ_Template, sizeof(kNUCLSQ_Template) / sizeof(TemplateLine),
			  std::regex(R"((NUCLSQ)(?: (\d+(?:\.\d+)?))?)"))
	{
	}

	virtual void fixup()
	{
		for (auto r : mDb["refine_hist"])
		{
			try
			{
				int p, n, h, s;
				cif::tie(p, n, h, s) = r.get("pdbx_number_atoms_protein", "pdbx_number_atoms_nucleic_acid", "pdbx_number_atoms_ligand", "number_atoms_solvent");
				r["number_atoms_total"] = p + n + h + s;
			}
			catch (...)
			{
			}
		}
	}
};

const TemplateLine kPROLSQ_Template[] = {
	/* 0 */ { R"(DATA USED IN REFINEMENT\.)", 1 },
	/* 1 */ { R"(RESOLUTION RANGE HIGH \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_high" } },
	/* 2 */ { R"(RESOLUTION RANGE LOW \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_low" } },
	/* 3 */ { R"(DATA CUTOFF \(SIGMA\(F\)\)\s*:\s*(.+))", 1, "refine", { "pdbx_ls_sigma_F" } },
	/* 4 */ { R"(COMPLETENESS FOR RANGE \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_obs" } },
	/* 5 */ { R"(NUMBER OF REFLECTIONS\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_obs" } },
	/* 6 */ { R"(FIT TO DATA USED IN REFINEMENT\.)", 1 },
	/* 7 */ { R"(CROSS-VALIDATION METHOD\s*:\s*(.+))", 1, "refine", { "pdbx_ls_cross_valid_method" } },
	/* 8 */ { R"(FREE R VALUE TEST SET SELECTION\s*:\s*(.+))", 1, "refine", { "pdbx_R_Free_selection_details" } },
	/* 9 */ { R"(R VALUE \(WORKING \+ TEST SET\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_obs" } },
	/* 10 */ { R"(R VALUE \(WORKING SET\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_work" } },
	/* 11 */ { R"(FREE R VALUE\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_free" } },
	/* 12 */ { R"(FREE R VALUE TEST SET SIZE \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_R_free" } },
	/* 13 */ { R"(FREE R VALUE TEST SET COUNT\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_R_free" } },
	/* 14 */ { R"(FIT/AGREEMENT OF MODEL WITH ALL DATA\.)", 1 },
	/* 15 */ { R"(R VALUE \(WORKING \+ TEST SET, NO CUTOFF\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_all" } },
	/* 16 */ { R"(R VALUE \(WORKING SET, NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "R_factor_obs_no_cutoff" } },
	/* 17 */ { R"(FREE R VALUE \(NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_factor_no_cutoff" } },
	/* 18 */ { R"(FREE R VALUE TEST SET SIZE \(%, NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_val_test_set_size_perc_no_cutoff" } },
	/* 19 */ { R"(FREE R VALUE TEST SET COUNT \(NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_val_test_set_ct_no_cutoff" } },
	/* 20 */ { R"(TOTAL NUMBER OF REFLECTIONS \(NO CUTOFF\)\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_all" } },
	/* 21 */ { R"(NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT\.)", 1 },
	/* 22 */ { R"(PROTEIN ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_protein" } },
	/* 23 */ { R"(NUCLEIC ACID ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_nucleic_acid" } },
	/* 24 */ { R"(HETEROGEN ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_ligand" } },
	/* 25 */ { R"(SOLVENT ATOMS\s*:\s*(.+))", 1, "refine_hist", { "number_atoms_solvent" } },
	/* 26 */ { R"(B VALUES\.)", 1 },
	/* 27 */ { R"(B VALUE TYPE\s*:\s*(.+))", 1, "refine", { "pdbx_TLS_residual_ADP_flag" } },
	/* 28 */ { R"(FROM WILSON PLOT \(A\*\*2\)\s*:\s*(.+))", 1, "reflns", { "B_iso_Wilson_estimate" } },
	/* 29 */ { R"(MEAN B VALUE \(OVERALL, A\*\*2\)\s*:\s*(.+))", 1, "refine", { "B_iso_mean" } },
	/* 30 */ { R"(OVERALL ANISOTROPIC B VALUE\.)", 1 },
	/* 31 */ { R"(B11 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][1]" } },
	/* 32 */ { R"(B22 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[2][2]" } },
	/* 33 */ { R"(B33 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[3][3]" } },
	/* 34 */ { R"(B12 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][2]" } },
	/* 35 */ { R"(B13 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][3]" } },
	/* 36 */ { R"(B23 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[2][3]" } },
	/* 37 */ { R"(ESTIMATED COORDINATE ERROR\.)", 1 },
	/* 38 */ { R"(ESD FROM LUZZATI PLOT \(A\)\s*:\s*(.+))", 1, "refine_analyze", { "Luzzati_coordinate_error_obs" } },
	/* 39 */ { R"(ESD FROM SIGMAA \(A\)\s*:\s*(.+))", 1, "refine_analyze", { "Luzzati_sigma_a_obs" } },
	/* 40 */ { R"(LOW RESOLUTION CUTOFF \(A\)\s*:\s*(.+))", 1, "refine_analyze", { "Luzzati_d_res_low_obs" } },
	/* 41 */ { R"(RMS DEVIATIONS FROM IDEAL VALUES\.)", 1 },
	/* 42 */ { R"(DISTANCE RESTRAINTS\. RMS SIGMA)", 1 },
	/* 43 */ { R"(BOND LENGTH \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_bond_d", false },
	/* 44 */ { R"(ANGLE DISTANCE \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_angle_d", false },
	/* 45 */ { R"(INTRAPLANAR 1-4 DISTANCE \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_planar_d", false },
	/* 46 */ { R"(H-BOND OR METAL COORDINATION \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_hb_or_metal_coord", false },
	/* 47 */ { R"(PLANE RESTRAINT \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_plane_restr", false },
	/* 48 */ { R"(CHIRAL-CENTER RESTRAINT \(A\*\*3\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_chiral_restr", false },
	/* 49 */ { R"(NON-BONDED CONTACT RESTRAINTS\.)", 1 },
	/* 50 */ { R"(SINGLE TORSION \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_singtor_nbd", false },
	/* 51 */ { R"(MULTIPLE TORSION \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_multtor_nbd", false },
	/* 52 */ { R"(H-BOND \(X\.\.\.Y\) \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_xyhbond_nbd", false },
	/* 53 */ { R"(H-BOND \(X-H\.\.\.Y\) \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_xhyhbond_nbd", false },
	/* 54 */ { R"(CONFORMATIONAL TORSION ANGLE RESTRAINTS\.)", 1 },
	/* 55 */ { R"(SPECIFIED \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_special_tor", false },
	/* 56 */ { R"(PLANAR \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_planar_tor", false },
	/* 57 */ { R"(STAGGERED \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_staggered_tor", false },
	/* 58 */ { R"(TRANSVERSE \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_transverse_tor", false },
	/* 59 */ { R"(ISOTROPIC THERMAL FACTOR RESTRAINTS\. RMS SIGMA)", 1 },
	/* 60 */ { R"(MAIN-CHAIN BOND \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_mcbond_it", false },
	/* 61 */ { R"(MAIN-CHAIN ANGLE \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_mcangle_it", false },
	/* 62 */ { R"(SIDE-CHAIN BOND \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_scbond_it", false },
	/* 63 */ { R"(SIDE-CHAIN ANGLE \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_scangle_it", false },
};

class PROLSQ_Remark3Parser : public Remark3Parser
{
  public:
	PROLSQ_Remark3Parser(const std::string &name, const std::string &expMethod, PDBRecord *r, cif::datablock &db)
		: Remark3Parser(name, expMethod, r, db, kPROLSQ_Template, sizeof(kPROLSQ_Template) / sizeof(TemplateLine),
			  std::regex(R"((PROLSQ)(?: (\d+(?:\.\d+)?))?)"))
	{
	}

	virtual void fixup()
	{
		for (auto r : mDb["refine_hist"])
		{
			try
			{
				int p, n, h, s;
				cif::tie(p, n, h, s) = r.get("pdbx_number_atoms_protein", "pdbx_number_atoms_nucleic_acid", "pdbx_number_atoms_ligand", "number_atoms_solvent");
				r["number_atoms_total"] = p + n + h + s;
			}
			catch (...)
			{
			}
		}
	}
};

const TemplateLine kREFMAC_Template[] = {
	/* 0 */ { "DATA USED IN REFINEMENT.", 1 },
	/* 1 */ { R"(RESOLUTION RANGE HIGH \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_high" } },
	/* 3 */ { R"(RESOLUTION RANGE LOW \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_low" } },
	/* 4 */ { R"(DATA CUTOFF \(SIGMA\(F\)\)\s*:\s*(.+))", 1, "refine", { "pdbx_ls_sigma_F" } },
	/* 5 */ { R"(COMPLETENESS FOR RANGE \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_obs" } },
	/* 6 */ { R"(NUMBER OF REFLECTIONS\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_obs" } },
	/* 7 */ { R"(FIT TO DATA USED IN REFINEMENT.)", 1 },
	/* 8 */ { R"(CROSS-VALIDATION METHOD\s*:\s*(.+))", 1, "refine", { "pdbx_ls_cross_valid_method" } },
	/* 9 */ { R"(FREE R VALUE TEST SET SELECTION\s*:\s*(.+))", 1, "refine", { "pdbx_R_Free_selection_details" } },
	/* 10 */ { R"(R VALUE \(WORKING \+ TEST SET\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_obs" } },
	/* 11 */ { R"(R VALUE \(WORKING SET\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_work" } },
	/* 12 */ { R"(FREE R VALUE\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_free" } },
	/* 13 */ { R"(FREE R VALUE TEST SET SIZE \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_R_free" } },
	/* 14 */ { R"(FREE R VALUE TEST SET COUNT\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_R_free" } },
	/* 15 */ { R"(NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.)", 1 },
	/* 16 */ { R"(PROTEIN ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_protein" } },
	/* 17 */ { R"(NUCLEIC ACID ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_nucleic_acid" } },
	/* 18 */ { R"(HETEROGEN ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_ligand" } },
	/* 19 */ { R"(SOLVENT ATOMS\s*:\s*(.+))", 1, "refine_hist", { "number_atoms_solvent" } },
	/* 20 */ { R"(ALL ATOMS\s*:\s*(.+))", 1, /* "refine_hist", "pdbx_number_atoms_protein" */ },
	/* 21 */ { R"(B VALUES\..*)", 1 },
	/* 22 */ { R"(B VALUE TYPE\s*:\s*(.+))", 1, "refine", { "pdbx_TLS_residual_ADP_flag" } },
	/* 23 */ { R"(FROM WILSON PLOT \(A\*\*2\)\s*:\s*(.+))", 1, "reflns", { "B_iso_Wilson_estimate" } },
	/* 24 */ { R"(MEAN B VALUE \(OVERALL, A\*\*2\)\s*:\s*(.+))", 1, "refine", { "B_iso_mean" } },
	/* 25 */ { R"(OVERALL ANISOTROPIC B VALUE.)", 1 },
	/* 26 */ { R"(B11 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][1]" } },
	/* 27 */ { R"(B22 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[2][2]" } },
	/* 28 */ { R"(B33 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[3][3]" } },
	/* 29 */ { R"(B12 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][2]" } },
	/* 30 */ { R"(B13 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][3]" } },
	/* 31 */ { R"(B23 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[2][3]" } },
	/* 32 */ { R"(ESTIMATED OVERALL COORDINATE ERROR.)", 1 },
	/* 33 */ { R"(ESU BASED ON R VALUE(?:\s*\(A\))?\s*:\s*(.+))", 1, "refine", { "pdbx_overall_ESU_R" } },
	/* 34 */ { R"(ESU BASED ON FREE R VALUE(?:\s*\(A\))?\s*:\s*(.+))", 1, "refine", { "pdbx_overall_ESU_R_Free" } },
	/* 35 */ { R"(ESU BASED ON MAXIMUM LIKELIHOOD(?:\s*\(A\))?\s*:\s*(.+))", 1, "refine", { "overall_SU_ML" } },
	/* 36 */ { R"(ESU FOR B VALUES BASED ON MAXIMUM LIKELIHOOD \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "overall_SU_B" } },
	/* 37 */ { R"(RMS DEVIATIONS FROM IDEAL VALUES.)", 1 },
	/* 38 */ { R"(DISTANCE RESTRAINTS. RMS SIGMA)", 1 },
	/* 39 */ { R"(BOND LENGTH \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_bond_d", false },
	/* 40 */ { R"(ANGLE DISTANCE \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_angle_d", false },
	/* 41 */ { R"(INTRAPLANAR 1-4 DISTANCE \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_planar_d", false },
	/* 42 */ { R"(H-BOND OR METAL COORDINATION \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_hb_or_metal_coord", false },
	/* 43 */ { R"(PLANE RESTRAINT \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_plane_restr", false },
	/* 44 */ { R"(CHIRAL-CENTER RESTRAINT \(A\*\*3\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_chiral_restr", false },
	/* 45 */ { R"(NON-BONDED CONTACT RESTRAINTS.)", 1 },
	/* 46 */ { R"(SINGLE TORSION \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_singtor_nbd", false },
	/* 47 */ { R"(MULTIPLE TORSION \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_multtor_nbd", false },
	/* 48 */ { R"(H-BOND \(X\.\..Y\) \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_xyhbond_nbd", false },
	/* 49 */ { R"(H-BOND \(X-H\.\.\.Y\) \(A\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_xhyhbond_nbd", false },
	/* 50 */ { R"(CONFORMATIONAL TORSION ANGLE RESTRAINTS.)", 1 },
	/* 51 */ { R"(SPECIFIED \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_special_tor", false },
	/* 52 */ { R"(PLANAR \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_planar_tor", false },
	/* 53 */ { R"(STAGGERED \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_staggered_tor", false },
	/* 54 */ { R"(TRANSVERSE \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_transverse_tor", false },
	/* 55 */ { R"(ISOTROPIC THERMAL FACTOR RESTRAINTS. RMS SIGMA)", 1 },
	/* 56 */ { R"(MAIN-CHAIN BOND \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_mcbond_it", false },
	/* 57 */ { R"(MAIN-CHAIN ANGLE \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_mcangle_it", false },
	/* 58 */ { R"(SIDE-CHAIN BOND \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_scbond_it", false },
	/* 59 */ { R"(SIDE-CHAIN ANGLE \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "p_scangle_it", false },
};

class REFMAC_Remark3Parser : public Remark3Parser
{
  public:
	REFMAC_Remark3Parser(const std::string &name, const std::string &expMethod, PDBRecord *r, cif::datablock &db)
		: Remark3Parser(name, expMethod, r, db, kREFMAC_Template, sizeof(kREFMAC_Template) / sizeof(TemplateLine),
			  std::regex(".+"))
	{
	}

	virtual std::string program() { return "REFMAC"; }
	virtual std::string version() { return ""; }
};

const TemplateLine kREFMAC5_Template[] = {
	/* 0 */ { R"(REFINEMENT TARGET\s*:\s*(.+))", 1, "refine", { "pdbx_stereochemistry_target_values" } },
	/* 1 */ { R"(DATA USED IN REFINEMENT\.)", 1 },
	/* 2 */ { R"(RESOLUTION RANGE HIGH \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_high" } },
	/* 3 */ { R"(RESOLUTION RANGE LOW \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_low" } },
	/* 4 */ { R"(DATA CUTOFF \(SIGMA\(F\)\)\s*:\s*(.+))", 1, "refine", { "pdbx_ls_sigma_F" } },
	/* 5 */ { R"(COMPLETENESS FOR RANGE \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_obs" } },
	/* 6 */ { R"(NUMBER OF REFLECTIONS\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_obs" } },
	/* 7 */ { R"(FIT TO DATA USED IN REFINEMENT.)", 1 },
	/* 8 */ { R"(CROSS-VALIDATION METHOD\s*:\s*(.+))", 1, "refine", { "pdbx_ls_cross_valid_method" } },
	/* 9 */ { R"(FREE R VALUE TEST SET SELECTION\s*:\s*(.+))", 1, "refine", { "pdbx_R_Free_selection_details" } },
	/* 10 */ { R"(R VALUE \(WORKING \+ TEST SET\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_obs" } },
	/* 11 */ { R"(R VALUE \(WORKING SET\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_work" } },
	/* 12 */ { R"(FREE R VALUE\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_free" } },
	/* 13 */ { R"(FREE R VALUE TEST SET SIZE \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_R_free" } },
	/* 14 */ { R"(FREE R VALUE TEST SET COUNT\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_R_free" } },
	/* 15 */ { R"(FIT IN THE HIGHEST RESOLUTION BIN.)", 1 },
	/* 16 */ { R"(TOTAL NUMBER OF BINS USED\s*:\s*(.+))", 1, "refine_ls_shell", { "pdbx_total_number_of_bins_used" } },
	/* 17 */ { R"(BIN RESOLUTION RANGE HIGH(?:\s*\(A\))?\s*:\s*(.+))", 1, "refine_ls_shell", { "d_res_high" } },
	/* 18 */ { R"(BIN RESOLUTION RANGE LOW(?:\s*\(A\))?\s*:\s*(.+))", 1, "refine_ls_shell", { "d_res_low" } },
	/* 19 */ { R"(REFLECTION IN BIN \(WORKING SET\)\s*:\s*(.+))", 1, "refine_ls_shell", { "number_reflns_R_work" } },
	/* 20 */ { R"(BIN COMPLETENESS \(WORKING\+TEST\) \(%\)\s*:\s*(.+))", 1, "refine_ls_shell", { "percent_reflns_obs" } },
	/* 21 */ { R"(BIN R VALUE \(WORKING SET\)\s*:\s*(.+))", 1, "refine_ls_shell", { "R_factor_R_work" } },
	/* 22 */ { R"(BIN FREE R VALUE SET COUNT\s*:\s*(.+))", 1, "refine_ls_shell", { "number_reflns_R_free" } },
	/* 23 */ { R"(BIN FREE R VALUE\s*:\s*(.+))", 1, "refine_ls_shell", { "R_factor_R_free" } },
	/* 24 */ { R"(NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.)", 1 },
	/* 25 */ { R"(PROTEIN ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_protein" } },
	/* 26 */ { R"(NUCLEIC ACID ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_nucleic_acid" } },
	/* 27 */ { R"(HETEROGEN ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_ligand" } },
	/* 28 */ { R"(SOLVENT ATOMS\s*:\s*(.+))", 1, "refine_hist", { "number_atoms_solvent" } },
	/* 29 */ { R"(ALL ATOMS\s*:\s*(.+))", 1, /* "refine_hist", { "pdbx_number_atoms_protein" } */ },
	/* 30 */ { R"(B VALUES\..*)", 1 },
	/* 31 */ { R"(B VALUE TYPE\s*:\s*(.+))", 1, "refine", { "pdbx_TLS_residual_ADP_flag" } },
	/* 32 */ { R"(FROM WILSON PLOT \(A\*\*2\)\s*:\s*(.+))", 1, "reflns", { "B_iso_Wilson_estimate" } },
	/* 33 */ { R"(MEAN B VALUE \(OVERALL, A\*\*2\)\s*:\s*(.+))", 1, "refine", { "B_iso_mean" } },
	/* 34 */ { R"(OVERALL ANISOTROPIC B VALUE.)", 1 },
	/* 35 */ { R"(B11 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][1]" } },
	/* 36 */ { R"(B22 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[2][2]" } },
	/* 37 */ { R"(B33 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[3][3]" } },
	/* 38 */ { R"(B12 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][2]" } },
	/* 39 */ { R"(B13 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[1][3]" } },
	/* 40 */ { R"(B23 \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "aniso_B[2][3]" } },
	/* 41 */ { R"(ESTIMATED OVERALL COORDINATE ERROR.)", 1 },
	/* 42 */ { R"(ESU BASED ON R VALUE(?:\s*\(A\))?\s*:\s*(.+))", 1, "refine", { "pdbx_overall_ESU_R" } },
	/* 43 */ { R"(ESU BASED ON FREE R VALUE(?:\s*\(A\))?\s*:\s*(.+))", 1, "refine", { "pdbx_overall_ESU_R_Free" } },
	/* 44 */ { R"(ESU BASED ON MAXIMUM LIKELIHOOD(?:\s*\(A\))?\s*:\s*(.+))", 1, "refine", { "overall_SU_ML" } },
	/* 45 */ { R"(ESU FOR B VALUES BASED ON MAXIMUM LIKELIHOOD \(A\*\*2\)\s*:\s*(.+))", 1, "refine", { "overall_SU_B" } },
	/* 46 */ { R"(CORRELATION COEFFICIENTS.)", 1 },
	/* 47 */ { R"(CORRELATION COEFFICIENT FO-FC\s*:\s*(.+))", 1, "refine", { "correlation_coeff_Fo_to_Fc" } },
	/* 48 */ { R"(CORRELATION COEFFICIENT FO-FC FREE\s*:\s*(.+))", 1, "refine", { "correlation_coeff_Fo_to_Fc_free" } },
	/* 49 */ { R"(RMS DEVIATIONS FROM IDEAL VALUES COUNT RMS WEIGHT)", 1 },
	/* 50 */ { R"(BOND LENGTHS REFINED ATOMS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_bond_refined_d", false },
	/* 51 */ { R"(BOND LENGTHS OTHERS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_bond_other_d", false },
	/* 52 */ { R"(BOND ANGLES REFINED ATOMS \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_angle_refined_deg", false },
	/* 53 */ { R"(BOND ANGLES OTHERS \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_angle_other_deg", false },
	/* 54 */ { R"(TORSION ANGLES, PERIOD 1 \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_dihedral_angle_1_deg", false },
	/* 55 */ { R"(TORSION ANGLES, PERIOD 2 \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_dihedral_angle_2_deg", false },
	/* 56 */ { R"(TORSION ANGLES, PERIOD 3 \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_dihedral_angle_3_deg", false },
	/* 57 */ { R"(TORSION ANGLES, PERIOD 4 \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_dihedral_angle_4_deg", false },
	/* 58 */ { R"(CHIRAL-CENTER RESTRAINTS \(A\*\*3\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_chiral_restr", false },
	/* 59 */ { R"(GENERAL PLANES REFINED ATOMS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_gen_planes_refined", false },
	/* 60 */ { R"(GENERAL PLANES OTHERS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_gen_planes_other", false },
	/* 61 */ { R"(NON-BONDED CONTACTS REFINED ATOMS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_nbd_refined", false },
	/* 62 */ { R"(NON-BONDED CONTACTS OTHERS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_nbd_other", false },
	/* 63 */ { R"(NON-BONDED TORSION REFINED ATOMS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_nbtor_refined", false },
	/* 64 */ { R"(NON-BONDED TORSION OTHERS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_nbtor_other", false },
	/* 65 */ { R"(H-BOND \(X...Y\) REFINED ATOMS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_xyhbond_nbd_refined", false },
	/* 66 */ { R"(H-BOND \(X...Y\) OTHERS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_xyhbond_nbd_other", false },
	/* 67 */ { R"(POTENTIAL METAL-ION REFINED ATOMS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_metal_ion_refined", false },
	/* 68 */ { R"(POTENTIAL METAL-ION OTHERS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_metal_ion_other", false },
	/* 69 */ { R"(SYMMETRY VDW REFINED ATOMS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_symmetry_vdw_refined", false },
	/* 70 */ { R"(SYMMETRY VDW OTHERS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_symmetry_vdw_other", false },
	/* 71 */ { R"(SYMMETRY H-BOND REFINED ATOMS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_symmetry_hbond_refined", false },
	/* 72 */ { R"(SYMMETRY H-BOND OTHERS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_symmetry_hbond_other", false },
	/* 73 */ { R"(SYMMETRY METAL-ION REFINED ATOMS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_symmetry_metal_ion_refined", false },
	/* 74 */ { R"(SYMMETRY METAL-ION OTHERS(?:\s*\(A\))?\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_symmetry_metal_ion_other", false },
	/* 75 */ { R"(ISOTROPIC THERMAL FACTOR RESTRAINTS. COUNT RMS WEIGHT)", 1 },
	/* 76 */ { R"(MAIN-CHAIN BOND REFINED ATOMS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_mcbond_it", false },
	/* 77 */ { R"(MAIN-CHAIN BOND OTHER ATOMS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_mcbond_other", false },
	/* 78 */ { R"(MAIN-CHAIN ANGLE REFINED ATOMS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_mcangle_it", false },
	/* 79 */ { R"(MAIN-CHAIN ANGLE OTHER ATOMS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_mcangle_other", false },
	/* 80 */ { R"(SIDE-CHAIN BOND REFINED ATOMS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_scbond_it", false },
	/* 81 */ { R"(SIDE-CHAIN BOND OTHER ATOMS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_scbond_other", false },
	/* 82 */ { R"(SIDE-CHAIN ANGLE REFINED ATOMS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_scangle_it", false },
	/* 83 */ { R"(SIDE-CHAIN ANGLE OTHER ATOMS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_scangle_other", false },
	/* 84 */ { R"(LONG RANGE B REFINED ATOMS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_long_range_B_refined", false },
	/* 85 */ { R"(LONG RANGE B OTHER ATOMS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_long_range_B_other", false },
	/* 86 */ { R"(ANISOTROPIC THERMAL FACTOR RESTRAINTS. COUNT RMS WEIGHT)", 1 },
	/* 87 */ { R"(RIGID-BOND RESTRAINTS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_rigid_bond_restr", false },
	/* 88 */ { R"(SPHERICITY; FREE ATOMS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_sphericity_free", false },
	/* 89 */ { R"(SPHERICITY; BONDED ATOMS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "number", "dev_ideal", "dev_ideal_target" }, "r_sphericity_bonded", false },
	//	Simply ignore NCS, you can ask Robbie why
	/* 90 */ { R"(NCS RESTRAINTS STATISTICS)", 1 },
	/* 91 */ { R"(NUMBER OF DIFFERENT NCS GROUPS\s*:\s*(.+))", 1 },
	/* 92 */ { R"(NCS GROUP NUMBER\s*:\s*(\d+))", 1, /*"struct_ncs_dom", { "pdbx_ens_id" }*/ },
	/* 93 */ { R"(CHAIN NAMES\s*:\s*(.+))", 1, /*"struct_ncs_dom", { "details" }*/ },
	/* 94 */ { R"(NUMBER OF COMPONENTS NCS GROUP\s*:\s*(\d+))", 1 },
	/* 95 */ { R"(COMPONENT C SSSEQI TO C SSSEQI CODE)", 1 },
	//// This sucks.... The following line is fixed format
	/* 97 */ { R"((\d+)\s+(.)\s+(\d+)(.)\s+(.)\s+(\d+)(.)\s+(.+))", 0 },                                                //, "struct_ncs_dom_lim", { "pdbx_component_id", "beg_auth_asym_id", "beg_auth_seq_id", "beg_auth_icode", "end_auth_asym_id", "end_auth_seq_id", "end_auth_icode", "pdbx_refine_code" }, {}, 1 },
	/* 98 */ { R"((\d+)\s+(.)\s+(\d+)\s+(.)\s+(\d+)\s+(.+))", 0 },                                                      //, "struct_ncs_dom_lim", { "pdbx_component_id", "beg_auth_asym_id", "beg_auth_seq_id", "end_auth_asym_id", "end_auth_seq_id", "pdbx_refine_code" }, {}, 1 },
	/* 96 */ { R"(GROUP CHAIN COUNT RMS WEIGHT)", 1 },                                                                  /*, "refine_ls_restr_ncs", { "pdbx_type", "dom_id", "pdbx_auth_asym_id", "pdbx_number", "rms_dev_position", "weight_position", }*/
	/* 99 */ { R"(TIGHT POSITIONAL\s+\d+\s+(.)\s+\(A\):\s+(\d+)\s*;\s*(\d+(?:\.\d*)?)\s*;\s*(\d+(?:\.\d*)?))", 0 },     // , "refine_ls_restr_ncs", {"pdbx_auth_asym_id", "pdbx_number", "rms_dev_position", "weight_position"}, { "pdbx_type", "tight positional"}, 1 },
	/* 100 */ { R"(MEDIUM POSITIONAL\s+\d+\s+(.)\s+\(A\):\s+(\d+)\s*;\s*(\d+(?:\.\d*)?)\s*;\s*(\d+(?:\.\d*)?))", 0 },   // , "refine_ls_restr_ncs", {"pdbx_auth_asym_id", "pdbx_number", "rms_dev_position", "weight_position"}, { "pdbx_type", "medium positional"}, 1 },
	/* 101 */ { R"(LOOSE POSITIONAL\s+\d+\s+(.)\s+\(A\):\s+(\d+)\s*;\s*(\d+(?:\.\d*)?)\s*;\s*(\d+(?:\.\d*)?))", 0 },    // , "refine_ls_restr_ncs", {"pdbx_auth_asym_id", "pdbx_number", "rms_dev_position", "weight_position"}, { "pdbx_type", "loose positional"}, 1 },
	/* 102 */ { R"(TIGHT THERMAL\s+\d+\s+(.)\s+\(A\*\*2\):\s+(\d+)\s*;\s*(\d+(?:\.\d*)?)\s*;\s*(\d+(?:\.\d*)?))", 0 },  // , "refine_ls_restr_ncs", {"pdbx_auth_asym_id", "pdbx_number", "rms_dev_position", "weight_position"}, { "pdbx_type", "tight thermal", }, 1 },
	/* 103 */ { R"(MEDIUM THERMAL\s+\d+\s+(.)\s+\(A\*\*2\):\s+(\d+)\s*;\s*(\d+(?:\.\d*)?)\s*;\s*(\d+(?:\.\d*)?))", 0 }, // , "refine_ls_restr_ncs", {"pdbx_auth_asym_id", "pdbx_number", "rms_dev_position", "weight_position"}, { "pdbx_type", "medium thermal", }, 1 },
	/* 104 */ { R"(LOOSE THERMAL\s+\d+\s+(.)\s+\(A\*\*2\):\s+(\d+)\s*;\s*(\d+(?:\.\d*)?)\s*;\s*(\d+(?:\.\d*)?))", 0 },  // , "refine_ls_restr_ncs", {"pdbx_auth_asym_id", "pdbx_number", "rms_dev_position", "weight_position"}, { "pdbx_type", "loose thermal", }, 10 },
	/* 105 */ { R"(NCS GROUP NUMBER\s*:\s*(\d+))", 93 - 105, /*"struct_ncs_dom", { "pdbx_ens_id" }*/ },
	/* 106 */ { R"(TWIN DETAILS)", 1 },
	/* 107 */ { R"(NUMBER OF TWIN DOMAINS\s*:\s*(\d*))", 1 },
	/* 108 */ { R"(TWIN DOMAIN\s*:\s*(.+))", 1, "pdbx_reflns_twin", { "domain_id" }, nullptr, true },
	/* 109 */ { R"(TWIN OPERATOR\s*:\s*(.+))", 1, "pdbx_reflns_twin", { "operator" } },
	/* 110 */ { R"(TWIN FRACTION\s*:\s*(.+))", 108 - 110, "pdbx_reflns_twin", { "fraction" } },
	/* 111 */ { R"(TLS DETAILS)", 1 },
	/* 112 */ { R"(NUMBER OF TLS GROUPS\s*:\s*(.+))", 1 },
	/* 113 */ { R"(TLS GROUP\s*:\s*(.+))", 1, "pdbx_refine_tls", { "id" }, nullptr, true },
	/* 114 */ { R"(NUMBER OF COMPONENTS GROUP\s*:\s*(.+))", 1 },
	/* 115 */ { R"(COMPONENTS C SSSEQI TO C SSSEQI)", 1 },
	/* 116 */ { R"(RESIDUE RANGE\s*:\s+(\S+)\s+(\d*\S)\s+(\S+)\s+(\d*\S))", 0, "pdbx_refine_tls_group", { "beg_auth_asym_id", "beg_auth_seq_id", "end_auth_asym_id", "end_auth_seq_id" }, nullptr, true },
	/* 117 */ { R"(ORIGIN FOR THE GROUP(?:\s*\(A\))?\s*:\s*([-+]?\d+(?:\.\d+)?)\s*([-+]?\d+(?:\.\d+)?)\s*([-+]?\d+(?:\.\d+)?))", 1, "pdbx_refine_tls", { "origin_x", "origin_y", "origin_z" } },
	/* 118 */ { R"(T TENSOR)", 1 },
	/* 119 */ { R"(T11\s*:\s*(.+) T22\s*:\s*(.+))", 1, "pdbx_refine_tls", { "T[1][1]", "T[2][2]" } },
	/* 120 */ { R"(T33\s*:\s*(.+) T12\s*:\s*(.+))", 1, "pdbx_refine_tls", { "T[3][3]", "T[1][2]" } },
	/* 121 */ { R"(T13\s*:\s*(.+) T23\s*:\s*(.+))", 1, "pdbx_refine_tls", { "T[1][3]", "T[2][3]" } },
	/* 122 */ { R"(L TENSOR)", 1 },
	/* 123 */ { R"(L11\s*:\s*(.+) L22\s*:\s*(.+))", 1, "pdbx_refine_tls", { "L[1][1]", "L[2][2]" } },
	/* 124 */ { R"(L33\s*:\s*(.+) L12\s*:\s*(.+))", 1, "pdbx_refine_tls", { "L[3][3]", "L[1][2]" } },
	/* 125 */ { R"(L13\s*:\s*(.+) L23\s*:\s*(.+))", 1, "pdbx_refine_tls", { "L[1][3]", "L[2][3]" } },
	/* 126 */ { R"(S TENSOR)", 1 },
	/* 127 */ { R"(S11\s*:\s*(.+) S12\s*:\s*(.+) S13\s*:\s*(.+))", 1, "pdbx_refine_tls", { "S[1][1]", "S[1][2]", "S[1][3]" } },
	/* 128 */ { R"(S21\s*:\s*(.+) S22\s*:\s*(.+) S23\s*:\s*(.+))", 1, "pdbx_refine_tls", { "S[2][1]", "S[2][2]", "S[2][3]" } },
	/* 129 */ { R"(S31\s*:\s*(.+) S32\s*:\s*(.+) S33\s*:\s*(.+))", 113 - 129, "pdbx_refine_tls", { "S[3][1]", "S[3][2]", "S[3][3]" } },
	/* 130 */ { R"(BULK SOLVENT MODELLING.)", 1 },
	/* 131 */ { R"(METHOD USED\s*:\s*(.+))", 1, "refine", { "solvent_model_details" } },
	/* 132 */ { R"(PARAMETERS FOR MASK CALCULATION)", 1 },
	/* 133 */ { R"(VDW PROBE RADIUS\s*:\s*(.+))", 1, "refine", { "pdbx_solvent_vdw_probe_radii" } },
	/* 134 */ { R"(ION PROBE RADIUS\s*:\s*(.+))", 1, "refine", { "pdbx_solvent_ion_probe_radii" } },
	/* 135 */ { R"(SHRINKAGE RADIUS\s*:\s*(.+))", 1, "refine", { "pdbx_solvent_shrinkage_radii" } },
};

class REFMAC5_Remark3Parser : public Remark3Parser
{
  public:
	REFMAC5_Remark3Parser(const std::string &name, const std::string &expMethod, PDBRecord *r, cif::datablock &db)
		: Remark3Parser(name, expMethod, r, db, kREFMAC5_Template, sizeof(kREFMAC5_Template) / sizeof(TemplateLine),
			  std::regex(R"((REFMAC)(?: (\d+(?:\..+)?))?)"))
	{
	}
};

const TemplateLine kSHELXL_Template[] = {
	/* 0 */ { R"(DATA USED IN REFINEMENT\.)", 1 },
	/* 1 */ { R"(RESOLUTION RANGE HIGH \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_high" } },
	/* 2 */ { R"(RESOLUTION RANGE LOW \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_low" } },
	/* 3 */ { R"(DATA CUTOFF \(SIGMA\(F\)\)\s*:\s*(.+))", 1, "refine", { "pdbx_ls_sigma_F" } },
	/* 4 */ { R"(COMPLETENESS FOR RANGE \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_obs" } },
	/* 5 */ { R"(CROSS-VALIDATION METHOD\s*:\s*(.+))", 1, "refine", { "pdbx_ls_cross_valid_method" } },
	/* 6 */ { R"(FREE R VALUE TEST SET SELECTION\s*:\s*(.+))", 1, "refine", { "pdbx_R_Free_selection_details" } },
	/* 7 */ { R"(FIT TO DATA USED IN REFINEMENT \(NO CUTOFF\)\.)", 1 },
	/* 8 */ { R"(R VALUE \(WORKING \+ TEST SET, NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "R_factor_all_no_cutoff" } },
	/* 9 */ { R"(R VALUE \(WORKING SET, NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "R_factor_obs_no_cutoff" } },
	/* 10 */ { R"(FREE R VALUE \(NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_factor_no_cutoff" } },
	/* 11 */ { R"(FREE R VALUE TEST SET SIZE \(%, NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_val_test_set_size_perc_no_cutoff" } },
	/* 12 */ { R"(FREE R VALUE TEST SET COUNT \(NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_val_test_set_ct_no_cutoff" } },
	/* 13 */ { R"(TOTAL NUMBER OF REFLECTIONS \(NO CUTOFF\)\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_all" } },
	/* 14 */ { R"(FIT/AGREEMENT OF MODEL FOR DATA WITH F>4SIG\(F\)\.)", 1 },
	/* 15 */ { R"(R VALUE \(WORKING \+ TEST SET, F>4SIG\(F\)\)\s*:\s*(.+))", 1, "pdbx_refine", { "R_factor_all_4sig_cutoff" } },
	/* 16 */ { R"(R VALUE \(WORKING SET, F>4SIG\(F\)\)\s*:\s*(.+))", 1, "pdbx_refine", { "R_factor_obs_4sig_cutoff" } },
	/* 17 */ { R"(FREE R VALUE \(F>4SIG\(F\)\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_factor_4sig_cutoff" } },
	/* 18 */ { R"(FREE R VALUE TEST SET SIZE \(%, F>4SIG\(F\)\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_val_test_set_size_perc_4sig_cutoff" } },
	/* 19 */ { R"(FREE R VALUE TEST SET COUNT \(F>4SIG\(F\)\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_val_test_set_ct_4sig_cutoff" } },
	/* 20 */ { R"(TOTAL NUMBER OF REFLECTIONS \(F>4SIG\(F\)\)\s*:\s*(.+))", 1, "pdbx_refine", { "number_reflns_obs_4sig_cutoff" } },
	/* 21 */ { R"(NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT\.)", 1 },
	/* 22 */ { R"(PROTEIN ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_protein" } },
	/* 23 */ { R"(NUCLEIC ACID ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_nucleic_acid" } },
	/* 24 */ { R"(HETEROGEN ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_ligand" } },
	/* 25 */ { R"(SOLVENT ATOMS\s*:\s*(.+))", 1, "refine_hist", { "number_atoms_solvent" } },
	/* 26 */ { R"(MODEL REFINEMENT\.)", 1 },
	/* 27 */ { R"(OCCUPANCY SUM OF NON-HYDROGEN ATOMS\s*:\s*(.+))", 1, "refine_analyze", { "occupancy_sum_non_hydrogen" } },
	/* 28 */ { R"(OCCUPANCY SUM OF HYDROGEN ATOMS\s*:\s*(.+))", 1, "refine_analyze", { "occupancy_sum_hydrogen" } },
	/* 29 */ { R"(NUMBER OF DISCRETELY DISORDERED RESIDUES\s*:\s*(.+))", 1, "refine_analyze", { "number_disordered_residues" } },
	/* 30 */ { R"(NUMBER OF LEAST-SQUARES PARAMETERS\s*:\s*(.+))", 1, "refine", { "ls_number_parameters" } },
	/* 31 */ { R"(NUMBER OF RESTRAINTS\s*:\s*(.+))", 1, "refine", { "ls_number_restraints" } },
	/* 32 */ { R"(RMS DEVIATIONS FROM RESTRAINT TARGET VALUES\.)", 1 },
	/* 33 */ { R"(BOND LENGTHS \(A\)\s*:\s*(.+))", 1, "refine_ls_restr", { "dev_ideal" }, "s_bond_d", false },
	/* 34 */ { R"(ANGLE DISTANCES \(A\)\s*:\s*(.+))", 1, "refine_ls_restr", { "dev_ideal" }, "s_angle_d", false },
	/* 35 */ { R"(SIMILAR DISTANCES \(NO TARGET VALUES\) \(A\)\s*:\s*(.+))", 1, "refine_ls_restr", { "dev_ideal" }, "s_similar_dist", false },
	/* 36 */ { R"(DISTANCES FROM RESTRAINT PLANES \(A\)\s*:\s*(.+))", 1, "refine_ls_restr", { "dev_ideal" }, "s_from_restr_planes", false },
	/* 37 */ { R"(ZERO CHIRAL VOLUMES \(A\*\*3\)\s*:\s*(.+))", 1, "refine_ls_restr", { "dev_ideal" }, "s_zero_chiral_vol", false },
	/* 38 */ { R"(NON-ZERO CHIRAL VOLUMES \(A\*\*3\)\s*:\s*(.+))", 1, "refine_ls_restr", { "dev_ideal" }, "s_non_zero_chiral_vol", false },
	/* 39 */ { R"(ANTI-BUMPING DISTANCE RESTRAINTS \(A\)\s*:\s*(.+))", 1, "refine_ls_restr", { "dev_ideal" }, "s_anti_bump_dis_restr", false },
	/* 40 */ { R"(RIGID-BOND ADP COMPONENTS \(A\*\*2\)\s*:\s*(.+))", 1, "refine_ls_restr", { "dev_ideal" }, "s_rigid_bond_adp_cmpnt", false },
	/* 41 */ { R"(SIMILAR ADP COMPONENTS \(A\*\*2\)\s*:\s*(.+))", 1, "refine_ls_restr", { "dev_ideal" }, "s_similar_adp_cmpnt", false },
	/* 42 */ { R"(APPROXIMATELY ISOTROPIC ADPS \(A\*\*2\)\s*:\s*(.+))", 1, "refine_ls_restr", { "dev_ideal" }, "s_approx_iso_adps", false },
	/* 43 */ { R"(BULK SOLVENT MODELING\.)", 1 },
	/* 44 */ { R"(METHOD USED\s*:\s*(.+))", 1, "refine", { "solvent_model_details" } },
	/* 45 */ { R"(STEREOCHEMISTRY TARGET VALUES\s*:\s*(.+))", 1, "refine", { "pdbx_stereochemistry_target_values" } },
	/* 46 */ { R"(SPECIAL CASE\s*:\s*(.+))", 1, "refine", { "pdbx_stereochem_target_val_spec_case" } },
};

class SHELXL_Remark3Parser : public Remark3Parser
{
  public:
	SHELXL_Remark3Parser(const std::string &name, const std::string &expMethod, PDBRecord *r, cif::datablock &db)
		: Remark3Parser(name, expMethod, r, db, kSHELXL_Template, sizeof(kSHELXL_Template) / sizeof(TemplateLine),
			  std::regex(R"((SHELXL)(?:-(\d+(?:\..+)?)))"))
	{
	}
};

const TemplateLine kTNT_Template[] = {
	/* 0 */ { R"(DATA USED IN REFINEMENT\.)", 1 },
	/* 1 */ { R"(RESOLUTION RANGE HIGH \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_high" } },
	/* 2 */ { R"(RESOLUTION RANGE LOW \(ANGSTROMS\)\s*:\s*(.+))", 1, "refine", { "ls_d_res_low" } },
	/* 3 */ { R"(DATA CUTOFF \(SIGMA\(F\)\)\s*:\s*(.+))", 1, "refine", { "pdbx_ls_sigma_F" } },
	/* 4 */ { R"(COMPLETENESS FOR RANGE \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_obs" } },
	/* 5 */ { R"(NUMBER OF REFLECTIONS\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_obs" } },
	/* 6 */ { R"(USING DATA ABOVE SIGMA CUTOFF\.)", 1 },
	/* 7 */ { R"(CROSS-VALIDATION METHOD\s*:\s*(.+))", 1, "refine", { "pdbx_ls_cross_valid_method" } },
	/* 8 */ { R"(FREE R VALUE TEST SET SELECTION\s*:\s*(.+))", 1, "refine", { "pdbx_R_Free_selection_details" } },
	/* 9 */ { R"(R VALUE \(WORKING \+ TEST SET\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_obs" } },
	/* 10 */ { R"(R VALUE \(WORKING SET\)\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_work" } },
	/* 11 */ { R"(FREE R VALUE\s*:\s*(.+))", 1, "refine", { "ls_R_factor_R_free" } },
	/* 12 */ { R"(FREE R VALUE TEST SET SIZE \(%\)\s*:\s*(.+))", 1, "refine", { "ls_percent_reflns_R_free" } },
	/* 13 */ { R"(FREE R VALUE TEST SET COUNT\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_R_free" } },
	/* 14 */ { R"(USING ALL DATA, NO SIGMA CUTOFF\.)", 1 },
	/* 15 */ { R"(R VALUE \(WORKING \+ TEST SET, NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "R_factor_all_no_cutoff" } },
	/* 16 */ { R"(R VALUE \(WORKING SET, NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "R_factor_obs_no_cutoff" } },
	/* 17 */ { R"(FREE R VALUE \(NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_factor_no_cutoff" } },
	/* 18 */ { R"(FREE R VALUE TEST SET SIZE \(%, NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_val_test_set_size_perc_no_cutoff" } },
	/* 19 */ { R"(FREE R VALUE TEST SET COUNT \(NO CUTOFF\)\s*:\s*(.+))", 1, "pdbx_refine", { "free_R_val_test_set_ct_no_cutoff" } },
	/* 20 */ { R"(TOTAL NUMBER OF REFLECTIONS \(NO CUTOFF\)\s*:\s*(.+))", 1, "refine", { "ls_number_reflns_all" } },
	/* 21 */ { R"(NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT\.)", 1 },
	/* 22 */ { R"(PROTEIN ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_protein" } },
	/* 23 */ { R"(NUCLEIC ACID ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_nucleic_acid" } },
	/* 24 */ { R"(HETEROGEN ATOMS\s*:\s*(.+))", 1, "refine_hist", { "pdbx_number_atoms_ligand" } },
	/* 25 */ { R"(SOLVENT ATOMS\s*:\s*(.+))", 1, "refine_hist", { "number_atoms_solvent" } },
	/* 26 */ { R"(WILSON B VALUE \(FROM FCALC, A\*\*2\)\s*:\s*(.+))", 1, "reflns", { "B_iso_Wilson_estimate" } },
	/* 27 */ { R"(RMS DEVIATIONS FROM IDEAL VALUES\. RMS WEIGHT COUNT)", 1 },
	/* 28 */ { R"(BOND LENGTHS \(A\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_bond_d", false },
	/* 29 */ { R"(BOND ANGLES \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_angle_deg", false },
	/* 30 */ { R"(TORSION ANGLES \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_dihedral_angle_d", false },
	/* 31 */ { R"(PSEUDOROTATION ANGLES \(DEGREES\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_pseud_angle", false },
	/* 32 */ { R"(TRIGONAL CARBON PLANES \(A\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_trig_c_planes", false },
	/* 33 */ { R"(GENERAL PLANES \(A\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_gen_planes", false },
	/* 34 */ { R"(ISOTROPIC THERMAL FACTORS \(A\*\*2\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_it", false },
	/* 35 */ { R"(NON-BONDED CONTACTS \(A\)\s*:\s*(.+)\s*;\s*(.+)\s*;\s*(.+))", 1, "refine_ls_restr", { "dev_ideal", "weight", "number" }, "t_nbd", false },
	/* 36 */ { R"(INCORRECT CHIRAL-CENTERS \(COUNT\)\s*:\s*(.+)\s*)", 1, "refine_ls_restr", { "number" }, "t_incorr_chiral_ct", false },
	/* 37 */ { R"(BULK SOLVENT MODELING\.)", 1 },
	/* 38 */ { R"(METHOD USED\s*:\s*(.+))", 1, "refine", { "solvent_model_details" } },
	/* 39 */ { R"(KSOL\s*:\s*(.+))", 1, "refine", { "solvent_model_param_ksol" } },
	/* 40 */ { R"(BSOL\s*:\s*(.+))", 1, "refine", { "solvent_model_param_bsol" } },
	/* 41 */ { R"(RESTRAINT LIBRARIES\.)", 1 },
	/* 42 */ { R"(STEREOCHEMISTRY\s*:\s*(.+))", 1, "refine", { "pdbx_stereochemistry_target_values" } },
	/* 43 */ { R"(ISOTROPIC THERMAL FACTOR RESTRAINTS\s*:\s*(.+))", 1, "refine", { "pdbx_isotropic_thermal_model" } },
};

class TNT_Remark3Parser : public Remark3Parser
{
  public:
	TNT_Remark3Parser(const std::string &name, const std::string &expMethod, PDBRecord *r, cif::datablock &db)
		: Remark3Parser(name, expMethod, r, db, kTNT_Template, sizeof(kTNT_Template) / sizeof(TemplateLine),
			  std::regex(R"((TNT)(?: V. (\d+.+)?)?)"))
	{
	}
};

const TemplateLine kXPLOR_Template[] = {
	/* 0 */ { R"(DATA USED IN REFINEMENT\.)", 1 },
	/* 1 */ { R"(RESOLUTION RANGE HIGH \(ANGSTROMS\) :\s+(.+))", 1, "refine", { "ls_d_res_high" } },
	/* 2 */ { R"(RESOLUTION RANGE LOW \(ANGSTROMS\) :\s+(.+))", 1, "refine", { "ls_d_res_low" } },
	/* 3 */ { R"(DATA CUTOFF \(SIGMA\(F\)\) :\s+(.+))", 1, "refine", { "pdbx_ls_sigma_F" } },
	/* 4 */ { R"(DATA CUTOFF HIGH \(ABS\(F\)\) :\s+(.+))", 1, "refine", { "pdbx_data_cutoff_high_absF" } },
	/* 5 */ { R"(DATA CUTOFF LOW \(ABS\(F\)\) :\s+(.+))", 1, "refine", { "pdbx_data_cutoff_low_absF" } },
	/* 6 */ { R"(COMPLETENESS \(WORKING\+TEST\) \(%\) :\s+(.+))", 1, "refine", { "ls_percent_reflns_obs" } },
	/* 7 */ { R"(NUMBER OF REFLECTIONS :\s+(.+))", 1, "refine", { "ls_number_reflns_obs" } },
	/* 8 */ { R"(FIT TO DATA USED IN REFINEMENT\.)", 1 },
	/* 9 */ { R"(CROSS-VALIDATION METHOD :\s+(.+))", 1, "refine", { "pdbx_ls_cross_valid_method" } },
	/* 10 */ { R"(FREE R VALUE TEST SET SELECTION :\s+(.+))", 1, "refine", { "pdbx_R_Free_selection_details" } },
	/* 11 */ { R"(R VALUE \(WORKING SET\) :\s+(.+))", 1, "refine", { "ls_R_factor_R_work" } },
	/* 12 */ { R"(FREE R VALUE :\s+(.+))", 1, "refine", { "ls_R_factor_R_free" } },
	/* 13 */ { R"(FREE R VALUE TEST SET SIZE \(%\) :\s+(.+))", 1, "refine", { "ls_percent_reflns_R_free" } },
	/* 14 */ { R"(FREE R VALUE TEST SET COUNT :\s+(.+))", 1, "refine", { "ls_number_reflns_R_free" } },
	/* 15 */ { R"(ESTIMATED ERROR OF FREE R VALUE :\s+(.+))", 1, "refine", { "ls_R_factor_R_free_error" } },
	/* 16 */ { R"(FIT IN THE HIGHEST RESOLUTION BIN\.)", 1 },
	/* 17 */ { R"(TOTAL NUMBER OF BINS USED :\s+(.+))", 1, "refine_ls_shell", { "pdbx_total_number_of_bins_used" } },
	/* 18 */ { R"(BIN RESOLUTION RANGE HIGH \(A\) :\s+(.+))", 1, "refine_ls_shell", { "d_res_high" } },
	/* 19 */ { R"(BIN RESOLUTION RANGE LOW \(A\) :\s+(.+))", 1, "refine_ls_shell", { "d_res_low" } },
	/* 20 */ { R"(BIN COMPLETENESS \(WORKING\+TEST\) \(%\) :\s+(.+))", 1, "refine_ls_shell", { "percent_reflns_obs" } },
	/* 21 */ { R"(REFLECTIONS IN BIN \(WORKING SET\) :\s+(.+))", 1, "refine_ls_shell", { "number_reflns_R_work" } },
	/* 22 */ { R"(BIN R VALUE \(WORKING SET\) :\s+(.+))", 1, "refine_ls_shell", { "R_factor_R_work" } },
	/* 23 */ { R"(BIN FREE R VALUE :\s+(.+))", 1, "refine_ls_shell", { "R_factor_R_free" } },
	/* 24 */ { R"(BIN FREE R VALUE TEST SET SIZE \(%\) :\s+(.+))", 1, "refine_ls_shell", { "percent_reflns_R_free" } },
	/* 25 */ { R"(BIN FREE R VALUE TEST SET COUNT :\s+(.+))", 1, "refine_ls_shell", { "number_reflns_R_free" } },
	/* 26 */ { R"(ESTIMATED ERROR OF BIN FREE R VALUE :\s+(.+))", 1, "refine_ls_shell", { "R_factor_R_free_error" } },
	/* 27 */ { R"(NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT\.)", 1 },
	/* 28 */ { R"(PROTEIN ATOMS :\s+(.+))", 1, "refine_hist", { "pdbx_number_atoms_protein" } },
	/* 29 */ { R"(NUCLEIC ACID ATOMS :\s+(.+))", 1, "refine_hist", { "pdbx_number_atoms_nucleic_acid" } },
	/* 30 */ { R"(HETEROGEN ATOMS :\s+(.+))", 1, "refine_hist", { "pdbx_number_atoms_ligand" } },
	/* 31 */ { R"(SOLVENT ATOMS :\s+(.+))", 1, "refine_hist", { "number_atoms_solvent" } },
	/* 32 */ { R"(B VALUES\.)", 1 },
	/* 33 */ { R"(B VALUE TYPE :\s+(.+))", 1, "refine", { "pdbx_TLS_residual_ADP_flag" } },
	/* 34 */ { R"(FROM WILSON PLOT \(A\*\*2\) :\s+(.+))", 1, "reflns", { "B_iso_Wilson_estimate" } },
	/* 35 */ { R"(MEAN B VALUE \(OVERALL, A\*\*2\) :\s+(.+))", 1, "refine", { "B_iso_mean" } },
	/* 36 */ { R"(OVERALL ANISOTROPIC B VALUE\.)", 1 },
	/* 37 */ { R"(B11 \(A\*\*2\) :\s+(.+))", 1, "refine", { "aniso_B[1][1]" } },
	/* 38 */ { R"(B22 \(A\*\*2\) :\s+(.+))", 1, "refine", { "aniso_B[2][2]" } },
	/* 39 */ { R"(B33 \(A\*\*2\) :\s+(.+))", 1, "refine", { "aniso_B[3][3]" } },
	/* 40 */ { R"(B12 \(A\*\*2\) :\s+(.+))", 1, "refine", { "aniso_B[1][2]" } },
	/* 41 */ { R"(B13 \(A\*\*2\) :\s+(.+))", 1, "refine", { "aniso_B[1][3]" } },
	/* 42 */ { R"(B23 \(A\*\*2\) :\s+(.+))", 1, "refine", { "aniso_B[2][3]" } },
	/* 43 */ { R"(ESTIMATED COORDINATE ERROR\.)", 1 },
	/* 44 */ { R"(ESD FROM LUZZATI PLOT \(A\) :\s+(.+))", 1, "refine_analyze", { "Luzzati_coordinate_error_obs" } },
	/* 45 */ { R"(ESD FROM SIGMAA \(A\) :\s+(.+))", 1, "refine_analyze", { "Luzzati_sigma_a_obs" } },
	/* 46 */ { R"(LOW RESOLUTION CUTOFF \(A\) :\s+(.+))", 1, "refine_analyze", { "Luzzati_d_res_low_obs" } },
	/* 47 */ { R"(CROSS-VALIDATED ESTIMATED COORDINATE ERROR\.)", 1 },
	/* 48 */ { R"(ESD FROM C-V LUZZATI PLOT \(A\) :\s+(.+))", 1, "refine_analyze", { "Luzzati_coordinate_error_free" } },
	/* 49 */ { R"(ESD FROM C-V SIGMAA \(A\) :\s+(.+))", 1, "refine_analyze", { "Luzzati_sigma_a_free" } },
	/* 50 */ { R"(RMS DEVIATIONS FROM IDEAL VALUES\..*)", 1 },
	/* 51 */ { R"(BOND LENGTHS \(A\) :\s+(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "x_bond_d", false },
	/* 52 */ { R"(BOND ANGLES \(DEGREES\) :\s+(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "x_angle_deg", false },
	/* 53 */ { R"(DIHEDRAL ANGLES \(DEGREES\) :\s+(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "x_dihedral_angle_d", false },
	/* 54 */ { R"(IMPROPER ANGLES \(DEGREES\) :\s+(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "x_improper_angle_d", false },
	/* 55 */ { R"(ISOTROPIC THERMAL MODEL :\s+(.+))", 1, "refine", { "pdbx_isotropic_thermal_model" } },
	/* 56 */ { R"(ISOTROPIC THERMAL FACTOR RESTRAINTS\. RMS SIGMA)", 1 },
	/* 57 */ { R"(MAIN-CHAIN BOND \(A\*\*2\) :\s+(.+?);\s+(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "x_mcbond_it", false },
	/* 58 */ { R"(MAIN-CHAIN ANGLE \(A\*\*2\) :\s+(.+?);\s+(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "x_mcangle_it", false },
	/* 59 */ { R"(SIDE-CHAIN BOND \(A\*\*2\) :\s+(.+?);\s+(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "x_scbond_it", false },
	/* 60 */ { R"(SIDE-CHAIN ANGLE \(A\*\*2\) :\s+(.+?);\s+(.+))", 1, "refine_ls_restr", { "dev_ideal", "dev_ideal_target" }, "x_scangle_it", false },
	/* 61 */ { R"(NCS MODEL :\s+(.+))", 1, /* "refine_ls_restr_ncs", { "ncs_model_details" } */ },
	/* 62 */ { R"(NCS RESTRAINTS\. RMS SIGMA/WEIGHT)", 1 },
	/* 63 */ { R"(GROUP (\d+) POSITIONAL \(A\) :\s+(.+?);\s+(.+))", 1, /* "refine_ls_restr_ncs", { ":dom_id", "rms_dev_position", "weight_position" } */ },
	/* 64 */ { R"(GROUP (\d+) B-FACTOR \(A\*\*2\) :\s+(.+?);\s+(.+))", 63 - 64, /* "refine_ls_restr_ncs", { ":dom_id", "rms_dev_B_iso", "weight_B_iso" } */ },
	/* 65 */ { R"(PARAMETER FILE (\d+) :\s+(.+))", 0, /* "pdbx_xplor_file", { "serial_no", "param_file" } */ },
	/* 66 */ { R"(TOPOLOGY FILE (\d+) :\s+(.+))", 0, /* "pdbx_xplor_file", { "serial_no", "topol_file" } */ },
};

class XPLOR_Remark3Parser : public Remark3Parser
{
  public:
	XPLOR_Remark3Parser(const std::string &name, const std::string &expMethod, PDBRecord *r, cif::datablock &db)
		: Remark3Parser(name, expMethod, r, db, kXPLOR_Template, sizeof(kXPLOR_Template) / sizeof(TemplateLine),
			  std::regex(R"((X-PLOR)(?: (\d+(?:\.\d+)?))?)"))
	{
	}
};

// --------------------------------------------------------------------

Remark3Parser::Remark3Parser(const std::string &name, const std::string &expMethod, PDBRecord *r, cif::datablock &db,
	const TemplateLine templatelines[], uint32_t templateLineCount, std::regex programversion)
	: mName(name)
	, mExpMethod(expMethod)
	, mRec(r)
	, mDb(db.name())
	, mTemplate(templatelines)
	, mTemplateCount(templateLineCount)
	, mProgramVersion(programversion)
{
	mDb.set_validator(db.get_validator());
}

std::string Remark3Parser::nextLine()
{
	mLine.clear();

	while (mRec != nullptr and mRec->is("REMARK   3"))
	{
		size_t valueIndent = 0;
		for (size_t i = 4; i < mRec->mVlen; ++i)
		{
			if (mRec->mValue[i] == ' ')
				continue;

			if (mRec->mValue[i] == ':')
			{
				valueIndent = i;
				while (valueIndent < mRec->mVlen and mRec->mValue[i] == ' ')
					++valueIndent;
				break;
			}
		}

		mLine = mRec->vS(12);
		mRec = mRec->mNext;

		if (mLine.empty())
			continue;

		// concatenate value that is wrapped over multiple lines (tricky code...)

		if (valueIndent > 4)
		{
			std::string indent(valueIndent - 4, ' ');

			while (mRec->is("REMARK   3") and mRec->mVlen > valueIndent)
			{
				std::string v(mRec->mValue + 4, mRec->mValue + mRec->mVlen);
				if (not cif::starts_with(v, indent))
					break;

				mLine += ' ';
				mLine.append(mRec->mValue + valueIndent, mRec->mValue + mRec->mVlen);

				mRec = mRec->mNext;
			}
		}

		// collapse multiple spaces
		bool space = false;
		auto i = mLine.begin(), j = i;

		while (i != mLine.end())
		{
			bool nspace = isspace(*i);

			if (nspace == false)
			{
				if (space)
					*j++ = ' ';
				*j++ = *i;
			}
			space = nspace;
			++i;
		}
		mLine.erase(j, mLine.end());

		break;
	}

	if (cif::VERBOSE >= 2)
		std::cerr << "RM3: " << mLine << '\n';

	return mLine;
}

bool Remark3Parser::match(const char *expr, int nextState)
{
	std::regex rx(expr);

	bool result = regex_match(mLine, mM, rx);

	if (result)
		mState = nextState;
	else if (cif::VERBOSE >= 3)
	{
		using namespace colour;

		std::cerr << coloured("No match:", white, red, bold) << " '" << expr << '\'' << '\n';
	}

	return result;
}

float Remark3Parser::parse()
{
	int lineCount = 0, dropped = 0;
	std::string remarks;
	mState = 0;

	while (mRec != nullptr)
	{
		nextLine();

		if (mLine.empty())
			break;

		++lineCount;

		// Skip over AUTHORS lines
		if (mState == 0 and match(R"(AUTHORS\s*:.+)", 0))
			continue;

		auto state = mState;
		for (state = mState; state < mTemplateCount; ++state)
		{
			const TemplateLine &tmpl = mTemplate[state];

			if (match(tmpl.rx, state + tmpl.nextStateOffset))
			{
				if (not(tmpl.category == nullptr or tmpl.items.size() == 0))
				{
					if (tmpl.lsRestrType == nullptr)
						storeCapture(tmpl.category, tmpl.items, tmpl.createNew);
					else if (tmpl.createNew)
						storeRefineLsRestr(tmpl.lsRestrType, tmpl.items);
					else
						updateRefineLsRestr(tmpl.lsRestrType, tmpl.items);
				}
				break;
			}
		}

		if (state < mTemplateCount)
			continue;

		if (state == mTemplateCount and match(R"(OTHER REFINEMENT REMARKS\s*:\s*(.*))", mTemplateCount + 1))
		{
			remarks = mM[1].str();
			continue;
		}

		if (state == mTemplateCount + 1)
		{
			remarks = remarks + '\n' + mLine;
			continue;
		}

		if (cif::VERBOSE >= 2)
		{
			using namespace colour;

			std::cerr << coloured("Dropping line:", white, red, bold) << " '" << mLine << '\'' << '\n';
		}

		++dropped;
	}

	if (not remarks.empty() and not iequals(remarks, "NULL"))
	{
		if (not mDb["refine"].empty())
			mDb["refine"].front()["details"] = remarks;
	}

	float score = float(lineCount - dropped) / lineCount;

	return score;
}

std::string Remark3Parser::program()
{
	std::string result = mName;

	std::smatch m;
	if (regex_match(mName, m, mProgramVersion))
		result = m[1].str();

	return result;
}

std::string Remark3Parser::version()
{
	std::string result;

	std::smatch m;
	if (regex_match(mName, m, mProgramVersion))
		result = m[2].str();

	return result;
}

void Remark3Parser::storeCapture(const char *category, std::initializer_list<const char *> items, bool createNew)
{
	int capture = 0;
	for (auto item : items)
	{
		++capture;

		std::string value = mM[capture].str();
		cif::trim(value);

		if (iequals(value, "NULL") or iequals(value, "NONE") or iequals(value, "Inf") or iequals(value, "+Inf") or iequals(value, std::string(value.length(), '*')))
			continue;

		if (cif::VERBOSE >= 3)
			std::cerr << "storing: '" << value << "' in _" << category << '.' << item << '\n';

		auto &cat = mDb[category];
		if (cat.empty() or createNew)
		{
			if (iequals(category, "refine"))
				cat.emplace({ { "pdbx_refine_id", mExpMethod },
					{ "entry_id", mDb.name() },
					//#warning("this diffrn-id is probably not correct?")
					{ "pdbx_diffrn_id", 1 } });
			else if (iequals(category, "refine_analyze") or iequals(category, "pdbx_refine"))
				cat.emplace({
					{ "pdbx_refine_id", mExpMethod },
					{ "entry_id", mDb.name() },
					//					{ "pdbx_diffrn_id", 1 }
				});
			else if (iequals(category, "refine_hist"))
			{
				std::string dResHigh, dResLow;
				for (auto r : mDb["refine"])
				{
					cif::tie(dResHigh, dResLow) = r.get("ls_d_res_high", "ls_d_res_low");
					break;
				}

				cat.emplace({ { "pdbx_refine_id", mExpMethod },
					{ "cycle_id", "LAST" },
					{ "d_res_high", dResHigh.empty() ? "." : dResHigh },
					{ "d_res_low", dResLow.empty() ? "." : dResLow } });
			}
			else if (iequals(category, "refine_ls_shell"))
			{
				cat.emplace({
					{ "pdbx_refine_id", mExpMethod },
				});
			}
			else if (iequals(category, "pdbx_refine_tls_group"))
			{
				std::string tlsID;
				if (not mDb["pdbx_refine_tls"].empty())
					tlsID = mDb["pdbx_refine_tls"].back()["id"].as<std::string>();
				std::string tlsGroupID = cat.get_unique_id("");

				cat.emplace({
					{ "pdbx_refine_id", mExpMethod },
					{ "id", tlsGroupID },
					{ "refine_tls_id", tlsID } });
			}
			else if (iequals(category, "pdbx_refine_tls"))
			{
				cat.emplace({ { "pdbx_refine_id", mExpMethod },
					{ "method", "refined" } });
			}
			//			else if (iequals(category, "struct_ncs_dom"))
			//			{
			//				size_t id = cat.size() + 1;
			//
			//				cat.emplace({
			//					{ "id", id }
			//				});
			//			}
			else if (iequals(category, "pdbx_reflns_twin"))
			{
				cat.emplace({ // #warning("crystal id, diffrn id, what should be put here?")
					{ "crystal_id", 1 },
					{ "diffrn_id", 1 },
					{ "operator", "" },
					{ "fraction", 0.f } });
			}
			else if (iequals(category, "reflns"))
				cat.emplace({ { "pdbx_ordinal", cat.size() + 1 },
					{ "entry_id", mDb.name() },
					{ "pdbx_diffrn_id", 1 } });
			else
				cat.emplace({});

			createNew = false;
		}

		cat.back()[item] = value;
	}
}

void Remark3Parser::storeRefineLsRestr(const char *type, std::initializer_list<const char *> items)
{
	cif::row_handle r;
	int capture = 0;

	for (auto item : items)
	{
		++capture;

		std::string value = mM[capture].str();
		cif::trim(value);
		if (value.empty() or iequals(value, "NULL") or iequals(value, "Inf") or iequals(value, "+Inf") or iequals(value, std::string(value.length(), '*')))
			continue;

		if (r.empty())
		{
			r = mDb["refine_ls_restr"].emplace({
				{"pdbx_refine_id", mExpMethod},
				{"type", type}
			});
		}

		r[item] = value;
	}
}

void Remark3Parser::updateRefineLsRestr(const char *type, std::initializer_list<const char *> items)
{
	auto rows = mDb["refine_ls_restr"].find(cif::key("type") == type and cif::key("pdbx_refine_id") == mExpMethod);
	if (rows.empty())
		storeRefineLsRestr(type, items);
	else
	{
		for (auto r : rows)
		{
			int capture = 0;
			for (auto item : items)
			{
				++capture;

				std::string value = mM[capture].str();
				cif::trim(value);
				if (iequals(value, "NULL") or iequals(value, std::string(value.length(), '*')))
					value.clear();

				r[item] = value;
			}

			break;
		}
	}
}

// --------------------------------------------------------------------

bool Remark3Parser::parse(const std::string &expMethod, PDBRecord *r, cif::datablock &db)
{
	// simple version, only for the first few lines
	auto getNextLine = [&]()
	{
		std::string result;

		while (result.empty() and r != nullptr and r->is("REMARK   3"))
		{
			result = r->vS(12);
			r = r->mNext;
		}

		return result;
	};

	// All remark 3 records should start with the same data.

	std::string line = getNextLine();

	if (line != "REFINEMENT.")
	{
		if (cif::VERBOSE > 0)
			std::cerr << "Unexpected data in REMARK 3\n";
		return false;
	}

	line = getNextLine();

	std::regex rxp(R"(^PROGRAM\s*:\s*(.+))");
	std::smatch m;

	if (not std::regex_match(line, m, rxp))
	{
		if (cif::VERBOSE > 0)
			std::cerr << "Expected valid PROGRAM line in REMARK 3\n";
		return false;
	}

	line = m[1].str();

	struct programScore
	{
		programScore(const std::string &program, Remark3Parser *parser, float score)
			: program(program)
			, parser(parser)
			, score(score)
		{
		}

		std::string program;
		std::unique_ptr<Remark3Parser> parser;
		float score;

		bool operator<(const programScore &rhs) const
		{
			return score > rhs.score;
		}
	};

	std::vector<programScore> scores;

	auto tryParser = [&](Remark3Parser *p)
	{
		std::unique_ptr<Remark3Parser> parser(p);
		float score;

		try
		{
			score = parser->parse();
		}
		catch (const std::exception &e)
		{
			if (cif::VERBOSE >= 0)
				std::cerr << "Error parsing REMARK 3 with " << parser->program() << '\n'
						  << e.what() << '\n';
			score = 0;
		}

		if (cif::VERBOSE >= 2)
			std::cerr << "Score for " << parser->program() << ": " << score << '\n';

		if (score > 0)
		{
			std::string program = parser->program();
			std::string version = parser->version();

			scores.emplace_back(program, parser.release(), score);
		}
	};

	for (auto program : cif::split<std::string>(line, ", ", true))
	{
		if (cif::starts_with(program, "BUSTER"))
			tryParser(new BUSTER_TNT_Remark3Parser(program, expMethod, r, db));
		else if (cif::starts_with(program, "CNS") or cif::starts_with(program, "CNX"))
			tryParser(new CNS_Remark3Parser(program, expMethod, r, db));
		else if (cif::starts_with(program, "PHENIX"))
			tryParser(new PHENIX_Remark3Parser(program, expMethod, r, db));
		else if (cif::starts_with(program, "NUCLSQ"))
			tryParser(new NUCLSQ_Remark3Parser(program, expMethod, r, db));
		else if (cif::starts_with(program, "PROLSQ"))
			tryParser(new PROLSQ_Remark3Parser(program, expMethod, r, db));
		else if (cif::starts_with(program, "REFMAC"))
		{
			// simply try both and take the best
			tryParser(new REFMAC_Remark3Parser(program, expMethod, r, db));
			tryParser(new REFMAC5_Remark3Parser(program, expMethod, r, db));
		}
		else if (cif::starts_with(program, "SHELXL"))
			tryParser(new SHELXL_Remark3Parser(program, expMethod, r, db));
		else if (cif::starts_with(program, "TNT"))
			tryParser(new TNT_Remark3Parser(program, expMethod, r, db));
		else if (cif::starts_with(program, "X-PLOR"))
			tryParser(new XPLOR_Remark3Parser(program, expMethod, r, db));
		else if (cif::VERBOSE > 0)
			std::cerr << "Skipping unknown program (" << program << ") in REMARK 3\n";
	}

	sort(scores.begin(), scores.end());

	bool guessProgram = scores.empty() or scores.front().score < 0.9f;
	if (guessProgram)
	{
		if (cif::VERBOSE > 0)
			std::cerr << "Unknown or untrusted program in REMARK 3, trying all parsers to see if there is a match\n";

		tryParser(new BUSTER_TNT_Remark3Parser("BUSTER-TNT", expMethod, r, db));
		tryParser(new CNS_Remark3Parser("CNS", expMethod, r, db));
		tryParser(new PHENIX_Remark3Parser("PHENIX", expMethod, r, db));
		tryParser(new NUCLSQ_Remark3Parser("NUCLSQ", expMethod, r, db));
		tryParser(new PROLSQ_Remark3Parser("PROLSQ", expMethod, r, db));
		tryParser(new REFMAC_Remark3Parser("REFMAC", expMethod, r, db));
		tryParser(new REFMAC5_Remark3Parser("REFMAC5", expMethod, r, db));
		tryParser(new SHELXL_Remark3Parser("SHELXL", expMethod, r, db));
		tryParser(new TNT_Remark3Parser("TNT", expMethod, r, db));
		tryParser(new XPLOR_Remark3Parser("X-PLOR", expMethod, r, db));
	}

	bool result = false;

	if (not scores.empty())
	{
		result = true;

		sort(scores.begin(), scores.end());

		auto &best = scores.front();

		if (cif::VERBOSE > 0)
			std::cerr << "Choosing " << best.parser->program() << " version '" << best.parser->version() << "' as refinement program. Score = " << best.score << '\n';

		auto &software = db["software"];
		std::string program = best.parser->program();
		std::string version = best.parser->version();

		software.emplace({ { "name", program },
			{ "classification", "refinement" },
			{ "version", version },
			{ "pdbx_ordinal", software.size() + 1 } });

		best.parser->fixup();

		for (auto &cat1 : best.parser->mDb)
		{
			auto &cat2 = db[cat1.name()];

			// copy only the values in the first row for the following categories
			if (cat1.name() == "reflns" or cat1.name() == "refine")
			{
				if (cat2.empty())
					cat2.emplace(cat1.front());
				else
				{
						
					auto r1 = cat1.front();
					auto r2 = cat2.front();

					for (auto item : cat1.key_items())
						r2[item] = r1[item].text();
				}
			}
			else
			{
				for (auto rs : cat1)
					cat2.emplace(rs);
			}
		}
	}

	return result;
}

} // namespace pdbx