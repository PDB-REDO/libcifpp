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

#include "cif++/model.hpp"

#include "cif++.hpp"
#include "cif++/point.hpp"
#include "cif++/utilities.hpp"

#include <algorithm>
#include <exception>
#include <initializer_list>
#include <numeric>
#include <stack>
#include <stdexcept>
#include <utility>

namespace fs = std::filesystem;

namespace cif::mm
{

// --------------------------------------------------------------------
// atom

void atom::atom_impl::moveTo(const point &p)
{
	if (m_symop != "1_555")
		throw std::runtime_error("Moving symmetry copy");

	auto r = row();

	r.assign("Cartn_x", std::format("{:.3f}", p.m_x), false, false);
	r.assign("Cartn_y", std::format("{:.3f}", p.m_y), false, false);
	r.assign("Cartn_z", std::format("{:.3f}", p.m_z), false, false);

	m_location = p;
}

// const compound *compound() const;

std::string atom::atom_impl::get_property(std::string_view name) const
{
	return row()[name].as<std::string>();
}

int atom::atom_impl::get_property_int(std::string_view name) const
{
	int result = 0;
	if (not row()[name].empty())
	{
		auto s = get_property(name);

		std::from_chars_result r = std::from_chars(s.data(), s.data() + s.length(), result);
		if (r.ec != std::errc{} and VERBOSE > 0)
			std::cerr << "Error converting " << s << " to number for property " << name << '\n';
	}
	return result;
}

float atom::atom_impl::get_property_float(std::string_view name) const
{
	float result = 0;
	if (not row()[name].empty())
	{
		auto s = get_property(name);

		std::from_chars_result r = cif::from_chars(s.data(), s.data() + s.length(), result);
		if (r.ec != std::errc{} and VERBOSE > 0)
			std::cerr << "Error converting " << s << " to number for property " << name << '\n';
	}
	return result;
}

void atom::atom_impl::set_property(const std::string_view name, const std::string &value)
{
	auto r = row();
	if (not r)
		throw std::runtime_error("Trying to modify a row that does not exist");
	r.assign(name, value, true, true);
}

int atom::atom_impl::compare(const atom_impl &b) const
{
	int d = get_property("label_asym_id").compare(b.get_property("label_asym_id"));
	if (d == 0)
		d = get_property_int("label_seq_id") - b.get_property_int("label_seq_id");
	if (d == 0)
		d = get_property_int("auth_seq_id") - b.get_property_int("auth_seq_id");
	if (d == 0)
		d = get_property("label_atom_id").compare(b.get_property("label_atom_id"));

	return d;
}

// bool atom::atom_impl::getAnisoU(float anisou[6]) const
// {
// 	bool result = false;

// 	auto cat = m_db.get("atom_site_anisotrop");
// 	if (cat)
// 	{
// 		for (auto r : cat->find(key("id") == m_id))
// 		{
// 			tie(anisou[0], anisou[1], anisou[2], anisou[3], anisou[4], anisou[5]) =
// 				r.get("U[1][1]", "U[1][2]", "U[1][3]", "U[2][2]", "U[2][3]", "U[3][3]");
// 			result = true;
// 			break;
// 		}
// 	}

// 	return result;
// }

int atom::atom_impl::get_charge() const
{
	auto formalCharge = row()["pdbx_formal_charge"].as<std::optional<int>>();

	if (not formalCharge.has_value())
	{
		auto c = cif::compound_factory::instance().create(get_property("label_comp_id"));

		if (c != nullptr and c->atoms().size() == 1)
			formalCharge = c->atoms().front().charge;
	}

	return formalCharge.value_or(0);
}

std::ostream &operator<<(std::ostream &os, const atom &atom)
{
	if (atom.is_water())
		os << atom.get_label_comp_id() << ' ' << atom.get_label_asym_id() << ':' << atom.get_auth_seq_id() << ' ' << atom.get_label_atom_id();
	else
	{
		os << atom.get_label_comp_id() << ' ' << atom.get_label_asym_id() << ':' << atom.get_label_seq_id() << ' ' << atom.get_label_atom_id();

		if (atom.is_alternate())
			os << '(' << atom.get_label_alt_id() << ')';
		if (atom.get_auth_asym_id() != atom.get_label_asym_id() or atom.get_auth_seq_id() != std::to_string(atom.get_label_seq_id()) or atom.get_pdb_ins_code().empty() == false)
			os << " [" << atom.get_auth_asym_id() << ':' << atom.get_auth_seq_id() << atom.get_pdb_ins_code() << ']';
	}

	return os;
}

// --------------------------------------------------------------------
// residue

residue::residue(structure &structure, const std::vector<atom> &atoms)
	: m_structure(&structure)
{
	if (atoms.empty())
		throw std::runtime_error("Empty list of atoms");

	auto &a = atoms.front();

	m_compound_id = a.get_label_comp_id();
	m_asym_id = a.get_label_asym_id();
	m_seq_id = a.get_label_seq_id();
	m_pdb_strand_id = a.get_auth_asym_id();
	m_pdb_seq_num = a.get_auth_seq_id();
	m_pdb_ins_code = a.get_pdb_ins_code();

	for (auto atom : atoms)
		m_atoms.push_back(atom);
}

std::string residue::get_entity_id() const
{
	std::string result;

	if (not m_atoms.empty())
		result = m_atoms.front().get_label_entity_id();
	else if (m_structure != nullptr and not m_asym_id.empty())
	{
		using namespace literals;

		auto &db = m_structure->get_datablock();
		result = db["struct_asym"].find1<std::string>("id"_key == m_asym_id, "entity_id");
	}

	return result;
}

EntityType residue::entity_type() const
{
	assert(m_structure);
	return m_structure->get_entity_type_for_entity_id(get_entity_id());
}

void residue::add_atom(atom &atom)
{
	// update atom since it is now part of this residue
	m_atoms.push_back(atom);
}

atom residue::create_new_atom(atom_type inType, const std::string &inAtomID, point inLocation)
{
	auto &db = m_structure->get_datablock();
	auto &atom_site = db["atom_site"];

	auto ai = atom_site.emplace({
		{ "group_PDB", "HETATM" },
		{ "id", atom_site.get_unique_id("") },
		{ "type_symbol", atom_type_traits(inType).symbol() },
		{ "label_entity_id", get_entity_id() },
		{ "label_atom_id", inAtomID },
		{ "label_asym_id", m_asym_id },
		{ "label_alt_id", "." },
		{ "label_comp_id", m_compound_id },
		{ "label_seq_id", m_seq_id },
		{ "auth_asym_id", m_pdb_strand_id },
		{ "auth_atom_id", inAtomID },
		{ "auth_comp_id", m_compound_id },
		{ "auth_seq_id", m_pdb_seq_num },
		{ "occupancy", 1.0f, 2 },
		{ "B_iso_or_equiv", 20.0f },
		{ "pdbx_PDB_model_num", m_structure->get_model_nr() },
	});

	atom a(db, *ai);

	m_atoms.push_back(a);

	a.set_location(inLocation);

	return a;
}

std::vector<atom> residue::unique_atoms() const
{
	std::vector<atom> result;
	std::string firstAlt;

	for (auto &atom : m_atoms)
	{
		auto alt = atom.get_label_alt_id();
		if (alt.empty())
		{
			result.push_back(atom);
			continue;
		}

		if (firstAlt.empty())
			firstAlt = alt;
		else if (alt != firstAlt)
		{
			if (VERBOSE > 0)
				std::cerr << "skipping alternate atom " << atom << '\n';
			continue;
		}

		result.push_back(atom);
	}

	return result;
}

std::set<std::string> residue::get_alternate_ids() const
{
	std::set<std::string> result;

	for (auto a : m_atoms)
	{
		auto alt = a.get_label_alt_id();
		if (not alt.empty())
			result.insert(alt);
	}

	return result;
}

atom residue::get_atom_by_atom_id(const std::string &atom_id) const
{
	atom result;

	for (auto &a : m_atoms)
	{
		if (a.get_label_atom_id() == atom_id)
		{
			result = a;
			break;
		}
	}

	if (not result and VERBOSE > 1)
		std::cerr << "atom with atom_id " << atom_id << " not found in residue " << m_asym_id << ':' << m_seq_id << '\n';

	return result;
}

atom residue::get_atom_by_atom_id(const std::string &atomID, const std::string &altID) const
{
	if (altID.empty())
		return get_atom_by_atom_id(atomID);

	atom result;

	for (auto &a : m_atoms)
	{
		if (auto a_alt_id = a.get_label_alt_id(); a.get_label_atom_id() == atomID and (a_alt_id.empty() or a_alt_id == altID))
		{
			result = a;
			break;
		}
	}

	if (not result and VERBOSE > 1)
		std::cerr << "atom with atom_id " << atomID << " and alt_id " << altID << " not found in residue " << m_asym_id << ':' << m_seq_id << '\n';

	return result;
}

// residue is a single entity if the atoms for the asym with m_asym_id is equal
// to the number of atoms in this residue...  hope this is correct....
bool residue::is_entity() const
{
	auto &db = m_structure->get_datablock();

	auto a1 = db["atom_site"].find(key("label_asym_id") == m_asym_id);
	//	auto a2 = atoms();
	auto &a2 = m_atoms;

	return a1.size() == a2.size();
}

std::tuple<point, float> residue::center_and_radius() const
{
	std::vector<point> pts;
	for (auto &a : m_atoms)
		pts.push_back(a.get_location());

	return smallest_sphere_around_points(pts);
}

bool residue::has_alternate_atoms() const
{
	return std::ranges::find_if(m_atoms, [](const atom &atom)
			   { return atom.is_alternate(); }) != m_atoms.end();
}

bool residue::has_alternate_atoms_for(const std::string &atomID) const
{
	return std::ranges::find_if(m_atoms, [atomID](const atom &atom)
			   { return atom.get_label_atom_id() == atomID and atom.is_alternate(); }) != m_atoms.end();
}

std::set<std::string> residue::get_atom_ids() const
{
	std::set<std::string> ids;
	for (auto a : m_atoms)
		ids.insert(a.get_label_atom_id());

	return ids;
}

std::vector<atom> residue::get_atoms_by_id(const std::string &atom_id) const
{
	std::vector<atom> atoms;
	for (auto a : m_atoms)
	{
		if (a.get_label_atom_id() == atom_id)
			atoms.push_back(a);
	}
	return atoms;
}

std::ostream &operator<<(std::ostream &os, const residue &res)
{
	os << res.get_compound_id() << ' ' << res.get_asym_id() << ':' << res.get_seq_id();

	if (res.get_pdb_strand_id() != res.get_asym_id() or res.get_pdb_seq_num() != std::to_string(res.get_seq_id()))
		os << " [" << res.get_pdb_strand_id() << ':' << res.get_pdb_seq_num() << ']';

	return os;
}

// --------------------------------------------------------------------
// monomer

monomer::monomer(const polymer &polymer, std::size_t index, int seqID, const std::string &authSeqID, const std::string &pdbInsCode, const std::string &compoundID)
	: residue(*polymer.get_structure(), compoundID, polymer.get_asym_id(), seqID, polymer.get_pdb_strand_id(), authSeqID, pdbInsCode)
	, m_polymer(&polymer)
	, m_index(index)
{
}

bool monomer::is_first_in_chain() const
{
	return m_index == 0;
}

bool monomer::is_last_in_chain() const
{
	return m_index + 1 == m_polymer->size();
}

const monomer &monomer::prev() const
{
	return m_polymer->at(m_index - 1);
}

const monomer &monomer::next() const
{
	return m_polymer->at(m_index + 1);
}

bool monomer::has_alpha() const
{
	return m_index >= 1 and m_index + 2 < m_polymer->size();
}

bool monomer::has_kappa() const
{
	return m_index >= 2 and m_index + 2 < m_polymer->size();
}

float monomer::phi() const
{
	float result = 360;

	if (m_index > 0)
	{
		auto &prev = m_polymer->at(m_index - 1);
		if (prev.m_seq_id + 1 == m_seq_id)
		{
			auto a1 = prev.C();
			auto a2 = N();
			auto a3 = CAlpha();
			auto a4 = C();

			if (a1 and a2 and a3 and a4)
				result = dihedral_angle(a1.get_location(), a2.get_location(), a3.get_location(), a4.get_location());
		}
	}

	return result;
}

float monomer::psi() const
{
	float result = 360;

	if (m_index + 1 < m_polymer->size())
	{
		auto &next = m_polymer->at(m_index + 1);
		if (m_seq_id + 1 == next.m_seq_id)
		{
			auto a1 = N();
			auto a2 = CAlpha();
			auto a3 = C();
			auto a4 = next.N();

			if (a1 and a2 and a3 and a4)
				result = dihedral_angle(a1.get_location(), a2.get_location(), a3.get_location(), a4.get_location());
		}
	}

	return result;
}

float monomer::alpha() const
{
	float result = 360;

	try
	{
		if (m_index >= 1 and m_index + 2 < m_polymer->size())
		{
			auto &prev = m_polymer->at(m_index - 1);
			auto &next = m_polymer->at(m_index + 1);
			auto &nextNext = m_polymer->at(m_index + 2);

			result = static_cast<float>(dihedral_angle(prev.CAlpha().get_location(), CAlpha().get_location(), next.CAlpha().get_location(), nextNext.CAlpha().get_location()));
		}
	}
	catch (const std::exception &ex)
	{
		if (VERBOSE > 0)
			std::cerr << ex.what() << '\n';
	}

	return result;
}

float monomer::kappa() const
{
	float result = 360;

	try
	{
		if (m_index >= 2 and m_index + 2 < m_polymer->size())
		{
			auto &prevPrev = m_polymer->at(m_index - 2);
			auto &nextNext = m_polymer->at(m_index + 2);

			if (prevPrev.m_seq_id + 4 == nextNext.m_seq_id)
			{
				double ckap = cosinus_angle(CAlpha().get_location(), prevPrev.CAlpha().get_location(), nextNext.CAlpha().get_location(), CAlpha().get_location());
				double skap = std::sqrt(1 - ckap * ckap);
				result = static_cast<float>(std::atan2(skap, ckap) * 180 / kPI);
			}
		}
	}
	catch (const std::exception &ex)
	{
		if (VERBOSE > 0)
			std::cerr << "When trying to calculate kappa for " << m_asym_id << ':' << m_seq_id << ": "
					  << ex.what() << '\n';
	}

	return result;
}

float monomer::tco() const
{
	float result = 0.0;

	try
	{
		if (m_index > 0)
		{
			auto &prev = m_polymer->at(m_index - 1);
			if (prev.m_seq_id + 1 == m_seq_id)
				result = static_cast<float>(cosinus_angle(C().get_location(), O().get_location(), prev.C().get_location(), prev.O().get_location()));
		}
	}
	catch (const std::exception &ex)
	{
		if (VERBOSE > 0)
			std::cerr << "When trying to calculate tco for " << get_asym_id() << ':' << get_seq_id() << ": "
					  << ex.what() << '\n';
	}

	return result;
}

float monomer::omega() const
{
	float result = 360;

	try
	{
		if (not is_last_in_chain())
			result = omega(*this, m_polymer->at(m_index + 1));
	}
	catch (const std::exception &ex)
	{
		if (VERBOSE > 0)
			std::cerr << "When trying to calculate omega for " << get_asym_id() << ':' << get_seq_id() << ": "
					  << ex.what() << '\n';
	}

	return result;
}

const std::map<std::string, std::vector<std::string>> kChiAtomsMap = { // NOLINT(bugprone-throwing-static-initialization,cert-err58-cpp)
	{ "ASP", { "CG", "OD1" } },
	{ "ASN", { "CG", "OD1" } },
	{ "ARG", { "CG", "CD", "NE", "CZ" } },
	{ "HIS", { "CG", "ND1" } },
	{ "GLN", { "CG", "CD", "OE1" } },
	{ "GLU", { "CG", "CD", "OE1" } },
	{ "SER", { "OG" } },
	{ "THR", { "OG1" } },
	{ "LYS", { "CG", "CD", "CE", "NZ" } },
	{ "TYR", { "CG", "CD1" } },
	{ "PHE", { "CG", "CD1" } },
	{ "LEU", { "CG", "CD1" } },
	{ "TRP", { "CG", "CD1" } },
	{ "CYS", { "SG" } },
	{ "ILE", { "CG1", "CD1" } },
	{ "MET", { "CG", "SD", "CE" } },
	{ "MSE", { "CG", "SE", "CE" } },
	{ "PRO", { "CG", "CD" } },
	{ "VAL", { "CG1" } }
};

std::size_t monomer::nr_of_chis() const
{
	std::size_t result = 0;

	auto i = kChiAtomsMap.find(m_compound_id);
	if (i != kChiAtomsMap.end())
		result = i->second.size();

	return result;
}

float monomer::chi(std::size_t nr) const
{
	float result = 0;

	try
	{
		auto i = kChiAtomsMap.find(m_compound_id);
		if (i != kChiAtomsMap.end() and nr < i->second.size())
		{
			std::vector<std::string> atoms{ "N", "CA", "CB" };

			atoms.insert(atoms.end(), i->second.begin(), i->second.end());

			// in case we have a positive chiral volume we need to swap atoms
			if (chiral_volume() > 0)
			{
				if (m_compound_id == "LEU")
					atoms.back() = "CD2";
				if (m_compound_id == "VAL")
					atoms.back() = "CG2";
			}

			auto atom_0 = get_atom_by_atom_id(atoms[nr + 0]);
			auto atom_1 = get_atom_by_atom_id(atoms[nr + 1]);
			auto atom_2 = get_atom_by_atom_id(atoms[nr + 2]);
			auto atom_3 = get_atom_by_atom_id(atoms[nr + 3]);

			if (atom_0 and atom_1 and atom_2 and atom_3)
				result = static_cast<float>(dihedral_angle(
					atom_0.get_location(),
					atom_1.get_location(),
					atom_2.get_location(),
					atom_3.get_location()));
		}
	}
	catch (const std::exception &e)
	{
		if (VERBOSE > 0)
			std::cerr << e.what() << '\n';
		result = 0;
	}

	return result;
}

bool monomer::is_cis() const
{
	bool result = false;

	if (m_index + 1 < m_polymer->size())
	{
		auto &next = m_polymer->at(m_index + 1);

		result = monomer::is_cis(*this, next);
	}

	return result;
}

bool monomer::is_complete() const
{
	int seen = 0;
	for (auto &a : m_atoms)
	{
		if (a.get_label_atom_id() == "CA")
			seen |= 1;
		else if (a.get_label_atom_id() == "C")
			seen |= 2;
		else if (a.get_label_atom_id() == "N")
			seen |= 4;
		else if (a.get_label_atom_id() == "O")
			seen |= 8;
		// else if (a.get_label_atom_id() == "OXT")		seen |= 16;
	}
	return seen == 15;
}

bool monomer::has_alternate_backbone_atoms() const
{
	bool result = false;

	for (auto &a : m_atoms)
	{
		if (not a.is_alternate())
			continue;

		auto atom_id = a.get_label_atom_id();
		if (atom_id == "CA" or atom_id == "C" or atom_id == "N" or atom_id == "O")
		{
			result = true;
			break;
		}
	}

	return result;
}

float monomer::chiral_volume() const
{
	float result = 0;

	if (m_compound_id == "LEU")
	{
		auto centre = get_atom_by_atom_id("CG");
		auto atom1 = get_atom_by_atom_id("CB");
		auto atom2 = get_atom_by_atom_id("CD1");
		auto atom3 = get_atom_by_atom_id("CD2");

		result = dot_product(atom1.get_location() - centre.get_location(),
			cross_product(atom2.get_location() - centre.get_location(), atom3.get_location() - centre.get_location()));
	}
	else if (m_compound_id == "VAL")
	{
		auto centre = get_atom_by_atom_id("CB");
		auto atom1 = get_atom_by_atom_id("CA");
		auto atom2 = get_atom_by_atom_id("CG1");
		auto atom3 = get_atom_by_atom_id("CG2");

		result = dot_product(atom1.get_location() - centre.get_location(),
			cross_product(atom2.get_location() - centre.get_location(), atom3.get_location() - centre.get_location()));
	}

	return result;
}

bool monomer::are_bonded(const monomer &a, const monomer &b, float errorMargin)
{
	bool result = false;

	try
	{
		point atoms[4] = {
			a.get_atom_by_atom_id("CA").get_location(),
			a.get_atom_by_atom_id("C").get_location(),
			b.get_atom_by_atom_id("N").get_location(),
			b.get_atom_by_atom_id("CA").get_location()
		};

		auto distanceCACA = distance(atoms[0], atoms[3]);
		double omega = dihedral_angle(atoms[0], atoms[1], atoms[2], atoms[3]);

		bool cis = std::abs(omega) <= 30.0;
		float maxCACADistance = cis ? 3.0f : 3.8f;

		result = std::abs(distanceCACA - maxCACADistance) < errorMargin;
	}
	catch (const std::exception &ex)
	{
		if (VERBOSE > 2)
			std::cerr << "missing atoms in monomer::are_bonded: " << ex.what() << '\n';
	}

	return result;
}

float monomer::omega(const monomer &a, const monomer &b)
{
	float result = 360;

	auto a1 = a.get_atom_by_atom_id("CA");
	auto a2 = a.get_atom_by_atom_id("C");
	auto a3 = b.get_atom_by_atom_id("N");
	auto a4 = b.get_atom_by_atom_id("CA");

	if (a1 and a2 and a3 and a4)
		result = static_cast<float>(dihedral_angle(
			a1.get_location(),
			a2.get_location(),
			a3.get_location(),
			a4.get_location()));

	return result;
}

bool monomer::is_cis(const monomer &a, const monomer &b)
{
	return std::abs(omega(a, b)) < 30.0f;
}

atom monomer::create_new_atom(atom_type inType, const std::string &inAtomID, point inLocation)
{
	atom a = residue::create_new_atom(inType, inAtomID, inLocation);

	a.set_property("group_PDB", "ATOM");

	return a;
}

// --------------------------------------------------------------------
// polymer

polymer::polymer(structure &s, std::string entityID, std::string asym_id, std::string auth_asym_id)
	: m_structure(const_cast<structure *>(&s))
	, m_entity_id(std::move(entityID))
	, m_asym_id(std::move(asym_id))
	, m_pdb_strand_id(std::move(auth_asym_id))
{
	using namespace cif::literals;

	std::map<std::size_t, std::size_t> ix;

	auto &poly_seq_scheme = s.get_datablock()["pdbx_poly_seq_scheme"];
	reserve(poly_seq_scheme.size());

	for (auto r : poly_seq_scheme.find("asym_id"_key == m_asym_id))
	{
		int seqID;
		std::string compoundID, pdbSeqNum, pdbInsCode;
		cif::tie(seqID, pdbSeqNum, compoundID, pdbInsCode) = r.get("seq_id", "pdb_seq_num", "mon_id", "pdb_ins_code");

		std::size_t index = size();

		// store only the first
		if (not ix.count(seqID))
		{
			ix[seqID] = index;
			emplace_back(*this, index, seqID, pdbSeqNum, pdbInsCode, compoundID);
		}
		else if (VERBOSE > 0)
		{
			monomer m{ *this, index, seqID, pdbSeqNum, pdbInsCode, compoundID };
			std::cerr << "Dropping alternate residue " << m << '\n';
		}
	}
}

// std::string polymer::chainID() const
// {
// 	return mPolySeq.front()["pdb_strand_id"].as<std::string>();
// }

// monomer &polymer::getBySeqID(int seqID)
// {
// 	for (auto &m : *this)
// 		if (m.get_seq_id() == seqID)
// 			return m;

// 	throw std::runtime_error("monomer with seqID " + std::to_string(seqID) + " not found in polymer " + m_asym_id);
// }

// const monomer &polymer::getBySeqID(int seqID) const
// {
// 	for (auto &m : *this)
// 		if (m.get_seq_id() == seqID)
// 			return m;

// 	throw std::runtime_error("monomer with seqID " + std::to_string(seqID) + " not found in polymer " + m_asym_id);
// }

// int polymer::Distance(const monomer &a, const monomer &b) const
// {
// 	int result = std::numeric_limits<int>::max();

// 	if (a.get_asym_id() == b.get_asym_id())
// 	{
// 		int ixa = std::numeric_limits<int>::max(), ixb = std::numeric_limits<int>::max();

// 		int ix = 0, f = 0;
// 		for (auto &m : *this)
// 		{
// 			if (m.get_seq_id() == a.get_seq_id())
// 				ixa = ix, ++f;
// 			if (m.get_seq_id() == b.get_seq_id())
// 				ixb = ix, ++f;
// 			if (f == 2)
// 			{
// 				result = std::abs(ixa - ixb);
// 				break;
// 			}
// 		}
// 	}

// 	return result;
// }

// --------------------------------------------------------------------

sugar::sugar(branch &branch, const std::string &compoundID,
	const std::string &asym_id, int authSeqID)
	: residue(branch.get_structure(), compoundID, asym_id, 0, asym_id, std::to_string(authSeqID), "")
	, m_branch(&branch)
{
}

// bool sugar::hasLinkedSugarAtLeavingO(int leavingO) const
// {
// 	return false;
// }

// sugar &sugar::operator[](int leavingO)
// {
// 	throw std::logic_error("not implemented");
// }

// const sugar &sugar::operator[](int leavingO) const
// {
// 	throw std::logic_error("not implemented");
// }

std::string sugar::name() const
{
	std::string result;

	if (m_compound_id == "MAN")
		result += "alpha-D-mannopyranose";
	else if (m_compound_id == "BMA")
		result += "beta-D-mannopyranose";
	else if (m_compound_id == "NAG")
		result += "2-acetamido-2-deoxy-beta-D-glucopyranose";
	else if (m_compound_id == "NDG")
		result += "2-acetamido-2-deoxy-alpha-D-glucopyranose";
	else if (m_compound_id == "FUC")
		result += "alpha-L-fucopyranose";
	else if (m_compound_id == "FUL")
		result += "beta-L-fucopyranose";
	else
	{
		auto compound = compound_factory::instance().create(m_compound_id);
		if (compound)
			result += compound->name();
		else
			result += m_compound_id;
	}

	return result;
}

cif::mm::atom sugar::add_atom(row_initializer atom_info)
{
	auto &db = m_structure->get_datablock();
	auto &atom_site = db["atom_site"];

	auto atom_id = atom_site.get_unique_id("");

	atom_info.set_value({ "group_PDB", "HETATM" });
	atom_info.set_value({ "id", atom_id });
	atom_info.set_value({ "label_entity_id", m_branch->get_entity_id() });
	atom_info.set_value({ "label_asym_id", m_branch->get_asym_id() });
	atom_info.set_value({ "label_comp_id", m_compound_id });
	atom_info.set_value({ "label_seq_id", "." });
	atom_info.set_value({ "label_alt_id", "." });
	atom_info.set_value({ "auth_asym_id", m_branch->get_asym_id() });
	atom_info.set_value({ "auth_comp_id", m_compound_id });
	atom_info.set_value({ "auth_seq_id", m_pdb_seq_num });
	atom_info.set_value({ "occupancy", 1.0, 2 });
	atom_info.set_value({ "B_iso_or_equiv", 30.0, 2 });
	atom_info.set_value({ "pdbx_PDB_model_num", 1 });

	auto row = atom_site.emplace(std::move(atom_info));
	auto result = m_structure->emplace_atom(db, row);

	residue::add_atom(result);

	return result;
}

branch::branch(structure &structure, std::string asym_id, std::string entity_id)
	: m_structure(&structure)
	, m_asym_id(std::move(asym_id))
	, m_entity_id(std::move(entity_id))
{
	using namespace literals;

	auto &db = structure.get_datablock();
	auto &struct_asym = db["struct_asym"];
	auto &branch_scheme = db["pdbx_branch_scheme"];
	auto &branch_link = db["pdbx_entity_branch_link"];

	for (const auto &asym_entity_id : struct_asym.find<std::string>("id"_key == m_asym_id, "entity_id"))
	{
		for (const auto &[comp_id, num] : branch_scheme.find<std::string, int>(
				 "asym_id"_key == m_asym_id, "mon_id", "pdb_seq_num"))
		{
			emplace_back(*this, comp_id, m_asym_id, num);
		}

		for (const auto &[num1, num2, atom1, atom2] : branch_link.find<std::size_t, std::size_t, std::string, std::string>(
				 "entity_id"_key == asym_entity_id, "entity_branch_list_num_1", "entity_branch_list_num_2", "atom_id_1", "atom_id_2"))
		{
			// if (not iequals(atom1, "c1"))
			// 	throw std::runtime_error("invalid pdbx_entity_branch_link");

			auto &s1 = at(num1 - 1);
			auto &s2 = at(num2 - 1);

			s1.set_link(s2.get_atom_by_atom_id(atom2));
		}

		break;
	}
}

void branch::link_atoms()
{
	if (not empty())
	{
		using namespace literals;

		auto &db = m_structure->get_datablock();
		auto &branch_link = db["pdbx_entity_branch_link"];

		auto entity_id = front().get_entity_id();

		for (const auto &[num1, num2, atom1, atom2] : branch_link.find<std::size_t, std::size_t, std::string, std::string>(
				 "entity_id"_key == entity_id, "entity_branch_list_num_1", "entity_branch_list_num_2", "atom_id_1", "atom_id_2"))
		{
			// if (not iequals(atom1, "c1"))
			// 	throw std::runtime_error("invalid pdbx_entity_branch_link");

			auto &s1 = at(num1 - 1);
			auto &s2 = at(num2 - 1);

			s1.set_link(s2.get_atom_by_atom_id(atom2));
		}
	}
}

sugar &branch::get_sugar_by_num(int nr)
{
	auto i = std::ranges::find_if(*this, [nr](const sugar &s)
		{ return s.num() == nr; });
	if (i == end())
		throw std::out_of_range("Sugar with num " + std::to_string(nr) + " not found in branch " + m_asym_id);

	return *i;
}

std::string branch::name() const
{
	return empty() ? "" : name(front());
}

sugar &branch::construct_sugar(const std::string &compound_id)
{
	auto &db = m_structure->get_datablock();

	auto compound = compound_factory::instance().create(compound_id);
	if (compound == nullptr)
		throw std::runtime_error("Trying to insert unknown compound " + compound_id + " (not found in CCD)");

	auto &chemComp = db["chem_comp"];
	auto r = chemComp.find(key("id") == compound_id);
	if (r.empty())
	{
		chemComp.emplace({ { "id", compound_id },
			{ "name", compound->name() },
			{ "formula", compound->formula() },
			{ "formula_weight", compound->formula_weight() },
			{ "type", compound->type() } });
	}

	sugar &result = emplace_back(*this, compound_id, m_asym_id, static_cast<int>(size() + 1));

	db["pdbx_branch_scheme"].emplace({ { "asym_id", result.get_asym_id() },
		{ "entity_id", result.get_entity_id() },
		{ "num", result.num() },
		{ "mon_id", result.get_compound_id() },

		{ "pdb_asym_id", result.get_asym_id() },
		{ "pdb_seq_num", result.get_pdb_seq_num() },
		{ "pdb_mon_id", result.get_compound_id() },

		{ "auth_asym_id", result.get_pdb_strand_id() },
		{ "auth_mon_id", result.get_compound_id() },
		{ "auth_seq_num", result.get_pdb_seq_num() },

		{ "hetero", "n" } });

	return result;
}

sugar &branch::construct_sugar(const std::string &compound_id, const std::string &atom_id,
	int linked_sugar_nr, const std::string &linked_atom_id)
{
	auto &result = construct_sugar(compound_id);

	auto &linked = get_sugar_by_num(linked_sugar_nr);
	result.set_link(linked.get_atom_by_atom_id(linked_atom_id));

	auto &db = m_structure->get_datablock();

	auto &pdbx_entity_branch_link = db["pdbx_entity_branch_link"];
	auto linkID = pdbx_entity_branch_link.get_unique_id("");

	db["pdbx_entity_branch_link"].emplace({ { "link_id", linkID },
		{ "entity_id", get_entity_id() },
		{ "entity_branch_list_num_1", result.num() },
		{ "comp_id_1", compound_id },
		{ "atom_id_1", atom_id },
		{ "leaving_atom_id_1", "O1" }, /// TODO: Need to fix this!
		{ "entity_branch_list_num_2", linked.num() },
		{ "comp_id_2", linked.get_compound_id() },
		{ "atom_id_2", linked_atom_id },
		{ "leaving_atom_id_2", "." },
		{ "value_order", "sing" } });

	return result;
}

std::string branch::name(const sugar &s) const
{
	using namespace literals;

	std::string result;

	for (auto &sn : *this)
	{
		if (not sn.get_link() or sn.get_link().get_auth_seq_id() != s.get_pdb_seq_num())
			continue;

		auto n = name(sn) + "-(1-" + sn.get_link().get_label_atom_id().substr(1) + ')';

		result = result.empty() ? n : result + "-[" + n + ']';
	}

	if (not result.empty() and result.back() != ']')
		result += '-';

	return result + s.name();
}

float branch::weight() const
{
	return std::accumulate(begin(), end(), 0.f, [](float sum, const sugar &s)
		{
		auto compound = compound_factory::instance().create(s.get_compound_id());
		if (compound)
			sum += compound->formula_weight();
		return sum; });
}

// --------------------------------------------------------------------
//	structure

structure::structure(file &p, std::size_t modelNr, structure_open_options options)
	: structure(p.front(), modelNr, options)
{
}

structure::structure(datablock &db, std::size_t modelNr, structure_open_options options)
	: m_db(db)
	, m_model_nr(modelNr)
{
	if (db.get_validator() == nullptr)
		db.load_dictionary();

	auto &atom_site = db["atom_site"];

	load_atoms_for_model(options);

	// Check to see if we should actually load another model?
	if (m_atoms.empty() and m_model_nr == 1)
	{
		std::optional<std::size_t> model_nr;
		cif::tie(model_nr) = atom_site.front().get("pdbx_PDB_model_num");
		if (model_nr and *model_nr != m_model_nr)
		{
			if (VERBOSE > 0)
				std::cerr << "No atoms loaded for model 1, trying model " << *model_nr << '\n';
			m_model_nr = *model_nr;
			load_atoms_for_model(options);
		}
	}

	if (m_atoms.empty())
	{
		if (VERBOSE >= 0)
			std::cerr << "Warning: no atoms loaded\n";
	}
	else
		load_data();
}

void structure::load_atoms_for_model(structure_open_options options)
{
	using namespace literals;

	auto &atom_site = m_db["atom_site"];

	condition c = "pdbx_PDB_model_num"_key == null or "pdbx_PDB_model_num"_key == m_model_nr;

	if (options.skip_hydrogen)
		c = std::move(c) and (cif::key("type_symbol") != "H" and cif::key("type_symbol") != "D");

	if (options.skip_water)
		c = std::move(c) and (cif::key("auth_comp_id") != "HOH" and cif::key("auth_comp_id") != "H20" and cif::key("auth_comp_id") != "WAT");

	if (options.skip_hetatom)
	{
		if (options.skip_water)
			c = std::move(c) and cif::key("group_PDB") != "HETATM";
		else
			c = std::move(c) and (cif::key("group_PDB") != "HETATM" or (cif::key("auth_comp_id") == "HOH" or cif::key("auth_comp_id") == "H20" or cif::key("auth_comp_id") == "WAT"));
	}

	if (options.min_b_factor.has_value())
		c = std::move(c) and cif::key("B_iso_or_equiv") >= *options.min_b_factor;

	if (options.max_b_factor.has_value())
		c = std::move(c) and cif::key("B_iso_or_equiv") <= *options.max_b_factor;

	if (not options.asyms.empty())
	{
		condition tmp_c;
		for (auto asym_id : options.asyms)
			tmp_c = std::move(tmp_c) or cif::key("label_asym_id") == asym_id;
		c = std::move(c) and std::move(tmp_c);
	}

	if (options.occupancy_mode == occupancy_policy::ALL)
	{
		for (auto id : atom_site.find<std::string>(std::move(c), "id"))
			emplace_atom(std::make_shared<atom::atom_impl>(m_db, id));
	}
	else if (options.occupancy_mode == occupancy_policy::UNOCCUPIED)
	{
		for (auto id : atom_site.find<std::string>(std::move(c), "id"))
		{
			auto a = std::make_shared<atom::atom_impl>(m_db, id);
			if (a->get_property_float("occupancy") > 0)
				continue;
			emplace_atom(a);
		}
	}
	else
	{
		std::vector<cif::mm::atom> atoms;
		std::map<std::tuple<std::string, int>, std::map<std::string, float>> alts;

		for (auto id : atom_site.find<std::string>(std::move(c), "id"))
		{
			auto a = atoms.emplace_back(std::make_shared<atom::atom_impl>(m_db, id));

			if (a.is_alternate())
			{
				auto key = std::make_tuple(a.get_label_asym_id(), a.get_label_seq_id());
				auto alt_id = a.get_label_alt_id();

				if (auto i = alts.find(key); i != alts.end())
					i->second[alt_id] += a.get_occupancy();
				else
					alts[key][alt_id] = a.get_occupancy();
			}
		}

		for (auto &&[key, value] : alts)
		{
			// const auto &[asym_id, seq_id] = key;

			// select highest occupancy for this residue's alternates
			std::string alt_id;
			float occupancy = options.occupancy_mode == occupancy_policy::MAX ? 0.f : std::numeric_limits<float>::max();
			for (const auto &[alt_key, alt_value] : value)
			{
				if (options.occupancy_mode == occupancy_policy::MAX)
				{
					if (occupancy < alt_value)
					{
						alt_id = alt_key;
						occupancy = alt_value;
					}
				}
				else
				{
					if (occupancy > alt_value)
					{
						alt_id = alt_key;
						occupancy = alt_value;
					}
				}
			}

			value.clear();
			value.emplace(alt_id, occupancy);
		}

		for (auto a : atoms)
		{
			if (a.is_alternate())
			{
				auto key = std::make_tuple(a.get_label_asym_id(), a.get_label_seq_id());

				if (alts[key].contains(a.get_label_alt_id()))
					emplace_atom(a);
			}
			else
				emplace_atom(a);
		}
	}
}

void structure::load_data()
{
	auto &polySeqScheme = m_db["pdbx_poly_seq_scheme"];

	for (const auto &[asym_id, auth_asym_id, entityID] : polySeqScheme.rows<std::string, std::string, std::string>("asym_id", "pdb_strand_id", "entity_id"))
	{
		if (m_polymers.empty() or m_polymers.back().get_asym_id() != asym_id)
			m_polymers.emplace_back(*this, entityID, asym_id, auth_asym_id);
	}

	auto &branchScheme = m_db["pdbx_branch_scheme"];

	for (const auto &[asym_id, entity_id] : branchScheme.rows<std::string, std::string>("asym_id", "entity_id"))
	{
		if (m_branches.empty() or m_branches.back().get_asym_id() != asym_id)
			m_branches.emplace_back(*this, asym_id, entity_id);
	}

	auto &nonPolyScheme = m_db["pdbx_nonpoly_scheme"];

	for (const auto &[asym_id, monID, pdbStrandID, pdbSeqNum, pdbInsCode] :
		nonPolyScheme.rows<std::string, std::string, std::string, std::string, std::string>("asym_id", "mon_id", "pdb_strand_id", "pdb_seq_num", "pdb_ins_code"))
		m_non_polymers.emplace_back(*this, monID, asym_id, 0, pdbStrandID, pdbSeqNum, pdbInsCode);

	// place atoms in residues

	using key_type = std::tuple<std::string, int, std::string>;
	std::map<key_type, residue *> resMap;

	for (auto &poly : m_polymers)
	{
		for (auto &res : poly)
			resMap[{ res.get_asym_id(), res.get_seq_id(), res.get_pdb_seq_num() }] = &res;
	}

	for (auto &res : m_non_polymers)
		resMap[{ res.get_asym_id(), res.get_seq_id(), res.get_pdb_seq_num() }] = &res;

	std::set<std::string> sugars;
	for (auto &branch : m_branches)
	{
		for (auto &sugar : branch)
		{
			resMap[{ sugar.get_asym_id(), sugar.get_seq_id(), sugar.get_pdb_seq_num() }] = &sugar;
			sugars.insert(sugar.get_compound_id());
		}
	}

	for (auto &atom : m_atoms)
	{
		key_type k(atom.get_label_asym_id(), atom.get_label_seq_id(), atom.get_auth_seq_id());
		auto ri = resMap.find(k);

		if (ri == resMap.end())
		{
			if (VERBOSE > 0)
				std::cerr << "Missing residue for atom " << atom << '\n';

			// see if it might match a non poly
			for (auto &res : m_non_polymers)
			{
				if (res.get_asym_id() != atom.get_label_asym_id())
					continue;

				res.add_atom(atom);
				break;
			}

			continue;
		}

		ri->second->add_atom(atom);
	}

	// what the ...
	std::erase_if(m_branches, [](const branch &b)
		{ return b.empty(); });

	for (auto &branch : m_branches)
		branch.link_atoms();
}

EntityType structure::get_entity_type_for_entity_id(const std::string entityID) const
{
	using namespace literals;

	auto &entity = m_db["entity"];
	auto entity_type = entity.find1<std::string>("id"_key == entityID, "type");

	EntityType result;

	if (iequals(entity_type, "polymer"))
		result = EntityType::Polymer;
	else if (iequals(entity_type, "non-polymer"))
		result = EntityType::NonPolymer;
	else if (iequals(entity_type, "macrolide"))
		result = EntityType::Macrolide;
	else if (iequals(entity_type, "water"))
		result = EntityType::Water;
	else if (iequals(entity_type, "branched"))
		result = EntityType::Branched;
	else
		throw std::runtime_error("Unknown entity type " + entity_type);

	return result;
}

EntityType structure::get_entity_type_for_asym_id(const std::string asym_id) const
{
	using namespace literals;

	auto &struct_asym = m_db["struct_asym"];
	auto entityID = struct_asym.find1<std::string>("id"_key == asym_id, "entity_id");

	return get_entity_type_for_entity_id(entityID);
}

bool structure::has_atom_id(const std::string &id) const
{
	assert(m_atoms.size() == m_atom_index.size());

	bool result = false;

	int L = 0, R = static_cast<int>(m_atoms.size() - 1);
	while (L <= R)
	{
		int i = (L + R) / 2;

		const atom &atom = m_atoms[m_atom_index[i]];

		int d = atom.id().compare(id);

		if (d == 0)
		{
			result = true;
			break;
		}

		if (d < 0)
			L = i + 1;
		else
			R = i - 1;
	}

	return result;
}

atom structure::get_atom_by_id(const std::string &id) const
{
	assert(m_atoms.size() == m_atom_index.size());

	int L = 0, R = static_cast<int>(m_atoms.size() - 1);
	while (L <= R)
	{
		int i = (L + R) / 2;

		const atom &atom = m_atoms[m_atom_index[i]];

		int d = atom.id().compare(id);

		if (d == 0)
			return atom;

		if (d < 0)
			L = i + 1;
		else
			R = i - 1;
	}

	throw std::out_of_range("Could not find atom with id " + id);
}

atom structure::get_atom_by_label(const std::string &atom_id, const std::string &asym_id, const std::string &compID, int seqID, const std::string &altID)
{
	for (auto &a : m_atoms)
	{
		if (a.get_label_atom_id() == atom_id and
			a.get_label_asym_id() == asym_id and
			a.get_label_comp_id() == compID and
			a.get_label_seq_id() == seqID and
			a.get_label_alt_id() == altID)
		{
			return a;
		}
	}

	throw std::out_of_range("Could not find atom with specified label");
}

atom structure::get_atom_by_position(point p) const
{
	double dist = std::numeric_limits<double>::max();
	std::size_t index = std::numeric_limits<std::size_t>::max();

	for (std::size_t i = 0; i < m_atoms.size(); ++i)
	{
		auto &a = m_atoms.at(i);

		auto d = distance(a.get_location(), p);
		if (d < dist)
		{
			dist = d;
			index = i;
		}
	}

	if (index < m_atoms.size())
		return m_atoms.at(index);

	return {};
}

atom structure::get_atom_by_position_and_type(point p, std::string_view type, std::string_view res_type) const
{
	double dist = std::numeric_limits<double>::max();
	std::size_t index = std::numeric_limits<std::size_t>::max();

	for (std::size_t i = 0; i < m_atoms.size(); ++i)
	{
		auto &a = m_atoms.at(i);

		if (a.get_label_comp_id() != res_type)
			continue;

		if (a.get_label_atom_id() != type)
			continue;

		auto d = distance(a.get_location(), p);
		if (dist > d)
		{
			dist = d;
			index = i;
		}
	}

	if (index < m_atoms.size())
		return m_atoms.at(index);

	return {};
}

polymer &structure::get_polymer_by_asym_id(const std::string &asym_id)
{
	for (auto &poly : m_polymers)
	{
		if (poly.get_asym_id() != asym_id)
			continue;

		return poly;
	}

	throw std::runtime_error("polymer with asym id " + asym_id + " not found");
}

residue &structure::create_residue(const std::vector<atom> &atoms)
{
	return m_non_polymers.emplace_back(*this, atoms);
}

residue &structure::get_residue(const std::string &asym_id, int seqID, const std::string &authSeqID)
{
	if (seqID == 0)
	{
		for (auto &res : m_non_polymers)
		{
			if (res.get_asym_id() == asym_id and (authSeqID.empty() or res.get_pdb_seq_num() == authSeqID))
				return res;
		}
	}

	for (auto &poly : m_polymers)
	{
		if (poly.get_asym_id() != asym_id)
			continue;

		for (auto &res : poly)
		{
			if (res.get_seq_id() == seqID)
				return res;
		}
	}

	for (auto &branch : m_branches)
	{
		if (branch.get_asym_id() != asym_id)
			continue;

		for (auto &sugar : branch)
		{
			if (sugar.get_asym_id() == asym_id and sugar.get_pdb_seq_num() == authSeqID)
				return sugar;
		}
	}

	std::string desc = asym_id;

	if (seqID != 0)
		desc += "/" + std::to_string(seqID);

	if (not authSeqID.empty())
		desc += "-" + authSeqID;

	throw std::out_of_range("Could not find residue " + desc);
}

residue &structure::get_residue(const std::string &asym_id, const std::string &compID, int seqID, const std::string &authSeqID)
{
	if (seqID == 0)
	{
		for (auto &res : m_non_polymers)
		{
			if (res.get_asym_id() == asym_id and res.get_pdb_seq_num() == authSeqID and res.get_compound_id() == compID)
				return res;
		}
	}

	for (auto &poly : m_polymers)
	{
		if (poly.get_asym_id() != asym_id)
			continue;

		for (auto &res : poly)
		{
			if (res.get_seq_id() == seqID and res.get_compound_id() == compID)
				return res;
		}
	}

	for (auto &branch : m_branches)
	{
		if (branch.get_asym_id() != asym_id)
			continue;

		for (auto &sugar : branch)
		{
			if (sugar.get_asym_id() == asym_id and sugar.get_pdb_seq_num() == authSeqID and sugar.get_compound_id() == compID)
				return sugar;
		}
	}

	std::string desc = asym_id;

	if (seqID != 0)
		desc += "/" + std::to_string(seqID);

	if (not authSeqID.empty())
		desc += "-" + authSeqID;

	throw std::out_of_range("Could not find residue " + desc + " of type " + compID);
}

branch &structure::get_branch_by_asym_id(const std::string &asym_id)
{
	for (auto &branch : m_branches)
	{
		if (branch.get_asym_id() == asym_id)
			return branch;
	}

	throw std::runtime_error("branch not found for asym id " + asym_id);
}

std::string structure::insert_compound(const std::string &compoundID, bool is_entity)
{
	using namespace literals;

	auto compound = compound_factory::instance().create(compoundID);
	if (compound == nullptr)
		throw std::runtime_error("Trying to insert unknown compound " + compoundID + " (not found in CCD)");

	auto &chemComp = m_db["chem_comp"];
	auto r = chemComp.find(key("id") == compoundID);
	if (r.empty())
	{
		chemComp.emplace({ { "id", compoundID },
			{ "name", compound->name() },
			{ "formula", compound->formula() },
			{ "formula_weight", compound->formula_weight() },
			{ "type", compound->type() } });
	}

	std::string entity_id;

	if (is_entity)
	{
		auto &pdbxEntityNonpoly = m_db["pdbx_entity_nonpoly"];

		entity_id = pdbxEntityNonpoly.find_first<std::string>("comp_id"_key == compoundID, "entity_id");

		if (entity_id.empty())
		{
			auto &entity = m_db["entity"];
			entity_id = entity.get_unique_id("");

			entity.emplace({ { "id", entity_id },
				{ "type", "non-polymer" },
				{ "pdbx_description", compound->name() },
				{ "formula_weight", compound->formula_weight() } });

			pdbxEntityNonpoly.emplace({ { "entity_id", entity_id },
				{ "name", compound->name() },
				{ "comp_id", compoundID } });
		}
	}

	return entity_id;
}

// --------------------------------------------------------------------

atom &structure::emplace_atom(atom &&atom)
{
	int L = 0, R = static_cast<int>(m_atom_index.size() - 1);
	while (L <= R)
	{
		int i = (L + R) / 2;

		auto &ai = m_atoms[m_atom_index[i]];

		int d = ai.id().compare(atom.id());

		if (d == 0)
			throw std::runtime_error("Duplicate atom ID " + atom.id());

		if (d < 0)
			L = i + 1;
		else
			R = i - 1;
	}

	if (R == -1) // msvc...
		m_atom_index.insert(m_atom_index.begin(), m_atoms.size());
	else
		m_atom_index.insert(m_atom_index.begin() + R + 1, m_atoms.size());

	// make sure the atom_type is known
	auto &atom_type = m_db["atom_type"];
	std::string symbol = atom.get_property("type_symbol");

	using namespace cif::literals;
	if (not atom_type.contains("symbol"_key == symbol))
		atom_type.emplace({ { "symbol", symbol } });

	return m_atoms.emplace_back(std::move(atom));
}

void structure::remove_atom(atom &a, bool removeFromResidue)
{
	using namespace literals;

	auto &atomSite = m_db["atom_site"];

	if (a.is_water())
	{
		auto ra = atomSite.find1("id"_key == a.id());
		if (ra)
		{
			auto &nps = m_db["pdbx_nonpoly_scheme"];
			for (auto rnp : atomSite.get_children(ra, nps))
				nps.erase(rnp);
		}
	}
	else if (removeFromResidue)
	{
		try
		{
			auto &res = get_residue(a);
			std::erase(res.m_atoms, a);
		}
		catch (const std::exception &ex)
		{
			if (VERBOSE > 0)
				std::cerr << "Error removing atom from residue: " << ex.what() << '\n';
		}
	}

	for (auto ri : atomSite.find("id"_key == a.id()))
	{
		// also remove struct_conn records for this atom
		auto &structConn = m_db["struct_conn"];

		condition cond;

		for (std::string prefix : { "ptnr1_", "ptnr2_", "pdbx_ptnr3_" })
		{
			if (a.get_label_seq_id() == 0)
				cond = std::move(cond) or (cif::key(prefix + "label_asym_id") == a.get_label_asym_id() and
											  cif::key(prefix + "label_seq_id") == null and
											  cif::key(prefix + "auth_seq_id") == a.get_auth_seq_id() and
											  cif::key(prefix + "label_atom_id") == a.get_label_atom_id());
			else
				cond = std::move(cond) or (cif::key(prefix + "label_asym_id") == a.get_label_asym_id() and
											  cif::key(prefix + "label_seq_id") == a.get_label_seq_id() and
											  cif::key(prefix + "auth_seq_id") == a.get_auth_seq_id() and
											  cif::key(prefix + "label_atom_id") == a.get_label_atom_id());
		}

		if (cond)
			structConn.erase(std::move(cond));

		atomSite.erase(ri);
		break;
	}

	assert(m_atom_index.size() == m_atoms.size());

#ifndef NDEBUG
	bool removed = false;
#endif

	int L = 0, R = static_cast<int>(m_atom_index.size() - 1);
	while (L <= R)
	{
		int i = (L + R) / 2;

		const atom &atom = m_atoms[m_atom_index[i]];

		int d = atom.id().compare(a.id());

		if (d == 0)
		{
			m_atoms.erase(m_atoms.begin() + m_atom_index[i]);

			auto ai = m_atom_index[i];
			m_atom_index.erase(m_atom_index.begin() + i);

			for (auto &j : m_atom_index)
			{
				if (j > ai)
					--j;
			}
#ifndef NDEBUG
			removed = true;
#endif
			break;
		}

		if (d < 0)
			L = i + 1;
		else
			R = i - 1;
	}
#ifndef NDEBUG
	assert(removed);
#endif
}

void structure::swap_atoms(atom a1, atom a2)
{
	auto &atomSites = m_db["atom_site"];

	try
	{
		auto r1 = atomSites.find1(key("id") == a1.id());
		auto r2 = atomSites.find1(key("id") == a2.id());

		for (std::string fld : std::initializer_list<std::string>{ "label_atom_id", "auth_atom_id", "type_symbol" })
		{
			auto l1 = r1[fld];
			auto l2 = r2[fld];
			l1.swap(l2);
		}
	}
	catch (const std::exception &ex)
	{
		std::throw_with_nested(std::runtime_error("Failed to swap atoms"));
	}
}

void structure::move_atom(atom a, point p)
{
	a.set_location(p);
}

void structure::change_residue(residue &res, const std::string &newCompound,
	const std::vector<std::tuple<std::string, std::string>> &remappedAtoms)
{
	using namespace literals;

	std::string asym_id = res.get_asym_id();

	const auto compound = compound_factory::instance().create(newCompound);
	if (not compound)
		throw std::runtime_error("Unknown compound " + newCompound);

	// First make sure the compound is already known or insert it.
	// And if the residue is an entity, we must make sure it exists

	std::string entityID;
	if (res.is_entity())
	{
		// create a copy of the entity first
		auto &entity = m_db["entity"];

		entityID = entity.find_first<std::string>("type"_key == "non-polymer" and "pdbx_description"_key == compound->name(), "id");

		if (entityID.empty())
		{
			entityID = entity.get_unique_id("");
			entity.emplace({ { "id", entityID },
				{ "type", "non-polymer" },
				{ "pdbx_description", compound->name() },
				{ "formula_weight", compound->formula_weight() } });

			auto &pdbxEntityNonpoly = m_db["pdbx_entity_nonpoly"];
			pdbxEntityNonpoly.emplace({ { "entity_id", entityID },
				{ "name", compound->name() },
				{ "comp_id", newCompound } });
		}

		auto &pdbxNonPolyScheme = m_db["pdbx_nonpoly_scheme"];
		for (auto nps : pdbxNonPolyScheme.find("asym_id"_key == asym_id))
		{
			nps.assign("mon_id", newCompound, true);
			nps.assign("pdb_mon_id", newCompound, true);
			nps.assign("auth_mon_id", newCompound, true);
			nps.assign("entity_id", entityID, true);
		}

		// create rest
		auto &chemComp = m_db["chem_comp"];
		if (not chemComp.contains(key("id") == newCompound))
		{
			chemComp.emplace({ { "id", newCompound },
				{ "name", compound->name() },
				{ "formula", compound->formula() },
				{ "formula_weight", compound->formula_weight() },
				{ "type", compound->type() } });
		}

		// update the struct_asym for the new entity
		m_db["struct_asym"].update_value("id"_key == asym_id, "entity_id", entityID);
	}
	else
		insert_compound(newCompound, false);

	res.set_compound_id(newCompound);

	auto &atomSites = m_db["atom_site"];
	auto atoms = res.atoms();

	for (const auto &[a1, a2] : remappedAtoms)
	{
		auto i = std::ranges::find_if(atoms, [id = a1](const atom &a)
			{ return a.get_label_atom_id() == id; });
		if (i == atoms.end())
		{
			if (VERBOSE >= 0)
				std::cerr << "Missing atom for atom ID " << a1 << '\n';
			continue;
		}

		auto r = atomSites.find(key("id") == i->id());

		if (r.size() != 1)
			continue;

		if (a2.empty() or a2 == ".")
		{
			i->set_property("label_comp_id", newCompound);
			remove_atom(*i);
		}
		else if (a1 != a2)
		{
			auto ra = r.front();
			ra["label_atom_id"] = a2;
			ra["auth_atom_id"] = a2;
			ra["type_symbol"] = atom_type_traits(compound->get_atom_by_atom_id(a2).type_symbol).symbol();
		}
	}

	for (auto a : atoms)
	{
		atomSites.update_value(key("id") == a.id(), "label_comp_id", newCompound);
		atomSites.update_value(key("id") == a.id(), "auth_comp_id", newCompound);
	}
}

void structure::remove_residue(const std::string &asym_id, int seq_id, const std::string &auth_seq_id)
{
	if (seq_id == 0)
	{
		for (auto &res : m_non_polymers)
		{
			if (res.get_asym_id() == asym_id and (auth_seq_id.empty() or res.get_pdb_seq_num() == auth_seq_id))
			{
				remove_residue(res);
				return;
			}
		}
	}

	for (auto &poly : m_polymers)
	{
		if (poly.get_asym_id() != asym_id)
			continue;

		for (auto &res : poly)
		{
			if (res.get_seq_id() == seq_id)
			{
				remove_residue(res);
				return;
			}
		}
	}

	for (auto &branch : m_branches)
	{
		if (branch.get_asym_id() != asym_id)
			continue;

		for (auto &sugar : branch)
		{
			if (sugar.get_asym_id() == asym_id and sugar.get_pdb_seq_num() == auth_seq_id)
			{
				remove_residue(sugar);
				return;
			}
		}
	}
}

void structure::remove_residue(residue &res)
{
	using namespace literals;

	auto atoms = res.atoms();

	switch (res.entity_type())
	{
		case EntityType::Polymer:
		{
			auto &m = dynamic_cast<monomer &>(res);

			m_db["pdbx_poly_seq_scheme"].erase(
				"asym_id"_key == res.get_asym_id() and
				"seq_id"_key == res.get_seq_id());

			for (auto &poly : m_polymers)
				std::erase(poly, m);
			break;
		}

		case EntityType::NonPolymer:
			m_db["pdbx_nonpoly_scheme"].erase("asym_id"_key == res.get_asym_id());
			m_db["struct_asym"].erase("id"_key == res.get_asym_id());
			std::erase(m_non_polymers, res);
			break;

		case EntityType::Water:
			m_db["pdbx_nonpoly_scheme"].erase("asym_id"_key == res.get_asym_id());
			std::erase(m_non_polymers, res);
			break;

		case EntityType::Branched:
		{
			auto &s = dynamic_cast<sugar &>(res);

			remove_sugar(s);

			atoms.clear();
			break;
		}

		case EntityType::Macrolide:
			// TODO: Fix this?
			throw std::runtime_error("no support for macrolides yet");
	}

	for (auto atom : atoms)
		remove_atom(atom, false);
}

void structure::remove_sugar(sugar &s)
{
	using namespace literals;

	std::string asym_id = s.get_asym_id();
	branch &branch = get_branch_by_asym_id(asym_id);
	auto si = std::ranges::find(branch, s);
	if (si == branch.end())
		throw std::runtime_error("sugar not part of branch");
	std::size_t six = si - branch.begin();

	if (six == 0) // first sugar, means the death of this branch
		remove_branch(branch);
	else
	{
		std::set<std::size_t> dix;
		std::stack<std::size_t> test;
		test.push(s.num());

		while (not test.empty())
		{
			auto tix = test.top();
			test.pop();

			if (dix.count(tix))
				continue;

			dix.insert(tix);

			for (auto &s2 : branch)
			{
				if (s2.get_link_nr() == tix)
					test.push(s2.num());
			}

			for (auto atom : branch[tix - 1].atoms())
				remove_atom(atom, false);
		}

		std::erase_if(branch, [dix](const sugar &s)
			{ return dix.count(s.num()); });

		auto entity_id = create_entity_for_branch(branch);

		// Update the entity id of the asym
		auto &struct_asym = m_db["struct_asym"];
		auto r = struct_asym.find1("id"_key == asym_id);
		r["entity_id"] = entity_id;

		for (auto &sugar : branch)
		{
			for (auto atom : sugar.atoms())
				atom.set_property("label_entity_id", entity_id);
		}

		auto &pdbx_branch_scheme = m_db["pdbx_branch_scheme"];
		pdbx_branch_scheme.erase("asym_id"_key == asym_id);

		for (auto &sugar : branch)
		{
			pdbx_branch_scheme.emplace({ { "asym_id", asym_id },
				{ "entity_id", entity_id },
				{ "num", sugar.num() },
				{ "mon_id", sugar.get_compound_id() },

				{ "pdb_asym_id", asym_id },
				{ "pdb_seq_num", sugar.num() },
				{ "pdb_mon_id", sugar.get_compound_id() },

				// TODO: need fix, collect from nag_atoms?
				{ "auth_asym_id", asym_id },
				{ "auth_mon_id", sugar.get_compound_id() },
				{ "auth_seq_num", sugar.get_pdb_seq_num() },

				{ "hetero", "n" } });
		}
	}
}

void structure::remove_branch(branch &branch)
{
	using namespace literals;

	for (auto &sugar : branch)
	{
		auto atoms = sugar.atoms();
		for (auto atom : atoms)
			remove_atom(atom);
	}

	m_db["pdbx_branch_scheme"].erase("asym_id"_key == branch.get_asym_id());
	m_db["struct_asym"].erase("id"_key == branch.get_asym_id());
	m_db["struct_conn"].erase("ptnr1_label_asym_id"_key == branch.get_asym_id() or "ptnr2_label_asym_id"_key == branch.get_asym_id());

	std::erase(m_branches, branch);
}

std::string structure::create_non_poly_entity(const std::string &comp_id)
{
	return insert_compound(comp_id, true);
}

std::string structure::create_non_poly(const std::string &entity_id, const std::vector<atom> &atoms)
{
	using namespace cif::literals;

	auto &struct_asym = m_db["struct_asym"];
	std::string asym_id = struct_asym.get_unique_id();

	struct_asym.emplace({ { "id", asym_id },
		{ "pdbx_blank_PDB_chainid_flag", "N" },
		{ "pdbx_modified", "N" },
		{ "entity_id", entity_id },
		{ "details", "?" } });

	auto comp_id = m_db["pdbx_entity_nonpoly"].find1<std::string>("entity_id"_key == entity_id, "comp_id");

	auto &atom_site = m_db["atom_site"];

	auto &res = m_non_polymers.emplace_back(*this, comp_id, asym_id, 0, asym_id, "1", "");

	for (auto &atom : atoms)
	{
		auto atom_id = atom_site.get_unique_id("");

		auto row = atom_site.emplace({ { "group_PDB", atom.get_property("group_PDB") },
			{ "id", atom_id },
			{ "type_symbol", atom.get_property("type_symbol") },
			{ "label_atom_id", atom.get_property("label_atom_id") },
			{ "label_alt_id", atom.get_property("label_alt_id") },
			{ "label_comp_id", comp_id },
			{ "label_asym_id", asym_id },
			{ "label_entity_id", entity_id },
			{ "label_seq_id", "." },
			{ "pdbx_PDB_ins_code", "" },
			{ "Cartn_x", atom.get_property("Cartn_x") },
			{ "Cartn_y", atom.get_property("Cartn_y") },
			{ "Cartn_z", atom.get_property("Cartn_z") },
			{ "occupancy", atom.get_property("occupancy") },
			{ "B_iso_or_equiv", atom.get_property("B_iso_or_equiv") },
			{ "pdbx_formal_charge", atom.get_property("pdbx_formal_charge") },
			{ "auth_seq_id", 1 },
			{ "auth_comp_id", comp_id },
			{ "auth_asym_id", asym_id },
			{ "auth_atom_id", atom.get_property("label_atom_id") },
			{ "pdbx_PDB_model_num", 1 } });

		auto &newAtom = emplace_atom(std::make_shared<atom::atom_impl>(m_db, atom_id));
		res.add_atom(newAtom);
	}

	auto &pdbx_nonpoly_scheme = m_db["pdbx_nonpoly_scheme"];
	std::size_t ndb_nr = pdbx_nonpoly_scheme.find("asym_id"_key == asym_id and "entity_id"_key == entity_id).size() + 1;
	pdbx_nonpoly_scheme.emplace({
		{ "asym_id", asym_id },
		{ "entity_id", entity_id },
		{ "mon_id", comp_id },
		{ "ndb_seq_num", ndb_nr },
		{ "pdb_seq_num", res.get_pdb_seq_num() },
		{ "auth_seq_num", res.get_pdb_seq_num() },
		{ "pdb_mon_id", comp_id },
		{ "auth_mon_id", comp_id },
		{ "pdb_strand_id", asym_id },
		{ "pdb_ins_code", "." },
	});

	return asym_id;
}

std::string structure::create_non_poly(const std::string &entity_id, std::vector<row_initializer> atoms)
{
	using namespace literals;

	auto &struct_asym = m_db["struct_asym"];
	std::string asym_id = struct_asym.get_unique_id();

	struct_asym.emplace({ { "id", asym_id },
		{ "pdbx_blank_PDB_chainid_flag", "N" },
		{ "pdbx_modified", "N" },
		{ "entity_id", entity_id },
		{ "details", "?" } });

	auto comp_id = m_db["pdbx_entity_nonpoly"].find1<std::string>("entity_id"_key == entity_id, "comp_id");

	auto &atom_site = m_db["atom_site"];

	auto &res = m_non_polymers.emplace_back(*this, comp_id, asym_id, 0, asym_id, "1", "");

	for (auto &atom : atoms)
	{
		auto atom_id = atom_site.get_unique_id("");

		atom.set_value("id", atom_id);
		atom.set_value("label_asym_id", asym_id);
		atom.set_value("auth_asym_id", asym_id);
		atom.set_value("label_entity_id", entity_id);

		atom.set_value_if_empty({ "group_PDB", "HETATM" });
		atom.set_value_if_empty({ "label_comp_id", comp_id });
		atom.set_value_if_empty({ "label_seq_id", "." });
		atom.set_value_if_empty({ "auth_comp_id", comp_id });
		atom.set_value_if_empty({ "auth_seq_id", 1 });
		atom.set_value_if_empty({ "pdbx_PDB_model_num", 1 });
		atom.set_value_if_empty({ "label_alt_id", "" });
		atom.set_value_if_empty({ "occupancy", 1.0, 2 });

		auto row = atom_site.emplace(atom.begin(), atom.end());

		auto &newAtom = emplace_atom(std::make_shared<atom::atom_impl>(m_db, atom_id));
		res.add_atom(newAtom);
	}

	auto &pdbx_nonpoly_scheme = m_db["pdbx_nonpoly_scheme"];
	std::size_t ndb_nr = pdbx_nonpoly_scheme.find("asym_id"_key == asym_id and "entity_id"_key == entity_id).size() + 1;
	pdbx_nonpoly_scheme.emplace({
		{ "asym_id", asym_id },
		{ "entity_id", entity_id },
		{ "mon_id", comp_id },
		{ "ndb_seq_num", ndb_nr },
		{ "pdb_seq_num", res.get_pdb_seq_num() },
		{ "auth_seq_num", res.get_pdb_seq_num() },
		{ "pdb_mon_id", comp_id },
		{ "auth_mon_id", comp_id },
		{ "pdb_strand_id", asym_id },
		{ "pdb_ins_code", "." },
	});

	return asym_id;
}

std::string structure::create_non_poly(const std::string &compound_id, bool skip_hydrogen)
{
	auto compound = cif::compound_factory::instance().create(compound_id);
	if (compound == nullptr)
		throw std::runtime_error(std::format("{} is not a known compound", compound_id));

	std::vector<cif::row_initializer> atoms;
	for (auto a : compound->atoms())
	{
		// We skip H-atoms, as fitting without H-atoms works better and we avoid conflicts in protonation states between CCD and MONLIB
		if (skip_hydrogen and cif::atom_type_traits(a.type_symbol).symbol() == "H")
			continue;

		auto ax = a.get_location().get_x();
		auto ay = a.get_location().get_y();
		auto az = a.get_location().get_z();

		atoms.emplace_back(cif::row_initializer{
			{ "type_symbol", cif::atom_type_traits(a.type_symbol).symbol() },
			{ "label_atom_id", a.id },
			{ "auth_atom_id", a.id },
			{ "Cartn_x", ax },
			{ "Cartn_y", ay },
			{ "Cartn_z", az },
			{ "B_iso_or_equiv", 30.00 } });
	}

	return create_non_poly(create_non_poly_entity(compound_id), atoms);
}

void structure::create_water(row_initializer atom)
{
	using namespace literals;

	auto entity_id = insert_compound("HOH", true);

	auto &struct_asym = m_db["struct_asym"];
	std::string asym_id;
	try
	{
		asym_id = struct_asym.find1<std::string>("entity_id"_key == entity_id, "id");
	}
	catch (const std::exception &)
	{
		asym_id = struct_asym.get_unique_id();

		struct_asym.emplace({ { "id", asym_id },
			{ "pdbx_blank_PDB_chainid_flag", "N" },
			{ "pdbx_modified", "N" },
			{ "entity_id", entity_id },
			{ "details", "?" } });
	}

	auto &atom_site = m_db["atom_site"];
	auto auth_seq_id = atom_site.find_max<int>("auth_seq_id", "label_entity_id"_key == entity_id) + 1;
	if (auth_seq_id < 0)
		auth_seq_id = 1;

	auto atom_id = atom_site.get_unique_id("");

	atom.set_value("id", atom_id);
	atom.set_value("label_asym_id", asym_id);
	atom.set_value("auth_asym_id", asym_id);
	atom.set_value("label_entity_id", entity_id);
	atom.set_value("auth_seq_id", std::to_string(auth_seq_id));

	atom.set_value_if_empty({ "group_PDB", "HETATM" });
	atom.set_value_if_empty({ "label_comp_id", "HOH" });
	atom.set_value_if_empty({ "label_seq_id", "." });
	atom.set_value_if_empty({ "auth_comp_id", "HOH" });
	atom.set_value_if_empty({ "pdbx_PDB_model_num", 1 });
	atom.set_value_if_empty({ "label_alt_id", "" });
	atom.set_value_if_empty({ "occupancy", 1.0, 2 });

	auto row = atom_site.emplace(atom.begin(), atom.end());

	emplace_atom(std::make_shared<atom::atom_impl>(m_db, atom_id));

	auto &pdbx_nonpoly_scheme = m_db["pdbx_nonpoly_scheme"];
	int ndb_nr = pdbx_nonpoly_scheme.find_max<int>("ndb_seq_num") + 1;
	pdbx_nonpoly_scheme.emplace({
		{ "asym_id", asym_id },
		{ "entity_id", entity_id },
		{ "mon_id", "HOH" },
		{ "ndb_seq_num", ndb_nr },
		{ "pdb_seq_num", auth_seq_id },
		{ "auth_seq_num", auth_seq_id },
		{ "pdb_mon_id", "HOH" },
		{ "auth_mon_id", "HOH" },
		{ "pdb_strand_id", asym_id },
		{ "pdb_ins_code", "." },
	});
}

std::string structure::create_link(atom a1, atom a2, const std::string &link_type, const std::string &role)
{
	using namespace literals;

	auto &struct_conn = m_db["struct_conn"];
	auto &struct_conn_type = m_db["struct_conn_type"];

	// This will validate link_type :-)
	if (not struct_conn_type.contains("id"_key == link_type))
		struct_conn_type.emplace({ { "id", link_type } });

	std::string link_id = struct_conn.get_unique_id(link_type + '_');

	item label_seq_id_1("ptnr1_label_seq_id");
	if (int nr = a1.get_label_seq_id(); nr != 0)
		label_seq_id_1.value(std::to_string(nr));

	item label_seq_id_2("ptnr2_label_seq_id");
	if (int nr = a2.get_label_seq_id(); nr != 0)
		label_seq_id_2.value(std::to_string(nr));

	struct_conn.emplace(
		{ //
			{ "id", link_id },
			{ "conn_type_id", link_type },
			{ "pdbx_leaving_atom_flag", "one" },

			{ "ptnr1_label_asym_id", a1.get_label_asym_id() },
			{ "ptnr1_label_comp_id", a1.get_label_comp_id() },
			label_seq_id_1,
			{ "ptnr1_label_atom_id", a1.get_label_atom_id() },
			{ "pdbx_ptnr1_label_alt_id", a1.get_label_alt_id() },
			{ "pdbx_ptnr1_PDB_ins_code", a1.get_pdb_ins_code() },
			{ "ptnr1_auth_asym_id", a1.get_auth_asym_id() },
			{ "ptnr1_auth_comp_id", a1.get_auth_comp_id() },
			{ "ptnr1_auth_seq_id", a1.get_auth_seq_id() },
			{ "ptnr1_symmetry", a1.symmetry() },

			{ "ptnr2_label_asym_id", a2.get_label_asym_id() },
			{ "ptnr2_label_comp_id", a2.get_label_comp_id() },
			label_seq_id_2,
			{ "ptnr2_label_atom_id", a2.get_label_atom_id() },
			{ "pdbx_ptnr2_label_alt_id", a2.get_label_alt_id() },
			{ "pdbx_ptnr2_PDB_ins_code", a2.get_pdb_ins_code() },
			{ "ptnr2_auth_asym_id", a2.get_auth_asym_id() },
			{ "ptnr2_auth_comp_id", a2.get_auth_comp_id() },
			{ "ptnr2_auth_seq_id", a2.get_auth_seq_id() },
			{ "ptnr2_symmetry", a2.symmetry() },

			{ "pdbx_dist_value", distance(a1.get_location(), a2.get_location()), 3 },
			{ "pdbx_role", role } });

	return link_id;
}

branch &structure::create_branch()
{
	auto &entity = m_db["entity"];
	auto &struct_asym = m_db["struct_asym"];

	auto entity_id = entity.get_unique_id("");
	auto asym_id = struct_asym.get_unique_id();

	entity.emplace({ { "id", entity_id },
		{ "type", "branched" } });

	struct_asym.emplace({ { "id", asym_id },
		{ "pdbx_blank_PDB_chainid_flag", "N" },
		{ "pdbx_modified", "N" },
		{ "entity_id", entity_id },
		{ "details", "?" } });

	return m_branches.emplace_back(*this, asym_id, entity_id);
}

// branch &structure::create_branch(std::vector<row_initializer> atoms)
// {
// // 	// sanity check
// // 	for (auto &nag_atom : atoms)
// // 	{
// // 		for (const auto &[name, value] : nag_atom)
// // 		{
// // 			if (name == "label_comp_id" and value != "NAG")
// // 				throw std::logic_error("The first sugar in a branch should be a NAG");
// // 		}
// // 	}

// // 	using namespace literals;

// // 	auto &branch = create_branch();
// // 	auto asym_id = branch.get_asym_id();
// // 	auto entity_id = branch.get_entity_id();

// // 	auto &sugar = branch.emplace_back(branch, "NAG", asym_id, 1);

// // 	auto &atom_site = m_db["atom_site"];

// // 	for (auto &atom : atoms)
// // 	{
// // 		auto atom_id = atom_site.get_unique_id("");

// // 		atom.set_value("id", atom_id);
// // 		atom.set_value("label_asym_id", asym_id);
// // 		atom.set_value("auth_asym_id", asym_id);
// // 		atom.set_value("label_entity_id", entity_id);
// // 		atom.set_value({ "auth_seq_id", 1 });

// // 		atom.set_value_if_empty({"group_PDB", "HETATM"});
// // 		atom.set_value_if_empty({"label_comp_id", "NAG"});
// // 		atom.set_value_if_empty({"label_seq_id", "."});
// // 		atom.set_value_if_empty({"auth_comp_id", "NAG"});
// // 		atom.set_value_if_empty({"pdbx_PDB_model_num", 1});
// // 		atom.set_value_if_empty({"label_alt_id", ""});

// // 		auto row = atom_site.emplace(atom.begin(), atom.end());

// // 		auto &newAtom = emplace_atom(std::make_shared<atom::atom_impl>(m_db, atom_id));
// // 		sugar.add_atom(newAtom);
// // 	}

// // 	// // now we can create the entity and get the real ID
// // 	// auto entity_id = create_entity_for_branch(branch);
// // 	// assert(not entity_id.empty());

// // 	for (auto &a : sugar.atoms())
// // 		a.set_property("label_entity_id", entity_id);

// // 	m_db["pdbx_branch_scheme"].emplace({
// // 		{"asym_id", asym_id},
// // 		{"entity_id", entity_id},
// // 		{"num", 1},
// // 		{"mon_id", "NAG"},

// // 		{"pdb_asym_id", asym_id},
// // 		{"pdb_seq_num", 1},
// // 		{"pdb_mon_id", "NAG"},

// // 		// TODO: need fix, collect from nag_atoms?
// // 		{"auth_asym_id", asym_id},
// // 		{"auth_mon_id", "NAG"},
// // 		{"auth_seq_num", 1},

// // 		{"hetero", "n"}
// // 	});

// // 	return branch;
// }

// branch &structure::extend_branch(const std::string &asym_id, std::vector<row_initializer> atom_info,
// 	int link_sugar, const std::string &link_atom)
// {
// // 	// sanity check
// // 	std::string compoundID;

// // 	for (auto &atom : atom_info)
// // 	{
// // 		for (const auto &[name, value] : atom)
// // 		{
// // 			if (name != "label_comp_id")
// // 				continue;

// // 			if (compoundID.empty())
// // 				compoundID = value;
// // 			else if (value != compoundID)
// // 				throw std::logic_error("All atoms should be of the same type");
// // 		}
// // 	}

// // 	using namespace literals;

// // 	// auto &branch = m_branches.emplace_back(*this, asym_id);
// // 	auto tmp_entity_id = m_db["entity"].get_unique_id("");

// // 	auto &atom_site = m_db["atom_site"];

// // 	auto bi = std::find_if(m_branches.begin(), m_branches.end(), [asym_id](branch &b)
// // 		{ return b.get_asym_id() == asym_id; });
// // 	if (bi == m_branches.end())
// // 		throw std::logic_error("Create a branch first!");

// // 	branch &branch = *bi;

// // 	int sugarNum = static_cast<int>(branch.size() + 1);

// // 	auto &sugar = branch.emplace_back(branch, compoundID, asym_id, sugarNum);

// // 	for (auto &atom : atom_info)
// // 	{
// // 		auto atom_id = atom_site.get_unique_id("");

// // 		atom.set_value("id", atom_id);
// // 		atom.set_value("label_asym_id", asym_id);
// // 		atom.set_value("auth_asym_id", asym_id);
// // 		atom.set_value("label_entity_id", tmp_entity_id);
// // 		atom.set_value({"auth_seq_id", sugarNum });

// // 		atom.set_value_if_empty({"group_PDB", "HETATM"});
// // 		atom.set_value_if_empty({"label_comp_id", compoundID});
// // 		atom.set_value_if_empty({"auth_comp_id", compoundID});
// // 		atom.set_value_if_empty({"pdbx_PDB_model_num", 1});
// // 		atom.set_value_if_empty({"label_alt_id", ""});

// // 		auto row = atom_site.emplace(atom.begin(), atom.end());

// // 		auto &newAtom = emplace_atom(std::make_shared<atom::atom_impl>(m_db, atom_id));
// // 		sugar.add_atom(newAtom);
// // 	}

// // 	sugar.set_link(branch.at(link_sugar - 1).get_atom_by_atom_id(link_atom));

// // 	auto entity_id = create_entity_for_branch(branch);

// // 	// Update the entity id of the asym
// // 	auto &struct_asym = m_db["struct_asym"];
// // 	auto r = struct_asym.find1("id"_key == asym_id);
// // 	r["entity_id"] = entity_id;

// // 	for (auto &s2 : branch)
// // 	{
// // 		for (auto atom : s2.atoms())
// // 			atom.set_property("label_entity_id", entity_id);
// // 	}

// // 	auto &pdbx_branch_scheme = m_db["pdbx_branch_scheme"];
// // 	pdbx_branch_scheme.erase("asym_id"_key == asym_id);

// // 	for (auto &s2 : branch)
// // 	{
// // 		pdbx_branch_scheme.emplace({
// // 			{"asym_id", asym_id},
// // 			{"entity_id", entity_id},
// // 			{"num", s2.num()},
// // 			{"mon_id", s2.get_compound_id()},

// // 			{"pdb_asym_id", asym_id},
// // 			{"pdb_seq_num", s2.num()},
// // 			{"pdb_mon_id", s2.get_compound_id()},

// // 			// TODO: need fix, collect from nag_atoms?
// // 			{"auth_asym_id", asym_id},
// // 			{"auth_mon_id", s2.get_compound_id()},
// // 			{"auth_seq_num", s2.get_auth_seq_id()},

// // 			{"hetero", "n"}
// // 		});
// // 	}

// // 	return branch;
// }

std::string structure::create_entity_for_branch(branch &branch)
{
	using namespace literals;

	std::string entityName = branch.name();

	auto &entity = m_db["entity"];

	auto entityID = entity.find_first<std::string>("type"_key == "branched" and "pdbx_description"_key == entityName, "id");

	if (entityID.empty())
	{
		entityID = entity.get_unique_id("");

		if (VERBOSE)
			std::cout << "Creating new entity " << entityID << " for branched sugar " << entityName << '\n';

		entity.emplace({ { "id", entityID },
			{ "type", "branched" },
			{ "src_method", "man" },
			{ "pdbx_description", entityName },
			{ "formula_weight", branch.weight() } });

		auto &pdbx_entity_branch_list = m_db["pdbx_entity_branch_list"];
		for (auto &sugar : branch)
		{
			pdbx_entity_branch_list.emplace({ { "entity_id", entityID },
				{ "comp_id", sugar.get_compound_id() },
				{ "num", sugar.num() },
				{ "hetero", "n" } });
		}

		auto &pdbx_entity_branch_link = m_db["pdbx_entity_branch_link"];
		for (auto &s1 : branch)
		{
			auto l2 = s1.get_link();

			if (not l2 or l2.get_auth_seq_id().empty())
				continue;

			auto &s2 = branch.at(stoi(l2.get_auth_seq_id()) - 1);
			auto l1 = s2.get_atom_by_atom_id("C1");

			pdbx_entity_branch_link.emplace({ { "link_id", pdbx_entity_branch_link.get_unique_id("") },
				{ "entity_id", entityID },
				{ "entity_branch_list_num_1", s1.get_pdb_seq_num() },
				{ "comp_id_1", s1.get_compound_id() },
				{ "atom_id_1", l1.get_label_atom_id() },
				{ "leaving_atom_id_1", "O1" },
				{ "entity_branch_list_num_2", s2.get_pdb_seq_num() },
				{ "comp_id_2", s2.get_compound_id() },
				{ "atom_id_2", l2.get_label_atom_id() },
				{ "leaving_atom_id_2", "H" + l2.get_label_atom_id() },
				{ "value_order", "sing" } });
		}
	}

	return entityID;
}

void structure::cleanup_empty_categories()
{

	using namespace literals;

	auto &atomSite = m_db["atom_site"];
	auto &pdbxPolySeqScheme = m_db["pdbx_poly_seq_scheme"];
	auto &entityPolySeq = m_db["entity_poly_seq"];

	// Remove chem_comp's for which there are no atoms at all
	auto &chem_comp = m_db["chem_comp"];
	std::vector<row_handle> obsoleteChemComps;

	for (auto chemComp : chem_comp)
	{
		auto compID = chemComp["id"].as<std::string>();
		if (atomSite.contains("label_comp_id"_key == compID or "auth_comp_id"_key == compID) or
			pdbxPolySeqScheme.contains("mon_id"_key == compID or "auth_mon_id"_key == compID or "pdb_mon_id"_key == compID) or
			entityPolySeq.contains("mon_id"_key == compID))
		{
			continue;
		}

		obsoleteChemComps.push_back(chemComp);
	}

	for (auto chemComp : obsoleteChemComps)
		chem_comp.erase(chemComp);

	// similarly, remove entities not referenced by any atom

	auto &entities = m_db["entity"];
	std::vector<row_handle> obsoleteEntities;

	for (auto entity : entities)
	{
		auto entityID = entity["id"].as<std::string>();
		if (atomSite.contains("label_entity_id"_key == entityID))
			continue;

		obsoleteEntities.push_back(entity);
	}

	auto validator = m_db.get_validator();

	for (auto entity : obsoleteEntities)
	{
		auto entityID = entity["id"].as<std::string>();
		if (validator)
		{
			for (auto linked : validator->get_links_for_parent("entity"))
			{
				if (auto cat = m_db.get(linked->m_child_category))
					cat->erase(cif::key(linked->m_child_keys.front()) == entityID);
			}
		}

		entities.erase(entity);
	}

	// the rest?

	for (const char *cat : { "pdbx_entity_nonpoly" })
	{
		auto &category = m_db[cat];

		std::vector<row_handle> empty;
		for (auto row : category)
		{
			if (not category.has_children(row) and not category.has_parents(row))
				empty.push_back(row);
		}

		for (auto row : empty)
			category.erase(row);
	}

	// count molecules
	for (auto entity : entities)
	{
		std::string type, id;
		cif::tie(type, id) = entity.get("type", "id");

		std::optional<std::size_t> count;
		if (type == "polymer")
			count = m_db["struct_asym"].find("entity_id"_key == id).size();
		else if (type == "non-polymer" or type == "water")
			count = m_db["pdbx_nonpoly_scheme"].find("entity_id"_key == id).size();
		else if (type == "branched")
		{
			// is this correct?
			std::set<std::string> asym_ids;
			for (const auto &asym_id : m_db["pdbx_branch_scheme"].find<std::string>("entity_id"_key == id, "asym_id"))
				asym_ids.insert(asym_id);
			count = asym_ids.size();
		}

		entity["pdbx_number_of_molecules"] = count.value_or(0);
	}
}

void structure::translate(point t)
{
	for (auto &a : m_atoms)
		a.translate(t);
}

void structure::rotate(quaternion q)
{
	for (auto &a : m_atoms)
		a.rotate(q);
}

void structure::translate_and_rotate(point t, quaternion q)
{
	for (auto &a : m_atoms)
		a.translate_and_rotate(t, q);
}

void structure::translate_rotate_and_translate(point t1, quaternion q, point t2)
{
	for (auto &a : m_atoms)
		a.translate_rotate_and_translate(t1, q, t2);
}

void structure::validate_atoms() const
{
	// validate order
	assert(m_atoms.size() == m_atom_index.size());
	for (std::size_t i = 0; i + i < m_atoms.size(); ++i)
		assert(m_atoms[m_atom_index[i]].id().compare(m_atoms[m_atom_index[i + 1]].id()) < 0);

	std::vector<atom> atoms = m_atoms;

	for (auto &poly : m_polymers)
	{
		for (auto &monomer : poly)
		{
			for (auto &atom : monomer.atoms())
				std::erase(atoms, atom);
		}
	}

	for (auto &branch : m_branches)
	{
		for (auto &sugar : branch)
		{
			for (auto &atom : sugar.atoms())
				std::erase(atoms, atom);
		}
	}

	for (auto &res : m_non_polymers)
	{
		for (auto &atom : res.atoms())
			std::erase(atoms, atom);
	}

	assert(atoms.empty());
}

static int compare_numbers(std::string_view a, std::string_view b)
{
	int result = 0;
	double da, db;

	using namespace cif;
	using namespace std;

	std::from_chars_result ra, rb;

	ra = from_chars(a.data(), a.data() + a.length(), da);
	rb = from_chars(b.data(), b.data() + b.length(), db);

	if (ra.ec == std::errc{} and rb.ec == std::errc{})
	{
		auto d = da - db;
		if (std::abs(d) > std::numeric_limits<double>::epsilon())
		{
			if (d > 0)
				result = 1;
			else if (d < 0)
				result = -1;
		}
	}
	else if (ra.ec != std::errc{})
		result = 1;
	else
		result = -1;

	return result;
}

int compare_cif_id(const std::string &a, const std::string &b)
{
	int d = static_cast<int>(a.length() - b.length());
	if (d == 0)
		d = a.compare(b);
	return d;
}

void structure::reorder_atoms()
{
	auto &atom_site = m_db["atom_site"];

	atom_site.sort([](row_handle a, row_handle b)
		{
			int d;

			// First by model number
			d = a.get<int>("pdbx_PDB_model_num") - b.get<int>("pdbx_PDB_model_num");
			if (d == 0)
				d = compare_cif_id(a.get<std::string>("label_asym_id"), b.get<std::string>("label_asym_id"));
			if (d == 0)
			{
				auto na = a.get<std::optional<int>>("label_seq_id");
				auto nb = b.get<std::optional<int>>("label_seq_id");

				if (na.has_value() and nb.has_value())
					d = *na - *nb;
				else if (na.has_value())
					d = 1;
				else if (nb.has_value())
					d = -1;
			}

			if (d == 0)
			{
				auto na = a.get<std::optional<int>>("auth_seq_id");
				auto nb = b.get<std::optional<int>>("auth_seq_id");

				if (na.has_value() and nb.has_value())
					d = *na - *nb;
				else if (na.has_value())
					d = 1;
				else if (nb.has_value())
					d = -1;
			}

			if (d == 0)
				d = compare_numbers(a.get<std::string>("id"), b.get<std::string>("id"));

			return d;
			//
		});

	// atom_site.set_validator(nullptr, m_db);

	// for (int nr = 1; auto r : atom_site)
	// 	r["id"] = nr++;

	// atom_site.set_validator(m_db.get_validator(), m_db);
}

} // namespace cif::mm
