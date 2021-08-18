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

#pragma once

/// \file This file contains the definition for the class Compound, encapsulating
/// the information found for compounds in the CCD.

#include <map>
#include <set>
#include <tuple>
#include <vector>

#include "cif++/AtomType.hpp"
#include "cif++/Cif++.hpp"

namespace mmcif
{

// --------------------------------------------------------------------

class Compound;
struct CompoundAtom;
class CompoundFactoryImpl;

/// \brief The bond type as defined in the CCD, possible values taken from the mmcif_pdbx_v50 file
enum class BondType
{
	sing, // 'single bond'
	doub, // 'double bond'
	trip, // 'triple bond'
	quad, // 'quadruple bond'
	arom, // 'aromatic bond'
	poly, // 'polymeric bond'
	delo, // 'delocalized double bond'
	pi,   // 'pi bond'
};

std::string to_string(BondType bondType);
BondType from_string(const std::string& bondType);

/// --------------------------------------------------------------------
/// \brief struct containing information about an atom in a chemical compound.
/// This is a subset of the available information. Contact the author if you need more fields.

struct CompoundAtom
{
	std::string id;
	AtomType typeSymbol;
	int charge = 0;
	bool aromatic = false;
	bool leavingAtom = false;
	bool stereoConfig = false;
	float x, y, z;
};

/// --------------------------------------------------------------------
/// \brief struct containing information about the bonds

struct CompoundBond
{
	std::string atomID[2];
	BondType type;
	bool aromatic = false, stereoConfig = false;
};

/// --------------------------------------------------------------------
/// \brief a class that contains information about a chemical compound.
/// This information is derived from the CDD by default.
///
/// To create compounds, you use the factory method. You can add your own
/// compound definitions by calling the addExtraComponents function and
/// pass it a valid CCD formatted file.

class Compound
{
  public:

	// accessors

	std::string id() const { return mID; }
	std::string name() const { return mName; }
	std::string type() const { return mType; }
	std::string formula() const { return mFormula; }
	float formulaWeight() const { return mFormulaWeight; }
	int formalCharge() const { return mFormalCharge; }

	const std::vector<CompoundAtom> &atoms() const { return mAtoms; }
	const std::vector<CompoundBond> &bonds() const { return mBonds; }

	CompoundAtom getAtomByID(const std::string &atomID) const;

	bool atomsBonded(const std::string &atomId_1, const std::string &atomId_2) const;
	// float atomBondValue(const std::string &atomId_1, const std::string &atomId_2) const;
	// float bondAngle(const std::string &atomId_1, const std::string &atomId_2, const std::string &atomId_3) const;
	// float chiralVolume(const std::string &centreID) const;

	bool isWater() const
	{
		return mID == "HOH" or mID == "H2O" or mID == "WAT";
	}

  private:

	friend class CompoundFactoryImpl;
	friend class CCDCompoundFactoryImpl;
	friend class CCP4CompoundFactoryImpl;

	Compound(cif::Datablock &db);
	Compound(cif::Datablock &db, const std::string &id, const std::string &name, const std::string &type);

	std::string mID;
	std::string mName;
	std::string mType;
	std::string mFormula;
	float mFormulaWeight = 0;
	int mFormalCharge = 0;
	std::vector<CompoundAtom> mAtoms;
	std::vector<CompoundBond> mBonds;
};

// --------------------------------------------------------------------
// Factory class for Compound and Link objects

CIFPP_EXPORT extern const std::map<std::string, char> kAAMap, kBaseMap;

class CompoundFactory
{
  public:

	/// \brief Initialise a singleton instance.
	///
	/// If you have a multithreaded application and want to have different
	/// compounds in each thread (e.g. a web service processing user requests
	/// with different sets of compounds) you can set the \a useThreadLocalInstanceOnly
	/// flag to true.

	static void init(bool useThreadLocalInstanceOnly);
	static CompoundFactory &instance();
	static void clear();

	void setDefaultDictionary(const std::string &inDictFile);
	void pushDictionary(const std::string &inDictFile);
	void popDictionary();

	bool isKnownPeptide(const std::string &res_name) const;
	bool isKnownBase(const std::string &res_name) const;

	/// \brief Create the Compound object for \a id
	///
	/// This will create the Compound instance for \a id if it doesn't exist already.
	/// The result is owned by this factory and should not be deleted by the user.
	/// \param id	The Compound ID, a three letter code usually
	/// \result		The compound, or nullptr if it could not be created (missing info)
	const Compound *create(std::string id);

	~CompoundFactory();

  private:
	CompoundFactory();

	CompoundFactory(const CompoundFactory &) = delete;
	CompoundFactory &operator=(const CompoundFactory &) = delete;

	static std::unique_ptr<CompoundFactory> sInstance;
	static thread_local std::unique_ptr<CompoundFactory> tlInstance;
	static bool sUseThreadLocalInstance;

	std::shared_ptr<CompoundFactoryImpl> mImpl;
};

} // namespace mmcif
