/* 
   Created by: Maarten L. Hekkelman
   Date: woensdag 06 juni, 2018
*/

// Calculate DSSP-like secondary structure information

#include "cif++/Config.h"

#include <numeric>

#include <boost/algorithm/string.hpp>

#include "cif++/Structure.h"
#include "cif++/Secondary.h"

using namespace std;
namespace ba = boost::algorithm;

// --------------------------------------------------------------------

namespace mmcif
{

struct Res;

enum ResidueType
{
	kUnknownResidue,
	
	//
	kAlanine,				// A	ala
	kArginine,				// R	arg
	kAsparagine,			// N	asn
	kAsparticAcid,			// D	asp
	kCysteine,				// C	cys
	kGlutamicAcid,			// E	glu
	kGlutamine,				// Q	gln
	kGlycine,				// G	gly
	kHistidine,				// H	his
	kIsoleucine,			// I	ile
	kLeucine,				// L	leu
	kLysine,				// K	lys
	kMethionine,			// M	met
	kPhenylalanine,			// F	phe
	kProline,				// P	pro
	kSerine,				// S	ser
	kThreonine,				// T	thr
	kTryptophan,			// W	trp
	kTyrosine,				// Y	tyr
	kValine,				// V	val
	
	kResidueTypeCount
};

struct ResidueInfo
{
	ResidueType		type;
	char			code;
	char			name[4];
};

const ResidueInfo kResidueInfo[] = {
	{ kUnknownResidue,	'X', "UNK" },
	{ kAlanine,			'A', "ALA" },
	{ kArginine,		'R', "ARG" },
	{ kAsparagine,		'N', "ASN" },
	{ kAsparticAcid,	'D', "ASP" },
	{ kCysteine,		'C', "CYS" },
	{ kGlutamicAcid,	'E', "GLU" },
	{ kGlutamine,		'Q', "GLN" },
	{ kGlycine,			'G', "GLY" },
	{ kHistidine,		'H', "HIS" },
	{ kIsoleucine,		'I', "ILE" },
	{ kLeucine,			'L', "LEU" },
	{ kLysine,			'K', "LYS" },
	{ kMethionine,		'M', "MET" },
	{ kPhenylalanine,	'F', "PHE" },
	{ kProline,			'P', "PRO" },
	{ kSerine,			'S', "SER" },
	{ kThreonine,		'T', "THR" },
	{ kTryptophan,		'W', "TRP" },
	{ kTyrosine,		'Y', "TYR" },
	{ kValine,			'V', "VAL" }
};

ResidueType MapResidue(string inName)
{
	ba::trim(inName);

	ResidueType result = kUnknownResidue;
	
	for (auto& ri: kResidueInfo)
	{
		if (inName == ri.name)
		{
			result = ri.type;
			break;
		}
	}
	
	return result;
}

struct HBond
{
	Res*	residue;
	double	energy;
};

enum BridgeType
{
	btNoBridge, btParallel, btAntiParallel
};

struct Bridge
{
	BridgeType		type;
	uint32			sheet, ladder;
	set<Bridge*>	link;
	deque<uint32>	i, j;
	string			chainI, chainJ;
	
	bool			operator<(const Bridge& b) const		{ return chainI < b.chainI or (chainI == b.chainI and i.front() < b.i.front()); }
};

struct BridgeParner
{
	Res*		residue;
	uint32		ladder;
	bool		parallel;
};

enum HelixFlag
{
	helixNone, helixStart, helixEnd, helixStartAndEnd, helixMiddle
};

// --------------------------------------------------------------------

const double
	kSSBridgeDistance = 3.0,
	kMinimalDistance = 0.5,
	kMinimalCADistance = 9.0,
	kMinHBondEnergy = -9.9,
	kMaxHBondEnergy = -0.5,
	kCouplingConstant = -27.888,	//	= -332 * 0.42 * 0.2
	kMaxPeptideBondLength = 2.5;

const uint32 kHistogramSize = 30;

struct Res
{
	Res(Monomer&& m)
		: mM(move(m))
		, mType(MapResidue(m.compoundID()))
	{
		for (auto& a: mM.atoms())
		{
			if (a.labelAtomId() == "CA")		mCAlpha = a.location();
			else if (a.labelAtomId() == "C")	mC = a.location();
			else if (a.labelAtomId() == "N")	mN = a.location();
			else if (a.labelAtomId() == "O")	mO = a.location();
//			else if (a.labelAtomId() == "H")	mH = a.location();
		}
	}

	void assignHydrogen()
	{
		// assign the Hydrogen
		mH = mN;
		
		if (mType != kProline and mPrev != nullptr)
		{
			auto pc = mPrev->mC;
			auto po = mPrev->mO;
			
			double CODistance = Distance(pc, po);
			
			mH.mX += (pc.mX - po.mX) / CODistance; 
			mH.mY += (pc.mY - po.mY) / CODistance; 
			mH.mZ += (pc.mZ - po.mZ) / CODistance; 
		}
	}

	void SetSecondaryStructure(SecondaryStructureType inSS)	{ mSecondaryStructure = inSS; }
	SecondaryStructureType GetSecondaryStructure() const	{ return mSecondaryStructure; }
	
	void SetBetaPartner(uint32 n, Res& inResidue, uint32 inLadder, bool inParallel)
	{
		assert(n == 0 or n == 1);
		
		mBetaPartner[n].residue = &inResidue;
		mBetaPartner[n].ladder = inLadder;
		mBetaPartner[n].parallel = inParallel;
	}
	
	BridgeParner GetBetaPartner(uint32 n) const
	{
		assert(n == 0 or n == 1);
		return mBetaPartner[n];
	}
						
	void SetSheet(uint32 inSheet)						{ mSheet = inSheet; }
	uint32 GetSheet() const								{ return mSheet; }
	
	bool IsBend() const									{ return mBend; }
	void SetBend(bool inBend)							{ mBend = inBend; }
	
	HelixFlag GetHelixFlag(uint32 inHelixStride) const
	{
		assert(inHelixStride == 3 or inHelixStride == 4 or inHelixStride == 5);
		return mHelixFlags[inHelixStride - 3];
	}

	bool IsHelixStart(uint32 inHelixStride) const
	{
		assert(inHelixStride == 3 or inHelixStride == 4 or inHelixStride == 5);
		return mHelixFlags[inHelixStride - 3] == helixStart or mHelixFlags[inHelixStride - 3] == helixStartAndEnd;
	}

	void SetHelixFlag(uint32 inHelixStride, HelixFlag inHelixFlag)
	{
		assert(inHelixStride == 3 or inHelixStride == 4 or inHelixStride == 5);
		mHelixFlags[inHelixStride - 3] = inHelixFlag;
	}

	void SetSSBridgeNr(uint8 inBridgeNr)
	{
		if (mType != kCysteine)
			throw runtime_error("Only cysteine residues can form sulphur bridges");
		mSSBridgeNr = inBridgeNr;
	}
	
	uint8 GetSSBridgeNr() const
	{
		if (mType != kCysteine)
			throw runtime_error("Only cysteine residues can form sulphur bridges");
		return mSSBridgeNr;
	}

	Res* mNext = nullptr;
	Res* mPrev = nullptr;
	Monomer	mM;

	Point mCAlpha, mC, mN, mO, mH;
	
	ResidueType mType;
	uint8 mSSBridgeNr = 0;
	SecondaryStructureType mSecondaryStructure = ssLoop;
	HBond mHBondDonor[2] = {}, mHBondAcceptor[2] = {};
	BridgeParner mBetaPartner[2] = {};
	uint32 mSheet = 0;
	HelixFlag mHelixFlags[3] = { helixNone, helixNone, helixNone };	//
	bool mBend = false;
};

// --------------------------------------------------------------------
// TODO: use the angle to improve bond energy calculation.

double CalculateHBondEnergy(Res& inDonor, Res& inAcceptor)
{
	double result = 0;
	
	if (inDonor.mType != kProline)
	{
		double distanceHO = Distance(inDonor.mH, inAcceptor.mO);
		double distanceHC = Distance(inDonor.mH, inAcceptor.mC);
		double distanceNC = Distance(inDonor.mN, inAcceptor.mC);
		double distanceNO = Distance(inDonor.mN, inAcceptor.mO);
		
		if (distanceHO < kMinimalDistance or distanceHC < kMinimalDistance or distanceNC < kMinimalDistance or distanceNO < kMinimalDistance)
			result = kMinHBondEnergy;
		else
			result = kCouplingConstant / distanceHO - kCouplingConstant / distanceHC + kCouplingConstant / distanceNC - kCouplingConstant / distanceNO;

		// DSSP compatibility mode:
		result = round(result * 1000) / 1000;

		if (result < kMinHBondEnergy)
			result = kMinHBondEnergy;
	}

	// update donor
	if (result < inDonor.mHBondAcceptor[0].energy)
	{
		inDonor.mHBondAcceptor[1] = inDonor.mHBondAcceptor[0];
		inDonor.mHBondAcceptor[0].residue = &inAcceptor;
		inDonor.mHBondAcceptor[0].energy = result;
	}
	else if (result < inDonor.mHBondAcceptor[1].energy)
	{
		inDonor.mHBondAcceptor[1].residue = &inAcceptor;
		inDonor.mHBondAcceptor[1].energy = result;
	}		

	// and acceptor
	if (result < inAcceptor.mHBondDonor[0].energy)
	{
		inAcceptor.mHBondDonor[1] = inAcceptor.mHBondDonor[0];
		inAcceptor.mHBondDonor[0].residue = &inDonor;
		inAcceptor.mHBondDonor[0].energy = result;
	}
	else if (result < inAcceptor.mHBondDonor[1].energy)
	{
		inAcceptor.mHBondDonor[1].residue = &inDonor;
		inAcceptor.mHBondDonor[1].energy = result;
	}		
	
	return result;
}


// --------------------------------------------------------------------

void CalculateHBondEnergies(vector<Res>& inResidues)
{
	// Calculate the HBond energies
	for (uint32 i = 0; i + 1 < inResidues.size(); ++i)
	{
		auto& ri = inResidues[i];
		
		for (uint32 j = i + 1; j < inResidues.size(); ++j)
		{
			auto& rj = inResidues[j];
			
			if (Distance(ri.mCAlpha, rj.mCAlpha) < kMinimalCADistance)
			{
				CalculateHBondEnergy(ri, rj);
				if (j != i + 1)
					CalculateHBondEnergy(rj, ri);
			}
		}
	}
}

// --------------------------------------------------------------------

bool NoChainBreak(const Res* a, const Res* b)
{
	return a->mM.asymID() == b->mM.asymID();
}

bool NoChainBreak(const Res& a, const Res& b)
{
	return a.mM.asymID() == b.mM.asymID();
}

// --------------------------------------------------------------------

bool TestBond(const Res* a, const Res* b)
{
	return
		(a->mHBondAcceptor[0].residue == b and a->mHBondAcceptor[0].energy < kMaxHBondEnergy) or
		(a->mHBondAcceptor[1].residue == b and a->mHBondAcceptor[1].energy < kMaxHBondEnergy);
}

// --------------------------------------------------------------------

BridgeType TestBridge(const Res& r1, const Res& r2)
{										// I.	a	d	II.	a	d		parallel    
	auto a = r1.mPrev;					//		  \			  /
	auto b = &r1;						//		b	e		b	e
	auto c = r1.mNext;					// 		  /			  \                      ..
	auto d = r2.mPrev;					//		c	f		c	f
	auto e = &r2;						//
	auto f = r2.mNext;					// III.	a <- f	IV. a	  f		antiparallel
										//		                                   
	BridgeType result = btNoBridge;		//		b	 e      b <-> e                  
										//                                          
										//		c -> d		c     d
										
	if (a and c and NoChainBreak(a, c) and d and f and NoChainBreak(d, f))
	{
		if ((TestBond(c, e) and TestBond(e, a)) or (TestBond(f, b) and TestBond(b, d)))
			result = btParallel;
		else if ((TestBond(c, d) and TestBond(f, a)) or (TestBond(e, b) and TestBond(b, e)))
			result = btAntiParallel;
	}
	
	return result;
}

// --------------------------------------------------------------------
// return true if any of the residues in bridge a is identical to any of the residues in bridge b
bool Linked(const Bridge& a, const Bridge& b)
{
	return
		find_first_of(a.i.begin(), a.i.end(), b.i.begin(), b.i.end()) != a.i.end() or
		find_first_of(a.i.begin(), a.i.end(), b.j.begin(), b.j.end()) != a.i.end() or
		find_first_of(a.j.begin(), a.j.end(), b.i.begin(), b.i.end()) != a.j.end() or
		find_first_of(a.j.begin(), a.j.end(), b.j.begin(), b.j.end()) != a.j.end();
}

// --------------------------------------------------------------------

void CalculateBetaSheets(vector<Res>& inResidues)
{
	// Calculate Bridges
	vector<Bridge> bridges;
	if (inResidues.size() > 4)
	{
		for (uint32 i = 1; i + 4 < inResidues.size(); ++i)
		{
			auto& ri = inResidues[i];
			
			for (uint32 j = i + 3; j + 1 < inResidues.size(); ++j)
			{
				auto& rj = inResidues[j];
				
				BridgeType type = TestBridge(ri, rj);
				if (type == btNoBridge)
					continue;
				
				bool found = false;
				for (Bridge& bridge : bridges)
				{
					if (type != bridge.type or i != bridge.i.back() + 1)
						continue;
					
					if (type == btParallel and bridge.j.back() + 1 == j)
					{
						bridge.i.push_back(i);
						bridge.j.push_back(j);
						found = true;
						break;
					}
	
					if (type == btAntiParallel and bridge.j.front() - 1 == j)
					{
						bridge.i.push_back(i);
						bridge.j.push_front(j);
						found = true;
						break;
					}
				}
				
				if (not found)
				{
					Bridge bridge = {};
					
					bridge.type = type;
					bridge.i.push_back(i);
					bridge.chainI = ri.mM.asymID();
					bridge.j.push_back(j);
					bridge.chainJ = rj.mM.asymID();
					
					bridges.push_back(bridge);
				}
			}
		}
	}

	// extend ladders
	sort(bridges.begin(), bridges.end());
	
	for (uint32 i = 0; i < bridges.size(); ++i)
	{
		for (uint32 j = i + 1; j < bridges.size(); ++j)
		{
			uint32 ibi = bridges[i].i.front();
			uint32 iei = bridges[i].i.back();
			uint32 jbi = bridges[i].j.front();
			uint32 jei = bridges[i].j.back();
			uint32 ibj = bridges[j].i.front();
			uint32 iej = bridges[j].i.back();
			uint32 jbj = bridges[j].j.front();
			uint32 jej = bridges[j].j.back();

			if (bridges[i].type != bridges[j].type or
				NoChainBreak(inResidues[min(ibi, ibj)], inResidues[max(iei, iej)]) == false or
				NoChainBreak(inResidues[min(jbi, jbj)], inResidues[max(jei, jej)]) == false or
				ibj - iei >= 6 or
				(iei >= ibj and ibi <= iej))
			{
				continue;
			}
			
			bool bulge;
			if (bridges[i].type == btParallel)
				bulge = ((jbj - jei < 6 and ibj - iei < 3) or (jbj - jei < 3));
			else
				bulge = ((jbi - jej < 6 and ibj - iei < 3) or (jbi - jej < 3));

			if (bulge)
			{
				bridges[i].i.insert(bridges[i].i.end(), bridges[j].i.begin(), bridges[j].i.end());
				if (bridges[i].type == btParallel)
					bridges[i].j.insert(bridges[i].j.end(), bridges[j].j.begin(), bridges[j].j.end());
				else
					bridges[i].j.insert(bridges[i].j.begin(), bridges[j].j.begin(), bridges[j].j.end());
				bridges.erase(bridges.begin() + j);
				--j;
			}
		}
	}

	// Sheet
	set<Bridge*> ladderset;
	for (Bridge& bridge : bridges)
	{
		ladderset.insert(&bridge);
		
//		uint32 n = bridge.i.size();
//		if (n > kHistogramSize)
//			n = kHistogramSize;
//		
//		if (bridge.type == btParallel)
//			mParallelBridgesPerLadderHistogram[n - 1] += 1;
//		else
//			mAntiparallelBridgesPerLadderHistogram[n - 1] += 1;
	}
	
	uint32 sheet = 1, ladder = 0;
	while (not ladderset.empty())
	{
		set<Bridge*> sheetset;
		sheetset.insert(*ladderset.begin());
		ladderset.erase(ladderset.begin());

		bool done = false;
		while (not done)
		{
			done = true;
			for (Bridge* a : sheetset)
			{
				for (Bridge* b : ladderset)
				{
					if (Linked(*a, *b))
					{
						sheetset.insert(b);
						ladderset.erase(b);
						done = false;
						break;
					}
				}
				if (not done)
					break;
			}
		}

		for (Bridge* bridge : sheetset)
		{
			bridge->ladder = ladder;
			bridge->sheet = sheet;
			bridge->link = sheetset;
			
			++ladder;
		}
		
//		uint32 nrOfLaddersPerSheet = sheetset.size();
//		if (nrOfLaddersPerSheet > kHistogramSize)
//			nrOfLaddersPerSheet = kHistogramSize;
//		if (nrOfLaddersPerSheet == 1 and (*sheetset.begin())->i.size() > 1)
//			mLaddersPerSheetHistogram[0] += 1;
//		else if (nrOfLaddersPerSheet > 1)
//			mLaddersPerSheetHistogram[nrOfLaddersPerSheet - 1] += 1;
		
		++sheet;
	}

	for (Bridge& bridge : bridges)
	{
		// find out if any of the i and j set members already have
		// a bridge assigned, if so, we're assigning bridge 2
		
		uint32 betai = 0, betaj = 0;
		
		for (uint32 l : bridge.i)
		{
			if (inResidues[l].GetBetaPartner(0).residue != nullptr)
			{
				betai = 1;
				break;
			}
		}

		for (uint32 l : bridge.j)
		{
			if (inResidues[l].GetBetaPartner(0).residue != nullptr)
			{
				betaj = 1;
				break;
			}
		}
		
		SecondaryStructureType ss = ssBetabridge;
		if (bridge.i.size() > 1)
			ss = ssStrand;
		
		if (bridge.type == btParallel)
		{
//			mNrOfHBondsInParallelBridges += bridge.i.back() - bridge.i.front() + 2;
			
			deque<uint32>::iterator j = bridge.j.begin();
			for (uint32 i : bridge.i)
				inResidues[i].SetBetaPartner(betai, inResidues[*j++], bridge.ladder, true);

			j = bridge.i.begin();
			for (uint32 i : bridge.j)
				inResidues[i].SetBetaPartner(betaj, inResidues[*j++], bridge.ladder, true);
		}
		else
		{
//			mNrOfHBondsInAntiparallelBridges += bridge.i.back() - bridge.i.front() + 2;

			deque<uint32>::reverse_iterator j = bridge.j.rbegin();
			for (uint32 i : bridge.i)
				inResidues[i].SetBetaPartner(betai, inResidues[*j++], bridge.ladder, false);

			j = bridge.i.rbegin();
			for (uint32 i : bridge.j)
				inResidues[i].SetBetaPartner(betaj, inResidues[*j++], bridge.ladder, false);
		}

		for (uint32 i = bridge.i.front(); i <= bridge.i.back(); ++i)
		{
			if (inResidues[i].GetSecondaryStructure() != ssStrand)
				inResidues[i].SetSecondaryStructure(ss);
			inResidues[i].SetSheet(bridge.sheet);
		}

		for (uint32 i = bridge.j.front(); i <= bridge.j.back(); ++i)
		{
			if (inResidues[i].GetSecondaryStructure() != ssStrand)
				inResidues[i].SetSecondaryStructure(ss);
			inResidues[i].SetSheet(bridge.sheet);
		}
	}
}

// --------------------------------------------------------------------
// TODO: improve alpha helix calculation by better recognizing pi-helices 

void CalculateAlphaHelices(vector<Res>& inResidues, bool inPreferPiHelices = true)
{
	// Helix and Turn
	for (uint32 stride = 3; stride <= 5; ++stride)
	{
		for (uint32 i = 0; i + stride < inResidues.size(); ++i)
		{
			if (NoChainBreak(inResidues[i], inResidues[i + stride]) and TestBond(&inResidues[i + stride], &inResidues[i]))
			{
				inResidues[i + stride].SetHelixFlag(stride, helixEnd);
				for (uint32 j = i + 1; j < i + stride; ++j)
				{
					if (inResidues[j].GetHelixFlag(stride) == helixNone)
						inResidues[j].SetHelixFlag(stride, helixMiddle);
				}
				
				if (inResidues[i].GetHelixFlag(stride) == helixEnd)
					inResidues[i].SetHelixFlag(stride, helixStartAndEnd);
				else
					inResidues[i].SetHelixFlag(stride, helixStart);
			}
		}
	}
	
	for (auto& r : inResidues)
	{
		double kappa = r.mM.kappa();
		r.SetBend(kappa != 360 and kappa > 70);
	}

	for (uint32 i = 1; i + 4 < inResidues.size(); ++i)
	{
		if (inResidues[i].IsHelixStart(4) and inResidues[i - 1].IsHelixStart(4))
		{
			for (uint32 j = i; j <= i + 3; ++j)
				inResidues[j].SetSecondaryStructure(ssAlphahelix);
		}
	}

	for (uint32 i = 1; i + 3 < inResidues.size(); ++i)
	{
		if (inResidues[i].IsHelixStart(3) and inResidues[i - 1].IsHelixStart(3))
		{
			bool empty = true;
			for (uint32 j = i; empty and j <= i + 2; ++j)
				empty = inResidues[j].GetSecondaryStructure() == ssLoop or inResidues[j].GetSecondaryStructure() == ssHelix_3;
			if (empty)
			{
				for (uint32 j = i; j <= i + 2; ++j)
					inResidues[j].SetSecondaryStructure(ssHelix_3);
			}
		}
	}

	for (uint32 i = 1; i + 5 < inResidues.size(); ++i)
	{
		if (inResidues[i].IsHelixStart(5) and inResidues[i - 1].IsHelixStart(5))
		{
			bool empty = true;
			for (uint32 j = i; empty and j <= i + 4; ++j)
				empty = inResidues[j].GetSecondaryStructure() == ssLoop or inResidues[j].GetSecondaryStructure() == ssHelix_5 or
							(inPreferPiHelices and inResidues[j].GetSecondaryStructure() == ssAlphahelix);
			if (empty)
			{
				for (uint32 j = i; j <= i + 4; ++j)
					inResidues[j].SetSecondaryStructure(ssHelix_5);
			}
		}
	}
			
	for (uint32 i = 1; i + 1 < inResidues.size(); ++i)
	{
		if (inResidues[i].GetSecondaryStructure() == ssLoop)
		{
			bool isTurn = false;
			for (uint32 stride = 3; stride <= 5 and not isTurn; ++stride)
			{
				for (uint32 k = 1; k < stride and not isTurn; ++k)
					isTurn = (i >= k) and inResidues[i - k].IsHelixStart(stride);
			}
			
			if (isTurn)
				inResidues[i].SetSecondaryStructure(ssTurn);
			else if (inResidues[i].IsBend())
				inResidues[i].SetSecondaryStructure(ssBend);
		}
	}
}

// --------------------------------------------------------------------

void CalculateSecondaryStructure(Structure& s)
{
	auto polymers = s.polymers();
	
	size_t nRes = accumulate(polymers.begin(), polymers.end(),
		0.0, [](double s, auto& p) { return s + p.size(); });
	
	vector<Res> residues;
	residues.reserve(nRes);
	
	for (auto& p: polymers)
	{
		for (auto m: p)
			residues.emplace_back(move(m));
	}
	
	for (size_t i = 0; i + 1 < residues.size(); ++i)
	{
		residues[i].mNext = &residues[i + 1];
		residues[i + 1].mPrev = &residues[i];
		
		residues[i + 1].assignHydrogen();
	}
	
	CalculateHBondEnergies(residues);
	CalculateBetaSheets(residues);
	CalculateAlphaHelices(residues);
	
	if (VERBOSE)
	{
		for (auto& r: residues)
		{
			auto& m = r.mM;
			
			cout << m.asymID() << ':' << m.seqID() << '/' << m.compoundID() << '\t'
				 << char(r.mSecondaryStructure)
				 << endl;
		}
	}
}

// --------------------------------------------------------------------

struct DSSPImpl
{
	DSSPImpl(const Structure& s);
	
	const Structure&	mStructure;
	vector<Polymer>		mPolymers;
	vector<Res>			mResidues;
};

DSSPImpl::DSSPImpl(const Structure& s)
	: mStructure(s)
	, mPolymers(mStructure.polymers())
{
	size_t nRes = accumulate(mPolymers.begin(), mPolymers.end(),
		0.0, [](double s, auto& p) { return s + p.size(); });
	
	mResidues.reserve(nRes);
	
	for (auto& p: mPolymers)
	{
		for (auto m: p)
			mResidues.emplace_back(move(m));
	}
	
	for (size_t i = 0; i + 1 < mResidues.size(); ++i)
	{
		mResidues[i].mNext = &mResidues[i + 1];
		mResidues[i + 1].mPrev = &mResidues[i];
		
		mResidues[i + 1].assignHydrogen();
	}
	
	CalculateHBondEnergies(mResidues);
	CalculateBetaSheets(mResidues);
	CalculateAlphaHelices(mResidues);
	
	if (VERBOSE)
	{
		for (auto& r: mResidues)
		{
			auto& m = r.mM;
			
			char helix[4] = { };
			for (size_t stride: { 3, 4, 5 })
			{
				switch (r.GetHelixFlag(stride))
				{
					case helixStart:		helix[stride - 3] = '>'; break;
					case helixMiddle:		helix[stride - 3] = '0' + stride; break;
					case helixStartAndEnd:	helix[stride - 3] = 'X'; break;
					case helixEnd:			helix[stride - 3] = '<'; break;
					case helixNone:			helix[stride - 3] = ' '; break;
				}
			}
			
			string id = m.asymID() + ':' + to_string(m.seqID()) + '/' + m.compoundID();
			
			cout << id << string(12 - id.length(), ' ')
				 << char(r.mSecondaryStructure) << ' '
				 << helix
				 << endl;
		}
	}
}

DSSP::DSSP(const Structure& s)
	: mImpl(new DSSPImpl(s))
{
}

DSSP::~DSSP()
{
	delete mImpl;
}

SecondaryStructureType DSSP::operator()(const string& inAsymID, int inSeqID) const
{
	SecondaryStructureType result = ssLoop;
	auto i = find_if(mImpl->mResidues.begin(), mImpl->mResidues.end(),
		[&](auto& r) { return r.mM.asymID() == inAsymID and r.mM.seqID() == inSeqID; });
	if (i != mImpl->mResidues.end())
		result = i->mSecondaryStructure;
	else if (VERBOSE)
		cerr << "Could not find secondary structure for " << inAsymID << ':' << inSeqID << endl;
	return result;
}

SecondaryStructureType DSSP::operator()(const Monomer& m) const
{
	return operator()(m.asymID(), m.seqID());
}

bool DSSP::isAlphaHelixEndBeforeStart(const Monomer& m) const
{
	return isAlphaHelixEndBeforeStart(m.asymID(), m.seqID());
}

bool DSSP::isAlphaHelixEndBeforeStart(const string& inAsymID, int inSeqID) const
{
	auto i = find_if(mImpl->mResidues.begin(), mImpl->mResidues.end(),
		[&](auto& r) { return r.mM.asymID() == inAsymID and r.mM.seqID() == inSeqID; });

	bool result = false;

	if (i != mImpl->mResidues.end() and i + 1 != mImpl->mResidues.end())
		result = i->GetHelixFlag(4) == helixEnd and (i + 1)->GetHelixFlag(4) == helixStart;
	else if (VERBOSE)
		cerr << "Could not find secondary structure for " << inAsymID << ':' << inSeqID << endl;

	return result;
}

}
