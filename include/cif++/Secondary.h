/* 
   Created by: Maarten L. Hekkelman
   Date: woensdag 06 juni, 2018
*/

// Calculate DSSP-like secondary structure information

#pragma once

namespace mmcif
{
	
class Structure;
class Monomer;

extern const double kCouplingConstant, kMinHBondEnergy, kMaxHBondEnergy;

enum SecondaryStructureType : char
{
	ssLoop			= ' ',
	ssAlphahelix	= 'H',
	ssBetabridge	= 'B',
	ssStrand		= 'E',
	ssHelix_3		= 'G',
	ssHelix_5		= 'I',
	ssTurn			= 'T',
	ssBend			= 'S'
};

//struct HBond
//{
//	std::string 				labelAsymID;
//	int							labelSeqID;
//	double						energy;
//};
//
//struct BridgePartner
//{
//	std::string					labelAsymID;
//	int							labelSeqID;
//	int							ladder;
//	bool						parallel;
//};

struct SecondaryStructure
{
	SecondaryStructureType		type;
//	HBond						donor[2], acceptor[2];
//	BridgePartner				beta[2];
//	int							sheet;
//	bool						bend;
};

//void CalculateSecondaryStructure(Structure& s);

class DSSP
{
  public:
	DSSP(const Structure& s);
	~DSSP();
	
	DSSP(const DSSP&) = delete;
	DSSP& operator=(const DSSP&) = delete;
	
	SecondaryStructureType operator()(const std::string& inAsymID, int inSeqID) const;
	SecondaryStructureType operator()(const Monomer& m) const;
	
  private:
	struct DSSPImpl* mImpl;
};


}
