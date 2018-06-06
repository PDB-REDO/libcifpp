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

void CalculateSecondaryStructure(Structure& s);

}
