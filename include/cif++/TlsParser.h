/* 
   Created by: Maarten L. Hekkelman
   Date: dinsdag 07 november, 2017

	Copyright 2017 NKI AVL
	
	Permission is hereby granted, free of charge, to any person obtaining 
	a copy of this software and associated documentation files (the 
	"Software"), to deal in the Software without restriction, including 
	without limitation the rights to use, copy, modify, merge, publish, 
	distribute, sublicense, and/or sell copies of the Software, and to 
	permit persons to whom the Software is furnished to do so, subject to 
	the following conditions:
	
	The above copyright notice and this permission notice shall be 
	included in all copies or substantial portions of the Software.
	
	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
	MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
	IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
	CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
	TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
	SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <vector>
#include <string>
#include <tuple>

#include "cif++/Cif++.h"

namespace cif
{
	
extern const int
	kResidueNrWildcard,
	kNoSeqNum;

struct TLSSelection;
typedef std::unique_ptr<TLSSelection> TLSSelectionPtr;

struct TLSResidue;

struct TLSSelection
{
	virtual ~TLSSelection() {}
	virtual void CollectResidues(cif::Datablock& db, std::vector<TLSResidue>& residues, int indentLevel = 0) const = 0;
	std::vector<std::tuple<std::string,int,int>> GetRanges(cif::Datablock& db, bool pdbNamespace) const;
};

// Low level: get the selections
TLSSelectionPtr ParseSelectionDetails(const std::string& program, const std::string& selection);

}
