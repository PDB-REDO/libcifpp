// --------------------------------------------------------------------

#pragma once

namespace mmcif
{

class Structure;

/// \brief Add missing backbone atoms
///
/// \param structure	The structure that should be fixed
/// \param simplified	Use a simplified algorithm

void CreateMissingBackboneAtoms(Structure& structure, bool simplified);


}