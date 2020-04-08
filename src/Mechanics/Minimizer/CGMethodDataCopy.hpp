#ifndef MEDYAN_Mechanics_Minimizer_CGMethodDataCopy_hpp
#define MEDYAN_Mechanics_Minimizer_CGMethodDataCopy_hpp

#include <algorithm>

#include "Mechanics/ForceField/Types.hpp"
#include "Mechanics/Minimizer/CGMethod.h"
#include "Structure/Bead.h"

// Copies all the system data to the CGMethod data vector
inline FFCoordinateStartingIndex initCGMethodData(
    CGMethod& cg,
    floatingpoint defaultGradTol
) {
    FFCoordinateStartingIndex si {};
    std::size_t curIdx = 0;
    cg.coord.clear();
    cg.forceTol.clear();

    //---------------------------------
    // Copy all the coordinate information here
    // Also initializes the force tolerance vector

    // Bead coord
    si.bead = curIdx;
    cg.coord.reserve(cg.coord.size() + 3 * Bead::getBeads().size());
    for(auto pb : Bead::getBeads()) {
        cg.coord.insert(cg.coord.end(), pb->coord.begin(), pb->coord.end());
        curIdx += 3;
    }
    cg.forceTol.resize(cg.coord.size(), defaultGradTol);

    // Vertex coord
    si.vertex = curIdx;
    // Add things for vertices here

    // Membrane 2d coord
    si.mem2d = curIdx;
    // Add things for membrane 2d coordinates here

    //---------------------------------
    // Initialize the rest of cg data
    const auto ndof = cg.coord.size();
    cg.coordLineSearch.assign(ndof, 0);
    cg.force.assign(ndof, 0);
    cg.forcePrev.assign(ndof, 0);
    cg.searchDir.assign(ndof, 0);

    //---------------------------------
    // Return the starting index information for vectorizing the force fields
    return si;
}

// Copies all the CGMethod data back to the system
//
// Note:
//   - The copying must be in the same order with the initCGMethodData
//     function.
inline void copyFromCGMethodData(const CGMethod& cg) {
    std::size_t curIdx = 0;

    // Bead
    for(auto pb : Bead::getBeads()) {
        std::copy(cg.coord.begin() + curIdx, cg.coord.begin() + curIdx + 3, pb->coord.begin());
        std::copy(cg.force.begin() + curIdx, cg.force.begin() + curIdx + 3, pb->force.begin());
        curIdx += 3;
    }

    // Vertex

    // Membrane 2d coord

    // Do not clear CGMethod data, because it might be useful for debug purposes.
}

#endif
