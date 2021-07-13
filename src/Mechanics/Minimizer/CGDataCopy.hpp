#ifndef MEDYAN_Mechanics_Minimizer_CGDataCopy_hpp
#define MEDYAN_Mechanics_Minimizer_CGDataCopy_hpp

// The file provides functions to do data copy between elements in the system
// and the vectorized data in the Conjugate Gradient method (CGMethod).
//
// Functions:
//
//   - initCGMethodData(...)
//
//     Copy the degree-of-freedom data from all system elements (such as the
//     bead coordinates) to the CGMethod coordinate array, and initialize other
//     CGMethod arrays like the forces.
//
//     Returns the starting indices of different types of elements, which is
//     useful when building interactions in force fields.
//
//   - copyFromCGMethodData(...)
//
//     Copy the vectorized coordinate and force data in the CGMethod to all the
//     element instances in the system.

#include "Structure/DofSerializer.hpp"
#include "Mechanics/Minimizer/CGMethod.h"

// Copies all the system data to the CGMethod data vector
inline FFCoordinateStartingIndex initCGMethodData(
    CGMethod& cg,
    floatingpoint defaultGradTol
) {
    const auto si = medyan::serializeDof(cg.coord);

    cg.numDof = si.ndof;
    cg.forcePrev.assign(cg.numDof, 0);
    cg.searchDir.assign(cg.numDof, 0);
    cg.forceTol.resize(cg.numDof, defaultGradTol); // Future: can set different default grad tol

    const auto nvar = cg.coord.size();
    cg.coordLineSearch.assign(nvar, 0);
    cg.force.assign(nvar, 0);

    return si;
}

// Copies all the CGMethod data back to the system
inline void copyFromCGMethodData(const CGMethod& cg) {
    medyan::deserializeDof(cg.coord, cg.force);
}

#endif
