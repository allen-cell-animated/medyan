#ifndef MEDYAN_Structure_DofSerializer_hpp
#define MEDYAN_Structure_DofSerializer_hpp

// The file provides functions to do data copy between elements in the system
// and the vectorized data in the mechanics energy/force calculations.
//
// Functions:
//
//   - serializeDof(...)
//
//     Copy the degree-of-freedom data from all system elements (such as the
//     bead coordinates) to the coordinate array.
//
//     Returns the starting indices of different types of elements, which is
//     useful when building interactions in force fields.
//
//   - deserializeDof(...)
//
//     Copy the vectorized coordinate and force data to all the element
//     instances in the system.

#include <algorithm>

#include "Mechanics/ForceField/Types.hpp"
#include "Structure/Bead.h"
#include "Structure/Bubble.h"
#include "Structure/SurfaceMesh/Vertex.hpp"

namespace medyan {

// Copies all the system data to the CGMethod data vector
inline FFCoordinateStartingIndex serializeDof(
    std::vector< floatingpoint >& coord
) {
    FFCoordinateStartingIndex si {};
    std::size_t curIdx = 0;
    coord.clear();

    //---------------------------------
    // Copy all the coordinate information here
    // Also initializes the force tolerance vector

    // Bead coord
    si.bead = curIdx;
    coord.reserve(coord.size() + 3 * Bead::getBeads().size());
    for(auto pb : Bead::getBeads()) {
        coord.insert(coord.end(), pb->coord.begin(), pb->coord.end());
        curIdx += 3;
    }

    // Bubble coord
    si.bubble = curIdx;
    coord.reserve(coord.size() + 3 * Bubble::getBubbles().size());
    for(auto pb : Bubble::getBubbles()) {
        coord.insert(coord.end(), pb->coord.begin(), pb->coord.end());
        curIdx += 3;
    }

    // Vertex coord
    si.vertex = curIdx;
    coord.reserve(coord.size() + 3 * Vertex::getVertices().size());
    for(auto pv : Vertex::getVertices()) {
        coord.insert(coord.end(), pv->coord.begin(), pv->coord.end());
        curIdx += 3;
    }

    // Membrane 2d coord
    si.mem2d = curIdx;
    // Add things for membrane 2d coordinates here

    //---------------------------------
    // The independent variables have all been assigned.
    si.ndof = curIdx;

    //---------------------------------
    // Return the starting index information for vectorizing the force fields
    return si;
}

// Copies all the CGMethod data back to the system
//
// Note:
//   - The copying must be in the same order with the initCGMethodData
//     function.
inline void deserializeDof(
    const std::vector< floatingpoint >& coord,
    const std::vector< floatingpoint >& force
) {
    std::size_t curIdx = 0;

    // Copy coord and force data to beads
    for(auto pb : Bead::getBeads()) {
        std::copy(coord.begin() + curIdx, coord.begin() + curIdx + 3, pb->coord.begin());
        std::copy(force.begin() + curIdx, force.begin() + curIdx + 3, pb->force.begin());
        curIdx += 3;
    }

    // Copy coord and force data to bubbles
    for(auto pb : Bubble::getBubbles()) {
        std::copy(coord.begin() + curIdx, coord.begin() + curIdx + 3, pb->coord.begin());
        std::copy(force.begin() + curIdx, force.begin() + curIdx + 3, pb->force.begin());
        curIdx += 3;
    }

    // Vertex
    for(auto pv : Vertex::getVertices()) {
        std::copy(coord.begin() + curIdx, coord.begin() + curIdx + 3, pv->coord.begin());
        std::copy(force.begin() + curIdx, force.begin() + curIdx + 3, pv->force.begin());
        curIdx += 3;
    }

    // Copy coord and force data to vertices

    // Copy coord and force data to Membrane 2d points

}


// Helper function to find the starting index of bead coordinate in the minimizer.
inline int findBeadCoordIndex(
    Bead&                            b,
    const FFCoordinateStartingIndex& si
) {
    // if(SysParams::Mechanics().globalFilamentModel == FilamentModel::bezier) {
    //     const auto& fil = *static_cast< Filament* >(b.getParent());
    //     return 
    //         si.bezierFilPos[fil.getIndex()].start
    //         + fil.getCoordIndexInKvecrWithBeadPosition(b.getPosition());
    // }
    // else if(SysParams::Mechanics().globalFilamentModel == FilamentModel::geodesic) {
    //     const auto& fil = *static_cast< Filament* >(b.getParent());
    //     const auto& gcInfo = si.gcFil[fil.getIndex()];
    //     const int   relPos = fil.getRelativePositionFromMinusEnd(b.getPosition());
    //     return relPos == 0
    //         // bead is the minus end bead
    //         ? gcInfo.rMinusStart
    //         // bead is not the minus end
    //         : gcInfo.rNotMinusStart + 3 * (relPos - 1);
    // }
    // else {
    {
        return b.getIndex() * 3 + si.bead;
    }
}

} // namespace medyan

#endif
