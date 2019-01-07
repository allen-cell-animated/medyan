#ifndef MEDYAN_GHalfEdge_h
#define MEDYAN_GHalfEdge_h

#include <array>

#include "MathFunctions.h"

struct GHalfEdge {

    double theta; // Angle formed by (this, opposite(next(this)))
    double sTheta; // Stretched theta
    double cotTheta;
    double sCotTheta;
    std::array<mathfunc::Vec3, 3> dCotTheta; // Indexed by [(source, target, target(next))]

    mathfunc::Vec3 dEdgeLength; // Derivative of length of edge on target. FIXME the source derivative wont exist if at boundary.
    mathfunc::Vec3 dTriangleArea; // Derivative of area of triangle on target

    mathfunc::Vec3 dNeighborArea; // Derivative of area of vcell of source on target
    mathfunc::Vec3 dNeighborCurv; // Derivative of curv of vcell of source on target
    
    // Auxilliary getters
    template< bool stretched > auto& getTheta() { return stretched ? sTheta : theta; }
    template< bool stretched > auto& getCotTheta() { return stretched ? sCotTheta : cotTheta; }
};


#endif
