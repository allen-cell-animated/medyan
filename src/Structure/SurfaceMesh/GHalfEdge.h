#ifndef MEDYAN_GHalfEdge_h
#define MEDYAN_GHalfEdge_h

#include <array>
#include <vector>

struct GHalfEdge {

    double theta; // Angle formed by (this, opposite(next(this)))
    std::array<std::array<double, 3>, 3> dTheta; // Indexed by [(source, target, target(next))][coord]
    double cotTheta;
    std::array<std::array<double, 3>, 3> dCotTheta;

    std::array<double, 3> dEdgeLength; // Derivative of length of edge on target. FIXME the source derivative wont exist if at boundary.
    std::array<double, 3> dTriangleArea; // Derivative of area of triangle on target

    std::array<double, 3> dNeighborArea; // Derivative of area of vcell of source on target
    std::array<double, 3> dNeighborCurv; // Derivative of curv of vcell of source on target

};


#endif
