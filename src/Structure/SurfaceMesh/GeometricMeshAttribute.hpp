#ifndef MEDYAN_GeometricMeshAttribute_hpp
#define MEDYAN_GeometricMeshAttribute_hpp

#include "MathFunctions.h"

struct GVertex {

    double area; // Current area
    mathfunc::Vec3 dArea; // Derivative of area on the central vertex, derivative on neighbors are stored in half edges
    double astar; // 1-ring area
    mathfunc::Vec3 dAstar; // Derivative of 1-ring area on the central vertex. Derivatives on neighbors are stored in half edges

    double curv; // Current mean curvature
    mathfunc::Vec3 dCurv;

    mathfunc::Vec3 dVolume; // Derivative of volume on this vertex

    mathfunc::Vec3 pseudoUnitNormal; // Pseudo unit normal around the vertex

};

struct GHalfEdge {

    double theta; // Angle formed by (this, opposite(next(this)))
    double cotTheta;
    std::array<mathfunc::Vec3, 3> dCotTheta; // Indexed by [(source, target, target(next))]

    mathfunc::Vec3 dEdgeLength; // Derivative of length of edge on target. FIXME the source derivative wont exist if at boundary.
    mathfunc::Vec3 dTriangleArea; // Derivative of area of triangle on target

    mathfunc::Vec3 dNeighborArea; // Derivative of area of vcell of source on target
    mathfunc::Vec3 dNeighborAstar; // Derivative (on target vertex) of 1-ring area of source vertex
    mathfunc::Vec3 dNeighborCurv; // Derivative of curv of vcell of source on target
    
};

struct GEdge {

    double length; // Current length

    mathfunc::Vec3 pseudoUnitNormal; // The pseudo unit normal vector at the edge pointing outward.

};

struct GTriangle {

    double area; // Current area

    mathfunc::Vec3 unitNormal; // The unit normal vector pointing outward (since the meshwork is orientable)

    double coneVolume; // Volume of the tetrahedral formed by this triangle and the origin (0, 0, 0)

};

#endif
