#ifndef MEDYAN_GeometricMeshAttribute_hpp
#define MEDYAN_GeometricMeshAttribute_hpp

#include "MathFunctions.h"

struct GVertex {

    double astar; // 1-ring area
    mathfunc::Vec3d dAstar; // Derivative of 1-ring area on the central vertex. Derivatives on neighbors are stored in half edges

    double curv; // Current mean curvature
    mathfunc::Vec3d dCurv;
    double curv2; // Square curvature
    mathfunc::Vec3d dCurv2;

    mathfunc::Vec3d dVolume; // Derivative of volume on this vertex

    mathfunc::Vec3d pseudoUnitNormal; // Pseudo unit normal around the vertex

};

struct GHalfEdge {

    double theta; // Angle formed by (this, opposite(next(this)))
    double cotTheta;
    std::array<mathfunc::Vec3, 3> dCotTheta; // Indexed by [(source, target, target(next))]

    mathfunc::Vec3d dTriangleArea; // Derivative of area of triangle on target

    mathfunc::Vec3d dNeighborAstar; // Derivative (on target vertex) of 1-ring area of source vertex
    mathfunc::Vec3d dNeighborCurv; // Derivative of curv of vertex of source on target
    mathfunc::Vec3d dNeighborCurv2; // Derivative of curv2 of vertex of source on target
    
};

struct GEdge {

    mathfunc::Vec3d pseudoUnitNormal; // The pseudo unit normal vector at the edge pointing outward.

};

struct GTriangle {

    double area; // Current area

    mathfunc::Vec3d unitNormal; // The unit normal vector pointing outward (since the meshwork is orientable)

    double coneVolume; // Volume of the tetrahedral formed by this triangle and the origin (0, 0, 0)

};

#endif
