#ifndef MEDYAN_GeometricMeshAttribute_hpp
#define MEDYAN_GeometricMeshAttribute_hpp

#include "MathFunctions.h"

struct GVertex {

    double area; // Current area
    mathfunc::Vec3 dArea; // Derivative of area on the central vertex, derivative on neighbors are stored in half edges
    double sArea; // Temporarily store the stretched area

    double curv; // Current mean curvature
    mathfunc::Vec3 dCurv;
    double sCurv; // Temporarily store the stretched mean curvature

    mathfunc::Vec3 dVolume; // Derivative of volume on this vertex

    mathfunc::Vec3 pseudoUnitNormal; // Pseudo unit normal around the vertex
    mathfunc::Vec3 sPseudoUnitNormal; // Pseudo unit normal under stretch

    // Auxilliary getters
    template< bool stretched > auto& getArea() { return stretched ? sArea : area; }
    template< bool stretched > auto& getCurv() { return stretched ? sCurv : curv; }
    template< bool stretched > auto& getPseudoUnitNormal() { return stretched ? sPseudoUnitNormal : pseudoUnitNormal; }

};

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

struct GEdge {

    double length; // Current length
    double sLength; // Temporarily store the stretched length

    mathfunc::Vec3 pseudoUnitNormal; // The pseudo unit normal vector at the edge pointing outward.
    mathfunc::Vec3 sPseudoUnitNormal; // The pseudo normal under stretch

    // Auxilliary getters
    template< bool stretched > auto& getLength() { return stretched ? sLength : length; }
    template< bool stretched > auto& getPseudoUnitNormal() { return stretched ? sPseudoUnitNormal : pseudoUnitNormal; }
    
};

struct GTriangle {

    double area; // Current area
    double sArea; // Temporarily store the stretched area

    mathfunc::Vec3 unitNormal; // The unit normal vector pointing outward (since the meshwork is orientable)
    mathfunc::Vec3 sUnitNormal; // Temporarily stores unit normal under stretched conditions.

    double coneVolume; // Volume of the tetrahedral formed by this triangle and the origin (0, 0, 0)
    double sConeVolume;

    // Auxilliary getters
    template< bool stretched > auto& getArea() { return stretched ? sArea : area; }
    template< bool stretched > auto& getUnitNormal() { return stretched ? sUnitNormal : unitNormal; }
    template< bool stretched > auto& getConeVolume() { return stretched ? sConeVolume : coneVolume; }

};

#endif
