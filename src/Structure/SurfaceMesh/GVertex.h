#ifndef MEDYAN_GVertex_h
#define MEDYAN_GVertex_h

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


#endif
