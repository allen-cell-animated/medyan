#ifndef MEDYAN_GTriangle_h
#define MEDYAN_GTriangle_h

#include "MathFunctions.h"

/******************************************************************************
Storing the geometric properties of the triangle patches.
******************************************************************************/

struct GTriangle {

    double area; // Current area
    double sArea; // Temporarily store the stretched area

    mathfunc::Vec3 unitNormal; // The unit normal vector pointing outward (since the meshwork is orientable)
    mathfunc::Vec3 sUnitNormal; // Temporarily stores unit normal under stretched conditions.

    double coneVolume; // Volume of the tetrahedral formed by this triangle and the origin (0, 0, 0)
    double sConeVolume;

    // Auxilliary getters
    template< bool stretched > double& getArea() { return stretched ? sArea : area; }
    template< bool stretched > mathfunc::Vec3& getUnitNormal() { return stretched ? sUnitNormal : unitNormal; }
    template< bool stretched > double& getConeVolume() { return stretched ? sConeVolume : coneVolume; }

};


#endif
