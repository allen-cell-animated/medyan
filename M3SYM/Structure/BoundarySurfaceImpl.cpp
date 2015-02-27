
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "BoundarySurfaceImpl.h"

#include "BoundaryElementImpl.h"

#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

Plane::Plane(vector<double> coords, vector<double> normal ) :
    BoundarySurface(3), _coords(coords), _normal(normal) {
    
    //Create a plane boundary element
    _boundaryElements.emplace_back(new PlaneBoundaryElement(coords, normal,
                                   SysParams::Boundaries().BoundaryK,
                                   SysParams::Boundaries().BScreenLength));
}

Sphere::Sphere(vector<double> coords, double radius)
    : BoundarySurface(3), _coords(coords) {
    
    //Create a sphere boundary element
    _boundaryElements.emplace_back(new SphereBoundaryElement(coords, radius,
                                   SysParams::Boundaries().BoundaryK,
                                   SysParams::Boundaries().BScreenLength));
    
}

CylinderZ::CylinderZ(vector<double> coords, double radius, double height)
    : BoundarySurface(3), _coords(coords) {
    
    //Create a cylindricalZ boundary element
    _boundaryElements.emplace_back(new CylindricalZBoundaryElement(coords, radius, height,
                                   SysParams::Boundaries().BoundaryK,
                                   SysParams::Boundaries().BScreenLength));
}

HalfSphereZ::HalfSphereZ(vector<double> coords, double radius, bool up)
    : BoundarySurface(3), _coords(coords) {
    
    //Create a half sphere Z boundary element
    _boundaryElements.emplace_back(new HalfSphereZBoundaryElement(coords, radius, up,
                                   SysParams::Boundaries().BoundaryK,
                                   SysParams::Boundaries().BScreenLength));
    
}