
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "BoundarySurfaceImpl.h"

#include "BoundaryElementImpl.h"

#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

Plane::Plane(vector<double> coords, vector<double> normal ) :
    _coords(coords), _normal(normal), BoundarySurface(3) {
    
    //Create a plane boundary element
    _boundaryElements.emplace_back(new PlaneBoundaryElement(coords, normal,
                                   SystemParameters::Boundaries().boundaryK,
                                   SystemParameters::Boundaries().screenLength));
}

Sphere::Sphere(vector<double> coords, double radius) : _coords(coords), _radius(radius), BoundarySurface(3){
    
    //Create a sphere boundary element
    _boundaryElements.emplace_back(new SphereBoundaryElement(coords, radius,
                                   SystemParameters::Boundaries().boundaryK,
                                   SystemParameters::Boundaries().screenLength));
    
}

CylinderZ::CylinderZ(vector<double> coords, double radius, double height)
                     : _coords(coords), _radius(radius), _height(height), BoundarySurface(3){
    
    //Create a cylindricalZ boundary element
    _boundaryElements.emplace_back(new CylindricalZBoundaryElement(coords, radius, height,
                                   SystemParameters::Boundaries().boundaryK,
                                   SystemParameters::Boundaries().screenLength));
}

HalfSphereZ::HalfSphereZ(vector<double> coords, double radius, bool up)
                          : _coords(coords), _radius(radius), _up(up), BoundarySurface(3){
    
    //Create a half sphere Z boundary element
    _boundaryElements.emplace_back(new HalfSphereZBoundaryElement(coords, radius, up,
                                   SystemParameters::Boundaries().boundaryK,
                                   SystemParameters::Boundaries().screenLength));
    
}