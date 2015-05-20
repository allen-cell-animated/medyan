
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

#include "SubSystem.h"

#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

Plane::Plane(SubSystem* s, vector<double> coords, vector<double> normal ) :
    BoundarySurface(s, 3), _coords(coords), _normal(normal) {
    
    //Create a plane boundary element
    _boundaryElements.emplace_back(s->addTrackable<PlaneBoundaryElement>
                                   (coords, normal,
                                   SysParams::Boundaries().BoundaryK,
                                   SysParams::Boundaries().BScreenLength));
}

Sphere::Sphere(SubSystem* s, vector<double> coords, double radius)
    : BoundarySurface(s, 3), _coords(coords) {
    
    //Create a sphere boundary element
    _boundaryElements.emplace_back(s->addTrackable<SphereBoundaryElement>
                                   (coords, radius,
                                   SysParams::Boundaries().BoundaryK,
                                   SysParams::Boundaries().BScreenLength));
    
}

CylinderZ::CylinderZ(SubSystem* s, vector<double> coords, double radius, double height)
    : BoundarySurface(s, 3), _coords(coords) {
    
    //Create a cylindricalZ boundary element
    _boundaryElements.emplace_back(s->addTrackable<CylindricalZBoundaryElement>
                                   (coords, radius, height,
                                   SysParams::Boundaries().BoundaryK,
                                   SysParams::Boundaries().BScreenLength));
}

HalfSphereZ::HalfSphereZ(SubSystem* s, vector<double> coords, double radius, bool up)
    : BoundarySurface(s, 3), _coords(coords) {
    
    //Create a half sphere Z boundary element
    _boundaryElements.emplace_back(s->addTrackable<HalfSphereZBoundaryElement>
                                   (coords, radius, up,
                                   SysParams::Boundaries().BoundaryK,
                                   SysParams::Boundaries().BScreenLength));
    
}