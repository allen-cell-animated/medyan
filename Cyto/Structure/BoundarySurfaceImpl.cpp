//
//  BoundarySurfaceImpl.cpp
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundarySurfaceImpl.h"

#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

Plane::Plane(vector<double> coords, vector<double> normal ) :
    _coords(coords), _normal(normal), BoundarySurface(3) {
    
    ///Create a plane boundary element (CHANGE REPULSION CONSTANT)
    _boundaryElements.emplace_back(BoundaryElementDB::instance(BEDBKey())->createPlaneBoundaryElement(coords, normal,
            SystemParameters::Boundaries().boundaryK, SystemParameters::Boundaries().screenLength));
}

Sphere::Sphere(vector<double> coords, double radius) : _coords(coords), _radius(radius), BoundarySurface(3){
    
    ///Create a sphere boundary element
    _boundaryElements.emplace_back(BoundaryElementDB::instance(BEDBKey())->createSphereBoundaryElement(coords, radius,
        SystemParameters::Boundaries().boundaryK, SystemParameters::Boundaries().screenLength));
    
}

CylinderZ::CylinderZ(vector<double> coords, double radius, double height)
                     : _coords(coords), _radius(radius), _height(height), BoundarySurface(3){
    
    ///Create a cylindricalZ boundary element
    _boundaryElements.emplace_back(BoundaryElementDB::instance(BEDBKey())->createCylindricalZBoundaryElement(coords, radius, height,
                                   SystemParameters::Boundaries().boundaryK, SystemParameters::Boundaries().screenLength));
    
    
}

HalfSphereZ::HalfSphereZ(vector<double> coords, double radius, bool up)
                          : _coords(coords), _radius(radius), _up(up), BoundarySurface(3){
    
    ///Create a cylindricalZ boundary element
    _boundaryElements.emplace_back(BoundaryElementDB::instance(BEDBKey())->createHalfSphereZBoundaryElement(coords, radius, up,
                                   SystemParameters::Boundaries().boundaryK, SystemParameters::Boundaries().screenLength));
    
}