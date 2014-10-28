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

Plane::Plane(std::vector<double> coords, std::vector<double> normal ) :
    _coords(coords), _normal(normal), BoundarySurface(3) {
    
    ///Create a plane boundary element (CHANGE REPULSION CONSTANT)
    _boundaryElements.emplace_back(BoundaryElementDB::Instance(BEDBKey())->CreatePlaneBoundaryElement(coords, normal,
            SystemParameters::Boundaries().boundaryK, SystemParameters::Boundaries().screenLength));
}

Sphere::Sphere(std::vector<double> coords, double radius) : _coords(coords), _radius(radius), BoundarySurface(3){
    
    ///Create a sphere boundary element
    _boundaryElements.emplace_back(BoundaryElementDB::Instance(BEDBKey())->CreateSphereBoundaryElement(coords, radius,
        SystemParameters::Boundaries().boundaryK, SystemParameters::Boundaries().screenLength));
    
}