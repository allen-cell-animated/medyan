//
//  BoundarySurfaceImpl.cpp
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundarySurfaceImpl.h"
#include "MathFunctions.h"

using namespace mathfunc;

Plane::Plane(std::vector<double> coords, std::vector<double> normal ) :
    _coords(coords), _normal(normal), BoundarySurface(3) {
    
    ///Create a plane boundary element
        BoundaryElementDB::Instance(BEDBKey())->CreatePlaneBoundaryElement(coords, normal);

}