//
//  BoundaryElementImpl.cpp
//  Cyto
//
//  Created by James Komianos on 9/15/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryElementImpl.h"

PlaneBoundaryElement::PlaneBoundaryElement(std::vector<double> coords, std::vector<double> normal) : BoundaryElement(coords) {
    
    ///set parameters
    _a = normal[0];
    _b = normal[1];
    _c = normal[2];
    
    _d = -_a * _coords[0] - _b * _coords[1] - _c * _coords[2];
}

double PlaneBoundaryElement::distance(const std::vector<double>& point) {
    
    return (_a * point[0] + _b * point[1] * _c * point[2] + _d) /
                    sqrt(pow(_a, 2) + pow(_b, 2) + pow(_c, 2));
}