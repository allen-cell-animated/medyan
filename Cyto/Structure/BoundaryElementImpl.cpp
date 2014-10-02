//
//  BoundaryElementImpl.cpp
//  Cyto
//
//  Created by James Komianos on 9/15/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryElementImpl.h"

PlaneBoundaryElement::PlaneBoundaryElement(std::vector<double> coords, std::vector<double> normal, double repulsConst, double sceenLength) : BoundaryElement(coords, normal), _k_rep(repulsConst), _r0(sceenLength) {
    
    ///set parameters
    _a = normal[0];
    _b = normal[1];
    _c = normal[2];
    
    _d = -_a * _coords[0] - _b * _coords[1] - _c * _coords[2];
}

double PlaneBoundaryElement::distance(const std::vector<double>& point) {
    
    return (_a * point[0] + _b * point[1] + _c * point[2] + _d) /
                            sqrt(pow(_a, 2) + pow(_b, 2) + pow(_c, 2));
}

double PlaneBoundaryElement::stretchedDistance(const std::vector<double>& point, const std::vector<double>& force, double d) {
    
    
    return (_a * (point[0] + d*force[0]) + _b * (point[1] + d*force[1]) + _c * (point[2] + d*force[2]) + _d) /
        sqrt(pow(_a, 2) + pow(_b, 2) + pow(_c, 2));
    
}

double PlaneBoundaryElement::getRepulsionConst(){
    
    return _k_rep;
}

double PlaneBoundaryElement::getScreeningLength(){
    
    return _r0;
}