//
//  BoundaryElementImpl.cpp
//  Cyto
//
//  Created by James Komianos on 9/15/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryElementImpl.h"

#include "MathFunctions.h"

using namespace mathfunc;

///PLANE BOUNDARY ELEMENT

PlaneBoundaryElement::PlaneBoundaryElement(vector<double> coords, vector<double> normal, double repulsConst, double sceenLength)
                                                                : BoundaryElement(coords), _k_rep(repulsConst), _r0(sceenLength) {
    
    ///set parameters
    _a = normal[0];
    _b = normal[1];
    _c = normal[2];
    
    _d = -_a * _coords[0] - _b * _coords[1] - _c * _coords[2];
}

double PlaneBoundaryElement::distance(const vector<double>& point) {
    
    return (_a * point[0] + _b * point[1] + _c * point[2] + _d) /
                            sqrt(pow(_a, 2) + pow(_b, 2) + pow(_c, 2));
}

double PlaneBoundaryElement::stretchedDistance(const vector<double>& point, const vector<double>& force, double d) {
    
    
    return (_a * (point[0] + d*force[0]) + _b * (point[1] + d*force[1]) + _c * (point[2] + d*force[2]) + _d) /
        sqrt(pow(_a, 2) + pow(_b, 2) + pow(_c, 2));
    
}

const vector<double> PlaneBoundaryElement::normal(const vector<double>& point) {
    return vector<double>{_a, _b, _c};
}

double PlaneBoundaryElement::getRepulsionConst(){ return _k_rep; }

double PlaneBoundaryElement::getScreeningLength(){ return _r0; }


///SPHERE BOUNDARY ELEMENT

SphereBoundaryElement::SphereBoundaryElement(vector<double> coords, double radius, double repulsConst, double sceenLength)
                                                    : BoundaryElement(coords), _k_rep(repulsConst), _r0(sceenLength), _radius(radius) {}

double SphereBoundaryElement::distance(const vector<double>& point) {
    
    return - (TwoPointDistance(_coords, point) - _radius);
}

double SphereBoundaryElement::stretchedDistance(const vector<double>& point, const vector<double>& force, double d) {
    
    vector<double> stretchedPoint{point[0] + d * force[0], point[1] + d * force[1], point[2] + d * force[2]};
    return - (TwoPointDistance(_coords, stretchedPoint) - _radius);
    
}

const vector<double> SphereBoundaryElement::normal(const vector<double>& point) {
    
    return TwoPointDirection(point, _coords);
}

double SphereBoundaryElement::getRepulsionConst(){ return _k_rep; }

double SphereBoundaryElement::getScreeningLength(){ return _r0; }



