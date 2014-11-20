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

PlaneBoundaryElement::PlaneBoundaryElement(vector<double> coords, vector<double> normal,
                                           double repulsConst, double sceenLength)
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

SphereBoundaryElement::SphereBoundaryElement(vector<double> coords, double radius,
                                             double repulsConst, double sceenLength)
                    : BoundaryElement(coords), _k_rep(repulsConst), _r0(sceenLength), _radius(radius) {}

double SphereBoundaryElement::distance(const vector<double>& point) {
    
    return _radius - TwoPointDistance(_coords, point);
}

double SphereBoundaryElement::stretchedDistance(const vector<double>& point, const vector<double>& force, double d) {
    
    vector<double> stretchedPoint{point[0] + d * force[0], point[1] + d * force[1], point[2] + d * force[2]};
    return _radius - TwoPointDistance(_coords, stretchedPoint);
    
}

const vector<double> SphereBoundaryElement::normal(const vector<double>& point) {
    
    return TwoPointDirection(point, _coords);
}

double SphereBoundaryElement::getRepulsionConst(){ return _k_rep; }
double SphereBoundaryElement::getScreeningLength(){ return _r0; }


///CYLINDRICAL Z BOUNDARY ELEMENT

CylindricalZBoundaryElement::CylindricalZBoundaryElement(vector<double> coords, double radius, double height,
                                                         double repulsConst, double sceenLength)
        : BoundaryElement(coords), _k_rep(repulsConst), _r0(sceenLength), _radius(radius), _height(height) {}

double CylindricalZBoundaryElement::distance(const vector<double>& point) {
    
    ///check z coordinate. If outside, return infinity
    if(point[2] > (_coords[2] + _height / 2) || point[2] < (_coords[2] - _height / 2))
        return numeric_limits<double>::infinity();
    
    return _radius - TwoPointDistance({_coords[0],_coords[1], 0}, {point[0],point[1],0});
}

double CylindricalZBoundaryElement::stretchedDistance(const vector<double>& point, const vector<double>& force, double d) {
    
    ///check z coordinate. If outside, return infinity
    if((point[2] + d * force[2]) > (_coords[2] + _height / 2) ||
       (point[2] + d * force[2]) < (_coords[2] - _height / 2)) return numeric_limits<double>::infinity();
    
    return _radius - TwoPointDistance({_coords[0],_coords[1], 0},
                     {(point[0] + d * force[0]), (point[1] + d * force[1]) ,0});
}

const vector<double> CylindricalZBoundaryElement::normal(const vector<double>& point) {
    
    return TwoPointDirection({point[0],point[1], 0}, {_coords[0],_coords[1], 0});
}

double CylindricalZBoundaryElement::getRepulsionConst(){ return _k_rep; }
double CylindricalZBoundaryElement::getScreeningLength(){ return _r0; }


///HALF SPHERE Z BOUNDARY ELEMENT

HalfSphereZBoundaryElement::HalfSphereZBoundaryElement(vector<double> coords, double radius,
                                                       bool up, double repulsConst, double sceenLength)
            : BoundaryElement(coords), _k_rep(repulsConst), _r0(sceenLength), _radius(radius), _up(up) {}

double HalfSphereZBoundaryElement::distance(const vector<double>& point) {
    
    ///check z coordinate. If outside, return infinity
    if(_up && (point[2] > _coords[2])) return numeric_limits<double>::infinity();
    if(!_up && (point[2] < _coords[2])) return numeric_limits<double>::infinity();
    
    return _radius - TwoPointDistance(_coords, point);
}

double HalfSphereZBoundaryElement::stretchedDistance(const vector<double>& point, const vector<double>& force, double d) {
    
    ///check z coordinate. If outside, return infinity
    if(_up && (point[2] + d * force[2] > _coords[2])) return numeric_limits<double>::infinity();
    if(!_up && (point[2] + d * force[2] < _coords[2])) return numeric_limits<double>::infinity();
    
    vector<double> stretchedPoint{point[0] + d * force[0], point[1] + d * force[1], point[2] + d * force[2]};
    return _radius - TwoPointDistance(_coords, stretchedPoint);
    
}

const vector<double> HalfSphereZBoundaryElement::normal(const vector<double>& point) {
    
    return TwoPointDirection(point, _coords);
}

double HalfSphereZBoundaryElement::getRepulsionConst(){ return _k_rep; }
double HalfSphereZBoundaryElement::getScreeningLength(){ return _r0; }




