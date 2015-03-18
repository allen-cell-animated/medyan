
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

#include "BoundaryElementImpl.h"

#include "MathFunctions.h"

using namespace mathfunc;

//PLANE

PlaneBoundaryElement::PlaneBoundaryElement(vector<double> coords, vector<double> normal,
                                           double repulsConst, double screenLength)
    : BoundaryElement(coords, repulsConst, screenLength) {
    
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

double PlaneBoundaryElement::stretchedDistance(const vector<double>& point,
                                               const vector<double>& force, double d) {
    
    
    return (_a * (point[0] + d*force[0]) +
            _b * (point[1] + d*force[1]) +
            _c * (point[2] + d*force[2]) + _d) /
        sqrt(pow(_a, 2) + pow(_b, 2) + pow(_c, 2));
    
}

const vector<double> PlaneBoundaryElement::normal(const vector<double>& point) {
    return vector<double>{_a, _b, _c};
}

//SPHERE

double SphereBoundaryElement::distance(const vector<double>& point) {
    
    return _radius - twoPointDistance(_coords, point);
}

double SphereBoundaryElement::stretchedDistance(const vector<double>& point,
                                                const vector<double>& force, double d) {
    
    vector<double> stretchedPoint{point[0] + d * force[0],
                                  point[1] + d * force[1],
                                  point[2] + d * force[2]};
    return _radius - twoPointDistance(_coords, stretchedPoint);
}

const vector<double> SphereBoundaryElement::normal(const vector<double>& point) {
    
    return twoPointDirection(point, _coords);
}

//CYLINDERZ

double CylindricalZBoundaryElement::distance(const vector<double>& point) {
    
    ///check z coordinate. If outside, return infinity
    if(point[2] > (_coords[2] + _height / 2) ||
       point[2] < (_coords[2] - _height / 2))
        return numeric_limits<double>::infinity();
    
    return _radius - twoPointDistance({_coords[0],_coords[1], 0},
                                      {point[0],point[1],0});
}

double CylindricalZBoundaryElement::stretchedDistance(const vector<double>& point,
                                                      const vector<double>& force,
                                                      double d) {
    
    // check z coordinate. If outside, return infinity
    if((point[2] + d * force[2]) > (_coords[2] + _height / 2) ||
       (point[2] + d * force[2]) < (_coords[2] - _height / 2))
        return numeric_limits<double>::infinity();
    
    return _radius - twoPointDistance({_coords[0],_coords[1], 0},
                     {(point[0] + d * force[0]), (point[1] + d * force[1]) ,0});
}

const vector<double> CylindricalZBoundaryElement::normal(const vector<double>& point) {
    
    return twoPointDirection({point[0],point[1], 0}, {_coords[0],_coords[1], 0});
}

//HALFSPHEREZ

double HalfSphereZBoundaryElement::distance(const vector<double>& point) {
    
    // check z coordinate. If outside, return infinity
    if(_up && (point[2] > _coords[2])) return numeric_limits<double>::infinity();
    if(!_up && (point[2] < _coords[2])) return numeric_limits<double>::infinity();
    
    return _radius - twoPointDistance(_coords, point);
}

double HalfSphereZBoundaryElement::stretchedDistance(const vector<double>& point,
                                                     const vector<double>& force,
                                                     double d) {
    
    // check z coordinate. If outside, return infinity
    if(_up && (point[2] + d * force[2] > _coords[2]))
        return numeric_limits<double>::infinity();
    if(!_up && (point[2] + d * force[2] < _coords[2]))
        return numeric_limits<double>::infinity();
    
    vector<double> stretchedPoint{point[0] + d * force[0],
                                  point[1] + d * force[1],
                                  point[2] + d * force[2]};
    
    return _radius - twoPointDistance(_coords, stretchedPoint);
    
}

const vector<double> HalfSphereZBoundaryElement::normal(const vector<double>& point) {
    
    return twoPointDirection(point, _coords);
}

