

//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BoundaryElementImpl_h
#define MEDYAN_BoundaryElementImpl_h

#include <cmath>

#include "common.h"

#include "BoundaryElement.h"

#include "MathFunctions.h"

using namespace mathfunc;

/// A plane implementation of a BoundaryElement.
class PlaneBoundaryElement : public BoundaryElement {
    
friend class BoundaryCubic;
    
private:
    /// Parameters of equation (ax + by + cz + d = 0)
    double _a, _b, _c, _d;

public:
    /// Constructor, sets parameters of equation
    PlaneBoundaryElement(vector<double> coords,
                         vector<double> normal,
                         double repulsConst,
                         double screenLength)
    
        : BoundaryElement(coords, repulsConst, screenLength) {
    	cout << "I am 0"<<endl;
        
        ///set parameters
        _a = normal[0]; _b = normal[1]; _c = normal[2];
        _d = -_a * _coords[0] - _b * _coords[1] - _c * _coords[2];
    }
    
    virtual double distance(const vector<double>& point) {
    
        return (_a * point[0] + _b * point[1] + _c * point[2] + _d) /
                sqrt(pow(_a, 2) + pow(_b, 2) + pow(_c, 2));
    }
    
    virtual double stretchedDistance(const vector<double>& point,
                                     const vector<double>& force,
                                     double d) {
        
        
        vector<double> movedPoint = {point[0] + d*force[0],
                                     point[1] + d*force[1],
                                     point[2] + d*force[2]};
        return distance(movedPoint);
        
    }
    
    virtual const vector<double> normal(const vector<double>& point) {
        
        return vector<double>{_a, _b, _c};
    }
    
    virtual void updateCoords(const vector<double> newCoords) {
        
        _coords = newCoords;
        
        //also update plane params
        _d = -_a * _coords[0] - _b * _coords[1] - _c * _coords[2];
    }
};

/// A spherical implementation of a BoundaryElement.
class SphereBoundaryElement : public BoundaryElement {
    
friend class BoundarySpherical;
    
private:
    double _radius; ///< Radius of sphere
    
public:
    /// Constructor, sets parameters of equation
    SphereBoundaryElement(vector<double> coords,
                          double radius,
                          double repulsConst,
                          double screenLength)
    
        : BoundaryElement(coords, repulsConst, screenLength),
          _radius(radius) {cout << "I am 1"<<endl;} //doesnt get called: pass test: fail
    

    ///update the radius of the spherical boundary element, added by jl135
    virtual void updateRads(const double newRads) {
    	_radius= newRads;
    	cout << "I am inside virtual void updateRads"<<endl;
    }

    virtual double distance(const vector<double>& point) {
        
        return _radius - twoPointDistance(_coords, point);
        cout << "I am inside virtual double distance"<<endl;
    }
    
    virtual double stretchedDistance(const vector<double>& point,
                                     const vector<double>& force,
                                     double d) {
        
        vector<double> movedPoint{point[0] + d * force[0],
                                  point[1] + d * force[1],
                                  point[2] + d * force[2]};
        
        return distance(movedPoint);
        
    }
    virtual const vector<double> normal(const vector<double>& point) {
        
        return twoPointDirection(point, _coords);
    }

    virtual void updateCoords(const vector<double> newCoords) {
    
        _coords = newCoords;
    }
};

/// A cylinder implementation of a BoundaryElement.
class CylindricalZBoundaryElement : public BoundaryElement {
    
friend class BoundaryCapsule;
    
private:
    double _radius; ///< Radius of cylinder
    double _height; ///< Height of cylinder
    
public:
    ///Constructor, sets parameters of equation
    CylindricalZBoundaryElement(vector<double> coords,
                                double radius,
                                double height,
                                double repulsConst,
                                double screenLength)
    
        : BoundaryElement(coords, repulsConst, screenLength),
          _radius(radius), _height(height) {}
    
    virtual double distance(const vector<double>& point) {
        
        
        ///check z coordinate. If outside, return infinity
        if(point[2] > (_coords[2] + _height / 2) ||
           point[2] < (_coords[2] - _height / 2))
            
            return numeric_limits<double>::infinity();
        
        return _radius - twoPointDistance({_coords[0],_coords[1], 0},
                                          {  point[0],  point[1], 0});
    }
    
    virtual double stretchedDistance(const vector<double>& point,
                                     const vector<double>& force,
                                     double d) {
        
        // check z coordinate. If outside, return infinity
        if((point[2] + d * force[2]) > (_coords[2] + _height / 2) ||
           (point[2] + d * force[2]) < (_coords[2] - _height / 2))
            
            return numeric_limits<double>::infinity();
        
        vector<double> movedPoint{point[0] + d * force[0],
                                  point[1] + d * force[1],
                                  point[2] + d * force[2]};
        
        return distance(movedPoint);
        
    }
    
    virtual const vector<double> normal(const vector<double>& point) {
        
        return twoPointDirection({point[0],  point[1], 0},
                               {_coords[0],_coords[1], 0});
    }
    
    virtual void updateCoords(const vector<double> newCoords) {
        
        _coords = newCoords;
    }
};

/// A half-sphere implementation of a BoundaryElement.
class HalfSphereZBoundaryElement : public BoundaryElement {
    
friend class BoundaryCapsule;
    
private:
    double _radius; ///< Radius of half sphere
    bool _up;       ///< Whether the half sphere faces up or down
    
public:
    /// Constructor, sets parameters of equation
    HalfSphereZBoundaryElement(vector<double> coords,
                               double radius, bool up,
                               double repulsConst,
                               double screenLength)
    
        : BoundaryElement(coords, repulsConst, screenLength),
          _radius(radius), _up(up){}

    virtual double distance(const vector<double>& point) {
        
        // check z coordinate. If outside, return infinity
        if((_up && (point[2] > _coords[2])) ||
          (!_up && (point[2] < _coords[2])))
            return numeric_limits<double>::infinity();
        
        return _radius - twoPointDistance(_coords, point);
    }
    
    virtual double stretchedDistance(const vector<double>& point,
                                     const vector<double>& force,
                                     double d) {
        
        vector<double> movedPoint{point[0] + d * force[0],
                                  point[1] + d * force[1],
                                  point[2] + d * force[2]};
        
        // check z coordinate. If outside, return infinity
        if((_up && (movedPoint[2] > _coords[2])) ||
          (!_up && (movedPoint[2] < _coords[2])))
            return numeric_limits<double>::infinity();
        
        return distance(movedPoint);
        
    }
    
    virtual const vector<double> normal(const vector<double>& point) {
        
        return twoPointDirection(point, _coords);
    }
    
    virtual void updateCoords(const vector<double> newCoords) {
        
        _coords = newCoords;
    }
};

#endif
