

//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
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
// IF USING CUDA VERSION, MAKE SURE TO ADD RELEVANT DISTANCE, STRETCHED DISTANCE AND NORMAL FUNCTIONS IN MATHFUNCTIONS.H
/// A plane implementation of a BoundaryElement.
class PlaneBoundaryElement : public BoundaryElement {
    
friend class BoundaryCubic;
    
private:
    /// Parameters of equation (ax + by + cz + d = 0)
    floatingpoint _a, _b, _c, _d;

public:
    /// Constructor, sets parameters of equation
    PlaneBoundaryElement(vector<floatingpoint> coords,
                         vector<floatingpoint> normal,
                         floatingpoint repulsConst,
                         floatingpoint screenLength)
    
        : BoundaryElement(coords, repulsConst, screenLength) {
        
        ///set parameters
        _a = normal[0]; _b = normal[1]; _c = normal[2];
        _d = -_a * _coords[0] - _b * _coords[1] - _c * _coords[2];
    }
    
    virtual floatingpoint distance(const vector<floatingpoint>& point) {
    
        return (_a * point[0] + _b * point[1] + _c * point[2] + _d) /
                sqrt(pow(_a, 2) + pow(_b, 2) + pow(_c, 2));
    }

    //Qin,  lower distance is the z axis
    virtual floatingpoint lowerdistance(const vector<floatingpoint>& point) {
        return point[2];
    }
    //Qin, side distance is either x or y axis
    virtual floatingpoint sidedistance(const vector<floatingpoint>& point) {
        if(point[0] > point[1]) {
            return point[1];
        }
        else
            return point[0];
    }

    virtual floatingpoint distance(floatingpoint const *point) {
        return (_a * point[0] + _b * point[1] + _c * point[2] + _d) /
        sqrt(pow(_a, 2) + pow(_b, 2) + pow(_c, 2));
    }

    virtual floatingpoint stretchedDistance(const vector<floatingpoint>& point,
                                     const vector<floatingpoint>& force,
                                     floatingpoint d) {
        
        
        vector<floatingpoint> movedPoint = {point[0] + d*force[0],
                                     point[1] + d*force[1],
                                     point[2] + d*force[2]};
        return distance(movedPoint);
        
    }
    virtual floatingpoint stretchedDistance(floatingpoint const *point,
                                     floatingpoint const *force,
                                     floatingpoint d) {


        vector<floatingpoint> movedPoint = {point[0] + d*force[0],
                                     point[1] + d*force[1],
                                     point[2] + d*force[2]};
        return distance(movedPoint);

    }

    virtual const vector<floatingpoint> normal(const vector<floatingpoint>& point) {
        
        return vector<floatingpoint>{_a, _b, _c};
    }
    virtual const vector<floatingpoint> normal(floatingpoint const *point) {

        return vector<floatingpoint>{_a, _b, _c};
    }

    virtual const void elementeqn(floatingpoint* var){
        var[0] = _a; var[1] = _b; var[2] = _c; var[3] = _d;
    }

    virtual void updateCoords(const vector<floatingpoint> newCoords) {
        
        _coords = newCoords;
        
        //also update plane params
        _d = -_a * _coords[0] - _b * _coords[1] - _c * _coords[2];
    }
};

/// A spherical implementation of a BoundaryElement.
class SphereBoundaryElement : public BoundaryElement {
    
friend class BoundarySpherical;
    
private:
    floatingpoint _radius; ///< Radius of sphere
    
public:
    /// Constructor, sets parameters of equation
    SphereBoundaryElement(vector<floatingpoint> coords,
                          floatingpoint radius,
                          floatingpoint repulsConst,
                          floatingpoint screenLength)
    
        : BoundaryElement(coords, repulsConst, screenLength),
          _radius(radius) {}
    
    virtual floatingpoint distance(const vector<floatingpoint>& point) {
        
        return _radius - twoPointDistance(_coords, point);
    }
    virtual floatingpoint distance(floatingpoint const *point) {

        return _radius - twoPointDistance(_coords, point);
    }

    virtual const void elementeqn(floatingpoint* var){
        var[0] = _radius;
    }
    //Qin, the same as distance
    virtual floatingpoint lowerdistance(const vector<floatingpoint>& point) {

        return _radius - twoPointDistance(_coords, point);
    }
    virtual floatingpoint sidedistance(const vector<floatingpoint>& point) {

        return _radius - twoPointDistance(_coords, point);
    }

    virtual floatingpoint stretchedDistance(const vector<floatingpoint>& point,
                                     const vector<floatingpoint>& force,
                                     floatingpoint d) {
        
        vector<floatingpoint> movedPoint{point[0] + d * force[0],
                                  point[1] + d * force[1],
                                  point[2] + d * force[2]};
        
        return distance(movedPoint);
        
    }
    virtual floatingpoint stretchedDistance(floatingpoint const *point,
                                     floatingpoint const *force,
                                     floatingpoint d) {

        vector<floatingpoint> movedPoint{point[0] + d * force[0],
                                  point[1] + d * force[1],
                                  point[2] + d * force[2]};

        return distance(movedPoint);

    }

    virtual const vector<floatingpoint> normal(const vector<floatingpoint>& point) {
        
        return twoPointDirection(point, _coords);
    }
    virtual const vector<floatingpoint> normal(floatingpoint const *point) {

        return twoPointDirection(point, _coords);
    }

    virtual void updateCoords(const vector<floatingpoint> newCoords) {
    
        _coords = newCoords;
    }
};

/// A cylinder implementation of a BoundaryElement.
class CylindricalZBoundaryElement : public BoundaryElement {
    
friend class BoundaryCapsule;
    
private:
    floatingpoint _radius; ///< Radius of cylinder
    floatingpoint _height; ///< Height of cylinder
    
public:
    ///Constructor, sets parameters of equation
    CylindricalZBoundaryElement(vector<floatingpoint> coords,
                                floatingpoint radius,
                                floatingpoint height,
                                floatingpoint repulsConst,
                                floatingpoint screenLength)
    
        : BoundaryElement(coords, repulsConst, screenLength),
          _radius(radius), _height(height) {}
    
    virtual floatingpoint distance(const vector<floatingpoint>& point) {
        
        
        ///check z coordinate. If outside, return infinity
        if(point[2] > (_coords[2] + _height / 2) ||
           point[2] < (_coords[2] - _height / 2))
            
            return numeric_limits<floatingpoint>::infinity();


	        return _radius - twoPointDistance({_coords[0],_coords[1], 0},
	                               {  point[0],  point[1], 0});
    }
    //TODO
    virtual floatingpoint distance(floatingpoint const *point) {
        cout<<"Distance not implemented for Boundary Element of CylinderZ type. "
                "Exiting"<<endl;
        exit(EXIT_FAILURE);
        return 0.0; }

    virtual const void elementeqn(floatingpoint* var){
        var[0] = _radius; var[1] = _height;
    }

    //Qin, find the distance for the lower boundary
    virtual floatingpoint lowerdistance(const vector<floatingpoint>& point) {


        ///check z coordinate. If outside, return infinity
        if(point[2] > (_coords[2] + _height / 2) ||
           point[2] < (_coords[2] - _height / 2)) {

            return numeric_limits<floatingpoint>::infinity();
        }
        else {
            return point[2];
        }
    }


    //Qin, find the distance for the side boundary
    virtual floatingpoint sidedistance(const vector<floatingpoint>& point) {


        ///check z coordinate. If outside, return infinity
        if(point[2] > (_coords[2] + _height / 2) ||
           point[2] < (_coords[2] - _height / 2)) {

            return numeric_limits<floatingpoint>::infinity();
        }
        else {
            auto dxy = _radius - twoPointDistance({_coords[0],_coords[1], 0},
                                                  {  point[0],  point[1], 0});
            return dxy;
        }
    }


    virtual floatingpoint stretchedDistance(const vector<floatingpoint>& point,
                                     const vector<floatingpoint>& force,
                                     floatingpoint d) {
        
        // check z coordinate. If outside, return infinity
        if((point[2] + d * force[2]) > (_coords[2] + _height / 2) ||
           (point[2] + d * force[2]) < (_coords[2] - _height / 2))
            
            return numeric_limits<floatingpoint>::infinity();
        
        vector<floatingpoint> movedPoint{point[0] + d * force[0],
                                  point[1] + d * force[1],
                                  point[2] + d * force[2]};
        
        return distance(movedPoint);
        
    }
    ///TODO
    virtual floatingpoint stretchedDistance(floatingpoint const *point,
                                     floatingpoint const *force,
                                     floatingpoint d) {return 0.0;}

    virtual const vector<floatingpoint> normal(const vector<floatingpoint>& point) {

      //Qin
      auto dxy = _radius - twoPointDistance({_coords[0],_coords[1], 0},
					    {  point[0],  point[1], 0});

      floatingpoint dzz = point[2];
      if((_coords[2] * 2 - point[2]) < dzz) {
	dzz = _coords[2]*2 - point[2];
      }
      else {
	dzz = point[2];
      }

      if(dxy > dzz) {
	return twoPointDirection({0,  0, point[2]},
				 {0,0, _coords[2]});
      }
      else {
	return twoPointDirection({point[0],  point[1], 0},
				 {_coords[0],_coords[1], 0});
      }
      //return twoPointDirection({point[0],  point[1], 0},
      //                         {_coords[0],_coords[1], 0});
    }
    
    ///TODO
    virtual const vector<floatingpoint> normal(floatingpoint const *point) {return vector<floatingpoint>{};};

    virtual void updateCoords(const vector<floatingpoint> newCoords) {
        
        _coords = newCoords;
    }
};

/// A half-sphere implementation of a BoundaryElement.
class HalfSphereZBoundaryElement : public BoundaryElement {
    
friend class BoundaryCapsule;
    
private:
    floatingpoint _radius; ///< Radius of half sphere
    bool _up;       ///< Whether the half sphere faces up or down
    
public:
    /// Constructor, sets parameters of equation
    HalfSphereZBoundaryElement(vector<floatingpoint> coords,
                               floatingpoint radius, bool up,
                               floatingpoint repulsConst,
                               floatingpoint screenLength)
    
        : BoundaryElement(coords, repulsConst, screenLength),
          _radius(radius), _up(up){}

    virtual floatingpoint distance(const vector<floatingpoint>& point) {
        
        // check z coordinate. If outside, return infinity
        if((_up && (point[2] > _coords[2])) ||
          (!_up && (point[2] < _coords[2])))
            return numeric_limits<floatingpoint>::infinity();
        
        return _radius - twoPointDistance(_coords, point);
    }
    
    ///TODO
    virtual floatingpoint distance(floatingpoint const *point) {
        cout<<"Distance not implemented for Boundary Element of HalfSphereZ type. "
                "Exiting"<<endl;
        exit(EXIT_FAILURE);
        return 0.0; }

    //Qin, the same as distance
    virtual floatingpoint lowerdistance(const vector<floatingpoint>& point) {

        // check z coordinate. If outside, return infinity
        if((_up && (point[2] > _coords[2])) ||
           (!_up && (point[2] < _coords[2])))
            return numeric_limits<floatingpoint>::infinity();

        return _radius - twoPointDistance(_coords, point);
    }
    virtual floatingpoint sidedistance(const vector<floatingpoint>& point) {

        // check z coordinate. If outside, return infinity
        if((_up && (point[2] > _coords[2])) ||
           (!_up && (point[2] < _coords[2])))
            return numeric_limits<floatingpoint>::infinity();

        return _radius - twoPointDistance(_coords, point);
    }

    virtual floatingpoint stretchedDistance(const vector<floatingpoint>& point,
                                     const vector<floatingpoint>& force,
                                     floatingpoint d) {
        
        vector<floatingpoint> movedPoint{point[0] + d * force[0],
                                  point[1] + d * force[1],
                                  point[2] + d * force[2]};
        
        // check z coordinate. If outside, return infinity
        if((_up && (movedPoint[2] > _coords[2])) ||
          (!_up && (movedPoint[2] < _coords[2])))
            return numeric_limits<floatingpoint>::infinity();
        
        return distance(movedPoint);
        
    }

    virtual const void elementeqn(floatingpoint* var){
        var[0] = _radius;
    }

    ///TODO
    virtual floatingpoint stretchedDistance(floatingpoint const *point,
                                     floatingpoint const *force,
                                     floatingpoint d) {return 0.0;}

    virtual const vector<floatingpoint> normal(const vector<floatingpoint>& point) {
        
        return twoPointDirection(point, _coords);
    }
    
    ///TODO
    virtual const vector<floatingpoint> normal(floatingpoint const *point) {return vector<floatingpoint>{};};

    virtual void updateCoords(const vector<floatingpoint> newCoords) {
        
        _coords = newCoords;
    }
};

//----------------------------------------------------------
/// A cylinder implementation of a BoundaryElement.
class CylindricalXYZBoundaryElement : public BoundaryElement {

    friend class BoundaryCylinder;

private:
    floatingpoint _radius; ///< Radius of cylinder
    floatingpoint _height; ///< Height of cylinder

public:
    ///Constructor, sets parameters of equation
    CylindricalXYZBoundaryElement(vector<floatingpoint> coords,
                                  floatingpoint radius,
                                  floatingpoint height,
                                  floatingpoint repulsConst,
                                  floatingpoint screenLength)

            : BoundaryElement(coords, repulsConst, screenLength),
              _radius(radius), _height(height) {}

    virtual floatingpoint distance(const vector<floatingpoint>& point) {


        ///check z coordinate. If outside, return infinity
        if(point[2] > (_coords[2] + _height / 2) ||
           point[2] < (_coords[2] - _height / 2))

            return numeric_limits<floatingpoint>::infinity();

        //Qin
        auto dxy = _radius - twoPointDistance({_coords[0],_coords[1], 0},
                                              {  point[0],  point[1], 0});

        floatingpoint dzz = point[2];
        if((_coords[2] * 2 - point[2]) < dzz) {
            dzz = _coords[2]*2 - point[2];
        }
        else {
            dzz = point[2];
        }


        if(dxy > dzz) {
            return dzz;
        }
        else {
            return dxy;
        }

        //        return _radius - twoPointDistance({_coords[0],_coords[1], 0},
        //                                {  point[0],  point[1], 0});
    }

    //Qin, find the distance for the lower boundary
    virtual floatingpoint lowerdistance(const vector<floatingpoint>& point) {


        ///check z coordinate. If outside, return infinity
        if(point[2] > (_coords[2] + _height / 2) ||
           point[2] < (_coords[2] - _height / 2)) {

            return numeric_limits<floatingpoint>::infinity();
        }
        else {
            return point[2];
        }
    }

    //Qin, find the distance for the side boundary
    virtual floatingpoint sidedistance(const vector<floatingpoint>& point) {


        ///check z coordinate. If outside, return infinity
        if(point[2] > (_coords[2] + _height / 2) ||
           point[2] < (_coords[2] - _height / 2)) {

            return numeric_limits<floatingpoint>::infinity();
        }
        else {
            auto dxy = _radius - twoPointDistance({_coords[0],_coords[1], 0},
                                                  {  point[0],  point[1], 0});
            return dxy;
        }
    }


    virtual floatingpoint stretchedDistance(const vector<floatingpoint>& point,
                                     const vector<floatingpoint>& force,
                                     floatingpoint d) {

        // check z coordinate. If outside, return infinity
        if((point[2] + d * force[2]) > (_coords[2] + _height / 2) ||
           (point[2] + d * force[2]) < (_coords[2] - _height / 2))

            return numeric_limits<floatingpoint>::infinity();

        vector<floatingpoint> movedPoint{point[0] + d * force[0],
                                  point[1] + d * force[1],
                                  point[2] + d * force[2]};

        return distance(movedPoint);

    }

    virtual const vector<floatingpoint> normal(const vector<floatingpoint>& point) {

        //Qin
        auto dxy = _radius - twoPointDistance({_coords[0],_coords[1], 0},
                                              {  point[0],  point[1], 0});

        floatingpoint dzz = point[2];
        if((_coords[2] * 2 - point[2]) < dzz) {
            dzz = _coords[2]*2 - point[2];
        }
        else {
            dzz = point[2];
        }

        if(dxy > dzz) {
            return twoPointDirection({0,  0, point[2]},
                                     {0,0, _coords[2]});
        }
        else {
            return twoPointDirection({point[0],  point[1], 0},
                                     {_coords[0],_coords[1], 0});
        }
        //return twoPointDirection({point[0],  point[1], 0},
        //                         {_coords[0],_coords[1], 0});
    }

    virtual void updateCoords(const vector<floatingpoint> newCoords) {

        _coords = newCoords;
    }
    //TODO needs implementation @{
    virtual floatingpoint distance(floatingpoint const *point) {
        cout<<"Distance not implemented for Boundary Element of CylinderXYZ type. "
                "Exiting"<<endl;
        exit(EXIT_FAILURE);
        return 0.0;
    }
    virtual floatingpoint stretchedDistance(floatingpoint const *point,
                                     floatingpoint const *force, floatingpoint d) {return 0.0;};
    virtual const vector<floatingpoint> normal(const floatingpoint *point) {vector<floatingpoint> a; return a;};
    virtual const void elementeqn(floatingpoint* var){cout<<" element eqn not implemented for "
                "CylindricalXYZ. Exiting."<<endl;
    exit(EXIT_FAILURE);} ;

    //@}


};
#endif
