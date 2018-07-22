

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
    double _a, _b, _c, _d;

public:
    /// Constructor, sets parameters of equation
    PlaneBoundaryElement(vector<double> coords,
                         vector<double> normal,
                         double repulsConst,
                         double screenLength)
    
        : BoundaryElement(coords, repulsConst, screenLength) {
        
        ///set parameters
        _a = normal[0]; _b = normal[1]; _c = normal[2];
        _d = -_a * _coords[0] - _b * _coords[1] - _c * _coords[2];
    }
    
    virtual double distance(const vector<double>& point) {
    
        return (_a * point[0] + _b * point[1] + _c * point[2] + _d) /
                sqrt(pow(_a, 2) + pow(_b, 2) + pow(_c, 2));
    }

    //Qin,  lower distance is the z axis
    virtual double lowerdistance(const vector<double>& point) {
        return point[2];
    }
    
    //Qin, side distance is either x or y axis
    virtual double sidedistance(const vector<double>& point) {
        if(point[0] > point[1]) {
            return point[1];
        }
        else
            return point[0];
    }

    virtual double distance(double const *point) {
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
    virtual double stretchedDistance(double const *point,
                                     double const *force,
                                     double d) {


        vector<double> movedPoint = {point[0] + d*force[0],
                                     point[1] + d*force[1],
                                     point[2] + d*force[2]};
        return distance(movedPoint);

    }

    virtual const vector<double> normal(const vector<double>& point) {
        
        return vector<double>{_a, _b, _c};
    }
    virtual const vector<double> normal(double const *point) {

        return vector<double>{_a, _b, _c};
    }

    virtual const void elementeqn(double* var){
        var[0] = _a; var[1] = _b; var[2] = _c; var[3] = _d;
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
          _radius(radius) {}
    
    virtual double distance(const vector<double>& point) {
        
        return _radius - twoPointDistance(_coords, point);
    }
    virtual double distance(double const *point) {

        return _radius - twoPointDistance(_coords, point);
    }

    virtual const void elementeqn(double* var){
        var[0] = _radius;
    }
    //Qin, the same as distance
    virtual double lowerdistance(const vector<double>& point) {

        return _radius - twoPointDistance(_coords, point);
    }
    virtual double sidedistance(const vector<double>& point) {

        return _radius - twoPointDistance(_coords, point);
    }

    virtual double stretchedDistance(const vector<double>& point,
                                     const vector<double>& force,
                                     double d) {
        
        vector<double> movedPoint{point[0] + d * force[0],
                                  point[1] + d * force[1],
                                  point[2] + d * force[2]};
        
        return distance(movedPoint);
        
    }
    virtual double stretchedDistance(double const *point,
                                     double const *force,
                                     double d) {

        vector<double> movedPoint{point[0] + d * force[0],
                                  point[1] + d * force[1],
                                  point[2] + d * force[2]};

        return distance(movedPoint);

    }

    virtual const vector<double> normal(const vector<double>& point) {
        
        return twoPointDirection(point, _coords);
    }
    virtual const vector<double> normal(double const *point) {

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
    //TODO
    virtual double distance(double const *point) {
        cout<<"Distance not implemented for Boundary Element of CylinderZ type. "
                "Exiting"<<endl;
        exit(EXIT_FAILURE);
        return 0.0; }

    virtual const void elementeqn(double* var){
        var[0] = _radius; var[1] = _height;
    }

    //Qin, find the distance for the lower boundary
    virtual double lowerdistance(const vector<double>& point) {


        ///check z coordinate. If outside, return infinity
        if(point[2] > (_coords[2] + _height / 2) ||
           point[2] < (_coords[2] - _height / 2)) {

            return numeric_limits<double>::infinity();
        }
        else {
            return point[2];
        }
    }


    //Qin, find the distance for the side boundary
    virtual double sidedistance(const vector<double>& point) {


        ///check z coordinate. If outside, return infinity
        if(point[2] > (_coords[2] + _height / 2) ||
           point[2] < (_coords[2] - _height / 2)) {

            return numeric_limits<double>::infinity();
        }
        else {
            auto dxy = _radius - twoPointDistance({_coords[0],_coords[1], 0},
                                                  {  point[0],  point[1], 0});
            return dxy;
        }
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
    ///TODO
    virtual double stretchedDistance(double const *point,
                                     double const *force,
                                     double d) {return 0.0;}

    virtual const vector<double> normal(const vector<double>& point) {

        //Qin
        auto dxy = _radius - twoPointDistance({_coords[0],_coords[1], 0},
                                              {  point[0],  point[1], 0});
        
        double dzz = point[2];
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
    virtual const vector<double> normal(double const *point) {return vector<double>{};};

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
    
    ///TODO
    virtual double distance(double const *point) {
        cout<<"Distance not implemented for Boundary Element of HalfSphereZ type. "
                "Exiting"<<endl;
        exit(EXIT_FAILURE);
        return 0.0; }

    //Qin, the same as distance
    virtual double lowerdistance(const vector<double>& point) {

        // check z coordinate. If outside, return infinity
        if((_up && (point[2] > _coords[2])) ||
           (!_up && (point[2] < _coords[2])))
            return numeric_limits<double>::infinity();

        return _radius - twoPointDistance(_coords, point);
    }
    virtual double sidedistance(const vector<double>& point) {

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

    virtual const void elementeqn(double* var){
        var[0] = _radius;
    }

    ///TODO
    virtual double stretchedDistance(double const *point,
                                     double const *force,
                                     double d) {return 0.0;}

    virtual const vector<double> normal(const vector<double>& point) {
        
        return twoPointDirection(point, _coords);
    }
    
    ///TODO
    virtual const vector<double> normal(double const *point) {return vector<double>{};};

    virtual void updateCoords(const vector<double> newCoords) {
        
        _coords = newCoords;
    }
};

//----------------------------------------------------------
/// A cylinder implementation of a BoundaryElement.
class CylindricalXYZBoundaryElement : public BoundaryElement {

    friend class BoundaryCylinder;

private:
    double _radius; ///< Radius of cylinder
    double _height; ///< Height of cylinder

public:
    ///Constructor, sets parameters of equation
    CylindricalXYZBoundaryElement(vector<double> coords,
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

        //Qin
        auto dxy = _radius - twoPointDistance({_coords[0],_coords[1], 0},
                                              {  point[0],  point[1], 0});

        double dzz = point[2];
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
    virtual double lowerdistance(const vector<double>& point) {


        ///check z coordinate. If outside, return infinity
        if(point[2] > (_coords[2] + _height / 2) ||
           point[2] < (_coords[2] - _height / 2)) {

            return numeric_limits<double>::infinity();
        }
        else {
            return point[2];
        }
    }

    //Qin, find the distance for the side boundary
    virtual double sidedistance(const vector<double>& point) {


        ///check z coordinate. If outside, return infinity
        if(point[2] > (_coords[2] + _height / 2) ||
           point[2] < (_coords[2] - _height / 2)) {

            return numeric_limits<double>::infinity();
        }
        else {
            auto dxy = _radius - twoPointDistance({_coords[0],_coords[1], 0},
                                                  {  point[0],  point[1], 0});
            return dxy;
        }
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

        //Qin
        auto dxy = _radius - twoPointDistance({_coords[0],_coords[1], 0},
                                              {  point[0],  point[1], 0});

        double dzz = point[2];
        if((_coords[2] * 2 - point[2]) < dzz) {
            dzz = _coords[2]*2 - point[2];
        }
        else {
            dzz = point[2];
        }

        if(dxy > dzz) {
            
            // when the Z coordinate is located at geometry center
            if(areEqual(point[2],_coords[2])) {
                return vector<double> {0.0, 0.0, 0.0};
            }
            else {
                return twoPointDirection({0,  0, point[2]},
                                         {0,0, _coords[2]});;
            }
            
        }
        else {
 
            // when the Z coordinate is located at geometry center
            if(areEqual(point[0],_coords[0]) && areEqual(point[1],_coords[1])) {
                return vector<double> {0.0, 0.0, 0.0};
            }
            else {
                return twoPointDirection({point[0],  point[1], 0},
                                         {_coords[0],_coords[1], 0});;
            }
        }
        //return twoPointDirection({point[0],  point[1], 0},
        //                         {_coords[0],_coords[1], 0});
    }

    virtual void updateCoords(const vector<double> newCoords) {

        _coords = newCoords;
    }

    virtual double distance(double const *point) {
        ///check z coordinate. If outside, return infinity
        if(point[2] > (_coords[2] + _height / 2) ||
           point[2] < (_coords[2] - _height / 2))
        
        return numeric_limits<double>::infinity();
        
        //Qin
        auto dxy = _radius - twoPointDistance({_coords[0],_coords[1], 0},
                                              {  point[0],  point[1], 0});
        
        double dzz = point[2];
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
    }
    virtual double stretchedDistance(double const *point,
                                     double const *force, double d) {
        
        // check z coordinate. If outside, return infinity
        if((point[2] + d * force[2]) > (_coords[2] + _height / 2) ||
           (point[2] + d * force[2]) < (_coords[2] - _height / 2))
        
        return numeric_limits<double>::infinity();
        
        vector<double> movedPoint{point[0] + d * force[0],
            point[1] + d * force[1],
            point[2] + d * force[2]};
        
        return distance(movedPoint);
        
    };
    
    virtual const vector<double> normal(const double *point) {
        auto dxy = _radius - twoPointDistance({_coords[0],_coords[1], 0},
                                              {  point[0],  point[1], 0});
        
        double dzz = point[2];
        if((_coords[2] * 2 - point[2]) < dzz) {
            dzz = _coords[2]*2 - point[2];
        }
        else {
            dzz = point[2];
        }
        
        if(dxy > dzz) {
           
            // when the Z coordinate is located at geometry center
            if(areEqual(point[2],_coords[2])) {
                return vector<double> {0.0, 0.0, 0.0};
            }
            else {
                return twoPointDirection({0,  0, point[2]},
                                          {0,0, _coords[2]});;
            }
            
        }
        else {
            
            // when the Z coordinate is located at geometry center
            if(areEqual(point[0],_coords[0]) && areEqual(point[1],_coords[1])) {
                return vector<double> {0.0, 0.0, 0.0};
            }
            else {
                return twoPointDirection({point[0],  point[1], 0},
                                         {_coords[0],_coords[1], 0});;
            }
        }
        
    };
    
    virtual const void elementeqn(double* var){
        cout<<"elementeqn not implemented for Boundary Element of CylindericalXYZ type. "
        "Exiting"<<endl;
        exit(EXIT_FAILURE);
        
    } ;

    //@}


};
#endif
