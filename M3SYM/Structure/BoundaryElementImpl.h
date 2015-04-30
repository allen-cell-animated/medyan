
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

#ifndef M3SYM_BoundaryElementImpl_h
#define M3SYM_BoundaryElementImpl_h

#include <cmath>

#include "common.h"

#include "BoundaryElement.h"

/// A plane implementation of a BoundaryElement.
class PlaneBoundaryElement : public BoundaryElement {
    
private:
    /// Parameters of equation (ax + by + cz + d = 0)
    double _a, _b, _c, _d;

public:
    /// Constructor, sets parameters of equation
    PlaneBoundaryElement(vector<double> coords, vector<double> normal,
                         double repulsConst, double screenLength);
    
    virtual double distance(const vector<double>& point);
    
    virtual double stretchedDistance(const vector<double>& point,
                                     const vector<double>& force, double d);
    
    virtual const vector<double> normal(const vector<double>& point);
};

/// A spherical implementation of a BoundaryElement.
class SphereBoundaryElement : public BoundaryElement {
    
private:
    double _radius; ///< Radius of sphere
    
public:
    /// Constructor, sets parameters of equation
    SphereBoundaryElement(vector<double> coords, double radius,
                          double repulsConst, double screenLength)
        : BoundaryElement(coords, repulsConst, screenLength), _radius(radius) {}
    
    virtual double distance(const vector<double>& point);
    
    virtual double stretchedDistance(const vector<double>& point,
                                     const vector<double>& force, double d);
    
    virtual const vector<double> normal(const vector<double>& point);
};

/// A cylinder implementation of a BoundaryElement.
class CylindricalZBoundaryElement : public BoundaryElement {
    
private:
    double _radius; ///< Radius of cylinder
    double _height; ///< Height of cylinder
    
public:
    ///Constructor, sets parameters of equation
    CylindricalZBoundaryElement(vector<double> coords, double radius,
                                double height, double repulsConst, double screenLength)
    
        : BoundaryElement(coords, repulsConst, screenLength),
          _radius(radius), _height(height) {}
    
    virtual double distance(const vector<double>& point);
    
    virtual double stretchedDistance(const vector<double>& point,
                                     const vector<double>& force, double d);
    
    virtual const vector<double> normal(const vector<double>& point);
};

/// A half-sphere implementation of a BoundaryElement.
class HalfSphereZBoundaryElement : public BoundaryElement {
    
private:
    double _radius; ///< Radius of half sphere
    bool _up; ///< Whether the half sphere faces up or down
    
public:
    /// Constructor, sets parameters of equation
    HalfSphereZBoundaryElement(vector<double> coords, double radius,
                               bool up, double repulsConst, double screenLength)
    
        : BoundaryElement(coords, repulsConst, screenLength),
          _radius(radius), _up(up){}

    virtual double distance(const vector<double>& point);
    
    virtual double stretchedDistance(const vector<double>& point,
                                     const vector<double>& force, double d);
    
    virtual const vector<double> normal(const vector<double>& point);
};

#endif
