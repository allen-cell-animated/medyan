
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_BoundaryElementImpl_h
#define M3SYM_BoundaryElementImpl_h

#include <cmath>

#include "common.h"

#include "BoundaryElement.h"

/// A plane implementation of a boundary element.
class PlaneBoundaryElement : public BoundaryElement {
    
private:
    /// Parameters of equation (ax + by + cz + d = 0)
    double _a, _b, _c, _d;
    double _k_rep;
    double _r0;

public:
    /// Constructor, sets parameters of equation
    PlaneBoundaryElement(vector<double> coords, vector<double> normal,
                         double repulsConst, double screenLength);
    
    virtual double distance(const vector<double>& point);
    virtual double stretchedDistance(const vector<double>& point, const vector<double>& force, double d);
    virtual const vector<double> normal(const vector<double>& point);
    
    virtual double getRepulsionConst();
    virtual double getScreeningLength();
};

/// A spherical implementation of a boundary element.
class SphereBoundaryElement : public BoundaryElement {
    
private:
    double _radius; ///< Radius of sphere
    double _k_rep;
    double _r0;
    
public:
    /// Constructor, sets parameters of equation
    SphereBoundaryElement(vector<double> coords, double radius,
                          double repulsConst, double screenLength);
    
    virtual double distance(const vector<double>& point);
    virtual double stretchedDistance(const vector<double>& point, const vector<double>& force, double d);
    virtual const vector<double> normal(const vector<double>& point);
    
    virtual double getRepulsionConst();
    virtual double getScreeningLength();
};

/// A cylinder implementation of a boundary element.
class CylindricalZBoundaryElement : public BoundaryElement {
    
private:
    double _radius; ///< Radius of cylinder
    double _height; ///< Height of cylinder
    double _k_rep;
    double _r0;
    
public:
    ///Constructor, sets parameters of equation
    CylindricalZBoundaryElement(vector<double> coords, double radius,
                                double height, double repulsConst, double screenLength);
    
    virtual double distance(const vector<double>& point);
    virtual double stretchedDistance(const vector<double>& point, const vector<double>& force, double d);
    virtual const vector<double> normal(const vector<double>& point);
    
    virtual double getRepulsionConst();
    virtual double getScreeningLength();
};

/// A half-sphere implementat of a boundary element.
class HalfSphereZBoundaryElement : public BoundaryElement {
    
private:
    double _radius; ///< Radius of half sphere
    bool _up; ///< whether the half sphere faces up or down
    double _k_rep;
    double _r0;
    
public:
    /// Constructor, sets parameters of equation
    HalfSphereZBoundaryElement(vector<double> coords, double radius,
                               bool up, double repulsConst, double screenLength);

    virtual double distance(const vector<double>& point);
    virtual double stretchedDistance(const vector<double>& point, const vector<double>& force, double d);
    virtual const vector<double> normal(const vector<double>& point);
    
    virtual double getRepulsionConst();
    virtual double getScreeningLength();
};

#endif
