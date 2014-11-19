//
//  BoundaryElementImpl.h
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundaryElementImpl__
#define __Cyto__BoundaryElementImpl__

#include <iostream>
#include <cmath>

#include "common.h"
#include "BoundaryElement.h"

///PlaneBoundaryElement is a plane implementation of a boundary element
class PlaneBoundaryElement : public BoundaryElement {
    
private:
    ///Parameters of equation (ax + by + cz + d = 0)
    double _a, _b, _c, _d;
    double _k_rep;
    double _r0;

public:
    ///Constructor, sets parameters of equation
    PlaneBoundaryElement(vector<double> coords, vector<double> normal, double repulsConst, double screenLength);
    
    ///Distance to this plane
    virtual double distance(const vector<double>& point);
    
    ///Stretched distance to plane
    virtual double stretchedDistance(const vector<double>& point, const vector<double>& force, double d);
    
    ///normal to plane
    virtual const vector<double> normal(const vector<double>& point);
    
    ///get the repulsion constant and screening length for this plane
    virtual double getRepulsionConst();
    virtual double getScreeningLength();
};

class SphereBoundaryElement : public BoundaryElement {
    
private:
    double _radius;
    double _k_rep;
    double _r0;
    
public:
    ///Constructor, sets parameters of equation
    SphereBoundaryElement(vector<double> coords, double radius, double repulsConst, double screenLength);
    
    ///Distance to this sphere
    virtual double distance(const vector<double>& point);
    
    ///Stretched distance to sphere
    virtual double stretchedDistance(const vector<double>& point, const vector<double>& force, double d);
    
    ///normal to sphere
    virtual const vector<double> normal(const vector<double>& point);
    
    ///get the repulsion constant and screening length for this plane
    virtual double getRepulsionConst();
    virtual double getScreeningLength();
};

class CylindricalZBoundaryElement : public BoundaryElement {
    
private:
    double _radius;
    double _height;
    double _k_rep;
    double _r0;
    
public:
    ///Constructor, sets parameters of equation
    CylindricalZBoundaryElement(vector<double> coords, double radius, double height, double repulsConst, double screenLength);
    
    ///Distance to this sphere
    virtual double distance(const vector<double>& point);
    
    ///Stretched distance to sphere
    virtual double stretchedDistance(const vector<double>& point, const vector<double>& force, double d);
    
    ///normal to sphere
    virtual const vector<double> normal(const vector<double>& point);
    
    ///get the repulsion constant and screening length for this plane
    virtual double getRepulsionConst();
    virtual double getScreeningLength();
};


class HalfSphereZBoundaryElement : public BoundaryElement {
    
private:
    double _radius;
    bool _up; ///< whether the half sphere faces up or down
    double _k_rep;
    double _r0;
    
public:
    ///Constructor, sets parameters of equation
    HalfSphereZBoundaryElement(vector<double> coords, double radius, bool up, double repulsConst, double screenLength);
    
    ///Distance to this sphere
    virtual double distance(const vector<double>& point);
    
    ///Stretched distance to sphere
    virtual double stretchedDistance(const vector<double>& point, const vector<double>& force, double d);
    
    ///normal to sphere
    virtual const vector<double> normal(const vector<double>& point);
    
    ///get the repulsion constant and screening length for this plane
    virtual double getRepulsionConst();
    virtual double getScreeningLength();
};



#endif /* defined(__Cyto__BoundaryElementImpl__) */
