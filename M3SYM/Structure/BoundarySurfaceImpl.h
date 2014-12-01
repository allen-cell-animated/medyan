
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

#ifndef M3SYM_BoundarySurfaceImpl_h
#define M3SYM_BoundarySurfaceImpl_h

#include <cmath>

#include "common.h"

#include "BoundarySurface.h"

/// BasicPlane is a simple implementation of the BoundarySurface class
class Plane: public BoundarySurface {
    
private:
    vector<double> _coords; ///< Coordinates of center
    vector<double> _normal; ///< Normal vector
    
public:
    
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of plane
    ///@param normal - normal vector to plane
    Plane(vector<double> coords, vector<double> normal);
};


/// Sphere is a simple implementation of the BoundarySurface class
/// @note this represents a full sphere, not just a half
class Sphere: public BoundarySurface {
    
private:
    vector<double> _coords; ///< Center of sphere
    double _radius; ///< Radius of sphere
    
public:
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of plane
    ///@param normal - normal vector to sphere
    Sphere(vector<double> coords, double radius);
};


/// CylinderZ is a simple implementation of the BoundarySurface class
class CylinderZ: public BoundarySurface {
    
private:
    vector<double> _coords; ///< Center of cylinder
    double _radius; ///< Radius of cylinder
    double _height; ///< Height of cylinder
    
public:
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of cylinder
    ///@param normal - normal vector to sphere
    CylinderZ(vector<double> coords, double radius, double height);
    
};

/// HalfSphereZ is a simple implementation of the BoundarySurface class
class HalfSphereZ: public BoundarySurface {
    
private:
    vector<double> _coords; ///< Center of half-sphere
    double _radius; ///< Radius of half-sphere
    bool _up; ///< Facing up or down
    
public:
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of half sphere
    ///@param normal - normal vector to sphere
    HalfSphereZ(vector<double> coords, double radius, bool up);
    
};

#endif