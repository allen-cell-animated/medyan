
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BoundarySurfaceImpl_h
#define MEDYAN_BoundarySurfaceImpl_h

#include <cmath>

#include "common.h"

#include "BoundarySurface.h"

/// A simple implementation of the BoundarySurface class.
class Plane: public BoundarySurface {
    
private:
    vector<double> _coords; ///< Coordinates of center
    vector<double> _normal; ///< Normal vector
    
public:
    
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of plane
    ///@param normal - normal vector to plane
    Plane(SubSystem* s, vector<double> coords, vector<double> normal);
};


/// A simple implementation of the BoundarySurface class.
/// @note this represents a full sphere, not just a half
class Sphere: public BoundarySurface {
    
private:
    vector<double> _coords; ///< Center of sphere
    
public:
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of plane
    ///@param normal - normal vector to sphere
    Sphere(SubSystem* s, vector<double> coords, double radius);
};


/// A simple implementation of the BoundarySurface class.
class CylinderZ: public BoundarySurface {
    
private:
    vector<double> _coords; ///< Center of cylinder
    
public:
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of cylinder
    ///@param normal - normal vector to sphere
    CylinderZ(SubSystem* s, vector<double> coords, double radius, double height);
    
};

/// A simple implementation of the BoundarySurface class.
class HalfSphereZ: public BoundarySurface {
    
private:
    vector<double> _coords; ///< Center of half-sphere
    
public:
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of half sphere
    ///@param normal - normal vector to sphere
    HalfSphereZ(SubSystem* s, vector<double> coords, double radius, bool up);
    
};

//Qin -----
/// A simple implementation of the BoundarySurface class.
class CylinderXYZ: public BoundarySurface {
    
private:
    vector<double> _coords; ///< Center of cylinder
    
public:
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of cylinder
    ///@param normal - normal vector to sphere
    CylinderXYZ(SubSystem* s, vector<double> coords, double radius, double height);
    
};

#endif
