//
//  BoundarySurfaceImpl.h
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundarySurfaceImpl__
#define __Cyto__BoundarySurfaceImpl__

#include "BoundarySurface.h"
#include <cmath>
#include <iostream>

///BasicPlane is a simple implementation of the BoundarySurface class
class BasicPlane: public BoundarySurface {
    
private:
    std::vector<std::vector<double>> _points; ///<the points that define this plane (4 for 3D, 2 for 2D)
    short _orientation; ///<Direction of normal vector to plane (for basic plane, can be unit vectors only)
                        ///@note 2 for XY, 1 for XZ, 0 for YZ
    
public:
    
    ///Constructor, creates boundary elements
    ///@note - points must be in order of LL, LR, UR, UL
    ///@param numDivisions - number of divisions per side for boundary elements
    BasicPlane(std::vector<std::vector<double>> points, std::vector<int> numDivisions, short orientation);

};


///BasicEllipsoid is a simple implementation of the BoundarySurface class
///@note this represents a half-ellipsoid. To create a full ellipsoid, one must define two of these objects.
class BasicEllipsoid: public BoundarySurface {
    
private:
    double _A, _B, _C; ///< parameters of equation for ellipsoid
    std::vector<double> _center; /// center of ellipsoid
    short _orientation; ///< Direction of normal vector to bottom/top of ellipsoid
                        ///@note +/- denotes upward or downward in that direction
                        /// +/- 0: X direction, +/- 1: Y, +/- 2: Z
    
//    
//    ///Returns a value based on the given ellipsoid equation
//    ///@note CASE 1: if in X direction, p1 and p2 must be Y, Z coordinates, returns X value
//    ///      CASE 2: if in Y direction, p1 and p2 must be X, Z coordinates, returns Y value
//    ///      CASE 3: if in Z direction, p1 and p2 must be X, Y coordinates, returns Z value
//    double value (double p1, double p2) {
//        
//        
//    }
//    
    
public:
    ///Constructor, creates boundary elements
    
    
};








#endif /* defined(__Cyto__BoundarySurfaceImpl__) */
