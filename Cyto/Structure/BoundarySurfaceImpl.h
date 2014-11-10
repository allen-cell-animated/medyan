//
//  BoundarySurfaceImpl.h
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundarySurfaceImpl__
#define __Cyto__BoundarySurfaceImpl__

#include <cmath>
#include <iostream>

#include "common.h"
#include "BoundarySurface.h"

///BasicPlane is a simple implementation of the BoundarySurface class
class Plane: public BoundarySurface {
    
private:
    vector<double> _coords; ///< coordinates of center
    vector<double> _normal; ///< normal vector
    
public:
    
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of plane
    ///@param normal - normal vector to plane
    Plane(vector<double> coords, vector<double> normal);
};


///Sphere is a simple implementation of the BoundarySurface class
///@note this represents a full sphere, not just a half
class Sphere: public BoundarySurface {
    
private:
    vector<double> _coords; ///< center of sphere
    double _radius; ///<radius of sphere
    
public:
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of plane
    Sphere(vector<double> coords, double radius);
};








#endif /* defined(__Cyto__BoundarySurfaceImpl__) */
