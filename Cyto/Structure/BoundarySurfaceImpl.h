//
//  BoundarySurfaceImpl.h
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundarySurfaceImpl__
#define __Cyto__BoundarySurfaceImpl__

#include "common.h"
#include "BoundarySurface.h"
#include <cmath>
#include <iostream>

///BasicPlane is a simple implementation of the BoundarySurface class
class Plane: public BoundarySurface {
    
private:
    std::vector<double> _coords; ///< coordinates of center
    std::vector<double> _normal; ///< normal vector
    
public:
    
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of plane
    ///@param normal - normal vector to plane
    Plane(std::vector<double> coords, std::vector<double> normal);
};


///Sphere is a simple implementation of the BoundarySurface class
///@note this represents a half-ellipsoid. To create a full ellipsoid, one must define two of these objects.
class Sphere: public BoundarySurface {
    
private:
    std::vector<double> _coords; ///< center of sphere
    double _radius; ///<radius of sphere
    
public:
    ///Constructor, creates boundary elements
    ///@param coords - coordinates of center of plane
    Sphere(std::vector<double> coords, double radius);
};








#endif /* defined(__Cyto__BoundarySurfaceImpl__) */
