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
#include <iostream>

///BasicPlane is a simple implementation of the BoundarySurface class
class BasicPlane: public BoundarySurface {
    
private:
    std::vector<std::vector<double>> _points; ///<the points that define this plane (4 for 3D, 2 for 2D)
    short _orientation; ///<Direction of normal vector to plane (for basic plane, can be unit vectors only)
    
public:
    
    ///Constructor, creates boundary elements
    BasicPlane(int nDim, std::vector<std::vector<double>> points, short orientation);

};




#endif /* defined(__Cyto__BoundarySurfaceImpl__) */
