//
//  BoundaryImpl.h
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundaryImpl__
#define __Cyto__BoundaryImpl__

#include "Boundary.h"
#include "BoundarySurfaceImpl.h"
#include <vector>
#include <iostream>


///Cubic boundary implementation
class BoundaryCubic: public Boundary {
    
public:
    ///Default constructor, this will create a cube with given corners
    ///@param points - 8 corners of cube, in the following order:
    /// lower plane - front left, front right, back right, back left
    /// upper plane - front left, front right, back right, back left
    ///@param numDivisions - number of divisions for each plane (into square boundary elements
    BoundaryCubic(int nDim, std::vector<std::vector<double>> points, std::vector<int> numDivisions);
};

///Cylindrical boundary implementation
class BoundaryCylindrical: public Boundary {
    ///not yet implemented
};

///Spherical boundary implementation
class BoundarySpherical: public Boundary {
    ///not yet implemented
};


#endif /* defined(__Cyto__BoundaryImpl__) */
