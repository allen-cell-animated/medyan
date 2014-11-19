//
//  BoundaryImpl.h
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundaryImpl__
#define __Cyto__BoundaryImpl__

#include <vector>
#include <iostream>

#include "common.h"
#include "Boundary.h"

///Cubic boundary implementation
class BoundaryCubic: public Boundary {
    
public:
    ///Default constructor, this will create a cube with given corners at edges of current grid
    BoundaryCubic();
    virtual bool within(const vector<double> coordinates);
};

///Spherical boundary implementation
class BoundarySpherical: public Boundary {
    
public:
    ///Default constructor, will create an sphere with axes equal in length to current grid
    BoundarySpherical();
    virtual bool within(const vector<double> coordinates);
};

//Capsule boundary implementation
class BoundaryCapsule: public Boundary {
    
public:
    ///Default constructor, will create a capsule with axes equal in length to current grid
    BoundaryCapsule(double diameter);
    virtual bool within(const vector<double> coordinates);
};


#endif /* defined(__Cyto__BoundaryImpl__) */
