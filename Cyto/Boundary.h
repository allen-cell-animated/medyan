//
//  Boundary.h
//  Cyto
//
//  Created by James Komianos on 8/6/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__Boundary__
#define __Cyto__Boundary__

#include <iostream>
class Bead;

class Boundary {

public:
    
    virtual double ComputeDistance(Bead* b) = 0;
};


#endif /* defined(__Cyto__Boundary__) */
