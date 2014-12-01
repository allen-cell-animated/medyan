//
//  Movable.h
//  Cyto
//
//  Created by James Komianos on 11/25/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__Movable__
#define __Cyto__Movable__

#include <iostream>

///Movable class is for a movable object in the subsystem
class Movable {
    
public:
    virtual void updatePosition() = 0;
};

#endif /* defined(__Cyto__Movable__) */
