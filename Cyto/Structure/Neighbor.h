//
//  Neighbor.h
//  Cyto
//
//  Created by James Komianos on 11/12/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__Neighbor__
#define __Cyto__Neighbor__

#include <stdio.h>

#include "common.h"

#include "NeighborListDB.h"

///Neighbor class is a class that can be added/removed from a neighborlist
///@note - any subclass that inherits from Neighbor MUST add itself to the neighbors list DB.
class Neighbor {
    
protected:
    ///Constructor
    Neighbor() {}
    
public:
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~Neighbor() noexcept {
        ///remove from neighbor lists
        NeighborListDB::instance()->removeNeighbor(this);
    }
    
};

#endif /* defined(__Cyto__Neighbor__) */
