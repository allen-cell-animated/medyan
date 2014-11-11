//
//  NeighborsManager.h
//  Cyto
//
//  Created by James Komianos on 11/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__NeighborsManager__
#define __Cyto__NeighborsManager__

#include <stdio.h>

#include "common.h"
#include "Graph.h"

struct CylinderProperties {
    Cylinder* _cylinder; ///< ptr to the cylinder that this vertex represents
};

struct PairProperties {
    ///Nothing for now, could store distances....?
};

typedef Graph<CylinderProperties, PairProperties> NeighborList;


class NeighborListManager {
    
private:
    std::vector<NeighborList> reactionNeighborLists;
    std::vector<NeighborList> forceNeighborLists;

public:
    
    
    
    
    
    
    
    
};



#endif /* defined(__Cyto__NeighborsManager__) */
