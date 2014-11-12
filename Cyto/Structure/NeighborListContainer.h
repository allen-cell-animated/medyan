//
//  NeighborListContainer.h
//  Cyto
//
//  Created by James Komianos on 11/12/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__NeighborListContainer__
#define __Cyto__NeighborListContainer__

#include <stdio.h>

#include "common.h"
#include "NeighborList.h"

class NeighborListContainer {
    
protected:
    NeighborList* _neighborList;
    
public:
    
    NeighborListContainer() {}
    virtual ~NeighborListContainer() noexcept {}
    
    void setNeighborList(NeighborList* neighborList) {_neighborList = neighborList;}
    NeighborList* getNeighborList() {return _neighborList;}
};


#endif /* defined(__Cyto__NeighborListContainer__) */
