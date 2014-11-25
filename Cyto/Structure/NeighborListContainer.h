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

#include "NeighborListDB.h"


/// CylinderNLContainer, holds a cylinder neighbor list.
class CylinderNLContainer {

private:
    CylinderNeighborList* _neighborList;
    
public:
    ///constructor, adds a cylinder neighbor list to the database
    CylinderNLContainer(float rMax = 0.0, float rMin = 0.0, bool crossFilamentOnly = false) {
        
        _neighborList = NeighborListDB::instance()->
                        createCylinderNeighborList(rMax, rMin, crossFilamentOnly);
    }
    ~CylinderNLContainer() { NeighborListDB::instance()->removeNeighborList(_neighborList); }
    
    CylinderNeighborList* getNeighborList() {return _neighborList;}
    
};

/// BoundaryElementNLContainer, holds a boundary element / bead neighbor list.
class BoundaryElementNLContainer {
    
private:
    BoundaryElementNeighborList* _neighborList;

public:
    ///constructor, adds a cylinder neighbor list to the database
    BoundaryElementNLContainer(float rMax = 0.0, float rMin = 0.0) {
        
        _neighborList = NeighborListDB::instance()->
                        createBoundaryElementNeighborList(rMax, rMin);
    }
    ~BoundaryElementNLContainer() { NeighborListDB::instance()->removeNeighborList(_neighborList); }
    
    BoundaryElementNeighborList* getNeighborList() {return _neighborList;}
};


#endif /* defined(__Cyto__NeighborListContainer__) */
