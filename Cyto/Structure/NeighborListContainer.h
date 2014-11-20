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
#include "NeighborListDB.h"


///NLContainer is an abstract class used for any interaction that needs a neighborlist
class NLContainer {
    
protected:
    NeighborList* _neighborList; ///< the neighborlist that this holds
    
    ///Default constructor
    NLContainer() {}
    
public:
    ///Destructor, removes this neighborlist from DB
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~NLContainer() noexcept {
        NeighborListDB::instance()->removeNeighborList(_neighborList);
    }
    
    ///Setter and getter for neighborlist
    NeighborList* getNeighborList() { return _neighborList; }
};


/// CylinderNLContainer, a concrete implementation of
/// neighbor list container. Holds a cylinder neighbor list.
class CylinderNLContainer : public NLContainer {
    
public:
    ///constructor, adds a cylinder neighbor list to the database
    CylinderNLContainer(float rMax = 0.0, float rMin = 0.0, bool crossFilamentOnly = false) {
        
        _neighborList = NeighborListDB::instance()->
                        createCylinderNeighborList(rMax, rMin, crossFilamentOnly);
    }
    ~CylinderNLContainer() {}
};

/// BoundaryElementNLContainer, a concrete implementation of
/// neighbor list container. Holds a boundary element / bead neighbor list.
class BoundaryElementNLContainer : public NLContainer {
    
public:
    ///constructor, adds a cylinder neighbor list to the database
    BoundaryElementNLContainer(float rMax = 0.0, float rMin = 0.0) {
        
        _neighborList = NeighborListDB::instance()->
                        createBoundaryElementNeighborList(rMax, rMin);
    }
    ~BoundaryElementNLContainer() {}
};


#endif /* defined(__Cyto__NeighborListContainer__) */
