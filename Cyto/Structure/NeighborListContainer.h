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


///NeighborListContainer is an abstract class used for any interaction that needs a neighborlist
class NeighborListContainer {
    
protected:
    NeighborList* _neighborList; ///< the neighborlist that this holds
    
    ///Default constructor
    NeighborListContainer() {}
    
public:
    ///Destructor, removes this neighborlist from DB
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~NeighborListContainer() noexcept {
        NeighborListDB::instance(NeighborListDBKey())->removeNeighborList(_neighborList);
    }
    
    ///Setter and getter for neighborlist
    NeighborList* getNeighborList() { return _neighborList; }
};


/// CylinderNeighborListContainer, a concrete implementation of
/// neighbor list container. Holds a cylinder neighbor list.
class CylinderNeighborListContainer : public NeighborListContainer {
    
public:
    ///constructor, adds a cylinder neighbor list to the database
    CylinderNeighborListContainer(float rMax = 0.0, float rMin = 0.0, bool crossFilamentOnly = false) {
        
        _neighborList = NeighborListDB::instance(NeighborListDBKey())->
                         createCylinderNeighborList(rMax, rMin, false);
    }
    ~CylinderNeighborListContainer() {}
};




#endif /* defined(__Cyto__NeighborListContainer__) */
