//
//  BoundarySurface.h
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundarySurface__
#define __Cyto__BoundarySurface__

#include "common.h"
#include "BoundaryElementDB.h"
#include <iostream>

/// BoundarySurface class is a boundary shape that holds boundary elements
/*!
 *  The BoundarySurface class is a basic boundary shape that is owned and
 *  controlled by the Boundary that created it. It holds a vector of
 *  BoundaryElements as well as any other geometrical information needed
 *  for the given implementation.
 */

class BoundarySurface {
    
protected:
    std::vector<std::unique_ptr<BoundaryElement>> _boundaryElements;
    ///<vector of boundary elements that make up this surface
    short _nDim; ///<dimensionality of surface
    
public:
    
    ///Default constructor and destructor
    BoundarySurface(int nDim) : _nDim(nDim) {};
    
    ~BoundarySurface() {
        //loop through boundary elements, remove from DB
        for (auto &b : _boundaryElements)
            BoundaryElementDB::Instance(BEDBKey())->RemoveBoundaryElement(b.get());
    };
    
    ///Access for all implementations of BoundarySurface to the DB key
    BoundaryElementDBKey BEDBKey() {return BoundaryElementDBKey();}
    ///Get boundary elements
    const std::vector<std::unique_ptr<BoundaryElement>>& boundaryElements() {return _boundaryElements;}
    
};

#endif /* defined(__Cyto__BoundarySurface__) */
