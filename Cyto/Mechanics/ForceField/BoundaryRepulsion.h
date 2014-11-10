//
//  BoundaryRepulsion.h
//  Cyto
//
//  Created by Konstantin Popov on 9/12/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundaryRepulsion__
#define __Cyto__BoundaryRepulsion__

#include <iostream>
#include <vector>

#include "common.h"
#include "BoundaryInteractions.h"

class BoundaryElement;


template <class BRepulsionInteractionType>
class BoundaryRepulsion : public BoundaryInteractions {
    
private:
    BRepulsionInteractionType _FFType;
    
public:
    virtual double computeEnergy( BoundaryElement*, double d);
    virtual void computeForces( BoundaryElement*);
    virtual void computeForcesAux( BoundaryElement*);
};
#endif /* defined(__Cyto__BoundaryRepulsion__) */
