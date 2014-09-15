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
#include "BoundaryElement.h"
#include "BoundaryInteractions.h"

template <class BRepulsionInteractionType>
class BoundaryRepulsion : public BoundaryInteractions
{
    
private:
    BRepulsionInteractionType _FFType;
    
public:
    virtual double ComputeEnergy( BoundaryElement*, double d);
    virtual void ComputeForces( BoundaryElement*);
    virtual void ComputeForcesAux( BoundaryElement*);
};
#endif /* defined(__Cyto__BoundaryRepulsion__) */
