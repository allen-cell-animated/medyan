//
//  BranchingStretching.h
//  M3SYM
//
//  Created by Konstantin Popov on 12/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __M3SYM__BranchingStretching__
#define __M3SYM__BranchingStretching__

#include "common.h"

#include "BranchingInteractions.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// Branching stratching class represents an interaction fixing a branching chain on the main chain
template <class BStretchingInteractionType>
class BranchingStretching : public BranchingInteractions {
    
private:
    BStretchingInteractionType _FFType;
    
public:
    virtual double computeEnergy(BranchingPoint*, double d);
    virtual void computeForces(BranchingPoint*);
    virtual void computeForcesAux(BranchingPoint*);
};



#endif /* defined(__M3SYM__BranchingStretching__) */
