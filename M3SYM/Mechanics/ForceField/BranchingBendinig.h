//
//  BranchingBendinig.h
//  M3SYM
//
//  Created by Konstantin Popov on 12/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __M3SYM__BranchingBendinig__
#define __M3SYM__BranchingBendinig__

#include "common.h"

#include "BranchingInteractions.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// Branching stratching class represents an interaction keeping branching angle at phi0 (~270 for Arp2/3)
template <class BBendingInteractionType>
class BranchingBending : public BranchingInteractions {
    
private:
    BBendingInteractionType _FFType;
    
public:
    virtual double computeEnergy(BranchingPoint*, double d);
    virtual void computeForces(BranchingPoint*);
    virtual void computeForcesAux(BranchingPoint*);
};

#endif /* defined(__M3SYM__BranchingBendinig__) */
