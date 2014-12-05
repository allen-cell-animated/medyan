//
//  BranchingDihedral.h
//  M3SYM
//
//  Created by Konstantin Popov on 12/4/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __M3SYM__BranchingDihedral__
#define __M3SYM__BranchingDihedral__


#include "common.h"

#include "BranchingInteractions.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// Branching Dihedral class represents an interaction keeping branching in plane
template <class BDihedralInteractionType>
class BranchingDihedral : public BranchingInteractions {
    
private:
    BDihedralInteractionType _FFType;
    
public:
    virtual double computeEnergy(BranchingPoint*, double d);
    virtual void computeForces(BranchingPoint*);
    virtual void computeForcesAux(BranchingPoint*);
};




#endif /* defined(__M3SYM__BranchingDihedral__) */
