//
//  BranchingFF.h
//  M3SYM
//
//  Created by Konstantin Popov on 12/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __M3SYM__BranchingFF__
#define __M3SYM__BranchingFF__

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class BranchingInteractions;

/// Branching FF is an implementation of the [ForceField](@ref ForceField) class that calculates [Branching] (@ref Branching)
/// position of branced chain, angle of branced chain , and dihedral for branchinf chain.
class BranchingFF : public ForceField {
    
private:
    vector<unique_ptr<BranchingInteractions>> _branchingInteractionVector; ///< Vector of initialized branching interactions
    
public:
    /// Constructor, intializes all interaction at the branching point
    BranchingFF(string& stretching, string& bending, string& dihedral);
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
};

#endif /* defined(__M3SYM__BranchingFF__) */
