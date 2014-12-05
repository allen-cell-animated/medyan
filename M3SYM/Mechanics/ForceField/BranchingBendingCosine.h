//
//  BranchingBendingCosine.h
//  M3SYM
//
//  Created by Konstantin Popov on 12/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __M3SYM__BranchingBendingCosine__
#define __M3SYM__BranchingBendingCosine__



#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// BranchingBendingCosine class is a cosine potential used by the [BrancingBending](@ref BrancingBending) template.
class BranchingBendingCosine {
    
public:
    double energy(Bead*, Bead*, Bead*, Bead*, double, double);
    double energy(Bead*, Bead*, Bead*, Bead*, double, double, double);
    void forces(Bead*, Bead*, Bead*, Bead*, double, double);
    void forcesAux(Bead*, Bead*, Bead*, Bead*, double, double);
};




#endif /* defined(__M3SYM__BranchingBendingCosine__) */
