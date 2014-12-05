//
//  BranchingDihedralCosine.h
//  M3SYM
//
//  Created by Konstantin Popov on 12/4/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __M3SYM__BranchingDihedralCosine__
#define __M3SYM__BranchingDihedralCosine__

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// BranchingDihedralCosine class is a cosine potential used by the [BrancingBending](@ref BrancingBending) template.
class BranchingDihedralCosine {
    
public:
    double energy(Bead*, Bead*, Bead*, Bead*, double, double);
    double energy(Bead*, Bead*, Bead*, Bead*, double, double, double);
    void forces(Bead*, Bead*, Bead*, Bead*, double, double);
    void forcesAux(Bead*, Bead*, Bead*, Bead*, double, double);
};


#endif /* defined(__M3SYM__BranchingDihedralCosine__) */
