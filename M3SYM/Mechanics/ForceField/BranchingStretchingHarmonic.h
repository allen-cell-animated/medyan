//
//  BranchingStretchingHarmonic.h
//  M3SYM
//
//  Created by Konstantin Popov on 12/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __M3SYM__BranchingStretchingHarmonic__
#define __M3SYM__BranchingStretchingHarmonic__

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// BranchingStretchingHarmonic class is a harmonic potential used by the [BranchingStretching](@ref BranchingStretching) template.
class BranchingStretchingHarmonic {
    
public:
    double energy(Bead*, Bead*, Bead*, double position, double kStr, double L0);
    double energy(Bead*, Bead*, Bead*, double position, double kStr, double L0, double d);
    void forces(Bead*, Bead*, Bead*, double position, double kStr, double L0);;
    void forcesAux(Bead*, Bead*, Bead*, double position, double kStr, double L0);;
    
};



#endif /* defined(__M3SYM__BranchingStretchingHarmonic__) */
