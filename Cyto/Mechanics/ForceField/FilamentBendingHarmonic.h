//
//  FilamentBendingHarmonic.h
//  Cyto
//
//  Created by Konstantin Popov on 8/27/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__FilamentBendingHarmonic__
#define __Cyto__FilamentBendingHarmonic__

#include <iostream>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

class FilamentBendingHarmonic {
    
public:
    double energy(Bead*, Bead*, Bead*, double);
    double energy(Bead*, Bead*, Bead*, double, double);
    void forces(Bead*, Bead*, Bead*, double);
    void forcesAux(Bead*, Bead*, Bead*, double);
};

#endif /* defined(__Cyto__FilamentBendingHarmonic__) */
