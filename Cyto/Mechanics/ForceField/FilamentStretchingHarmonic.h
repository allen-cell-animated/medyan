//
//  FilamentStretchingHarmonic.h
//  Cyto
//
//  Created by Konstantin Popov on 8/27/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_FilamentStretchingHarmonic_h
#define Cyto_FilamentStretchingHarmonic_h

#include "common.h"

///FORWARD DECLARATIONS
class Bead;

class FilamentStretchingHarmonic {
    
public:
    double energy(Bead*, Bead*, double, double);
    double energy(Bead*, Bead*, double, double, double);
    void forces(Bead*, Bead*, double, double);
    void forcesAux(Bead*, Bead*, double, double);
};

#endif /* defined(__Cyto__FilamentStretchingHarmonic__) */
