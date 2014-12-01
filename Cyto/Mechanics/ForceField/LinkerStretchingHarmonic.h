//
//  LinkerStretchingHarmonic.h
//  Cyto
//
//  Created by Konstantin Popov on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__LinkerStretchingHarmonic__
#define __Cyto__LinkerStretchingHarmonic__

#include <iostream>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

class LinkerStretchingHarmonic {
    
public:
    double energy(Bead*, Bead*, Bead*, Bead*, double position1, double position2, double kStr, double L);
    double energy(Bead*, Bead*, Bead*, Bead*, double position1, double position2, double kStr, double L, double d);
    void forces(Bead*, Bead*, Bead*, Bead*, double position1, double position2, double kStr, double L);
    void forcesAux(Bead*, Bead*, Bead*, Bead*, double position1, double position2, double kStr, double L );
};



#endif /* defined(__Cyto__LinkerStretchingHarmonic__) */
