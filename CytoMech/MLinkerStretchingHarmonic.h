//
//  MLinkerStretchingHarmonic.h
//  Cyto
//
//  Created by Konstantin Popov on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MLinkerStretchingHarmonic__
#define __Cyto__MLinkerStretchingHarmonic__

#include <iostream>
#include "MBead.h"


class LinkerStretchingHarmonic
{
    
public:
    double Energy(Bead*, Bead*, Bead*, Bead*, double position1, double position2, double kStr, double L);
    double Energy(Bead*, Bead*, Bead*, Bead*, double position1, double position2, double kStr, double L, double d);
    void Forces(Bead*, Bead*, Bead*, Bead*, double position1, double position2, double kStr, double L);
    void ForcesAux(Bead*, Bead*, Bead*, Bead*, double position1, double position2, double kStr, double L );
    
};



#endif /* defined(__Cyto__MLinkerStretchingHarmonic__) */
