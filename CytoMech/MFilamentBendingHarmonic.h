//
//  MFilamentBendingHarmonic.h
//  Cyto
//
//  Created by Konstantin Popov on 8/27/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MFilamentBendingHarmonic__
#define __Cyto__MFilamentBendingHarmonic__

#include <iostream>
#include "MBead.h"


class FilamentBendingHarmonic
{
    
public:
    double Energy(Bead*, Bead*, Bead*, double);
    double Energy(Bead*, Bead*, Bead*, double, double);
    void Forces(Bead*, Bead*, Bead*, double);
    void ForcesAux(Bead*, Bead*, Bead*, double);
    
};



#endif /* defined(__Cyto__MFilamentBendingHarmonic__) */
