//
//  MFilamentHarmonicStretching.h
//  Cyto
//
//  Created by Konstantin Popov on 8/27/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_MFilamentStretchingHarmonic_h
#define Cyto_FilamentStretchingHarmonic_h

#include "MBead.h"


class FilamentStretchingHarmonic
{
    
public:
    double Energy(Bead*, Bead*, double, double);
    double Energy(Bead*, Bead*, double, double, double);
    void Forces(Bead*, Bead*, double, double);
    void ForcesAux(Bead*, Bead*, double, double);
    
};



#endif
