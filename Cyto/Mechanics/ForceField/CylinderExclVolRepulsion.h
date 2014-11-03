//
//  CylindricalVolRepulsion.h
//  Cyto
//
//  Created by Konstantin Popov on 10/29/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__CylinderExclVolRepulsion__
#define __Cyto__CylinderExclVolRepulsion__

#include <stdio.h>
#include "common.h"

class Bead;

class CylinderExclVolRepulsion {
    
public:
    double Energy(Bead*, Bead*, Bead*, Bead*, double Krepuls);
    double Energy(Bead*, Bead*, Bead*, Bead*, double Krepuls, double d);
    void Forces(Bead*, Bead*, Bead*, Bead*, double Krepuls);
    void ForcesAux(Bead*, Bead*, Bead*, Bead*, double Krepuls );
};


#endif /* defined(__Cyto__CylinderExclVolRepulsion__) */
