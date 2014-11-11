//
//  BoundaryRepulsionExp.h
//  Cyto
//
//  Created by Konstantin Popov on 10/1/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundaryRepulsionExp__
#define __Cyto__BoundaryRepulsionExp__

#include <iostream>
#include <vector>

#include "common.h"

///FORWARD DECLARATIONS
class Bead;

class BoundaryRepulsionExp {
    
public:
    double computeEnergy(Bead*, double, double, double);
    void computeForces(Bead*, double, vector<double>& norm, double, double);
    void computeForcesAux(Bead*, double, vector<double>& norm, double, double);
};

#endif /* defined(__Cyto__BoundaryRepulsionExp__) */
