//
//  BoundaryRepulsionLJ.h
//  Cyto
//
//  Created by Konstantin Popov on 9/16/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundaryRepulsionLJ__
#define __Cyto__BoundaryRepulsionLJ__

#include <iostream>
#include <vector>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

class BoundaryRepulsionLJ {

public:
    double computeEnergy(Bead*, double, double, double);
    void computeForces(Bead*, double, vector<double>& norm, double, double);
    void computeForcesAux(Bead*, double, vector<double>& norm, double, double);
};
#endif /* defined(__Cyto__BoundaryRepulsionLJ__) */
