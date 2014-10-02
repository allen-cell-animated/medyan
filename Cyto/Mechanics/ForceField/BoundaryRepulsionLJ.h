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

class Bead;

class BoundaryRepulsionLJ {

public:
    double ComputeEnergy(Bead*, double, double, double);
    void ComputeForces(Bead*, double, std::vector<double>& norm, double, double);
    void ComputeForcesAux(Bead*, double, std::vector<double>& norm, double, double);
};
#endif /* defined(__Cyto__BoundaryRepulsionLJ__) */
