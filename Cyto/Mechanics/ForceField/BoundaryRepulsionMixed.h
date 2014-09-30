//
//  BoundaryRepulsionMixed.h
//  Cyto
//
//  Created by James Komianos on 9/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundaryRepulsionMixed__
#define __Cyto__BoundaryRepulsionMixed__

#include <iostream>
#include <vector>

#include "common.h"

class Bead;

class BoundaryRepulsionMixed {
    
public:
    double ComputeEnergy(Bead*, double, double);
    void ComputeForces(Bead*, double, std::vector<double>& norm, double );
    void ComputeForcesAux(Bead*, double, std::vector<double>& norm, double );
};

#endif /* defined(__Cyto__BoundaryRepulsionMixed__) */
