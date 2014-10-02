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

class Bead;

class BoundaryRepulsionExp {
    
public:
    double ComputeEnergy(Bead*, double, double, double);
    void ComputeForces(Bead*, double, std::vector<double>& norm, double, double );
    void ComputeForcesAux(Bead*, double, std::vector<double>& norm, double, double );
};

#endif /* defined(__Cyto__BoundaryRepulsionExp__) */
