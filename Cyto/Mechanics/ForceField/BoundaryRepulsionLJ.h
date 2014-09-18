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
#include "Bead.h"

class BoundaryRepulsionLJ{

public:
    double ComputeEnergy(Bead*, double);
    double ComputeEnergy(Bead*, double, double);
    void ComputeForces(Bead*, double, std::vector<double> norm);
    void ComputeForcesAux(Bead*, double, double, std::vector<double> norm);
};
#endif /* defined(__Cyto__BoundaryRepulsionLJ__) */
