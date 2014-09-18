//
//  BoundaryRepulsionLJ.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/16/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryRepulsionLJ.h"

double BoundaryRepulsionLJ::ComputeEnergy(Bead* pb, double r, double k_rep)
{
    double inv_r4 = 1/r * 1/r * 1/r * 1/r;
    return k_rep * inv_r4 * inv_r4 * inv_r4;
    
}

void BoundaryRepulsionLJ::ComputeForces(Bead* pb, double r,  std::vector<double> norm ){
    
    double inv_r4 = 1/r * 1/r * 1/r * 1/r;
    
    pb1->force[0] +=   k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[0];
    
    pb1->force[1] +=   k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[1];
    
    pb1->force[2] +=   k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[2];
    
}

void BoundaryRepulsionLJ::ComputeForcesAux(Bead* pb, double r, std::vector<double> norm ){
    
    double inv_r4 = 1/r * 1/r * 1/r * 1/r;
    
    pb1->forceAux[0] +=   k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[0];
    
    pb1->forceAux[1] +=   k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[1];
    
    pb1->forceAux[2] +=   k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[2];
    
}