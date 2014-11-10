//
//  BoundaryRepulsionLJ.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/16/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryRepulsionLJ.h"
#include "Bead.h"

double BoundaryRepulsionLJ::computeEnergy(Bead* b, double r, double k_rep, double screenLength)
{
    assert(r != 0 && "Boundary repulsion cannot be calculated, distance from boundary is zero.");
    
    double inv_r4 = 1/r * 1/r * 1/r * 1/r;
    return k_rep * inv_r4;
}

void BoundaryRepulsionLJ::computeForces(Bead* b, double r, vector<double>& norm, double k_rep, double screenLength){
    
    assert(r != 0 && "Boundary repulsion cannot be calculated, distance from boundary is zero.");
    
    double inv_r4 = 1/r * 1/r * 1/r * 1/r;
    
    b->force[0] +=  k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[0];
    b->force[1] +=  k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[1];
    b->force[2] +=  k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[2];
}

void BoundaryRepulsionLJ::computeForcesAux(Bead* b, double r, vector<double>& norm,  double k_rep, double screenLength){
    
    assert(r != 0 && "Boundary repulsion cannot be calculated, distance from boundary is zero.");
    
    double inv_r4 = 1/r * 1/r * 1/r * 1/r;
    
    b->forceAux[0] +=  k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[0];
    b->forceAux[1] +=  k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[1];
    b->forceAux[2] +=  k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[2];

}