//
//  BoundaryRepulsionLJ.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/16/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryRepulsionLJ.h"
#include "Bead.h"

double BoundaryRepulsionLJ::ComputeEnergy(Bead* pb, double r, double k_rep)
{
    assert(r != 0 && "Boundary repulsion cannot be calculated, distance from boundary is zero.");
    
    double inv_r4 = 1/r;// * 1/r * 1/r * 1/r;
    return k_rep * inv_r4;
    
    std::cout << "r = " << r << std::endl;
}

void BoundaryRepulsionLJ::ComputeForces(Bead* pb, double r, std::vector<double>& norm, double k_rep){
    
    assert(r != 0 && "Boundary repulsion cannot be calculated, distance from boundary is zero.");
    
    double inv_r4 = 1/r;// * 1/r * 1/r * 1/r;
    
    //pb->force[0] +=  k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[0];
    
    //pb->force[1] +=  k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[1];
    
    //pb->force[2] +=  k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[2];
    
    pb->force[0] +=  k_rep * inv_r4 * inv_r4 * inv_r4 *norm[0];
    pb->force[1] +=  k_rep *inv_r4* inv_r4 *inv_r4 *norm[1];
    pb->force[2] +=  k_rep *inv_r4* inv_r4 *inv_r4 *norm[2];
    
    std::cout << "r = " << r << std::endl;
    
}

void BoundaryRepulsionLJ::ComputeForcesAux(Bead* pb, double r, std::vector<double>& norm,  double k_rep){
    
    assert(r != 0 && "Boundary repulsion cannot be calculated, distance from boundary is zero.");
    
    double inv_r4 = 1/r;// * 1/r * 1/r * 1/r;
    
    //pb->forceAux[0] +=  k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[0];
    
    //pb->forceAux[1] +=  k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[1];
    
    //pb->forceAux[2] +=  k_rep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[2];
    
    pb->forceAux[0] +=  k_rep * inv_r4 * inv_r4 *inv_r4 *norm[0];
    pb->forceAux[1] +=  k_rep *inv_r4* inv_r4 *inv_r4 *norm[1];
    pb->forceAux[2] +=  k_rep *inv_r4* inv_r4 *inv_r4 *norm[2];
    
    std::cout << "r = " << r << std::endl;
    
}