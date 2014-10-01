//
//  BoundaryRepulsionMixed.cpp
//  Cyto
//
//  Created by James Komianos on 9/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryRepulsionMixed.h"
#include "Bead.h"

double BoundaryRepulsionMixed::ComputeEnergy(Bead* pb, double r, double k_rep)
{
//    double r0;
//    
//    if (r > r0){
//        double inv_r4 = 1/r;// * 1/r * 1/r * 1/r;
//        return k_rep * inv_r4;
//    }
//    
//    else if (r >= 0 && r <= r0){
//        return 3*k_rep;
//    }
//    
//    if (r < 0){
//        
//        return 0.5* k_rep * r*r;
//    }
//    
    
    assert(r != 0 && "Boundary repulsion cannot be calculated, distance from boundary is zero.");
    
    double inv_r4 = 1/r;// * 1/r * 1/r * 1/r;
    return k_rep * inv_r4;
    
    std::cout << "r = " << r << std::endl;
}

void BoundaryRepulsionMixed::ComputeForces(Bead* pb, double r, std::vector<double>& norm, double k_rep){
    
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

void BoundaryRepulsionMixed::ComputeForcesAux(Bead* pb, double r, std::vector<double>& norm,  double k_rep){
    
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