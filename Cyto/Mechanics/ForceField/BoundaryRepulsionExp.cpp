//
//  BoundaryRepulsionExp.cpp
//  Cyto
//
//  Created by Konstantin Popov on 10/1/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryRepulsionExp.h"
#include "Bead.h"
#include <math.h>

double BoundaryRepulsionExp::ComputeEnergy(Bead* pb, double r, double k_rep, double screenLength)
{
    double R = -r/screenLength;
    return k_rep * exp(R);
    
//    std::cout << "r = " << r << std::endl;
}

void BoundaryRepulsionExp::ComputeForces(Bead* pb, double r, std::vector<double>& norm, double k_rep, double screenLength){
    
    double R = -r/screenLength;
    double f0 = k_rep * exp(R)/screenLength;
    
    pb->force[0] +=  f0 *norm[0];
    pb->force[1] +=  f0 *norm[1];
    pb->force[2] +=  f0 *norm[2];
    
//    std::cout << "r = " << r << std::endl;
    
}

void BoundaryRepulsionExp::ComputeForcesAux(Bead* pb, double r, std::vector<double>& norm,  double k_rep, double screenLength){
    
    double R = -r/screenLength;
    double f0 = k_rep * exp(R)/screenLength;
    
    pb->forceAux[0] +=  f0 *norm[0];
    pb->forceAux[1] +=  f0 *norm[1];
    pb->forceAux[2] +=  f0 *norm[2];
    
//    std::cout << "r = " << r << std::endl;
    
}