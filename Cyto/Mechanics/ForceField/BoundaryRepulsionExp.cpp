//
//  BoundaryRepulsionExp.cpp
//  Cyto
//
//  Created by Konstantin Popov on 10/1/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryRepulsionExp.h"

#include "Bead.h"

double BoundaryRepulsionExp::computeEnergy(Bead* b, double r, double k_rep, double screenLength)
{
    double R = -r/screenLength;
    return k_rep * exp(R);
}

void BoundaryRepulsionExp::computeForces(Bead* b, double r, vector<double>& norm, double k_rep, double screenLength){
    
    double R = -r/screenLength;
    double f0 = k_rep * exp(R)/screenLength;
    
    b->force[0] += f0 *norm[0];
    b->force[1] += f0 *norm[1];
    b->force[2] += f0 *norm[2];
    
}

void BoundaryRepulsionExp::computeForcesAux(Bead* b, double r, vector<double>& norm,  double k_rep, double screenLength){
    
    double R = -r/screenLength;
    double f0 = k_rep * exp(R)/screenLength;
    
    b->forceAux[0] += f0 *norm[0];
    b->forceAux[1] += f0 *norm[1];
    b->forceAux[2] += f0 *norm[2];
    
}