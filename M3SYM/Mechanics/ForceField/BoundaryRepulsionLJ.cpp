
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "BoundaryRepulsionLJ.h"

#include "Bead.h"

double BoundaryRepulsionLJ::computeEnergy(Bead* b, double r, double kRep,
                                          double screenLength)
{
    double inv_r4 = 1/r * 1/r * 1/r * 1/r;
    return kRep * inv_r4;
}

void BoundaryRepulsionLJ::computeForces(Bead* b, double r, vector<double>& norm,
                                        double kRep, double screenLength) {
    
    double inv_r4 = 1/r * 1/r * 1/r * 1/r;
    
    b->force[0] +=  kRep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[0];
    b->force[1] +=  kRep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[1];
    b->force[2] +=  kRep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[2];
}

void BoundaryRepulsionLJ::computeForcesAux(Bead* b, double r, vector<double>& norm,
                                           double kRep, double screenLength) {
    
    double inv_r4 = 1/r * 1/r * 1/r * 1/r;
    
    b->forceAux[0] +=  kRep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[0];
    b->forceAux[1] +=  kRep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[1];
    b->forceAux[2] +=  kRep * inv_r4 * inv_r4 * inv_r4 * 1/r *norm[2];

}