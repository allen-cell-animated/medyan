
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

#include "BoundaryRepulsionExp.h"

#include "Bead.h"

double BoundaryRepulsionExp::computeEnergy(Bead* b, double r,
                                           double kRep, double screenLength) {
    double R = -r/screenLength;
    return kRep * exp(R);
}

void BoundaryRepulsionExp::computeForces(Bead* b, double r, vector<double>& norm,
                                         double kRep, double screenLength) {
    
    double R = -r/screenLength;
    double f0 = kRep * exp(R)/screenLength;

    //update the load force of the bead
    b->loadForce = f0;
    
    b->force[0] += f0 *norm[0];
    b->force[1] += f0 *norm[1];
    b->force[2] += f0 *norm[2];
    
}

void BoundaryRepulsionExp::computeForcesAux(Bead* b, double r, vector<double>& norm,
                                            double kRep, double screenLength) {
    
    double R = -r/screenLength;
    double f0 = kRep * exp(R)/screenLength;
    
    //update the load force of the bead
    b->loadForce = f0;
    
    b->forceAux[0] += f0 *norm[0];
    b->forceAux[1] += f0 *norm[1];
    b->forceAux[2] += f0 *norm[2];
    
}