
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "BoundaryRepulsionExp.h"

#include "Bead.h"

#include "SysParams.h"

double BoundaryRepulsionExp::computeEnergy(Bead* b, double r,
                                           double kRep, double screenLength) {
    //ceiling to avoid blowups
    double R = -r/screenLength;
    return min(kRep * exp(R), SysParams::Boundaries().BCeiling * kRep);
}

void BoundaryRepulsionExp::computeForces(Bead* b, double r, vector<double>& norm,
                                         double kRep, double screenLength) {
    
    //ceiling to avoid blowups
    double R = -r/screenLength;
    double f0 = min(kRep * exp(R)/screenLength, kRep * SysParams::Boundaries().BCeiling  / screenLength);

    //update the load force of the bead
    b->loadForce = f0;
    
    b->force[0] += f0 *norm[0];
    b->force[1] += f0 *norm[1];
    b->force[2] += f0 *norm[2];
    
}

void BoundaryRepulsionExp::computeForcesAux(Bead* b, double r, vector<double>& norm,
                                            double kRep, double screenLength) {
    
    //ceiling to avoid blowups
    double R = -r/screenLength;
    double f0 = min(kRep * exp(R)/screenLength, kRep * SysParams::Boundaries().BCeiling  / screenLength);
    
    //update the load force of the bead
    b->loadForce = f0;
    
    b->forceAux[0] += f0 *norm[0];
    b->forceAux[1] += f0 *norm[1];
    b->forceAux[2] += f0 *norm[2];
    
}