
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BoundaryCylinderRepulsionExp.h"

#include "Bead.h"

double BoundaryCylinderRepulsionExp::energy(Bead* b, double r,
                                            double kRep, double screenLength) {
    double R = -r/screenLength;
    return kRep * exp(R);
}

void BoundaryCylinderRepulsionExp::forces(Bead* b, double r, vector<double>& norm,
                                          double kRep, double screenLength) {
    
    double R = -r/screenLength;
    double f0 = kRep * exp(R)/screenLength;
    
    b->force[0] += f0 *norm[0];
    b->force[1] += f0 *norm[1];
    b->force[2] += f0 *norm[2];
    
}

void BoundaryCylinderRepulsionExp::forcesAux(Bead* b, double r, vector<double>& norm,
                                             double kRep, double screenLength) {
    
    double R = -r/screenLength;
    double f0 = kRep * exp(R)/screenLength;
    
    b->forceAux[0] += f0 *norm[0];
    b->forceAux[1] += f0 *norm[1];
    b->forceAux[2] += f0 *norm[2];
    
}

double BoundaryCylinderRepulsionExp::loadForces(double r, double kRep, double screenLength) {
    
    double R = -r/screenLength;
    return kRep * exp(R)/screenLength;

}
