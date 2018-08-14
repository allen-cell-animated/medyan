
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BoundaryCylinderRepulsionExpIn.h"

#include "Bead.h"

double BoundaryCylinderRepulsionExpIn::energy(Bead* b, double r,
                                            double kRep, double screenLength) {
    //Boundary repulsion occurs 100 nm below the actual boundary
    double R = -r/screenLength + 100/screenLength;
    return kRep * exp(R);
}

void BoundaryCylinderRepulsionExpIn::forces(Bead* b, double r, vector<double>& norm,
                                          double kRep, double screenLength) {
    
    double R = -r/screenLength + 100/screenLength;
    double f0 = kRep * exp(R)/screenLength;
    
    b->force[0] += f0 *norm[0];
    b->force[1] += f0 *norm[1];
    b->force[2] += f0 *norm[2];
    // Qin add brforce
    b->brforce[0] = f0 *norm[0];
    b->brforce[1] = f0 *norm[1];
    b->brforce[2] = f0 *norm[2];
    
}

void BoundaryCylinderRepulsionExpIn::forcesAux(Bead* b, double r, vector<double>& norm,
                                             double kRep, double screenLength) {
    
    double R = -r/screenLength + 100/screenLength;
    double f0 = kRep * exp(R)/screenLength;
    
    b->forceAux[0] += f0 *norm[0];
    b->forceAux[1] += f0 *norm[1];
    b->forceAux[2] += f0 *norm[2];
    
}

double BoundaryCylinderRepulsionExpIn::loadForces(double r, double kRep, double screenLength) {
    
    double R = -r/screenLength + 100/screenLength;
    return kRep * exp(R)/screenLength;
    
}
