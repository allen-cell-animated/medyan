
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

#include "BoundaryBubbleRepulsionExp.h"

#include "Bead.h"

double BoundaryBubbleRepulsionExp::energy(Bead* b, double r, double r0,
                                          double kRep, double screenLength) {
    double R = -(r - r0) / screenLength;
    return kRep * exp(R);
}

void BoundaryBubbleRepulsionExp::forces(Bead* b, double r, double r0,
                                        vector<double>& norm, double kRep,
                                        double screenLength) {
    
    double R = -(r - r0) / screenLength;
    double f0 = kRep * exp(R) / screenLength;
    
    b->force[0] += f0 *norm[0];
    b->force[1] += f0 *norm[1];
    b->force[2] += f0 *norm[2];
    
}

void BoundaryBubbleRepulsionExp::forcesAux(Bead* b, double r, double r0,
                                           vector<double>& norm, double kRep,
                                           double screenLength) {
    
    double R = -(r - r0) / screenLength;
    double f0 = kRep * exp(R) / screenLength;
    
    b->forceAux[0] += f0 *norm[0];
    b->forceAux[1] += f0 *norm[1];
    b->forceAux[2] += f0 *norm[2];
    
}