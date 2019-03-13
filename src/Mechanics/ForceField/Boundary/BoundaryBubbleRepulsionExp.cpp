
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BoundaryBubbleRepulsionExp.h"

#include "Bead.h"

floatingpoint BoundaryBubbleRepulsionExp::energy(Bead* b, floatingpoint r, floatingpoint r0,
                                          floatingpoint kRep, floatingpoint screenLength) {
    floatingpoint R = -(r - r0) / screenLength;
    return kRep * exp(R);
}

void BoundaryBubbleRepulsionExp::forces(Bead* b, floatingpoint r, floatingpoint r0,
                                        vector<floatingpoint>& norm, floatingpoint kRep,
                                        floatingpoint screenLength) {
    
    floatingpoint R = -(r - r0) / screenLength;
    floatingpoint f0 = kRep * exp(R) / screenLength;
    
    b->force[0] += f0 *norm[0];
    b->force[1] += f0 *norm[1];
    b->force[2] += f0 *norm[2];
    
}

void BoundaryBubbleRepulsionExp::forcesAux(Bead* b, floatingpoint r, floatingpoint r0,
                                           vector<floatingpoint>& norm, floatingpoint kRep,
                                           floatingpoint screenLength) {
    
    floatingpoint R = -(r - r0) / screenLength;
    floatingpoint f0 = kRep * exp(R) / screenLength;
    
    b->forceAux[0] += f0 *norm[0];
    b->forceAux[1] += f0 *norm[1];
    b->forceAux[2] += f0 *norm[2];
    
}
