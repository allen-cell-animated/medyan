
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

#include "BubbleCylinderRepulsionExp.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

floatingpoint BubbleCylinderRepulsionExp::energy(Bead* b1, Bead* b2, floatingpoint radius,
                                          floatingpoint kRep, floatingpoint screenLength) {
    
    floatingpoint dist = twoPointDistance(b1->coordinate, b2->coordinate);
    
    floatingpoint effd = dist - radius;
    
    floatingpoint R = -effd / screenLength;
    return kRep * exp(R);
}

floatingpoint BubbleCylinderRepulsionExp::energy(Bead* b1, Bead* b2, floatingpoint radius,
                                          floatingpoint kRep, floatingpoint screenLength, floatingpoint d) {
    
    floatingpoint dist = twoPointDistanceStretched(b1->coordinate, b1->force,
                                            b2->coordinate, b2->force, d);
    floatingpoint effd = dist - radius;
    
    floatingpoint R = -effd / screenLength;
    return kRep * exp(R);
}

void BubbleCylinderRepulsionExp::forces(Bead* b1, Bead* b2, floatingpoint radius,
                                        floatingpoint kRep, floatingpoint screenLength) {
    
    //get dist
    floatingpoint dist = twoPointDistance(b1->coordinate, b2->coordinate);
    
    floatingpoint effd = dist - radius;
    
    floatingpoint R = -effd / screenLength;
    floatingpoint f0 = kRep * exp(R) / screenLength;
    
    //get norm
    auto norm = normalizeVector(twoPointDirection(b1->coordinate, b2->coordinate));
    
    b1->force[0] += - f0 *norm[0];
    b1->force[1] += - f0 *norm[1];
    b1->force[2] += - f0 *norm[2];
    
    b2->force[0] += f0 *norm[0];
    b2->force[1] += f0 *norm[1];
    b2->force[2] += f0 *norm[2];
}

void BubbleCylinderRepulsionExp::forcesAux(Bead* b1, Bead* b2, floatingpoint radius,
                                           floatingpoint kRep, floatingpoint screenLength) {
    
    //get dist
    floatingpoint dist = twoPointDistance(b1->coordinate, b2->coordinate);
    
    floatingpoint effd = dist - radius;
    
    floatingpoint R = -effd / screenLength;
    floatingpoint f0 = kRep * exp(R) / screenLength;
    
    //get norm
    auto norm = normalizeVector(twoPointDirection(b1->coordinate, b2->coordinate));
    
    b1->force[0] += - f0 *norm[0];
    b1->force[1] += - f0 *norm[1];
    b1->force[2] += - f0 *norm[2];
    
    b2->force[0] += f0 *norm[0];
    b2->force[1] += f0 *norm[1];
    b2->force[2] += f0 *norm[2];
    
}

floatingpoint BubbleCylinderRepulsionExp::loadForces(Bead* b1, Bead* b2, floatingpoint radius,
                                              floatingpoint kRep, floatingpoint screenLength) {
    
    //get dist
    floatingpoint dist = twoPointDistance(b1->coordinate, b2->coordinate);
    
    floatingpoint effd = dist - radius;
    
    floatingpoint R = -effd / screenLength;
    return kRep * exp(R) / screenLength;
}
