
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

#include "BubbleBubbleRepulsionExp.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

floatingpoint BubbleBubbleRepulsionExp::energy(Bead* b1, Bead* b2, floatingpoint r1, floatingpoint r2,
                                        floatingpoint kRep, floatingpoint screenLength) {
    
    floatingpoint dist = twoPointDistance(b1->vcoordinate(), b2->vcoordinate());
    
    floatingpoint effd = dist - r1 - r2;
    
    floatingpoint R = -effd / screenLength;
    return kRep * exp(R);
}

floatingpoint BubbleBubbleRepulsionExp::energy(Bead* b1, Bead* b2, floatingpoint r1, floatingpoint r2,
                                        floatingpoint kRep, floatingpoint screenLength, floatingpoint d) {
    
    floatingpoint dist = twoPointDistanceStretched(b1->vcoordinate(), b1->vforce(),
                                            b2->vcoordinate(), b2->vforce(), d);
    floatingpoint effd = dist - r1 - r2;
    
    floatingpoint R = -effd / screenLength;
    return kRep * exp(R);
}

void BubbleBubbleRepulsionExp::forces(Bead* b1, Bead* b2, floatingpoint r1, floatingpoint r2,
                                      floatingpoint kRep, floatingpoint screenLength) {
    
    //get dist
    floatingpoint dist = twoPointDistance(b1->vcoordinate(), b2->vcoordinate());
    
    floatingpoint effd = dist - r1 - r2;
    
    floatingpoint R = -effd / screenLength;
    floatingpoint f0 = kRep * exp(R) / screenLength;
    
    //get norm
    auto norm = normalizeVector(twoPointDirection(b1->vcoordinate(), b2->vcoordinate()));

    b1->force()[0] += - f0 *norm[0];
    b1->force()[1] += - f0 *norm[1];
    b1->force()[2] += - f0 *norm[2];
    
    b2->force()[0] += f0 *norm[0];
    b2->force()[1] += f0 *norm[1];
    b2->force()[2] += f0 *norm[2];
}

void BubbleBubbleRepulsionExp::forcesAux(Bead* b1, Bead* b2, floatingpoint r1, floatingpoint r2,
                                         floatingpoint kRep, floatingpoint screenLength) {
    
    //get dist
    floatingpoint dist = twoPointDistance(b1->vcoordinate(), b2->vcoordinate());
    
    floatingpoint effd = dist - r1 - r2;
    
    floatingpoint R = -effd / screenLength;
    floatingpoint f0 = kRep * exp(R) / screenLength;
    
    //get norm
    auto norm = normalizeVector(twoPointDirection(b1->vcoordinate(), b2->vcoordinate()));
    
    b1->force()[0] += - f0 *norm[0];
    b1->force()[1] += - f0 *norm[1];
    b1->force()[2] += - f0 *norm[2];
    
    b2->force()[0] += f0 *norm[0];
    b2->force()[1] += f0 *norm[1];
    b2->force()[2] += f0 *norm[2];
    
}
