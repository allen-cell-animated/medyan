
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

#include "BoundaryCylinderAttachmentHarmonic.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

double BoundaryCylinderAttachmentHarmonic::energy(Bead* b, double kAttr, bool stretched) {
    
    double dist = twoPointDistance(stretched ? b->getCoordinate<true>() : b->getCoordinate<false>(), b->pinnedPosition);
    return 0.5 * kAttr * dist * dist;
}

void BoundaryCylinderAttachmentHarmonic::forces(Bead* b, double kAttr) {
    
    
    double dist = twoPointDistance(b->coordinate, b->pinnedPosition);
    if(areEqual(dist, 0.0)) return;
    
    auto dir = normalizedVector(twoPointDirection(b->coordinate, b->pinnedPosition));
    double f0 = kAttr * dist;
    
    b->force[0] += f0 * dir[0];
    b->force[1] += f0 * dir[1];
    b->force[2] += f0 * dir[2];
    
    //Qin, add pinforce
    b->pinforce[0] += f0 * dir[0];
    b->pinforce[1] += f0 * dir[1];
    b->pinforce[2] += f0 * dir[2];
}

void BoundaryCylinderAttachmentHarmonic::forcesAux(Bead* b, double kAttr) {
    
    double dist = twoPointDistance(b->coordinate, b->pinnedPosition);
    if(areEqual(dist, 0.0)) return;
    
    auto dir = normalizedVector(twoPointDirection(b->coordinate, b->pinnedPosition));
    double f0 = kAttr * dist;
    
    b->forceAux[0] += f0 * dir[0];
    b->forceAux[1] += f0 * dir[1];
    b->forceAux[2] += f0 * dir[2];
}
