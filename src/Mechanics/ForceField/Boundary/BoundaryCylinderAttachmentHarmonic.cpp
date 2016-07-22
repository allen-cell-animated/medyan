
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
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

double BoundaryCylinderAttachmentHarmonic::energy(Bead* b, double kAttr) {
    
    double dist = twoPointDistance(b->coordinate, b->pinnedPosition);
    return 0.5 * kAttr * dist * dist;
}

double BoundaryCylinderAttachmentHarmonic::energy(Bead* b, double kAttr, double d) {
    
    vector<double> zeros{0,0,0};
    
    double dist = twoPointDistanceStretched(b->coordinate, b->force, b->pinnedPosition, zeros, d);
    return 0.5 * kAttr * dist * dist;
}

void BoundaryCylinderAttachmentHarmonic::forces(Bead* b, double kAttr) {
    
    auto dir = normalizedVector(twoPointDirection(b->coordinate, b->pinnedPosition));
    
    double f0 = kAttr * twoPointDistance(b->coordinate, b->pinnedPosition);
    
    cout << "Boundary attachment force = " << f0 << endl;
    
    
    b->force[0] += f0 * dir[0];
    b->force[1] += f0 * dir[1];
    b->force[2] += f0 * dir[2];
}

void BoundaryCylinderAttachmentHarmonic::forcesAux(Bead* b, double kAttr) {
    
    auto dir = normalizedVector(twoPointDirection(b->coordinate, b->pinnedPosition));
    
    double f0 = kAttr * twoPointDistance(b->coordinate, b->pinnedPosition);
    
    b->forceAux[0] += f0 * dir[0];
    b->forceAux[1] += f0 * dir[1];
    b->forceAux[2] += f0 * dir[2];
}
