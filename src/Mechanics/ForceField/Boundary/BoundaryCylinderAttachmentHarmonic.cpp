
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
    
    double dist = twoPointDistance(b->coordinate, b->_pinnedPosition);
    return 0.5 * kAttr * dist * dist;
}

double BoundaryCylinderAttachmentHarmonic::energy(Bead* b, double kAttr, double d) {
    
    vector<double> zeros{0,0,0};
    
    double dist = twoPointDistanceStretched(b->coordinate, b->force, b->_pinnedPosition, zeros, d);
    return 0.5 * kAttr * dist * dist;
}

void BoundaryCylinderAttachmentHarmonic::forces(Bead* b, double kAttr) {
    
    
    
}

void BoundaryCylinderAttachmentHarmonic::forcesAux(Bead* b, double kAttr) {
    
    
}
