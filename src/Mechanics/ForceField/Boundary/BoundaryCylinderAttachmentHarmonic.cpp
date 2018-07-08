
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

#include "BoundaryCylinderAttachmentHarmonic.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

double BoundaryCylinderAttachmentHarmonic::energy(Bead* b, double kAttr) {
    
    auto x0 = b->pinnedPosition[0] - 2000;
    auto y0 = b->pinnedPosition[1] - 2000;
    auto r0 = sqrt(x0 * x0 + y0 * y0);
    
    vector<double> pinboundary{0,0,0};
    pinboundary[0] = x0 * 2000 / r0 + 2000;
    pinboundary[1] = x0 * 2000 / r0 + 2000;
    pinboundary[2] = b->pinnedPosition[2];
    
    auto eqL = twoPointDistance(b->pinnedPosition,pinboundary);
    double dist2 = twoPointDistance(b->coordinate, pinboundary);
    
    //double dist = twoPointDistance(b->coordinate, b->pinnedPosition);
    //return 0.5 * kAttr * dist * dist;
    return 0.5 * kAttr * (dist2 - eqL) * (dist2 - eqL);
}

double BoundaryCylinderAttachmentHarmonic::energy(Bead* b, double kAttr, double d) {
    
    vector<double> zeros{0,0,0};
    
    auto x0 = b->pinnedPosition[0] - 2000;
    auto y0 = b->pinnedPosition[1] - 2000;
    auto r0 = sqrt(x0 * x0 + y0 * y0);
    
    vector<double> pinboundary{0,0,0};
    pinboundary[0] = x0 * 2000 / r0 + 2000;
    pinboundary[1] = x0 * 2000 / r0 + 2000;
    pinboundary[2] = b->pinnedPosition[2];
    
    auto eqL = twoPointDistance(b->pinnedPosition,pinboundary);
    double dist2 = twoPointDistanceStretched(b->coordinate, b->force, pinboundary, zeros, d);
    
    //double dist = twoPointDistanceStretched(b->coordinate, b->force, b->pinnedPosition, zeros, d);
    //return 0.5 * kAttr * dist * dist;
    return 0.5 * kAttr * (dist2 - eqL) * (dist2 - eqL);
}

void BoundaryCylinderAttachmentHarmonic::forces(Bead* b, double kAttr) {
    
    
    double dist = twoPointDistance(b->coordinate, b->pinnedPosition);
    if(areEqual(dist, 0.0)) return;
    
    
    auto x0 = b->pinnedPosition[0] - 2000;
    auto y0 = b->pinnedPosition[1] - 2000;
    auto r0 = sqrt(x0 * x0 + y0 * y0);
    
    vector<double> pinboundary{0,0,0};
    pinboundary[0] = x0 * 2000 / r0 + 2000;
    pinboundary[1] = x0 * 2000 / r0 + 2000;
    pinboundary[2] = b->pinnedPosition[2];
    
    auto eqL = twoPointDistance(b->pinnedPosition,pinboundary);
    double dist2 = twoPointDistance(b->coordinate, pinboundary);
    
    auto dir = normalizedVector(twoPointDirection(b->coordinate, pinboundary));
    double f0 = kAttr * (dist2 - eqL);
    //auto dir = normalizedVector(twoPointDirection(b->coordinate, b->pinnedPosition));
    //double f0 = kAttr * dist;
    
    
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
    
    auto x0 = b->pinnedPosition[0] - 2000;
    auto y0 = b->pinnedPosition[1] - 2000;
    auto r0 = sqrt(x0 * x0 + y0 * y0);
    
    vector<double> pinboundary{0,0,0};
    pinboundary[0] = x0 * 2000 / r0 + 2000;
    pinboundary[1] = x0 * 2000 / r0 + 2000;
    pinboundary[2] = b->pinnedPosition[2];
    
    auto eqL = twoPointDistance(b->pinnedPosition,pinboundary);
    double dist2 = twoPointDistance(b->coordinate, pinboundary);

    auto dir = normalizedVector(twoPointDirection(b->coordinate, pinboundary));
    double f0 = kAttr * (dist2 - eqL);
    
    //auto dir = normalizedVector(twoPointDirection(b->coordinate, b->pinnedPosition));
    //double f0 = kAttr * dist;
    
    b->forceAux[0] += f0 * dir[0];
    b->forceAux[1] += f0 * dir[1];
    b->forceAux[2] += f0 * dir[2];
}
