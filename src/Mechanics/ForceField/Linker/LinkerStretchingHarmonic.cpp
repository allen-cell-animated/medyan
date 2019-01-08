
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

#include "LinkerStretchingHarmonic.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

double LinkerStretchingHarmonic::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                        double position1, double position2,
                                        double kStretch, double eqLength, bool stretched) {

    const auto& c1 = stretched ? b1->getCoordinate<true>() : b1->getCoordinate<false>();
    const auto& c2 = stretched ? b2->getCoordinate<true>() : b2->getCoordinate<false>();
    const auto& c3 = stretched ? b3->getCoordinate<true>() : b3->getCoordinate<false>();
    const auto& c4 = stretched ? b4->getCoordinate<true>() : b4->getCoordinate<false>();

    auto v1 = midPointCoordinate(c1, c2, position1);
    auto v2 = midPointCoordinate(c3, c4, position2);
    
    double dist = twoPointDistance(v1, v2) - eqLength;
    return 0.5 * kStretch * dist * dist;
}

double LinkerStretchingHarmonic::forces(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                        double position1, double position2,
                                        double kStretch, double eqLength){
    
    auto v1 = midPointCoordinate(b1->coordinate, b2->coordinate, position1);
    auto v2 = midPointCoordinate(b3->coordinate, b4->coordinate, position2);
    
    double dist = twoPointDistance( v1, v2);
    
    double invL = 1 / dist;
    
    double f0 = kStretch * ( dist - eqLength ) * invL;

    //force on i
    b1->force[0] +=   -f0 * ( v1[0] - v2[0] ) * (1 - position1);
    b1->force[1] +=   -f0 * ( v1[1] - v2[1] ) * (1 - position1);
    b1->force[2] +=   -f0 * ( v1[2] - v2[2] ) * (1 - position1);
    
    // force i+1
    b2->force[0] +=   -f0 * ( v1[0] - v2[0] ) * (position1);
    b2->force[1] +=   -f0 * ( v1[1] - v2[1] ) * (position1);
    b2->force[2] +=   -f0 * ( v1[2] - v2[2] ) * (position1);

    //force on j
    b3->force[0] +=   f0 * ( v1[0] - v2[0] ) * (1 - position2);
    b3->force[1] +=   f0 * ( v1[1] - v2[1] ) * (1 - position2);
    b3->force[2] +=   f0 * ( v1[2] - v2[2] ) * (1 - position2);

    // force j+1
    b4->force[0] +=   f0 * ( v1[0] - v2[0] ) * (position2);
    b4->force[1] +=   f0 * ( v1[1] - v2[1] ) * (position2);
    b4->force[2] +=   f0 * ( v1[2] - v2[2] ) * (position2);
    
    return f0 / invL;
}

double LinkerStretchingHarmonic::forcesAux(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                           double position1, double position2,
                                           double kStretch, double eqLength){
    
    auto v1 = midPointCoordinate(b1->coordinate, b2->coordinate, position1);
    auto v2 = midPointCoordinate(b3->coordinate, b4->coordinate, position2);
    
    double dist = twoPointDistance( v1, v2);
    double invL = 1 / dist;
    double f0 = kStretch * ( dist - eqLength ) * invL;
    
    //force on i
    b1->forceAux[0] +=   -f0 * ( v1[0] - v2[0] ) * (1 - position1);
    b1->forceAux[1] +=   -f0 * ( v1[1] - v2[1] ) * (1 - position1);
    b1->forceAux[2] +=   -f0 * ( v1[2] - v2[2] ) * (1 - position1);
    
    // force i+1
    b2->forceAux[0] +=   -f0 * ( v1[0] - v2[0] ) * (position1);
    b2->forceAux[1] +=   -f0 * ( v1[1] - v2[1] ) * (position1);
    b2->forceAux[2] +=   -f0 * ( v1[2] - v2[2] ) * (position1);
    
    //force on j
    b3->forceAux[0] +=   f0 * ( v1[0] - v2[0] ) * (1 - position2);
    b3->forceAux[1] +=   f0 * ( v1[1] - v2[1] ) * (1 - position2);
    b3->forceAux[2] +=   f0 * ( v1[2] - v2[2] ) * (1 - position2);

    // force j+1
    b4->forceAux[0] +=   f0 * ( v1[0] - v2[0] ) * (position2);
    b4->forceAux[1] +=   f0 * ( v1[1] - v2[1] ) * (position2);
    b4->forceAux[2] +=   f0 * ( v1[2] - v2[2] ) * (position2);
    
    return f0 / invL;
}

