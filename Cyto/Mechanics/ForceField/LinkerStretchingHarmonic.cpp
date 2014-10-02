//
//  LinkerStretchingHarmonic.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "LinkerStretchingHarmonic.h"

#include "MathFunctions.h"
#include "Bead.h"

using namespace std;
using namespace mathfunc;


// Energy calculation methods:
double LinkerStretchingHarmonic::Energy(Bead* pb1, Bead* pb2, Bead* pb3, Bead* pb4, double position1, double position2, double kStr, double L) {
    
    auto v1 = MidPointCoordinate(pb1->coordinate, pb2->coordinate, position1);
    auto v2 = MidPointCoordinate(pb3->coordinate, pb4->coordinate, position2);
    
    double dist = TwoPointDistance(v1, v2) - L;
    
    return 0.5 * kStr * dist * dist;
}

double LinkerStretchingHarmonic::Energy(Bead* pb1, Bead* pb2, Bead* pb3, Bead* pb4, double position1, double position2, double kStr, double L, double d ){
    
    auto v1 = MidPointCoordinateStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, position1, d);
    auto v2 = MidPointCoordinateStretched(pb3->coordinate, pb3->force, pb4->coordinate, pb4->force, position2, d);
    
    double distStretched = TwoPointDistance(v1, v2) - L;
    
    return 0.5 * kStr * distStretched * distStretched ;

}
// Force calculation methods:
void LinkerStretchingHarmonic::Forces(Bead* pb1, Bead* pb2, Bead* pb3, Bead* pb4, double position1, double position2, double kStr, double L ){
    
    auto v1 = MidPointCoordinate(pb1->coordinate, pb2->coordinate, position1);
    auto v2 = MidPointCoordinate(pb3->coordinate, pb4->coordinate, position2);
    
    double dist = TwoPointDistance( v1, v2);
    
    double invL = 1 / dist;
    
    double f0 = kStr * ( dist - L ) * invL;
    
    
    //force on i
    
    pb1->force[0] +=   -f0 * ( v1[0] - v2[0] ) * (1 - position1);
    
    pb1->force[1] +=   -f0 * ( v1[1] - v2[1] ) * (1 - position1);
    
    pb1->force[2] +=   -f0 * ( v1[2] - v2[2] ) * (1 - position1);
    
    
    // force i+1
    
    
    pb2->force[0] +=   -f0 * ( v1[0] - v2[0] ) * (position1);
    
    pb2->force[1] +=   -f0 * ( v1[1] - v2[1] ) * (position1);
    
    pb2->force[2] +=   -f0 * ( v1[2] - v2[2] ) * (position1);
    
    
    //force on j
    
    pb3->force[0] +=   f0 * ( v1[0] - v2[0] ) * (1 - position2);
    
    pb3->force[1] +=   f0 * ( v1[1] - v2[1] ) * (1 - position2);
    
    pb3->force[2] +=   f0 * ( v1[2] - v2[2] ) * (1 - position2);
    
    
    // force j+1
    
    
    pb4->force[0] +=   f0 * ( v1[0] - v2[0] ) * (position2);
    
    pb4->force[1] +=   f0 * ( v1[1] - v2[1] ) * (position2);
    
    pb4->force[2] +=   f0 * ( v1[2] - v2[2] ) * (position2);
    
    
}

void LinkerStretchingHarmonic::ForcesAux(Bead* pb1, Bead* pb2, Bead* pb3, Bead* pb4, double position1, double position2, double kStr, double L ){
    
    auto v1 = MidPointCoordinate(pb1->coordinate, pb2->coordinate, position1);
    auto v2 = MidPointCoordinate(pb3->coordinate, pb4->coordinate, position2);
    
    double dist = TwoPointDistance( v1, v2);
    
    double invL = 1 / dist;
    
    double f0 = kStr * ( dist - L ) * invL;
    
    
    //force on i
    
    pb1->forceAux[0] +=   -f0 * ( v1[0] - v2[0] ) * (1 - position1);
    
    pb1->forceAux[1] +=   -f0 * ( v1[1] - v2[1] ) * (1 - position1);
    
    pb1->forceAux[2] +=   -f0 * ( v1[2] - v2[2] ) * (1 - position1);
    
    
    // force i+1
    
    
    pb2->forceAux[0] +=   -f0 * ( v1[0] - v2[0] ) * (position1);
    
    pb2->forceAux[1] +=   -f0 * ( v1[1] - v2[1] ) * (position1);
    
    pb2->forceAux[2] +=   -f0 * ( v1[2] - v2[2] ) * (position1);
    
    
    //force on j
    
    pb3->forceAux[0] +=   f0 * ( v1[0] - v2[0] ) * (1 - position2);
    
    pb3->forceAux[1] +=   f0 * ( v1[1] - v2[1] ) * (1 - position2);
    
    pb3->forceAux[2] +=   f0 * ( v1[2] - v2[2] ) * (1 - position2);
    
    
    // force j+1
    
    
    pb4->forceAux[0] +=   f0 * ( v1[0] - v2[0] ) * (position2);
    
    pb4->forceAux[1] +=   f0 * ( v1[1] - v2[1] ) * (position2);
    
    pb4->forceAux[2] +=   f0 * ( v1[2] - v2[2] ) * (position2);
    
}

