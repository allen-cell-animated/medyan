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

using namespace mathfunc;

double LinkerStretchingHarmonic::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                        double position1, double position2, double kStr, double L) {
    
    auto v1 = MidPointCoordinate(b1->coordinate, b2->coordinate, position1);
    auto v2 = MidPointCoordinate(b3->coordinate, b4->coordinate, position2);
    
    double dist = TwoPointDistance(v1, v2) - L;
    
    return 0.5 * kStr * dist * dist;
}

double LinkerStretchingHarmonic::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                        double position1, double position2, double kStr, double L, double d ){
    
    auto v1 = MidPointCoordinateStretched(b1->coordinate, b1->force, b2->coordinate, b2->force, position1, d);
    auto v2 = MidPointCoordinateStretched(b3->coordinate, b3->force, b4->coordinate, b4->force, position2, d);
    
    double distStretched = TwoPointDistance(v1, v2) - L;
    
    return 0.5 * kStr * distStretched * distStretched ;

}
void LinkerStretchingHarmonic::forces(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                      double position1, double position2, double kStr, double L ){
    
    auto v1 = MidPointCoordinate(b1->coordinate, b2->coordinate, position1);
    auto v2 = MidPointCoordinate(b3->coordinate, b4->coordinate, position2);
    
    double dist = TwoPointDistance( v1, v2);
    
    double invL = 1 / dist;
    
    double f0 = kStr * ( dist - L ) * invL;
    
    
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
}

void LinkerStretchingHarmonic::forcesAux(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                         double position1, double position2, double kStr, double L ){
    
    auto v1 = MidPointCoordinate(b1->coordinateAux, b2->coordinateAux, position1);
    auto v2 = MidPointCoordinate(b3->coordinateAux, b4->coordinateAux, position2);
    
    double dist = TwoPointDistance( v1, v2);
    
    double invL = 1 / dist;
    
    double f0 = kStr * ( dist - L ) * invL;
    
    
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
    
}

