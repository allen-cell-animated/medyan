//
//  BranchingStretchingHarmonic.cpp
//  M3SYM
//
//  Created by Konstantin Popov on 12/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BranchingStretchingHarmonic.h"
#include "MathFunctions.h"
#include "Bead.h"

using namespace mathfunc;

double BranchingStretchingHarmonic::energy(Bead* b1, Bead* b2, Bead* b3,
                                            double position, double kStr, double L){
    
    auto v1 = MidPointCoordinate(b1->coordinate, b2->coordinate, position);
    
    
    double dist = TwoPointDistance(v1, b3->coordinate) - L;
    return 0.5 * kStr * dist * dist ;
}

double BranchingStretchingHarmonic::energy(Bead* b1, Bead* b2, Bead* b3,
                                            double position, double kStr, double L, double d ){
    
    auto v1 = MidPointCoordinateStretched(b1->coordinate, b1->force, b2->coordinate, b2->force, position, d);
    
    
    double dist = TwoPointDistance(v1, b3->coordinate) - L;
    return 0.5 * kStr * dist * dist;
}
// Force calculation methods:
void BranchingStretchingHarmonic::forces(Bead* b1, Bead* b2, Bead* b3,
                                          double position, double kStr, double L ){
    
    auto v1 = MidPointCoordinate(b1->coordinate, b2->coordinate, position);
    
    
    double dist = TwoPointDistance( v1, b3->coordinate);
    
    double invL = 1 / dist;
    double f0 = kStr * ( dist - L ) * invL;
    
    //force on i
    b1->force[0] +=   -f0 * ( b3->coordinate[0] - v1[0] ) * (position - 1);
    b1->force[1] +=   -f0 * ( b3->coordinate[1] - v1[1] ) * (position - 1);
    b1->force[2] +=   -f0 * ( b3->coordinate[2] - v1[2] ) * (position - 1);
    
    // force i+1
    b2->force[0] +=   f0 * ( b3->coordinate[0] - v1[0] ) * position;
    b2->force[1] +=   f0 * ( b3->coordinate[1] - v1[1] ) * position;
    b2->force[2] +=   f0 * ( b3->coordinate[2] - v1[2] ) * position;
    
    //force on j
    b3->force[0] +=   -f0 * ( b3->coordinate[0] - v1[0] );
    b3->force[1] +=   -f0 * ( b3->coordinate[1] - v1[1] );
    b3->force[2] +=   -f0 * ( b3->coordinate[2] - v1[2] );
    
   }
void BranchingStretchingHarmonic::forcesAux(Bead* b1, Bead* b2, Bead* b3,
                                         double position, double kStr, double L ){
    
    auto v1 = MidPointCoordinate(b1->coordinateAux, b2->coordinateAux, position);
    
    
    double dist = TwoPointDistance( v1, b3->coordinateAux);
    
    double invL = 1 / dist;
    double f0 = kStr * ( dist - L ) * invL;
    
    //force on i
    b1->forceAux[0] +=   -f0 * ( b3->coordinateAux[0] - v1[0] ) * (position - 1);
    b1->forceAux[1] +=   -f0 * ( b3->coordinateAux[1] - v1[1] ) * (position - 1);
    b1->forceAux[2] +=   -f0 * ( b3->coordinateAux[2] - v1[2] ) * (position - 1);
    
    // force i+1
    b2->forceAux[0] +=   f0 * ( b3->coordinateAux[0] - v1[0] ) * position;
    b2->forceAux[1] +=   f0 * ( b3->coordinateAux[1] - v1[1] ) * position;
    b2->forceAux[2] +=   f0 * ( b3->coordinateAux[2] - v1[2] ) * position;
    
    //force on j
    b3->forceAux[0] +=   -f0 * ( b3->coordinateAux[0] - v1[0] );
    b3->forceAux[1] +=   -f0 * ( b3->coordinateAux[1] - v1[1] );
    b3->forceAux[2] +=   -f0 * ( b3->coordinateAux[2] - v1[2] );
    
}
