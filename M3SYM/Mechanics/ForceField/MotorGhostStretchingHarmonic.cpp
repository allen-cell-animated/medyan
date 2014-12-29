
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "MotorGhostStretchingHarmonic.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

double MotorGhostStretchingHarmonic::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                            double position1, double position2,
                                            double kStretch, double eqLength){
    
    auto v1 = midPointCoordinate(b1->coordinate, b2->coordinate, position1);
    auto v2 = midPointCoordinate(b3->coordinate, b4->coordinate, position2);
    
    double dist = twoPointDistance(v1, v2) - eqLength;
    return 0.5 * kStretch * dist * dist ;
}

double MotorGhostStretchingHarmonic::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                            double position1, double position2,
                                            double kStretch, double eqLength, double d){
    
    auto v1 = midPointCoordinateStretched(b1->coordinate, b1->force,
                                          b2->coordinate, b2->force, position1, d);
    auto v2 = midPointCoordinateStretched(b3->coordinate, b3->force,
                                          b4->coordinate, b4->force, position2, d);
    
    double dist = twoPointDistance(v1, v2) - eqLength;
    return 0.5 * kStretch * dist * dist;
}
// Force calculation methods:
double MotorGhostStretchingHarmonic::forces(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
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
    
    return f0;
}

double MotorGhostStretchingHarmonic::forcesAux(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                               double position1, double position2,
                                               double kStretch, double eqLength){
    
    auto v1 = midPointCoordinate(b1->coordinateAux, b2->coordinateAux, position1);
    auto v2 = midPointCoordinate(b3->coordinateAux, b4->coordinateAux, position2);
    
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
    
    return f0;
}