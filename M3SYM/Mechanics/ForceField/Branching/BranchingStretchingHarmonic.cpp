
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "BranchingStretchingHarmonic.h"

#include "MathFunctions.h"
#include "Bead.h"

using namespace mathfunc;

double BranchingStretchingHarmonic::energy(Bead* b1, Bead* b2, Bead* b3,
                                           double position, double kStretch,
                                           double eqLength){
    
    auto v1 = midPointCoordinate(b1->coordinate, b2->coordinate, position);
    
    double dist = twoPointDistance(v1, b3->coordinate) - eqLength;
    
    double energy =  0.5 * kStretch * dist * dist ;

    return energy;
}

double BranchingStretchingHarmonic::energy(Bead* b1, Bead* b2, Bead* b3,
                                           double position, double kStretch,
                                           double eqLength, double d){
    
    auto v1 = midPointCoordinateStretched(b1->coordinate, b1->force,
                                          b2->coordinate, b2->force, position, d);
    
    double dist = twoPointDistance(v1, b3->coordinate) - eqLength;
    
    double energy =  0.5 * kStretch * dist * dist ;

    return energy;
}

void BranchingStretchingHarmonic::forces(Bead* b1, Bead* b2, Bead* b3,
                                         double position, double kStretch,
                                         double eqLength){
    
    auto v1 = midPointCoordinate(b1->coordinate, b2->coordinate, position);
    
    
    double dist = twoPointDistance( v1, b3->coordinate);
    
    double invL = 1 / dist;
    double f0 = kStretch * ( dist - eqLength ) * invL;
    
    //force on i
    b1->force[0] +=  -f0 * ( b3->coordinate[0] - v1[0] ) * (position - 1);
    b1->force[1] +=  -f0 * ( b3->coordinate[1] - v1[1] ) * (position - 1);
    b1->force[2] +=  -f0 * ( b3->coordinate[2] - v1[2] ) * (position - 1);
    
    // force i+1
    b2->force[0] +=  f0 * ( b3->coordinate[0] - v1[0] ) * position;
    b2->force[1] +=  f0 * ( b3->coordinate[1] - v1[1] ) * position;
    b2->force[2] +=  f0 * ( b3->coordinate[2] - v1[2] ) * position;
    
    //force on j
    b3->force[0] +=  -f0 * ( b3->coordinate[0] - v1[0] );
    b3->force[1] +=  -f0 * ( b3->coordinate[1] - v1[1] );
    b3->force[2] +=  -f0 * ( b3->coordinate[2] - v1[2] );
    
}

void BranchingStretchingHarmonic::forcesAux(Bead* b1, Bead* b2, Bead* b3,
                                            double position, double kStretch,
                                            double eqLength){
    
    auto v1 = midPointCoordinate(b1->coordinate, b2->coordinate, position);
    
    double dist = twoPointDistance( v1, b3->coordinate);
    
    double invL = 1 / dist;
    double f0 = kStretch * ( dist - eqLength ) * invL;
    
    //force on i
    b1->forceAux[0] +=  -f0 * ( b3->coordinate[0] - v1[0] ) * (position - 1);
    b1->forceAux[1] +=  -f0 * ( b3->coordinate[1] - v1[1] ) * (position - 1);
    b1->forceAux[2] +=  -f0 * ( b3->coordinate[2] - v1[2] ) * (position - 1);
    
    // force i+1
    b2->forceAux[0] +=  f0 * ( b3->coordinate[0] - v1[0] ) * position;
    b2->forceAux[1] +=  f0 * ( b3->coordinate[1] - v1[1] ) * position;
    b2->forceAux[2] +=  f0 * ( b3->coordinate[2] - v1[2] ) * position;
    
    //force on j
    b3->forceAux[0] +=  -f0 * ( b3->coordinate[0] - v1[0] );
    b3->forceAux[1] +=  -f0 * ( b3->coordinate[1] - v1[1] );
    b3->forceAux[2] +=  -f0 * ( b3->coordinate[2] - v1[2] );
    
}
