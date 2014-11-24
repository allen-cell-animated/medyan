//
//  FilamentStretchingHarmonic.cpp
//  Cyto
//
//  Created by Konstantin Popov on 8/19/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "FilamentStretchingHarmonic.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

double FilamentStretchingHarmonic::energy(Bead* b1, Bead* b2, double k_str, double L){

    double dist = TwoPointDistance( b1->coordinate, b2->coordinate) - L;
    return 0.5 * k_str* dist * dist;
    
}


double FilamentStretchingHarmonic::energy(Bead* b1, Bead* b2, double k_str, double L, double d){

    double distStretched = TwoPointDistanceStretched(b1->coordinate, b1->force, b2->coordinate, b2->force, d) - L;
    return 0.5 * k_str * distStretched * distStretched;
    

}

void FilamentStretchingHarmonic::forces(Bead* b1, Bead* b2, double k_str, double L ){
    
    double dist = TwoPointDistance( b1->coordinate, b2->coordinate);
    double invL = 1 / dist;
    
    double f0 = k_str * ( dist - L ) * invL;
    
    //force on i
    b2->force[0] +=   f0 * ( b1->coordinate[0] - b2->coordinate[0] );
    b2->force[1] +=   f0 * ( b1->coordinate[1] - b2->coordinate[1] );
    b2->force[2] +=   f0 * ( b1->coordinate[2] - b2->coordinate[2] );
    
    // force i-1
    b1->force[0] +=   f0 * ( b2->coordinate[0] - b1->coordinate[0] );
    b1->force[1] +=   f0 * ( b2->coordinate[1] - b1->coordinate[1] );
    b1->force[2] +=  f0 * ( b2->coordinate[2] - b1->coordinate[2] );
   
  
    
}

void FilamentStretchingHarmonic::forcesAux(Bead* b1, Bead* b2, double k_str, double L ){
    
    double dist = TwoPointDistance( b1->coordinateAux, b2->coordinateAux);
    double invL = 1 / dist;
    double f0 = k_str * ( dist - L ) * invL;
    
    
    b2->forceAux[0] +=   f0 * ( b1->coordinateAux[0] - b2->coordinateAux[0] );
    b2->forceAux[1] +=   f0 * ( b1->coordinateAux[1] - b2->coordinateAux[1] );
    b2->forceAux[2] +=   f0 * ( b1->coordinateAux[2] - b2->coordinateAux[2] );
    
    
    // force i-1
    b1->forceAux[0] +=   f0 * ( b2->coordinateAux[0] - b1->coordinateAux[0] );
    b1->forceAux[1] +=   f0 * ( b2->coordinateAux[1] - b1->coordinateAux[1] );
    b1->forceAux[2] +=  f0 * ( b2->coordinateAux[2] - b1->coordinateAux[2] );
    

}

