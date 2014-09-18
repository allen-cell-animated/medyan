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

double FilamentStretchingHarmonic::Energy(Bead* pb1, Bead* pb2, double k_str, double L){

    return 0.5 * k_str* ( TwoPointDistance( pb1->coordinate, pb2->coordinate) - L ) * ( TwoPointDistance( pb1->coordinate, pb2->coordinate) - L );

}


double FilamentStretchingHarmonic::Energy(Bead* pb1, Bead* pb2, double k_str, double L, double d){
   
    
    return 0.5 * k_str * ( TwoPointDistanceStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, d) - L ) * ( TwoPointDistanceStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, d) - L );

}

void FilamentStretchingHarmonic::Forces(Bead* pb1, Bead* pb2, double k_str, double L ){
    
    
    double invL = 1/TwoPointDistance( pb1->coordinate, pb2->coordinate);
    
    double f0 = k_str * ( TwoPointDistance( pb1->coordinate, pb2->coordinate) - L ) * invL;
    
    
    //force on i
    
    pb2->force[0] +=   -f0 * ( pb1->coordinate[0] - pb2->coordinate[0] );
    
    pb2->force[1] +=   -f0 * ( pb1->coordinate[1] - pb2->coordinate[1] );
    
    pb2->force[2] +=   -f0 * ( pb1->coordinate[2] - pb2->coordinate[2] );
    
    
    // force i-1
    
    
    pb1->force[0] +=   -f0 * ( pb2->coordinate[0] - pb1->coordinate[0] );
    
    pb1->force[1] +=   -f0 * ( pb2->coordinate[1] - pb1->coordinate[1] );
    
    pb1->force[2] +=  -f0 * ( pb2->coordinate[2] - pb1->coordinate[2] );
   
  
    
}

void FilamentStretchingHarmonic::ForcesAux(Bead* pb1, Bead* pb2, double k_str, double L ){
    
    
    double invL = 1/TwoPointDistance( pb1->coordinate, pb2->coordinate);
    
    double f0 = k_str * ( TwoPointDistance( pb1->coordinate, pb2->coordinate) - L ) * invL;
    
    
    pb2->forceAux[0] +=   -f0 * ( pb1->coordinate[0] - pb2->coordinate[0] );
    
    pb2->forceAux[1] +=   -f0 * ( pb1->coordinate[1] - pb2->coordinate[1] );
    
    pb2->forceAux[2] +=   -f0 * ( pb1->coordinate[2] - pb2->coordinate[2] );
    
    
    // force i-1
    
    
    pb1->forceAux[0] +=   -f0 * ( pb2->coordinate[0] - pb1->coordinate[0] );
    
    pb1->forceAux[1] +=   -f0 * ( pb2->coordinate[1] - pb1->coordinate[1] );
    
    pb1->forceAux[2] +=  -f0 * ( pb2->coordinate[2] - pb1->coordinate[2] );
    

}

