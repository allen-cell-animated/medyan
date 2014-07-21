//
//  MLinkerStretching.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MLinker.h"
#include "MBead.h"

using namespace std;
using namespace mathfunc;


double Linker::EnergyHarmonicStretching(Bead* pb1, Bead* pb2 ){
    
    double u = 0.5 * _kStretch * ( TwoPointDistance( pb1->coordinate, pb2->coordinate) - _eqLength ) * ( TwoPointDistance( pb1->coordinate, pb2->coordinate) - _eqLength );
    
    return u;
}


double Linker::EnergyHarmonicStretching(Bead* pb1, Bead* pb2, double d){
    
    double u = 0.5 * _kStretch * ( TwoPointDistanceStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, d) - _eqLength ) * ( TwoPointDistanceStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, d) - _eqLength );
    return u;
}


void Linker::ForceHarmonicStretching(Bead* pb1, Bead* pb2 ){
    
    
    double invL = 1/TwoPointDistance( pb1->coordinate, pb2->coordinate);
    
    double f0 = _kStretch * ( TwoPointDistance( pb1->coordinate, pb2->coordinate) - _eqLength ) * invL;
    
    cout<< "f0= "<<f0<<endl;
    
    //force on b2
    
    pb2->force[0] +=   -f0 * ( pb1->coordinate[0] - pb2->coordinate[0] );
    
    pb2->force[1] +=   -f0 * ( pb1->coordinate[1] - pb2->coordinate[1] );
    
    pb2->force[2] +=   -f0 * ( pb1->coordinate[2] - pb2->coordinate[2] );
    
    
    // force b1
    
    
    pb1->force[0] +=   -f0 * ( pb2->coordinate[0] - pb1->coordinate[0] );
    
    pb1->force[1] +=   -f0 * ( pb2->coordinate[1] - pb1->coordinate[1] );
    
    pb1->force[2] +=  -f0 * ( pb2->coordinate[2] - pb1->coordinate[2] );
    
}

void Linker::ForceHarmonicStretchingAux(Bead* pb1, Bead* pb2 ){
    
    
    double invL = 1/TwoPointDistance( pb1->coordinate, pb2->coordinate);
    
    double f0 = _kStretch * ( TwoPointDistance( pb1->coordinate, pb2->coordinate) - _eqLength) * invL;
    
    cout<< "f0= "<<f0<<endl;
    
    //force on b2
    
    pb2->forceAux[0] +=   -f0 * ( pb1->coordinate[0] - pb2->coordinate[0] );
    
    pb2->forceAux[1] +=   -f0 * ( pb1->coordinate[1] - pb2->coordinate[1] );
    
    pb2->forceAux[2] +=   -f0 * ( pb1->coordinate[2] - pb2->coordinate[2] );
    
    
    // force b1
    
    
    pb1->forceAux[0] +=   -f0 * ( pb2->coordinate[0] - pb1->coordinate[0] );
    
    pb1->forceAux[1] +=   -f0 * ( pb2->coordinate[1] - pb1->coordinate[1] );
    
    pb1->forceAux[2] +=  -f0 * ( pb2->coordinate[2] - pb1->coordinate[2] );
    
}
