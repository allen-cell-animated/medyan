//
//  MFilamentStretching.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MFilament.h"
#include "MBead.h"
#include "MCylinder.h"

using namespace std;
using namespace mathfunc;

double Filament::EnergyHarmonicStretching(Cylinder* pc){
    
    Bead* pb1 = pc->GetFirstBead();
    Bead* pb2 = pc->GetSecondBead();
    
    double u = 0.5 * pc->GetStretchingConst() * ( TwoPointDistance( pb1->coordinate, pb2->coordinate) - pc->GetEqLength() ) * ( TwoPointDistance( pb1->coordinate, pb2->coordinate) - pc->GetEqLength() );
    
    pb1 = NULL;
    pb2 = NULL;
    
    return u;
}


double Filament::EnergyHarmonicStretching(Cylinder* pc, double d){
    
    Bead* pb1 = pc->GetFirstBead();
    Bead* pb2 = pc->GetSecondBead();
    
    double u = 0.5 * pc->GetStretchingConst() * ( TwoPointDistanceStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, d) - pc->GetEqLength() ) * ( TwoPointDistanceStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, d) - pc->GetEqLength() );
    
    pb1 = NULL;
    pb2 = NULL;
    
    return u;
}

void Filament::ForceHarmonicStretching(Cylinder* pc ){
    
    Bead* pb1 = pc->GetFirstBead();
    Bead* pb2 = pc->GetSecondBead();
    
    double invL = 1/TwoPointDistance( pb1->coordinate, pb2->coordinate);
    
    double f0 = pc->GetStretchingConst() * ( TwoPointDistance( pb1->coordinate, pb2->coordinate) - pc->GetEqLength() ) * invL;
    
    cout<< "f0= "<<f0<<endl;
    
    //force on i
    
    pb2->force[0] +=   -f0 * ( pb1->coordinate[0] - pb2->coordinate[0] );
    
    pb2->force[1] +=   -f0 * ( pb1->coordinate[1] - pb2->coordinate[1] );
    
    pb2->force[2] +=   -f0 * ( pb1->coordinate[2] - pb2->coordinate[2] );
    
    
    // force i-1
    
    
    pb1->force[0] +=   -f0 * ( pb2->coordinate[0] - pb1->coordinate[0] );
    
    pb1->force[1] +=   -f0 * ( pb2->coordinate[1] - pb1->coordinate[1] );
    
    pb1->force[2] +=  -f0 * ( pb2->coordinate[2] - pb1->coordinate[2] );
    
    pb1 = NULL;
    pb2 = NULL;
    
}

void Filament::ForceHarmonicStretchingAux(Cylinder* pc ){
    
    
    
    Bead* pb1 = pc->GetFirstBead();
    Bead* pb2 = pc->GetSecondBead();
    
    double invL = 1/TwoPointDistance( pb1->coordinate, pb2->coordinate);
    
    double f0 = pc->GetStretchingConst() * ( TwoPointDistance( pb1->coordinate, pb2->coordinate) - pc->GetEqLength() ) * invL;
    
    cout<< "f0= "<<f0<<endl;
    
    pb2->forceAux[0] +=   -f0 * ( pb1->coordinate[0] - pb2->coordinate[0] );
    
    pb2->forceAux[1] +=   -f0 * ( pb1->coordinate[1] - pb2->coordinate[1] );
    
    pb2->forceAux[2] +=   -f0 * ( pb1->coordinate[2] - pb2->coordinate[2] );
    
    
    // force i-1
    
    
    pb1->forceAux[0] +=   -f0 * ( pb2->coordinate[0] - pb1->coordinate[0] );
    
    pb1->forceAux[1] +=   -f0 * ( pb2->coordinate[1] - pb1->coordinate[1] );
    
    pb1->forceAux[2] +=  -f0 * ( pb2->coordinate[2] - pb1->coordinate[2] );
    
    pb1 = NULL;
    pb2 = NULL;
}

