//
//  BoundaryRepulsion.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/12/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryRepulsion.h"
//#include "BoundaryRepulsionOneOverTwelve.h"

template <class BRepulsionInteractionType>
double BoundaryRepulsion<class BRepulsionInteractionType>::ComputeEnergy(BoundaryElement* pbe, double d)
{
    double U = 0.0;
    
    for (auto pbe-> )
    
    if (d == 0.0){
                }
    

    }
    else {
        
        }
    }
    return U;
}

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::ComputeForces(Filament* pf)
{
    for(auto it : pf->getCylinderVector()){
        
        Bead* pb1 = it->getMCylinder()->GetFirstBead();
        Bead* pb2 = it->getMCylinder()->GetSecondBead();
        
        
        
        double kStr =it->getMCylinder()->GetStretchingConst();
        double L = it->getMCylinder()->GetEqLength();
        _FFType.Forces(pb1, pb2, kStr, L);
    }
}


template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::ComputeForcesAux(Filament* pf) /// Needed for Conjugated Gradient minimization;
{
    for(auto it : pf->getCylinderVector()){
        
        Bead* pb1 = it->getMCylinder()->GetFirstBead();
        Bead* pb2 = it->getMCylinder()->GetSecondBead();
        double kStr =it->getMCylinder()->GetStretchingConst();
        double L = it->getMCylinder()->GetEqLength();
        
        _FFType.ForcesAux(pb1, pb2, kStr, L);
    }
}

///Template specializations
template double FilamentStretching<FilamentStretchingHarmonic>::ComputeEnergy(Filament* pf, double d);
template void  FilamentStretching<FilamentStretchingHarmonic>::ComputeForces(Filament* pf);
template void  FilamentStretching<FilamentStretchingHarmonic>::ComputeForcesAux(Filament* pf);