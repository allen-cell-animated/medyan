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
    
    if (d == 0.0){
        
        for (auto pb: pbe->beads()){
           U+= _FFType.ComputeEnergy(pb, pbe->distance(pb->coordinates));
        }
    }
    

    
    else {
        for (auto pb: pbe->beads()){
            U+=_FFType.ComputeEnergy(pb, pbe->distance(pb->coordinates), d);
        }
    }

    return U;
}

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::ComputeForces(Filament* pf)
{
    for (auto pb: pbe->beads()){
        _FFType.ComputeForces(pb, pbe->distance(pb->coordinates), pbe->normal);
    }
}




template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::ComputeForcesAux(Filament* pf) /// Needed for Conjugated Gradient minimization;
{
    for (auto pb: pbe->beads()){
        _FFType.ComputeForcesAuX(pb, pbe->distance(pb->coordinates), pbe->normal);
    }
}

///Template specializations
template double BoundaryRepulsion<BoundaryRepulsionLJ>::ComputeEnergy(BoundaryElement* pbe, double d);
template void BoundaryRepulsion<BoundaryRepulsionLJ>::ComputeForces(BoundaryElement* pbe);
template void BoundaryRepulsion<BoundaryRepulsionLJ>::ComputeForcesAux(BoundaryElement* pbe);