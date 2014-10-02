//
//  FilamentStretching.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "FilamentStretching.h"
#include "FilamentStretchingHarmonic.h"
#include "Filament.h"

template <class FStretchingInteractionType>
double FilamentStretching<FStretchingInteractionType>::ComputeEnergy(Filament* pf, double d) {
    double U = 0.0;
    
    if (d == 0.0){
        for(auto it : pf->getCylinderVector()){
            
            Bead* pb1 = it->GetFirstBead();
            Bead* pb2 = it->GetSecondBead();
            double kStr = it->getMCylinder()->GetStretchingConst();
            double L = it->getMCylinder()->GetEqLength();
            U += _FFType.Energy(pb1, pb2, kStr, L);
        }
    }
    else {
        for(auto it : pf->getCylinderVector()){
            Bead* pb1 = it->GetFirstBead();
            Bead* pb2 = it->GetSecondBead();
            double kStr =it->getMCylinder()->GetStretchingConst();
            double L = it->getMCylinder()->GetEqLength();
            
            U += _FFType.Energy(pb1, pb2, kStr, L, d);   ///This type of function needed for conjugated gradient minimisation only;
        }
    }
    return U;
}

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::ComputeForces(Filament* pf) {
   for(auto it : pf->getCylinderVector()){
       
       Bead* pb1 = it->GetFirstBead();
       Bead* pb2 = it->GetSecondBead();
       
       double kStr =it->getMCylinder()->GetStretchingConst();
       double L = it->getMCylinder()->GetEqLength();
       _FFType.Forces(pb1, pb2, kStr, L);
   }
}


template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::ComputeForcesAux(Filament* pf) {/// Needed for Conjugated Gradient minimization;
    for(auto it : pf->getCylinderVector()){
        
        Bead* pb1 = it->GetFirstBead();
        Bead* pb2 = it->GetSecondBead();
        double kStr =it->getMCylinder()->GetStretchingConst();
        double L = it->getMCylinder()->GetEqLength();
        
        _FFType.ForcesAux(pb1, pb2, kStr, L);
    }
}

///Template specializations
template double FilamentStretching<FilamentStretchingHarmonic>::ComputeEnergy(Filament* pf, double d);
template void  FilamentStretching<FilamentStretchingHarmonic>::ComputeForces(Filament* pf);
template void  FilamentStretching<FilamentStretchingHarmonic>::ComputeForcesAux(Filament* pf);
