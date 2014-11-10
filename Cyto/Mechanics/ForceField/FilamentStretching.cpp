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
double FilamentStretching<FStretchingInteractionType>::computeEnergy(Filament* f, double d) {
    double U = 0.0;
    
    if (d == 0.0){
        for(auto it : f->getCylinderVector()){
            
            Bead* b1 = it->getFirstBead();
            Bead* b2 = it->getSecondBead();
            double kStr = it->getMCylinder()->getStretchingConst();
            double L = it->getMCylinder()->getEqLength();
            U += _FFType.energy(b1, b2, kStr, L);
        }
    }
    else {
        for(auto it : f->getCylinderVector()){
            Bead* b1 = it->getFirstBead();
            Bead* b2 = it->getSecondBead();
            double kStr =it->getMCylinder()->getStretchingConst();
            double L = it->getMCylinder()->getEqLength();
            
            U += _FFType.energy(b1, b2, kStr, L, d);   ///This type of function needed for conjugated gradient minimisation only;
        }
    }
    
    //cout << "Stretching Energy = " << U << endl;
    
    return U;
}

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::computeForces(Filament* f) {
   for(auto it : f->getCylinderVector()){
       
       Bead* b1 = it->getFirstBead();
       Bead* b2 = it->getSecondBead();
       
       double kStr =it->getMCylinder()->getStretchingConst();
       double L = it->getMCylinder()->getEqLength();
       _FFType.forces(b1, b2, kStr, L);
   }
}


template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::computeForcesAux(Filament* f) {/// Needed for Conjugated Gradient minimization;
    for(auto it : f->getCylinderVector()){
        
        Bead* b1 = it->getFirstBead();
        Bead* b2 = it->getSecondBead();
        double kStr =it->getMCylinder()->getStretchingConst();
        double L = it->getMCylinder()->getEqLength();
        
        _FFType.forcesAux(b1, b2, kStr, L);
    }
}

///Template specializations
template double FilamentStretching<FilamentStretchingHarmonic>::computeEnergy(Filament* f, double d);
template void  FilamentStretching<FilamentStretchingHarmonic>::computeForces(Filament* f);
template void  FilamentStretching<FilamentStretchingHarmonic>::computeForcesAux(Filament* f);
