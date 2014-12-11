
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

#include "FilamentStretching.h"

#include "FilamentStretchingHarmonic.h"
#include "Filament.h"
#include "Cylinder.h"

template <class FStretchingInteractionType>
double FilamentStretching<FStretchingInteractionType>::computeEnergy(Filament* f, double d) {
    double U = 0.0;
    
    if (d == 0.0){
        for(auto it : f->getCylinderVector()){
            
            Bead* b1 = it->getFirstBead();
            Bead* b2 = it->getSecondBead();
            double kStretch = it->getMCylinder()->getStretchingConst();
            double eqLength = it->getMCylinder()->getEqLength();

            U += _FFType.energy(b1, b2, kStretch, eqLength);
        }
    }
    else {
        for(auto it : f->getCylinderVector()){
            Bead* b1 = it->getFirstBead();
            Bead* b2 = it->getSecondBead();
            double kStretch =it->getMCylinder()->getStretchingConst();
            double eqLength = it->getMCylinder()->getEqLength();

            U += _FFType.energy(b1, b2, kStretch, eqLength, d);
        }
    }
    
    return U;
}

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::computeForces(Filament* f) {
   for(auto it : f->getCylinderVector()){
       
       Bead* b1 = it->getFirstBead();
       Bead* b2 = it->getSecondBead();
       double kStretch =it->getMCylinder()->getStretchingConst();
       double eqLength = it->getMCylinder()->getEqLength();
       
       _FFType.forces(b1, b2, kStretch, eqLength);
   }
}


template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::computeForcesAux(Filament* f) {
    for(auto it : f->getCylinderVector()){
        
        Bead* b1 = it->getFirstBead();
        Bead* b2 = it->getSecondBead();
        double kStretch =it->getMCylinder()->getStretchingConst();
        double eqLength = it->getMCylinder()->getEqLength();
        
        _FFType.forcesAux(b1, b2, kStretch, eqLength);
    }
}

///Template specializations
template double FilamentStretching<FilamentStretchingHarmonic>::computeEnergy(Filament* f, double d);
template void  FilamentStretching<FilamentStretchingHarmonic>::computeForces(Filament* f);
template void  FilamentStretching<FilamentStretchingHarmonic>::computeForcesAux(Filament* f);
