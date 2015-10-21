
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

#include "FilamentStretching.h"

#include "FilamentStretchingHarmonic.h"
#include "Filament.h"
#include "Cylinder.h"

template <class FStretchingInteractionType>
double FilamentStretching<FStretchingInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto f: Filament::getFilaments()) {
        
        U_i = 0;
        
        if (d == 0.0){
            for(auto it : f->getCylinderVector()){
                
                Bead* b1 = it->getFirstBead();
                Bead* b2 = it->getSecondBead();
                double kStretch = it->getMCylinder()->getStretchingConst();
                double eqLength = it->getMCylinder()->getEqLength();
                
                U_i += _FFType.energy(b1, b2, kStretch, eqLength);
            }
        }
        else {
            for(auto it : f->getCylinderVector()){
                Bead* b1 = it->getFirstBead();
                Bead* b2 = it->getSecondBead();
                double kStretch =it->getMCylinder()->getStretchingConst();
                double eqLength = it->getMCylinder()->getEqLength();
                
                U_i += _FFType.energy(b1, b2, kStretch, eqLength, d);
            }
        }
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            _filamentCulprit = f;
            
            return -1;
        }
        else
            U += U_i;
    }
    
    return U;
}

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::computeForces() {
    
    for (auto f: Filament::getFilaments()) {
    
        for(auto it : f->getCylinderVector()){
            
            Bead* b1 = it->getFirstBead();
            Bead* b2 = it->getSecondBead();
            double kStretch =it->getMCylinder()->getStretchingConst();
            double eqLength = it->getMCylinder()->getEqLength();
           
            _FFType.forces(b1, b2, kStretch, eqLength);
        }
    }
}


template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::computeForcesAux() {
    
    for (auto f: Filament::getFilaments()) {
        
        for(auto it : f->getCylinderVector()){
            
            Bead* b1 = it->getFirstBead();
            Bead* b2 = it->getSecondBead();
            double kStretch =it->getMCylinder()->getStretchingConst();
            double eqLength = it->getMCylinder()->getEqLength();
            
            _FFType.forcesAux(b1, b2, kStretch, eqLength);
        }
    }
}

///Template specializations
template double FilamentStretching<FilamentStretchingHarmonic>::computeEnergy(double d);
template void FilamentStretching<FilamentStretchingHarmonic>::computeForces();
template void FilamentStretching<FilamentStretchingHarmonic>::computeForcesAux();
