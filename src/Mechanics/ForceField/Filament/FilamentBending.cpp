
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------
#include "FilamentBending.h"

#include "FilamentBendingHarmonic.h"
#include "FilamentBendingCosine.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "cross_check.h"
template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::vectorize() {
    
    int numInteractions = Bead::getBeads().size() - 2 * Filament::getFilaments().size();
    
    beadSet = new int[n * numInteractions];
    kbend = new double[numInteractions];
    eqt = new double[numInteractions];
    
    int i = 0;
    
    for (auto f: Filament::getFilaments()) {
        
        if (f->getCylinderVector().size() > 1){
        
            for (auto it = f->getCylinderVector().begin()+1;
                     it != f->getCylinderVector().end(); it++){
        
                auto it2 = it - 1;
                beadSet[n * i] = (*it2)->getFirstBead()->_dbIndex;
                beadSet[n * i + 1] = (*it)->getFirstBead()->_dbIndex;;
                beadSet[n * i + 2] = (*it)->getSecondBead()->_dbIndex;;
                
                kbend[i] = (*it)->getMCylinder()->getBendingConst();
                eqt[i]  = (*it)->getMCylinder()->getEqTheta();
                
                i++;
            }
        }
    }
}

template<class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::deallocate() {
    
    delete beadSet;
    delete kbend;
    delete eqt;
}


template <class FBendingInteractionType>
double FilamentBending<FBendingInteractionType>::computeEnergy(double *coord, double *f, double d){
    
    double U_i;
    
    if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, kbend, eqt);
    else
        U_i = _FFType.energy(coord, f, beadSet, kbend, eqt, d);
    
    return U_i;
    
}

template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::computeForces(double *coord, double *f) {
    
    _FFType.forces(coord, f, beadSet, kbend, eqt);
#ifdef CROSSCHECK
    for (auto f: Filament::getFilaments()) {
        
        if (f->getCylinderVector().size()>1){
            for (auto it = f->getCylinderVector().begin()+1;
                 it != f->getCylinderVector().end(); it++){
                
                auto it2 = it - 1;
                Bead* b1 = (*it2)->getFirstBead();
                Bead* b2 = (*it)->getFirstBead();
                Bead* b3 = (*it)->getSecondBead();
                double kBend = (*it)->getMCylinder()->getBendingConst();
                double eqTheta = (*it)->getMCylinder()->getEqTheta();
                
                _FFType.forces(b1, b2, b3, kBend, eqTheta);
            }
        }
    }
    auto state=cross_check::crosscheckforces(f);
    std::cout<<"F S+B YES "<<state<<endl;
#endif
}

///Template specializations
template double FilamentBending<FilamentBendingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void FilamentBending<FilamentBendingHarmonic>::computeForces(double *coord, double *f);
template void FilamentBending<FilamentBendingHarmonic>::vectorize();
template void FilamentBending<FilamentBendingHarmonic>::deallocate();


template double FilamentBending<FilamentBendingCosine>::computeEnergy(double *coord, double *f, double d);
template void FilamentBending<FilamentBendingCosine>::computeForces(double *coord, double *f);
template void FilamentBending<FilamentBendingCosine>::vectorize();
template void FilamentBending<FilamentBendingCosine>::deallocate();
