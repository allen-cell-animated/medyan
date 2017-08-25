
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

#include "FilamentStretching.h"

#include "FilamentStretchingHarmonic.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "cross_check.h"
template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::vectorize() {
    
    beadSet = new int[n * Cylinder::getCylinders().size()];
    kstr = new double[Cylinder::getCylinders().size()];
    eql = new double[Cylinder::getCylinders().size()];
    
    int i = 0;
    
    for (auto c: Cylinder::getCylinders()) {
        beadSet[n * i] = c->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = c->getSecondBead()->_dbIndex;
//        std::cout<<beadSet[n*i]<<" "<<beadSet[n*i+1]<<endl;
        kstr[i] = c->getMCylinder()->getStretchingConst();
        eql[i] = c->getMCylinder()->getEqLength();
        
        i++;
    }
}

template<class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::deallocate() {
    
    delete beadSet;
    delete kstr;
    delete eql;
}


template <class FStretchingInteractionType>
double FilamentStretching<FStretchingInteractionType>::computeEnergy(double* coord, double *f, double d){
    
    double U_i;
    
    if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, kstr, eql);
    else
        U_i = _FFType.energy(coord, f, beadSet, kstr, eql, d);
    
    return U_i;
}

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::computeForces(double *coord, double *f) {
    
    _FFType.forces(coord, f, beadSet, kstr, eql);
#ifdef CROSSCHECK
    for (auto f: Filament::getFilaments()) {
        
        for(auto it : f->getCylinderVector()){
            
            Bead* b1 = it->getFirstBead();
            Bead* b2 = it->getSecondBead();
            double kStretch =it->getMCylinder()->getStretchingConst();
            double eqLength = it->getMCylinder()->getEqLength();
            
            _FFType.forcesAux(b1, b2, kStretch, eqLength);
        }
    }
    auto state=cross_check::crosscheckforces(f);
    std::cout<<"F S YES "<<state<<endl;
#endif

}


///Temlate specializations
template double FilamentStretching<FilamentStretchingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void FilamentStretching<FilamentStretchingHarmonic>::computeForces(double *coord, double *f);
template void FilamentStretching<FilamentStretchingHarmonic>::vectorize();
template void FilamentStretching<FilamentStretchingHarmonic>::deallocate();

