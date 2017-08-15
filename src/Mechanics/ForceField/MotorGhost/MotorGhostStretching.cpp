
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

#include "MotorGhostStretching.h"

#include "MotorGhostStretchingHarmonic.h"

#include "MotorGhost.h"
#include "Cylinder.h"
#include "Bead.h"

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::vectorize() {
    
    beadSet = new int[n * MotorGhost::getMotorGhosts().size()];
    kstr = new double[MotorGhost::getMotorGhosts().size()];
    eql = new double[MotorGhost::getMotorGhosts().size()];
    pos1 = new double[MotorGhost::getMotorGhosts().size()];
    pos2 = new double[MotorGhost::getMotorGhosts().size()];
    
    int i = 0;
    
    for (auto m: MotorGhost::getMotorGhosts()) {
        
        beadSet[n * i] = m->getFirstCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = m->getFirstCylinder()->getSecondBead()->_dbIndex;
        beadSet[n * i + 2] = m->getSecondCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 3] = m->getSecondCylinder()->getSecondBead()->_dbIndex;
        
        kstr[i] = m->getMMotorGhost()->getStretchingConstant();
        eql[i] = m->getMMotorGhost()->getEqLength();
        pos1[i] = m->getFirstPosition();
        pos2[i] = m->getSecondPosition();
        
        i++;
    }
}

template<class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::deallocate() {
    
    delete beadSet;
    delete kstr;
    delete eql;
    delete pos1;
    delete pos2;
}


template <class MStretchingInteractionType>
double MotorGhostStretching<MStretchingInteractionType>::computeEnergy(double* coord, double *f, double d){
    
    double U_i;
    
    if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, kstr, eql, pos1, pos2);
    else
        U_i = _FFType.energy(coord, f, beadSet, kstr, eql, pos1, pos2, d);
    
    return U_i;
}

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::computeForces(double *coord, double *f) {
    
    _FFType.forces(coord, f, beadSet, kstr, eql, pos1, pos2);
}


///Temlate specializations
template double MotorGhostStretching<MotorGhostStretchingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::computeForces(double *coord, double *f);
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::vectorize();
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::deallocate();


