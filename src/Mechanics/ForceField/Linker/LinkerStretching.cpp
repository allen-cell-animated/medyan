
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

#include "LinkerStretching.h"

#include "LinkerStretchingHarmonic.h"

#include "Cylinder.h"
#include "Linker.h"
#include "Bead.h"

template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::vectorizeInteractions() {
    
    int i = 0;
    for(auto b: Bead::getBeads()) {
        b->_dbIndex = i;
        i++;
    }
    i = 0;
    
    for (auto l: Linker::getLinkers()) {
        
        beadSet[n * i] = l->getFirstCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = l->getFirstCylinder()->getSecondBead()->_dbIndex;
        beadSet[n * i + 2] = l->getSecondCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 3] = l->getSecondCylinder()->getSecondBead()->_dbIndex;
        
        kstr[i] = l->getMLinker()->getStretchingConstant();
        eql[i] = l->getMLinker()->getEqLength();
        pos1[i] = l->getFirstPosition();
        pos2[i] = l->getSecondPosition();
        
        i++;
    }
}



template <class LStretchingInteractionType>
double LinkerStretching<LStretchingInteractionType>::computeEnergy(double* coord, double *f, double d){
    
    double U_i;
    
    if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, kstr, eql, pos1, pos2);
    else
        U_i = _FFType.energy(coord, f, beadSet, kstr, eql, pos1, pos2, d);
    
    return U_i;
}

template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::computeForces(double *coord, double *f) {
    
    _FFType.forces(coord, f, beadSet, kstr, eql, pos1, pos2);
}


///Temlate specializations
template double LinkerStretching<LinkerStretchingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void LinkerStretching<LinkerStretchingHarmonic>::computeForces(double *coord, double *f);
