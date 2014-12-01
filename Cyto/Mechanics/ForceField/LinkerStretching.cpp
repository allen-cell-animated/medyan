
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

#include "LinkerStretching.h"

#include "LinkerStretchingHarmonic.h"

#include "Cylinder.h"
#include "Linker.h"
#include "Bead.h"

template <class LStretchingInteractionType>
double LinkerStretching<LStretchingInteractionType>::computeEnergy(Linker* l, double d) {
    
    Bead* b1 = l->getFirstCylinder()->getFirstBead();
    Bead* b2 = l->getFirstCylinder()->getSecondBead();
    Bead* b3 = l->getSecondCylinder()->getFirstBead();
    Bead* b4 = l->getSecondCylinder()->getSecondBead();
    double kStretch = l->getMLinker()->getStretchingConstant();
    double L = l->getMLinker()->getEqLength();
    double pos1 = l->getFirstPosition();
    double pos2 = l->getSecondPosition();
    
    if (d == 0.0)
        return _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, L);
    else
        return _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, L, d);
    
}

template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::computeForces(Linker* l) {
    
    Bead* b1 = l->getFirstCylinder()->getFirstBead();
    Bead* b2 = l->getFirstCylinder()->getSecondBead();
    Bead* b3 = l->getSecondCylinder()->getFirstBead();
    Bead* b4 = l->getSecondCylinder()->getSecondBead();
    double kStretch = l->getMLinker()->getStretchingConstant();
    double L = l->getMLinker()->getEqLength();
    
    double pos1 = l->getFirstPosition();
    double pos2 = l->getSecondPosition();
    
    _FFType.forces(b1, b2, b3, b4, pos1, pos2, kStretch, L);
}


template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::computeForcesAux(Linker* l) {

    Bead* b1 = l->getFirstCylinder()->getFirstBead();
    Bead* b2 = l->getFirstCylinder()->getSecondBead();
    Bead* b3 = l->getSecondCylinder()->getFirstBead();
    Bead* b4 = l->getSecondCylinder()->getSecondBead();
    double kStretch = l->getMLinker()->getStretchingConstant();
    double L = l->getMLinker()->getEqLength();
    
    double pos1 = l->getFirstPosition();
    double pos2 = l->getSecondPosition();
    
    _FFType.forcesAux(b1, b2, b3, b4, pos1, pos2, kStretch, L);
}

///Temlate specializations
template double LinkerStretching<LinkerStretchingHarmonic>::computeEnergy(Linker* l, double d);
template void  LinkerStretching<LinkerStretchingHarmonic>::computeForces(Linker* l);
template void  LinkerStretching<LinkerStretchingHarmonic>::computeForcesAux(Linker* l);
