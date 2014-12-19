
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

#include "MotorGhostStretching.h"

#include "MotorGhostStretchingHarmonic.h"

#include "MotorGhost.h"
#include "Cylinder.h"
#include "Bead.h"

template <class MStretchingInteractionType>
double MotorGhostStretching<MStretchingInteractionType>::computeEnergy(
                                                MotorGhost* m, double d) {
    
    Bead* b1 = m->getFirstCylinder()->getFirstBead();
    Bead* b2 = m->getFirstCylinder()->getSecondBead();
    Bead* b3 = m->getSecondCylinder()->getFirstBead();
    Bead* b4 = m->getSecondCylinder()->getSecondBead();
    double kStretch = m->getMMotorGhost()->getStretchingConstant();
    double eqLength = m->getMMotorGhost()->getEqLength();
    
    double pos1 = m->getFirstPosition();
    double pos2 = m->getSecondPosition();
    
    if (d == 0.0)
        return _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
    else
        return _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength, d);
    
}

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::computeForces(MotorGhost* m) {
    
    Bead* b1 = m->getFirstCylinder()->getFirstBead();
    Bead* b2 = m->getFirstCylinder()->getSecondBead();
    Bead* b3 = m->getSecondCylinder()->getFirstBead();
    Bead* b4 = m->getSecondCylinder()->getSecondBead();
    double kStretch = m->getMMotorGhost()->getStretchingConstant();
    double eqLength = m->getMMotorGhost()->getEqLength();
    
    double pos1 = m->getFirstPosition();
    double pos2 = m->getSecondPosition();
    
    _FFType.forces(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
    
}


template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::computeForcesAux(MotorGhost* m) {

    Bead* b1 = m->getFirstCylinder()->getFirstBead();
    Bead* b2 = m->getFirstCylinder()->getSecondBead();
    Bead* b3 = m->getSecondCylinder()->getFirstBead();
    Bead* b4 = m->getSecondCylinder()->getSecondBead();
    double kStretch = m->getMMotorGhost()->getStretchingConstant();
    double eqLength = m->getMMotorGhost()->getEqLength();
    
    double pos1 = m->getFirstPosition();
    double pos2 = m->getSecondPosition();
    
    _FFType.forcesAux(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
}

///Template specializations
template double
MotorGhostStretching<MotorGhostStretchingHarmonic>::computeEnergy(MotorGhost* m, double d);
template void
MotorGhostStretching<MotorGhostStretchingHarmonic>::computeForces(MotorGhost* m);
template void
MotorGhostStretching<MotorGhostStretchingHarmonic>::computeForcesAux(MotorGhost* m);

