
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

#include "CylinderExclVolume.h"

#include "CylinderExclVolRepulsion.h"

#include "Cylinder.h"
#include "Bead.h"

template <class CVolumeInteractionType>
double CylinderExclVolume<CVolumeInteractionType>::computeEnergy(
                             Cylinder* c1, Cylinder* c2, double d) {
    
    Bead* b1 = c1->getFirstBead();
    Bead* b2 = c1->getSecondBead();
    Bead* b3 = c2->getFirstBead();
    Bead* b4 = c2->getSecondBead();
    double kRepuls = c1->getMCylinder()->getExVolConst();
    
    
    if (d == 0.0)
        return _FFType.energy(b1, b2, b3, b4, kRepuls);
    else
        return _FFType.energy(b1, b2, b3, b4, kRepuls, d);
    
}

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::computeForces(
                                     Cylinder* c1, Cylinder* c2) {

    Bead* b1 = c1->getFirstBead();
    Bead* b2 = c1->getSecondBead();
    Bead* b3 = c2->getFirstBead();
    Bead* b4 = c2->getSecondBead();
    double kRepuls = c1->getMCylinder()->getExVolConst();

    
    _FFType.forces(b1, b2, b3, b4, kRepuls );
    
}


template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::computeForcesAux(
                                        Cylinder* c1, Cylinder* c2) {
    Bead* b1 = c1->getFirstBead();
    Bead* b2 = c1->getSecondBead();
    Bead* b3 = c2->getFirstBead();
    Bead* b4 = c2->getSecondBead();
    double kRepuls = c1->getMCylinder()->getExVolConst();
    
    
    _FFType.forcesAux(b1, b2, b3, b4, kRepuls );
}

///Template specializations
template double
CylinderExclVolume<CylinderExclVolRepulsion>::computeEnergy(Cylinder* c1, Cylinder* c2, double d);
template void
CylinderExclVolume<CylinderExclVolRepulsion>::computeForces(Cylinder* c1, Cylinder* c2);
template void
CylinderExclVolume<CylinderExclVolRepulsion>::computeForcesAux(Cylinder* c1, Cylinder* c2);

