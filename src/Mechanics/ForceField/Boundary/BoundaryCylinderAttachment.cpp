
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BoundaryCylinderAttachment.h"

#include "BoundaryCylinderAttachmentHarmonic.h"
#include "BoundaryElement.h"

#include "Bubble.h"
#include "Bead.h"

template <class BAttachmentInteractionType>
void BoundaryCylinderAttachment<BAttachmentInteractionType>::vectorize() {
    
    //first coord in beadset is bead, then pin position
    beadSet = new int[Bead::getPinnedBeads().size()];
    kattr = new double[Bead::getPinnedBeads().size()];
    
    int i = 0;
    for(auto b : Bead::getPinnedBeads()) {

        beadSet[n * i] = b->getDbIndex();
        kattr[n * i] = SysParams::Mechanics().pinK;

        pins[3 * (n * i)] = b->getPinPosition()[0];
        pins[3 * (n * i) + 1] = b->getPinPosition()[1];
        pins[3 * (n * i) + 2] = b->getPinPosition()[2];

        i++;
    }
}

template <class BAttachmentInteractionType>
void BoundaryCylinderAttachment<BAttachmentInteractionType>::deallocate() {

    delete beadSet;
    delete kattr;
    delete pins;
}

template <class BAttachmentInteractionType>
double BoundaryCylinderAttachment<BAttachmentInteractionType>::computeEnergy(double *coord, double *f, double d) {


    double U;

    if (d == 0.0)
        U =  _FFType.energy(coord, f, beadSet, kattr, pins);
    else
        U =  _FFType.energy(coord, f, beadSet, kattr, pins, d);

    return U;
}

template <class BAttachmentInteractionType>
void BoundaryCylinderAttachment<BAttachmentInteractionType>::computeForces(double *coord, double *f) {

    _FFType.forces(coord, f, beadSet, kattr, pins);
}

///Template specializations
template double BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::computeEnergy(double *coord, double *f, double d);
template void BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::computeForces(double *coord, double *f);
template void BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::vectorize();
template void BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::deallocate();
