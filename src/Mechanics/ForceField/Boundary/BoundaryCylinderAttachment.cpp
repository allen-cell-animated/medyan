
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
    kattr = new floatingpoint[Bead::getPinnedBeads().size()];
    
    int i = 0;
    for(auto b : Bead::getPinnedBeads()) {

        beadSet[n * i] = b->getIndex();
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
floatingpoint BoundaryCylinderAttachment<BAttachmentInteractionType>::computeEnergy(floatingpoint *coord) {

    return _FFType.energy(coord, beadSet, kattr, pins);

}

template <class BAttachmentInteractionType>
void BoundaryCylinderAttachment<BAttachmentInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {

    _FFType.forces(coord, f, beadSet, kattr, pins);
}

///Template specializations
template floatingpoint BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::computeEnergy(floatingpoint *coord);
template void BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::computeForces(floatingpoint *coord, floatingpoint *f);
template void BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::vectorize();
template void BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::deallocate();
