
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

#include "AFMAttachment.h"

#include "AFMAttachmentHarmonic.h"

#include "AFM.h"
#include "Bubble.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

template <class AFMInteractionType>
void AFMAttachment<AFMInteractionType>::vectorize(const FFCoordinateStartingIndex& si) {

    numInteractions_ = 0;
    for(auto afm : AFM::getAFMs()) numInteractions_ += afm->getFilaments().size();

    beadSet_.resize(2 * numInteractions_);
    radii_.resize(numInteractions_);
    kstr_.resize(numInteractions_);

    int ci = 0;
    for(auto afm : AFM::getAFMs()) {

        for (int fIndex = 0; fIndex < afm->getFilaments().size(); fIndex++) {
            const auto f = afm->getFilaments()[fIndex];

            beadSet_[2*ci    ] = afm->getBubble()->getIndex() * 3 + si.bubble;
            beadSet_[2*ci + 1] = f->getMinusEndCylinder()->getFirstBead()->getIndex() * 3 + si.bead;

            radii_[ci] = afm->getBubble()->getRadius();

            kstr_[ci] = f->getMinusEndCylinder()->getMCylinder()->getStretchingConst();

            ci++;
        }
    }
}

template <class AFMInteractionType>
floatingpoint AFMAttachment<AFMInteractionType>::computeEnergy(floatingpoint *coord, bool stretched) {
    return _FFType.energy(coord, numInteractions_, beadSet_.data(), kstr_.data(), radii_.data());
}

template <class AFMInteractionType>
void AFMAttachment<AFMInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
    _FFType.forces(coord, f, numInteractions_, beadSet_.data(), kstr_.data(), radii_.data());
}


///Template specializations
template floatingpoint AFMAttachment<AFMAttachmentHarmonic>::computeEnergy(floatingpoint *coord, bool stretched);
template void AFMAttachment<AFMAttachmentHarmonic>::computeForces(floatingpoint *coord, floatingpoint *f);
//template void AFMAttachment<AFMAttachmentHarmonic>::computeForcesAux(double *coord, double *f);
template void AFMAttachment<AFMAttachmentHarmonic>::vectorize(const FFCoordinateStartingIndex&);
template void AFMAttachment<AFMAttachmentHarmonic>::deallocate();


