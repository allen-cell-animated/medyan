
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "MTOCAttachment.h"

#include "MTOCAttachmentHarmonic.h"

#include "MTOC.h"
#include "Bubble.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"



template <class MTOCInteractionType>
void MTOCAttachment<MTOCInteractionType>::vectorize(const FFCoordinateStartingIndex& si) {


	//get total number of interactions
	uint nint = 0;
	for(auto mtoc : MTOC::getMTOCs()) {
		nint += mtoc->getFilaments().size();
	}
	//create vectors
	beadSet = new int[n * nint];
    beadStartIndex_   = si.bead;
    bubbleStartIndex_ = si.bubble;
	kstr = new floatingpoint[nint];
	radiusvec = new floatingpoint[nint];
	//Get the interactions
	uint interaction_counter = 0;

    if(MTOC::getMTOCs().size() > 1) {
        cout << "Should not have more than 1 MTOC" << endl;
        exit(EXIT_FAILURE);
    }
    

    for(auto mtoc : MTOC::getMTOCs()) {

        for (int fIndex = 0; fIndex < mtoc->getFilaments().size(); fIndex++) {
            Filament *f = mtoc->getFilaments()[fIndex];
			//get mtoc bead
	        beadSet[n*interaction_counter] = mtoc->getBubble()->getIndex() * 3 + si.bubble;
            //get filament bead
            beadSet[n*interaction_counter + 1] = f->getMinusEndCylinder()->getFirstBead()
            		->getIndex() * 3 + si.bead;
			//The MTOC attachment constant is the same as stretching constant
            kstr[interaction_counter] = f->getMinusEndCylinder()->getMCylinder()->getStretchingConst();
            radiusvec[interaction_counter] = mtoc->getBubble()->getRadius();
            
	        interaction_counter++;
        }
    }
    numInteractions = interaction_counter;
}

template <class MTOCInteractionType>
void MTOCAttachment<MTOCInteractionType>::deallocate() {

    delete [] beadSet;
    delete [] kstr;
    delete [] radiusvec;
}


template <class MTOCInteractionType>
floatingpoint MTOCAttachment<MTOCInteractionType>::computeEnergy(floatingpoint* coord, bool stretched) {

    floatingpoint U = 0.0;
    floatingpoint U_i=0.0;

    U_i = _FFType.energy(coord, beadSet, beadStartIndex_, bubbleStartIndex_, kstr, radiusvec);

    return U_i;
}

template <class MTOCInteractionType>
void MTOCAttachment<MTOCInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
        _FFType.forces(coord, f, beadSet, kstr, radiusvec);

}



///Template specializations
template <class MTOCInteractionType>
int MTOCAttachment<MTOCInteractionType>::numInteractions;
template floatingpoint MTOCAttachment<MTOCAttachmentHarmonic>::computeEnergy(floatingpoint *coord, bool stretched);
template void MTOCAttachment<MTOCAttachmentHarmonic>::computeForces(floatingpoint *coord, floatingpoint *f);
//template void MTOCAttachment<MTOCAttachmentHarmonic>::computeForcesAux(double *coord, double *f);
template void MTOCAttachment<MTOCAttachmentHarmonic>::vectorize(const FFCoordinateStartingIndex&);
template void MTOCAttachment<MTOCAttachmentHarmonic>::deallocate();

