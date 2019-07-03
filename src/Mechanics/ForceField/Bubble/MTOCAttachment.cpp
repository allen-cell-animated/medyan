
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
void MTOCAttachment<MTOCInteractionType>::vectorize() {

	//get total number of interactions
	uint nint = 0;
	for(auto mtoc : MTOC::getMTOCs()) {
		nint += mtoc->getFilaments().size();
	}
	//create vectors
	beadSet = new int[n * nint];
	kstr = new floatingpoint[nint];
	radiusvec = new floatingpoint[nint];
	//Get the interactions
	uint interaction_counter = 0;
    for(auto mtoc : MTOC::getMTOCs()) {

        for (int fIndex = 0; fIndex < mtoc->getFilaments().size(); fIndex++) {
            Filament *f = mtoc->getFilaments()[fIndex];
			//get mtoc bead
	        beadSet[n*interaction_counter] = mtoc->getBubble()->getBead()->_dbIndex;
            //get filament bead
            beadSet[n*interaction_counter + 1] = f->getMinusEndCylinder()->getFirstBead()
            		->_dbIndex;
			//The MTOC attachment constant is the same as stretching constant
            kstr[interaction_counter] = f->getMinusEndCylinder()->getMCylinder()->getStretchingConst();

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
floatingpoint MTOCAttachment<MTOCInteractionType>::computeEnergy(floatingpoint* coord,
		floatingpoint *f, floatingpoint d) {

	floatingpoint U = 0.0;
	floatingpoint U_i=0.0;
	if (d == 0.0)
		U_i = _FFType.energy(coord, f, beadSet, kstr, radiusvec);
	else
		U_i = _FFType.energy(coord, f, beadSet, kstr, radiusvec, d);
	return U_i;
}

template <class MTOCInteractionType>
void MTOCAttachment<MTOCInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
        _FFType.forces(coord, f, beadSet, kstr, radiusvec);

}


//template <class MTOCInteractionType>
//void MTOCAttachment<MTOCInteractionType>::computeForcesAux(double *coord, double *f) {
//    cout << "MTOCAttachment<MTOCInteractionType>::computeForcesAux should not be called in vectorized version." << endl;
//
//    for(auto mtoc : MTOC::getMTOCs()) {
//
//        Bead* b1 = mtoc->getBubble()->getBead();
//
//        for(int fIndex = 0; fIndex < mtoc->getFilaments().size(); fIndex++) {
//
//            Filament *f = mtoc->getFilaments()[fIndex];
//
//            Cylinder* c = f->getMinusEndCylinder();
//
//            Bead* b2 = c->getFirstBead();
//            double kStretch = c->getMCylinder()->getStretchingConst();
//            double radius = mtoc->getBubble()->getRadius();
//
//            _FFType.forcesAux(coord, f, beadSet, kstr);
//        }
//    }
//}

///Template specializations
template <class MTOCInteractionType>
int MTOCAttachment<MTOCInteractionType>::numInteractions;
template floatingpoint MTOCAttachment<MTOCAttachmentHarmonic>::computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d);
template void MTOCAttachment<MTOCAttachmentHarmonic>::computeForces(floatingpoint *coord, floatingpoint *f);
//template void MTOCAttachment<MTOCAttachmentHarmonic>::computeForcesAux(double *coord, double *f);
template void MTOCAttachment<MTOCAttachmentHarmonic>::vectorize();
template void MTOCAttachment<MTOCAttachmentHarmonic>::deallocate();

