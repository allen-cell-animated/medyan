
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

#include "MTOCAttachment.h"

#include "MTOCAttachmentHarmonic.h"

#include "MTOC.h"
#include "Bubble.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

template <class MTOCInteractionType>
void MTOCAttachment<MTOCInteractionType>::vectorize() {
    
    for(auto mtoc : MTOC::getMTOCs()) {
        beadSet = new int[n * mtoc->getFilaments().size() + 1];
        kstr = new floatingpoint[n * Cylinder::getCylinders().size() + 1];

        beadSet[0] = mtoc->getBubble()->getBead()->getIndex();
        kstr[0] = 0;

        int i = 1;
        
        for (int fIndex = 0; fIndex < mtoc->getFilaments().size(); fIndex++) {
            Filament *f = mtoc->getFilaments()[fIndex];
            
            beadSet[n * i] = f->getMinusEndCylinder()->getFirstBead()->getIndex();

            kstr[n * i] = f->getMinusEndCylinder()->getMCylinder()->getStretchingConst();
            
            i++;
        }
    }
}

template <class MTOCInteractionType>
void MTOCAttachment<MTOCInteractionType>::deallocate() {

    delete [] beadSet;
    delete [] kstr;
}


template <class MTOCInteractionType>
floatingpoint MTOCAttachment<MTOCInteractionType>::computeEnergy(floatingpoint* coord, bool stretched) {
    
    floatingpoint U = 0.0;
    floatingpoint U_i=0.0;
    
    //TO DO, for loop may be removed

    for(auto mtoc : MTOC::getMTOCs()) {

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
        floatingpoint radius = mtoc->getBubble()->getRadius();
        
        U_i = _FFType.energy(coord, beadSet, kstr, radius);
    }

    return U_i;

    //            if(fabs(U_i) == numeric_limits<double>::infinity()
    //               || U_i != U_i || U_i < -1.0) {
    //
    //                //set culprits and return
    //                _otherCulprit = f;
    //                _bubbleCulprit = mtoc->getBubble();
    //
    //                return -1;
    //            }
    //            else
    //                U += U_i;
    //        }
    //    }
    //    return U;

}

template <class MTOCInteractionType>
void MTOCAttachment<MTOCInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
    
    for(auto mtoc : MTOC::getMTOCs()) {
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
	    floatingpoint radius = mtoc->getBubble()->getRadius();
        _FFType.forces(coord, f, beadSet, kstr, radius);
        //        }
    }

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
template floatingpoint MTOCAttachment<MTOCAttachmentHarmonic>::computeEnergy(floatingpoint *coord, bool stretched);
template void MTOCAttachment<MTOCAttachmentHarmonic>::computeForces(floatingpoint *coord, floatingpoint *f);
//template void MTOCAttachment<MTOCAttachmentHarmonic>::computeForcesAux(double *coord, double *f);
template void MTOCAttachment<MTOCAttachmentHarmonic>::vectorize();
template void MTOCAttachment<MTOCAttachmentHarmonic>::deallocate();

