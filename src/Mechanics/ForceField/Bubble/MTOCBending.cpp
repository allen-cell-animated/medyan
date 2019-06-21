
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

#include "MTOCBending.h"

#include "MTOCBendingCosine.h"

#include "MTOC.h"
#include "Bubble.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

template <class MTOCInteractionType>
void MTOCBending<MTOCInteractionType>::vectorize() {
    //Watch out! Only one MTOC is allowed
    for(auto mtoc : MTOC::getMTOCs()) {
        beadSet = new int[n *  mtoc->getFilaments().size() + 1];
        kbend = new floatingpoint[n *  mtoc->getFilaments().size() + 1];
        
        beadSet[0] = mtoc->getBubble()->getBead()->getStableIndex();
        kbend[0] = 0.0;
        
        int i = 1;
        
        for (int fIndex = 0; fIndex < mtoc->getFilaments().size(); fIndex++) {
            Filament *f = mtoc->getFilaments()[fIndex];
            
            beadSet[n * i - 1] = f->getMinusEndCylinder()->getFirstBead()->getStableIndex();
            beadSet[n * i] = f->getMinusEndCylinder()->getSecondBead()->getStableIndex();
            
            floatingpoint kk = mtoc->getBubble()->getMTOCBendingK();
            kbend[n * i - 1] = kk;
            kbend[n * i] = kk;
            
//            kbend[2 * n * i - 1] = f->getMinusEndCylinder()->getMCylinder()->getBendingConst();
//            kbend[2 * n * i] = f->getMinusEndCylinder()->getMCylinder()->getBendingConst();
            i++;
        }
    }
}

template <class MTOCInteractionType>
void MTOCBending<MTOCInteractionType>::deallocate() {
    
    delete [] beadSet;
    delete [] kbend;
}


template <class MTOCInteractionType>
floatingpoint MTOCBending<MTOCInteractionType>::computeEnergy(floatingpoint* coord, floatingpoint *f, floatingpoint d) {
    
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
        //            floatingpoint kStretch = c->getMCylinder()->getStretchingConst();
        floatingpoint radius = mtoc->getBubble()->getRadius();
        
        if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, kbend, radius);
        else
        U_i = _FFType.energy(coord, f, beadSet, kbend, radius, d);
    }
    
    return U_i;
    
    //            if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
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
void MTOCBending<MTOCInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
    
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
        //            floatingpoint kStretch = c->getMCylinder()->getStretchingConst();
        floatingpoint radius = mtoc->getBubble()->getRadius();
        _FFType.forces(coord, f, beadSet, kbend, radius);
        //        }
    }
    
}


//template <class MTOCInteractionType>
//void MTOCAttachment<MTOCInteractionType>::computeForcesAux(floatingpoint *coord, floatingpoint *f) {
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
//            floatingpoint kStretch = c->getMCylinder()->getStretchingConst();
//            floatingpoint radius = mtoc->getBubble()->getRadius();
//
//            _FFType.forcesAux(coord, f, beadSet, kstr);
//        }
//    }
//}

///Template specializations
template floatingpoint MTOCBending<MTOCBendingCosine>::computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d);
template void MTOCBending<MTOCBendingCosine>::computeForces(floatingpoint *coord, floatingpoint *f);
//template void MTOCAttachment<MTOCAttachmentHarmonic>::computeForcesAux(floatingpoint *coord, floatingpoint *f);
template void MTOCBending<MTOCBendingCosine>::vectorize();
template void MTOCBending<MTOCBendingCosine>::deallocate();


