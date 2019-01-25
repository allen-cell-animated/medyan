
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
    
    for(auto mtoc : MTOC::getMTOCs()) {
        beadSet = new int[n *  mtoc->getFilaments().size() + 1];
        kbend = new double[Cylinder::getCylinders().size() + 1];
        
        beadSet[0] = mtoc->getBubble()->getBead()->_dbIndex;
        kbend[0] = 0;
        
        int i = 1;
        
        for (int fIndex = 0; fIndex < mtoc->getFilaments().size(); fIndex++) {
            Filament *f = mtoc->getFilaments()[fIndex];
            
            beadSet[n * i] = f->getMinusEndCylinder()->getFirstBead()->_dbIndex;
            
            kbend[i] = f->getMinusEndCylinder()->getMCylinder()->getStretchingConst();
            
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
double MTOCBending<MTOCInteractionType>::computeEnergy(double* coord, double *f, double d) {
    
    double U = 0.0;
    double U_i=0.0;
    
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
        double radius = mtoc->getBubble()->getRadius();
        
        if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, kbend, radius);
        else
        U_i = _FFType.energy(coord, f, beadSet, kbend, radius, d);
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
void MTOCBending<MTOCInteractionType>::computeForces(double *coord, double *f) {
    
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
        double radius = mtoc->getBubble()->getRadius();
        _FFType.forces(coord, f, beadSet, kbend, radius);
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
template double MTOCBending<MTOCBendingCosine>::computeEnergy(double *coord, double *f, double d);
template void MTOCBending<MTOCBendingCosine>::computeForces(double *coord, double *f);
//template void MTOCAttachment<MTOCAttachmentHarmonic>::computeForcesAux(double *coord, double *f);
template void MTOCBending<MTOCBendingCosine>::vectorize();
template void MTOCBending<MTOCBendingCosine>::deallocate();


