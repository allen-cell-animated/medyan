
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
void AFMAttachment<AFMInteractionType>::vectorize() {
    
    if(AFM::getAFMs().size() > 1) {
        cout << "Should not have more than 1 AFM" << endl;
        exit(EXIT_FAILURE);
    }
    
    for(auto afm : AFM::getAFMs()) {
        beadSet = new int[n * afm->getFilaments().size() + 1];
        kstr = new double[n * Cylinder::getCylinders().size() + 1];
        
        beadSet[0] = afm->getBubble()->getBead()->_dbIndex;
        kstr[0] = 0;
        
        int i = 1;
        
        for (int fIndex = 0; fIndex < afm->getFilaments().size(); fIndex++) {
            Filament *f = afm->getFilaments()[fIndex];
            
            beadSet[n * i] = f->getMinusEndCylinder()->getFirstBead()->_dbIndex;
            
            kstr[n * i] = f->getMinusEndCylinder()->getMCylinder()->getStretchingConst();
            
            i++;
        }
    }
}

template <class AFMInteractionType>
void AFMAttachment<AFMInteractionType>::deallocate() {
    
    delete [] beadSet;
    delete [] kstr;
}


template <class AFMInteractionType>
double AFMAttachment<AFMInteractionType>::computeEnergy(double* coord, double *f, double d) {
    
    double U = 0.0;
    double U_i=0.0;
    
    //TO DO, for loop may be removed
    
    for(auto afm : AFM::getAFMs()) {
        
        //        Bead* b1 = afm->getBubble()->getBead();
        //
        //        for(int fIndex = 0; fIndex < afm->getFilaments().size(); fIndex++) {
        //
        //            Filament *f = afm->getFilaments()[fIndex];
        //
        //            Cylinder* c = f->getMinusEndCylinder();
        //
        //            Bead* b2 = c->getFirstBead();
        //            double kStretch = c->getMCylinder()->getStretchingConst();
        double radius = afm->getBubble()->getRadius();
        
        if (d == 0.0)
            U_i = _FFType.energy(coord, f, beadSet, kstr, radius);
        else
            U_i = _FFType.energy(coord, f, beadSet, kstr, radius, d);
    }
    
    return U_i;
    
    //            if(fabs(U_i) == numeric_limits<double>::infinity()
    //               || U_i != U_i || U_i < -1.0) {
    //
    //                //set culprits and return
    //                _otherCulprit = f;
    //                _bubbleCulprit = afm->getBubble();
    //
    //                return -1;
    //            }
    //            else
    //                U += U_i;
    //        }
    //    }
    //    return U;
    
}

template <class AFMInteractionType>
void AFMAttachment<AFMInteractionType>::computeForces(double *coord, double *f) {
    
    for(auto afm : AFM::getAFMs()) {
        //
        //        Bead* b1 = afm->getBubble()->getBead();
        //
        //        for(int fIndex = 0; fIndex < afm->getFilaments().size(); fIndex++) {
        //
        //            Filament *f = afm->getFilaments()[fIndex];
        //
        //            Cylinder* c = f->getMinusEndCylinder();
        //
        //            Bead* b2 = c->getFirstBead();
        //            double kStretch = c->getMCylinder()->getStretchingConst();
        double radius = afm->getBubble()->getRadius();
        _FFType.forces(coord, f, beadSet, kstr, radius);
        //        }
    }
    
}


//template <class AFMInteractionType>
//void AFMAttachment<AFMInteractionType>::computeForcesAux(double *coord, double *f) {
//    cout << "AFMAttachment<AFMInteractionType>::computeForcesAux should not be called in vectorized version." << endl;
//
//    for(auto afm : AFM::getAFMs()) {
//
//        Bead* b1 = afm->getBubble()->getBead();
//
//        for(int fIndex = 0; fIndex < afm->getFilaments().size(); fIndex++) {
//
//            Filament *f = afm->getFilaments()[fIndex];
//
//            Cylinder* c = f->getMinusEndCylinder();
//
//            Bead* b2 = c->getFirstBead();
//            double kStretch = c->getMCylinder()->getStretchingConst();
//            double radius = afm->getBubble()->getRadius();
//
//            _FFType.forcesAux(coord, f, beadSet, kstr);
//        }
//    }
//}

///Template specializations
template double AFMAttachment<AFMAttachmentHarmonic>::computeEnergy(double *coord, double *f, double d);
template void AFMAttachment<AFMAttachmentHarmonic>::computeForces(double *coord, double *f);
//template void AFMAttachment<AFMAttachmentHarmonic>::computeForcesAux(double *coord, double *f);
template void AFMAttachment<AFMAttachmentHarmonic>::vectorize();
template void AFMAttachment<AFMAttachmentHarmonic>::deallocate();


