
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
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
double MTOCAttachment<MTOCInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for(auto mtoc : MTOC::getMTOCs()) {
        
        Bead* b1 = mtoc->getBubble()->getBead();
        
        for(int fIndex = 0; fIndex < mtoc->getFilaments().size(); fIndex++) {
            
            Filament *f = mtoc->getFilaments()[fIndex];
            
            Cylinder* c = f->getMinusEndCylinder();
            
            Bead* b2 = c->getFirstBead();
            double kStretch = c->getMCylinder()->getStretchingConst();
            double radius = mtoc->getBubble()->getRadius();
            
            if (d == 0.0)
                U_i = _FFType.energy(b1, b2, kStretch, radius);
            else
                U_i = _FFType.energy(b1, b2, kStretch, radius, d);
            
            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {
                
                //set culprits and return
                _otherCulprit = f;
                _bubbleCulprit = mtoc->getBubble();
                
                return -1;
            }
            else
                U += U_i;
        }
    }
    return U;
}

template <class MTOCInteractionType>
void MTOCAttachment<MTOCInteractionType>::computeForces() {
    
    for(auto mtoc : MTOC::getMTOCs()) {
        
        Bead* b1 = mtoc->getBubble()->getBead();
        
        for(int fIndex = 0; fIndex < mtoc->getFilaments().size(); fIndex++) {
            
            Filament *f = mtoc->getFilaments()[fIndex];
            
            Cylinder* c = f->getMinusEndCylinder();
            
            Bead* b2 = c->getFirstBead();
            double kStretch = c->getMCylinder()->getStretchingConst();
            double radius = mtoc->getBubble()->getRadius();

            _FFType.forces(b1, b2, kStretch, radius);
        }
    }
}


template <class MTOCInteractionType>
void MTOCAttachment<MTOCInteractionType>::computeForcesAux() {
    
    for(auto mtoc : MTOC::getMTOCs()) {
    
        Bead* b1 = mtoc->getBubble()->getBead();
        
        for(int fIndex = 0; fIndex < mtoc->getFilaments().size(); fIndex++) {
            
            Filament *f = mtoc->getFilaments()[fIndex];
            
            Cylinder* c = f->getMinusEndCylinder();
            
            Bead* b2 = c->getFirstBead();
            double kStretch = c->getMCylinder()->getStretchingConst();
            double radius = mtoc->getBubble()->getRadius();
            
            _FFType.forcesAux(b1, b2, kStretch, radius);
        }
    }
}

///Template specializations
template double MTOCAttachment<MTOCAttachmentHarmonic>::computeEnergy(double d);
template void MTOCAttachment<MTOCAttachmentHarmonic>::computeForces();
template void MTOCAttachment<MTOCAttachmentHarmonic>::computeForcesAux();
