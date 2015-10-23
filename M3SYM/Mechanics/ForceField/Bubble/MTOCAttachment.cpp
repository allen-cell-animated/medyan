
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "MTOCAttachment.h"

#include "FilamentStretchingHarmonic.h"

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
        
        for(auto &f : mtoc->getFilaments()) {
            
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
        
        for(auto &f : mtoc->getFilaments()) {
            
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
        
        for(auto &f : mtoc->getFilaments()) {
            
            Cylinder* c = f->getMinusEndCylinder();
            
            Bead* b2 = c->getFirstBead();
            double kStretch = c->getMCylinder()->getStretchingConst();
            
            double radius = mtoc->getBubble()->getRadius();
            
            _FFType.forcesAux(b1, b2, kStretch, radius);
        }
    }
}

///Template specializations
template double MTOCAttachment<FilamentStretchingHarmonic>::computeEnergy(double d);
template void MTOCAttachment<FilamentStretchingHarmonic>::computeForces();
template void MTOCAttachment<FilamentStretchingHarmonic>::computeForcesAux();
