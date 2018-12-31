
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
double BoundaryCylinderAttachment<BAttachmentInteractionType>::computeEnergy(bool stretched) {
    
    double U = 0.0;
    double U_i=0.0;
    
    
    for(auto b : Bead::getPinnedBeads()) {
    
        double kAttr = SysParams::Mechanics().pinK;
            
        U_i =  _FFType.energy(b, kAttr, stretched);
            
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
                
            //set culprits and return
            _otherCulprit = b;
                
            return -1;
        }
        else
        U += U_i;
    }
    
    return U;
}

template <class BAttachmentInteractionType>
void BoundaryCylinderAttachment<BAttachmentInteractionType>::computeForces() {
    //Qin
    for(auto b : Bead::getPinnedBeads()) {
        double kAttr = SysParams::Mechanics().pinK;
        _FFType.forces(b, kAttr);
    }
}


template <class BAttachmentInteractionType>
void BoundaryCylinderAttachment<BAttachmentInteractionType>::computeForcesAux() {
    
    for(auto b : Bead::getPinnedBeads()) {
        
        double kAttr = SysParams::Mechanics().pinK;
        _FFType.forcesAux(b, kAttr);
    }
}

///Template specializations
template double BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::computeEnergy(bool stretched);
template void BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::computeForces();
template void BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::computeForcesAux();
