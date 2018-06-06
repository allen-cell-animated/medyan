
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

#include "BoundaryCylinderAttachment.h"

#include "BoundaryCylinderAttachmentHarmonic.h"
#include "BoundaryElement.h"

#include "Bubble.h"
#include "Bead.h"

template <class BAttachmentInteractionType>
double BoundaryCylinderAttachment<BAttachmentInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    
    for(auto b : Bead::getPinnedBeads()) {
    
        double kAttr = SysParams::Mechanics().pinK;
            
        if (d == 0.0)
            U_i =  _FFType.energy(b, kAttr);
        else
            U_i =  _FFType.energy(b, kAttr, d);
            
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
template double BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::computeEnergy(double d);
template void BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::computeForces();
template void BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::computeForcesAux();
