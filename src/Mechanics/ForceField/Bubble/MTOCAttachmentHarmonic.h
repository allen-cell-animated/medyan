
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

#ifndef MEDYAN_MTOCAttachmentHarmonic_h
#define MEDYAN_MTOCAttachmentHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the MTOCAttachment template.
class MTOCAttachmentHarmonic {
    
public:
    floatingpoint energy(Bead*, Bead*, floatingpoint, floatingpoint);
    floatingpoint energy(Bead*, Bead*, floatingpoint, floatingpoint, floatingpoint);
    
    void forces(Bead*, Bead*, floatingpoint, floatingpoint);
    void forcesAux(Bead*, Bead*, floatingpoint, floatingpoint);
};

#endif
