
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef M3SYM_MTOCAttachmentHarmonic_h
#define M3SYM_MTOCAttachmentHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the MTOCAttachment template.
class MTOCAttachmentHarmonic {
    
public:
    double energy(Bead*, Bead*, double, double);
    double energy(Bead*, Bead*, double, double, double);
    
    void forces(Bead*, Bead*, double, double);
    void forcesAux(Bead*, Bead*, double, double);
};

#endif