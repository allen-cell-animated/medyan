
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

#ifndef MEDYAN_AFMAttachmentHarmonic_h
#define MEDYAN_AFMAttachmentHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the AFMAttachment template.
class AFMAttachmentHarmonic {
    
public:
    floatingpoint energy(
        floatingpoint *coord,
        int numInteractions, int *beadSet, floatingpoint *kstr, const floatingpoint* radii
    ) const;

    void forces(
        floatingpoint *coord, floatingpoint* f,
        int numInteractions, int *beadSet, floatingpoint *kstr, const floatingpoint* radii
    ) const;

    //void forcesAux(Bead*, Bead*, floatingpoint, floatingpoint);
};

#endif

