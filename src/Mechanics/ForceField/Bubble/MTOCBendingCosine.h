
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

#ifndef MEDYAN_MTOCBendingCosine_h
#define MEDYAN_MTOCBendingCosine_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the MTOCAttachment template.
class MTOCBendingCosine {
    
public:
    floatingpoint energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                  floatingpoint *kbend, floatingpoint);
    floatingpoint energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                  floatingpoint *kbend, floatingpoint, floatingpoint);
    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                floatingpoint *kbend, floatingpoint);
    
    //void forcesAux(Bead*, Bead*, floatingpoint, floatingpoint);
};

#endif


