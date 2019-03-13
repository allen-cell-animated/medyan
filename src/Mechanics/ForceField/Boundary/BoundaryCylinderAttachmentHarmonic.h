
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

#ifndef MEDYAN_BoundaryCylinderAttachmentHarmonic_h
#define MEDYAN_BoundaryCylinderAttachmentHarmonic_h

#include <vector>
#include <cmath>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic attractive potential used by the BoundaryCylinderAttachment template.
class BoundaryCylinderAttachmentHarmonic {
    
public:
    floatingpoint energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                  floatingpoint *kattr, floatingpoint *pins);
    floatingpoint energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                  floatingpoint *kattr, floatingpoint *pins, floatingpoint d);
    
    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                floatingpoint *kattr, floatingpoint *pins);
};

#endif
