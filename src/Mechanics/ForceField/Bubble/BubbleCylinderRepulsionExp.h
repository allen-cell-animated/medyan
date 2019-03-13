
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

#ifndef MEDYAN_BubbleCylinderRepulsionExp_h
#define MEDYAN_BubbleCylinderRepulsionExp_h

#include <vector>
#include <cmath>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A exponential repulsive potential used by the BubbleCylinderRepulsion template.
class BubbleCylinderRepulsionExp {
    
public:
    floatingpoint energy(Bead*, Bead*, floatingpoint, floatingpoint, floatingpoint);
    floatingpoint energy(Bead*, Bead*, floatingpoint, floatingpoint, floatingpoint, floatingpoint);
    
    void forces(Bead*, Bead*, floatingpoint, floatingpoint, floatingpoint);
    void forcesAux(Bead*, Bead*, floatingpoint, floatingpoint, floatingpoint);
    
    floatingpoint loadForces(Bead*, Bead*, floatingpoint, floatingpoint, floatingpoint);
};

#endif
