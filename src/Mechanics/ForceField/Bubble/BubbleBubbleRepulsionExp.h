
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

#ifndef MEDYAN_BubbleBubbleRepulsionExp_h
#define MEDYAN_BubbleBubbleRepulsionExp_h

#include <vector>
#include <cmath>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A exponential repulsive potential used by the BubbleBubbleRepulsion template.
class BubbleBubbleRepulsionExp {
    
public:
    floatingpoint energy(Bead*, Bead*, floatingpoint, floatingpoint, floatingpoint, floatingpoint);
    floatingpoint energy(Bead*, Bead*, floatingpoint, floatingpoint, floatingpoint, floatingpoint, floatingpoint);
    
    void forces(Bead*, Bead*, floatingpoint, floatingpoint, floatingpoint, floatingpoint);
    void forcesAux(Bead*, Bead*, floatingpoint, floatingpoint, floatingpoint, floatingpoint);
};

#endif
