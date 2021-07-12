
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
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
    floatingpoint energy(
        const floatingpoint* coord,
        std::size_t i1, std::size_t i2, floatingpoint, floatingpoint, floatingpoint, floatingpoint);
    
    void forces(
        const floatingpoint* coord, floatingpoint* force,
        std::size_t i1, std::size_t i2, floatingpoint, floatingpoint, floatingpoint, floatingpoint);
};

#endif
