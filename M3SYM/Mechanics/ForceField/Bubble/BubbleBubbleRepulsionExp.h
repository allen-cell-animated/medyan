
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

#ifndef M3SYM_BubbleBubbleRepulsionExp_h
#define M3SYM_BubbleBubbleRepulsionExp_h

#include <vector>
#include <cmath>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A exponential repulsive potential used by the BubbleBubbleRepulsion template.
class BubbleBubbleRepulsionExp {
    
public:
    double energy(Bead*, Bead*, double, double, double, double);
    double energy(Bead*, Bead*, double, double, double, double, double);
    
    void forces(Bead*, Bead*, double, double, double, double);
    void forcesAux(Bead*, Bead*, double, double, double, double);
};

#endif
