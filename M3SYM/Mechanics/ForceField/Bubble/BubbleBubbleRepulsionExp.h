
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
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
    double computeEnergy(Bead*, Bead*, double, double, double, double);
    double computeEnergy(Bead*, Bead*, double, double, double, double, double);
    
    void computeForces(Bead*, Bead*, double, double, double, double);
    void computeForcesAux(Bead*, Bead*, double, double, double, double);
};

#endif
