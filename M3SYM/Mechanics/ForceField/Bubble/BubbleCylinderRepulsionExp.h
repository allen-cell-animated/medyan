
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

#ifndef M3SYM_BubbleCylinderRepulsionExp_h
#define M3SYM_BubbleCylinderRepulsionExp_h

#include <vector>
#include <cmath>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A exponential repulsive potential used by the BubbleCylinderRepulsion template.
class BubbleCylinderRepulsionExp {
    
public:
    double energy(Bead*, Bead*, double, double, double);
    double energy(Bead*, Bead*, double, double, double, double);
    
    void forces(Bead*, Bead*, double, double, double);
    void forcesAux(Bead*, Bead*, double, double, double);
};

#endif
