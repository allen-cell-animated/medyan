
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

#ifndef M3SYM_BoundaryBubbleRepulsionExp_h
#define M3SYM_BoundaryBubbleRepulsionExp_h

#include <vector>
#include <cmath>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A exponential repulsive potential used by the BoundaryBubbleRepulsion template.
class BoundaryBubbleRepulsionExp {
    
public:
    double energy(Bead*, double, double, double, double);
    void forces(Bead*, double, double, vector<double>& norm, double, double);
    void forcesAux(Bead*, double, double, vector<double>& norm, double, double);
};

#endif
