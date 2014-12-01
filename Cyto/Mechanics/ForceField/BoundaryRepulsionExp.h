
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_BoundaryRepulsionExp_h
#define M3SYM_BoundaryRepulsionExp_h

#include <vector>
#include <cmath>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// BoundaryRepulsionLJ class is a exponential repulsive potential.
class BoundaryRepulsionExp {
    
public:
    double computeEnergy(Bead*, double, double, double);
    void computeForces(Bead*, double, vector<double>& norm, double, double);
    void computeForcesAux(Bead*, double, vector<double>& norm, double, double);
};

#endif
