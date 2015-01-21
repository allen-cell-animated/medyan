
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

#ifndef M3SYM_BoundaryRepulsionLJ_h
#define M3SYM_BoundaryRepulsionLJ_h

#include <vector>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A Lennard-Jones repulsive potential used by the BoundaryRepulsion template.
class BoundaryRepulsionLJ {

public:
    double computeEnergy(Bead*, double, double, double);
    void computeForces(Bead*, double, vector<double>& norm, double, double);
    void computeForcesAux(Bead*, double, vector<double>& norm, double, double);
};
#endif
