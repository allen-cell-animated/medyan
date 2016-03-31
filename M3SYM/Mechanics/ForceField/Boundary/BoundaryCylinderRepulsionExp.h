
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

#ifndef M3SYM_BoundaryCylinderRepulsionExp_h
#define M3SYM_BoundaryCylinderRepulsionExp_h

#include <vector>
#include <cmath>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A exponential repulsive potential used by the BoundaryCylinderRepulsion template.
class BoundaryCylinderRepulsionExp {
    
public:
    double energy(Bead*, double, double, double);
    void forces(Bead*, double, vector<double>& norm, double, double);
    void forcesAux(Bead*, double, vector<double>& norm, double, double);
};

#endif
