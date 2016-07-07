
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
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
    double energy(Bead*, Bead*, double, double, double);
    double energy(Bead*, Bead*, double, double, double, double);
    
    void forces(Bead*, Bead*, double, double, double);
    void forcesAux(Bead*, Bead*, double, double, double);
    
    double loadForces(Bead*, Bead*, double, double, double);
};

#endif
