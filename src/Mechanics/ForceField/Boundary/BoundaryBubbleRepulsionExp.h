
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BoundaryBubbleRepulsionExp_h
#define MEDYAN_BoundaryBubbleRepulsionExp_h

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
