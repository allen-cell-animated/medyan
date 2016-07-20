
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

#ifndef MEDYAN_BoundaryCylinderRepulsionExp_h
#define MEDYAN_BoundaryCylinderRepulsionExp_h

#include <vector>
#include <cmath>

#include <BoundaryElement.h>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A exponential repulsive potential used by the BoundaryCylinderRepulsion template.
class BoundaryCylinderRepulsionExp {
    
public:
    double energy(Bead*, double, double, double);
    void forces(Bead*, double, vector<double>& norm, double, double, BoundaryElement*);
    void forcesAux(Bead*, double, vector<double>& norm, double, double, BoundaryElement*); //edited by jl135
    double loadForces(double, double, double);
};

#endif
