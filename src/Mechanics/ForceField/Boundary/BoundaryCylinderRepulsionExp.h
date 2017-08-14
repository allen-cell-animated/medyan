
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

#ifndef MEDYAN_BoundaryCylinderRepulsionExp_h
#define MEDYAN_BoundaryCylinderRepulsionExp_h

#include <vector>
#include <cmath>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A exponential repulsive potential used by the BoundaryCylinderRepulsion template.
class BoundaryCylinderRepulsionExp {
    
public:
    inline double energy(double *coord, double *f, int *beadSet,
                         double *krep, double *slen, int *nneighbors);
    
    inline double energy(double *coord, double *f, int *beadSet,
                         double *krep, double *slen, int *nnneighbors, double d);
    
    inline void forces(double *coord, double *f, int *beadSet,
                       double *krep, double *slen, int *nneighbors);
    
    inline double loadForces(double r, double krep , double slen);
};

#endif
