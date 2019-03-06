
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
    double energy(double *coord, double *f, int *beadSet, int *bubbleSet,
                  double *krep, double *slen, double *radius, int *nneighbors);
    
    double energy(double *coord, double *f, int *beadSet, int *bubbleSet,
                  double *krep, double *slen, double *radius, int *nnneighbors, double d);
    
    void forces(double *coord, double *f, int *beadSet, int *bubbleSet,
                double *krep, double *slen, double *radius, int *nneighbors);
    void forcesAux(Bead*, Bead*, double, double, double);
    
    double loadForces(Bead*, Bead*, double, double, double);
};

#endif
