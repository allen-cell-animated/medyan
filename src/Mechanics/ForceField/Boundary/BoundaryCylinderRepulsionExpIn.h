
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

#ifndef MEDYAN_BoundaryCylinderRepulsionExpIn_h
#define MEDYAN_BoundaryCylinderRepulsionExpIn_h

#include <vector>
#include <cmath>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A exponential repulsive potential used by the BoundaryCylinderRepulsion template.
class BoundaryCylinderRepulsionExpIn {
    
public:
    //TODO needs implementation @{
    double energy(double *coord, double *f, int *beadSet,
                  double *krep, double *slen, int *nneighbors){return 0.0;};

    double energy(double *coord, double *f, int *beadSet,
                  double *krep, double *slen, int *nnneighbors, double d){return 0.0;};

    void forces(double *coord, double *f, int *beadSet,
                double *krep, double *slen, int *nneighbors){};
    //@}
    double energy(Bead*, double, double, double);
    void forces(Bead*, double, vector<double>& norm, double, double);
    void forcesAux(Bead*, double, vector<double>& norm, double, double);
    double loadForces(double, double, double);
};

#endif