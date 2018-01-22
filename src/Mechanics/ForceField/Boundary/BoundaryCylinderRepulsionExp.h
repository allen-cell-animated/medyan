
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
    double energy(double *coord, double *f, int *beadSet,
                  double *krep, double *slen, int *nneighbors);
    
    double energy(double *coord, double *f, int *beadSet,
                  double *krep, double *slen, int *nnneighbors, double d);
    
    void forces(double *coord, double *f, int *beadSet,
                double *krep, double *slen, int *nneighbors);
    
    double loadForces(double r, double krep , double slen);
#ifdef CROSSCHECK
    double energy(Bead*, double, double, double);
    void forces(Bead*, double, vector<double>& norm, double, double);
    void forcesAux(Bead*, double, vector<double>& norm, double, double);
#endif
private:
#ifdef CUDAACCL
//    double *F_i;
//    double *forcecopy;
#endif
};

#endif
