
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

#ifndef MEDYAN_BranchingPositionCosine_h
#define MEDYAN_BranchingPositionCosine_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the BranchingPosition template.
class BranchingPositionCosine {
    
public:
    double energy(double *coord, double *f, int *beadSet,
                  double *kpos, double *pos);
    
    double energy(double *coord, double *f, int *beadSet,
                  double *kpos, double *pos, double d);
    
    void forces(double *coord, double *f, int *beadSet,
                double *kpos, double *pos);
};

#endif
