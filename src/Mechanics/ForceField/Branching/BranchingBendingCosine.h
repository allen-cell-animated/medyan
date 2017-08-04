
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

#ifndef MEDYAN_BranchingBendingCosine_h
#define MEDYAN_BranchingBendingCosine_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the BranchingBending template.
class BranchingBendingCosine {
    
public:
    inline double energy(double *coord, double *f, int *beadSet,
                         double *kbend, double *eqt);
    
    inline double energy(double *coord, double *f, int *beadSet,
                         double *kbend, double *eqt, double d);
    
    inline void forces(double *coord, double *f, int *beadSet,
                       double *kbend, double *eqt);
};

#endif
