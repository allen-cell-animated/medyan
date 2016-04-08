
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

#ifndef MEDYAN_BranchingBendingCosine_h
#define MEDYAN_BranchingBendingCosine_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the BranchingBending template.
class BranchingBendingCosine {
    
public:
    double energy(Bead*, Bead*, Bead*, Bead*, double, double);
    double energy(Bead*, Bead*, Bead*, Bead*, double, double, double);
    
    void forces(Bead*, Bead*, Bead*, Bead*, double, double);
    void forcesAux(Bead*, Bead*, Bead*, Bead*, double, double);
};

#endif
