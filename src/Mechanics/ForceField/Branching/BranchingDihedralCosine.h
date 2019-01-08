
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BranchingDihedralCosine_h
#define MEDYAN_BranchingDihedralCosine_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the BranchingDihedraltemplate.
class BranchingDihedralCosine {
    
public:
    double energy(Bead*, Bead*, Bead*, Bead*, double, double, bool stretched);
    
    double forces(Bead*, Bead*, Bead*, Bead*, double, double);
    void forcesAux(Bead*, Bead*, Bead*, Bead*, double, double);
};

#endif
