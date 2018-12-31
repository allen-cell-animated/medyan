
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

#ifndef MEDYAN_FilamentBendingCosine_h
#define MEDYAN_FilamentBendingCosine_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the [FilamentBending](@ref FilamentBending) template.
class FilamentBendingCosine {
    
public:
    double energy(Bead*, Bead*, Bead*, double, double, bool stretched);
    [[deprecated]] double energy(Bead*, Bead*, Bead*, double, double, double);
    
    void forces(Bead*, Bead*, Bead*, double, double);
    void forcesAux(Bead*, Bead*, Bead*, double, double);
};

#endif
