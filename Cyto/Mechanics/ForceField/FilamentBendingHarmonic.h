
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_FilamentBendingHarmonic_h
#define M3SYM_FilamentBendingHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// FilamentBendingHarmonic class is a harmonic potential used by the [FilamentBending](@ref FilamentBending) template.
class FilamentBendingHarmonic {
    
public:
    double energy(Bead*, Bead*, Bead*, double);
    double energy(Bead*, Bead*, Bead*, double, double);
    void forces(Bead*, Bead*, Bead*, double);
    void forcesAux(Bead*, Bead*, Bead*, double);
};

#endif
