
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

#ifndef M3SYM_FilamentStretchingHarmonic_h
#define M3SYM_FilamentStretchingHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the FilamentStretching template.
class FilamentStretchingHarmonic {
    
public:
    double energy(Bead*, Bead*, double, double);
    double energy(Bead*, Bead*, double, double, double);
    
    void forces(Bead*, Bead*, double, double);
    void forcesAux(Bead*, Bead*, double, double);
};

#endif
