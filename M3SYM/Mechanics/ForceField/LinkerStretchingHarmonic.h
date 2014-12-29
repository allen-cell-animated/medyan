
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

#ifndef M3SYM_LinkerStretchingHarmonic_h
#define M3SYM_LinkerStretchingHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the LinkerStretching template.
class LinkerStretchingHarmonic {
    
public:
    double energy(Bead*, Bead*, Bead*, Bead*,
                  double position1, double position2,
                  double kStretch, double eqLength);
    double energy(Bead*, Bead*, Bead*, Bead*,
                  double position1, double position2,
                  double kStretch, double eqLength, double d);
    
    double forces(Bead*, Bead*, Bead*, Bead*,
                double position1, double position2,
                double kStretch, double eqLength);
    double forcesAux(Bead*, Bead*, Bead*, Bead*,
                double position1, double position2,
                double kStretch, double eqLength);
};

#endif 
