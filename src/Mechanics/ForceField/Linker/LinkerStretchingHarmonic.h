
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

#ifndef MEDYAN_LinkerStretchingHarmonic_h
#define MEDYAN_LinkerStretchingHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the LinkerStretching template.
class LinkerStretchingHarmonic {
    
public:
    double energy(Bead*, Bead*, Bead*, Bead*,
                  double position1, double position2,
                  double kStretch, double eqLength, bool stretched);
    [[deprecated]] double energy(Bead*, Bead*, Bead*, Bead*,
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
