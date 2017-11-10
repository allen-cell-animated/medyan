#ifndef MEDYAN_MembraneStretchingHarmonic_h
#define MEDYAN_MembraneStretchingHarmonic_h

#include <array>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the MembraneStretching
class MembraneStretchingHarmonic {
    
public:
    double energy(const std::array<Bead*, 3>, double, double);
    double energy(Bead*, Bead*, double, double, double);
    
    void forces(Bead*, Bead*, double, double);
    void forcesAux(Bead*, Bead*, double, double);
};

#endif
