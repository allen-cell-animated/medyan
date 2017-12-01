#ifndef MEDYAN_MembraneStretchingVoronoiHarmonic_h
#define MEDYAN_MembraneStretchingVoronoiHarmonic_h

#include <array>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the MembraneStretching
class MembraneStretchingVoronoiHarmonic {
    
public:
    double energy(const std::array<Bead*, 3>&, double, double);
    double energy(const std::array<Bead*, 3>&, double, double, double);
    
    void forces(const std::array<Bead*, 3>&, double, double);
    void forcesAux(const std::array<Bead*, 3>&, double, double);
};

#endif
