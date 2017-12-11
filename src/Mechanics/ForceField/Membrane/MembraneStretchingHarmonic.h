#ifndef MEDYAN_MembraneStretchingHarmonic_h
#define MEDYAN_MembraneStretchingHarmonic_h

#include <array>

#include "common.h"

//FORWARD DECLARATIONS
class Vertex;
class Triangle;

/// A harmonic potential used by the MembraneStretching
class MembraneStretchingHarmonic {
    
public:
    double energy(double, double, double);
    double energy(double, double, double, double);
    
    void forces(const std::array<Vertex*, 3>&, double, const std::array<std::array<double, 3>, 3>&, double, double);
    void forcesAux(const std::array<Vertex*, 3>&, double, const std::array<std::array<double, 3>, 3>&, double, double);
};

#endif
