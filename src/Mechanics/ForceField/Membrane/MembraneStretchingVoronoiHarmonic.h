#ifndef MEDYAN_MembraneStretchingVoronoiHarmonic_h
#define MEDYAN_MembraneStretchingVoronoiHarmonic_h

#include <array>
#include <vector>

#include "common.h"

//FORWARD DECLARATIONS
class Vertex;

/// A harmonic potential used by the MembraneStretching
class MembraneStretchingVoronoiHarmonic {
    
public:
    double energy(double, double, double);
    double energy(double, double, double, double);
    
    void forces(Vertex*, const std::vector<Vertex*>&, double, const std::array<double, 3>&, const std::vector<std::array<double, 3>>&, double, double);
    void forcesAux(Vertex*, const std::vector<Vertex*>&, double, const std::array<double, 3>&, const std::vector<std::array<double, 3>>&, double, double);
};

#endif
