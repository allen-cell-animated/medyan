#ifndef MEDYAN_VolumeConservationMembraneHarmonic_h
#define MEDYAN_VolumeConservationMembraneHarmonic_h

#include <array>
#include <vector>

// Forward declaration
class Vertex;

class VolumeConservationMembraneHarmonic {

public:
    double energy(double, double, double);
    double energy(double, double, double, double);
    
    void forces(const std::vector<Vertex*>&, double, const std::vector<std::array<double, 3>>&, double, double);
    void forcesAux(const std::vector<Vertex*>&, double, const std::vector<std::array<double, 3>>&, double, double);
};


#endif
