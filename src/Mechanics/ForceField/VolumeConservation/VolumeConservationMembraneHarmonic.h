#ifndef MEDYAN_VolumeConservationMembraneHarmonic_h
#define MEDYAN_VolumeConservationMembraneHarmonic_h

#include "MathFunctions.h"

// Forward declaration
class Vertex;

class VolumeConservationMembraneHarmonic {

public:
    double energy(double, double, double);
    
    void forces(floatingpoint* force, double, const mathfunc::Vec3&, double, double);
};


#endif
