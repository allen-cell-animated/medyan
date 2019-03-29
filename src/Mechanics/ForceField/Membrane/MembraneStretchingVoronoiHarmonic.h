#ifndef MEDYAN_MembraneStretchingVoronoiHarmonic_h
#define MEDYAN_MembraneStretchingVoronoiHarmonic_h

#include <array>
#include <vector>

#include "common.h"
#include "MathFunctions.h"

//FORWARD DECLARATIONS
class Vertex;

/// A harmonic potential used by the MembraneStretching
class MembraneStretchingVoronoiHarmonic {
    
public:
    double energy(double, double, double);
    
    void forces(double* force, double area, const mathfunc::Vec3& dArea, double kElastic, double eqArea);
    void forcesAux(Vertex* v, double area, const mathfunc::Vec3& dArea, double kElastic, double eqArea);
};

#endif
