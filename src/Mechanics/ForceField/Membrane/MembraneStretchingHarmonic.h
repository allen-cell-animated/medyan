#ifndef MEDYAN_MembraneStretchingHarmonic_h
#define MEDYAN_MembraneStretchingHarmonic_h

#include <array>

#include "common.h"
#include "MathFunctions.h"

//FORWARD DECLARATIONS
class Vertex;
class Triangle;

/// A harmonic potential used by the MembraneStretching
class MembraneStretchingHarmonic {
    
public:
    double energy(double, double, double);
    
    void forces(double* force, double area, const mathfunc::Vec3& dArea, double kElastic, double eqArea);
    void forcesAux(Vertex* v, double area, const mathfunc::Vec3& dArea, double kElastic, double eqArea);
};

#endif
