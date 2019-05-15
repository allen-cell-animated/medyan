#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneStretchingHarmonic_Hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneStretchingHarmonic_Hpp

#include "MathFunctions.h"

/// A harmonic potential used by the MembraneStretching
class MembraneStretchingHarmonic {
    
public:
    double energy(double, double, double);
    
    void forces(double* force, double area, const mathfunc::Vec3& dArea, double kElastic, double eqArea);
};

#endif
