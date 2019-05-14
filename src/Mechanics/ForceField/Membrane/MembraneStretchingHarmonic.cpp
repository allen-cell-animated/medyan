#include "Mechanics/ForceField/Membrane/MembraneStretchingHarmonic.hpp"

using namespace mathfunc;

double MembraneStretchingHarmonic::energy(double area,
                                          double kElastic, double eqArea){
    // kElastic is the elastic modulus, which is independent of the actual eqArea

    double dist = area - eqArea;
    
    return 0.5 * kElastic * dist * dist / eqArea;
    
}

void MembraneStretchingHarmonic::forces(
    double* force, double area, const mathfunc::Vec3& dArea, double kElastic, double eqArea
) {
    // F_i = -grad_i U = -k / A_0 * (A - A_0) * grad_i A
    // A(rea) and grad_i A(rea) are obtained as function parameters

    const auto deltaF = (- kElastic * (area - eqArea) / eqArea) * dArea;

    for(size_t i = 0; i < 3; ++i) force[i] += deltaF[i]; // v->force += deltaF;
}
