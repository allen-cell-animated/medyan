#include <array>
#include <vector>

#include "MembraneStretchingVoronoiHarmonic.h"

#include "Vertex.h"

#include "MathFunctions.h"

using namespace mathfunc;

double MembraneStretchingVoronoiHarmonic::energy(double area,
                                                 double kElastic, double eqArea){
    // kElastic is the elastic modulus, which is independent of the actual eqArea

    double dist = area - eqArea;
    
    return 0.5 * kElastic * dist * dist / eqArea;
    
}

double MembraneStretchingVoronoiHarmonic::energy(double areaStretched,
                                                 double kElastic, double eqArea, double d){
    // In fact, d is a dummy variable here, as areaStretched is already dependent on d.

    double distStretched = areaStretched - eqArea;
    return 0.5 * kElastic * distStretched * distStretched / eqArea;
}

void MembraneStretchingVoronoiHarmonic::forces(
    Vertex* v, double area, const mathfunc::Vec3& dArea, double kElastic, double eqArea
) {
    // F_i = -grad_i U = -k / A_0 * (A - A_0) * grad_i A
    // A(rea) and grad_i A(rea) are obtained as function parameters

    const auto deltaF = (- kElastic * (area - eqArea) / eqArea) * dArea;

    for(size_t i = 0; i < 3; ++i) v->force[i] += deltaF[i]; // v->force += deltaF;
}

void MembraneStretchingVoronoiHarmonic::forcesAux(
    Vertex* v, double area, const mathfunc::Vec3& dArea, double kElastic, double eqArea
) {
    // Same as force calculation

    const auto deltaF = (- kElastic * (area - eqArea) / eqArea) * dArea;

    for(size_t i = 0; i < 3; ++i) v->forceAux[i] += deltaF[i]; // v->force += deltaF;
}
