#include "MembraneStretchingHarmonic.h"

#include "Vertex.h"
#include "Triangle.h"

#include "MathFunctions.h"

using namespace mathfunc;

double MembraneStretchingHarmonic::energy(double area,
                                          double kElastic, double eqArea){
    // kElastic is the elastic modulus, which is independent of the actual eqArea

    double dist = area - eqArea;
    
    return 0.5 * kElastic * dist * dist / eqArea;
    
}

double MembraneStretchingHarmonic::energy(double areaStretched,
                                          double kElastic, double eqArea, double d){
    // In fact, d is a dummy variable here, as areaStretched is already dependent on d.

    double distStretched = areaStretched - eqArea;
    return 0.5 * kElastic * distStretched * distStretched / eqArea;
}

void MembraneStretchingHarmonic::forces(const std::array<Vertex*, 3>& v,
                                        double area, const std::array<std::array<double, 3>, 3>& dArea,
                                        double kElastic, double eqArea) {
    // F_i = -grad_i U = -k / A_0 * (A - A_0) * grad_i A
    // A(rea) and grad_i A(rea) are obtained as function parameters

    double coeff = -kElastic / eqArea * (area - eqArea);
    for(size_t vtxIdx = 0; vtxIdx < 3; ++vtxIdx) {
        for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
            v[vtxIdx]->force[coordIdx] += coeff * dArea[vtxIdx][coordIdx];
        }
    }
}

void MembraneStretchingHarmonic::forcesAux(const std::array<Vertex*, 3>& v,
                                           double area, const std::array<std::array<double, 3>, 3>& dArea,
                                           double kElastic, double eqArea) {
    // Same as force calculation

    double coeff = -kElastic / eqArea * (area - eqArea);
    for(size_t vtxIdx = 0; vtxIdx < 3; ++vtxIdx) {
        for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
            v[vtxIdx]->forceAux[coordIdx] += coeff * dArea[vtxIdx][coordIdx];
        }
    }
}
