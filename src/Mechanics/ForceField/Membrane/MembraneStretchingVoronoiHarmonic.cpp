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

void MembraneStretchingVoronoiHarmonic::forces(Vertex* vCenter, const std::vector<Vertex*>& v,
                                               double area, const std::array<double, 3>& dArea, const std::vector<std::array<double, 3>>& dNeighborArea,
                                               double kElastic, double eqArea) {
    // F_i = -grad_i U = -k / A_0 * (A - A_0) * grad_i A
    // A(rea) and grad_i A(rea) are obtained as function parameters

    size_t n = v.size();

    double coeff = -kElastic / eqArea * (area - eqArea);

    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        vCenter->force[coordIdx] += coeff * dArea[coordIdx];
        for(size_t nIdx = 0; nIdx < n; ++nIdx) {
            v[nIdx]->force[coordIdx] += coeff * dNeighborArea[nIdx][coordIdx];
        }
    }
}

void MembraneStretchingVoronoiHarmonic::forcesAux(Vertex* vCenter, const std::vector<Vertex*>& v,
                                                  double area, const std::array<double, 3>& dArea, const std::vector<std::array<double, 3>>& dNeighborArea,
                                                  double kElastic, double eqArea) {
    // Same as force calculation

    size_t n = v.size();

    double coeff = -kElastic / eqArea * (area - eqArea);

    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        vCenter->forceAux[coordIdx] += coeff * dArea[coordIdx];
        for(size_t nIdx = 0; nIdx < n; ++nIdx) {
            v[nIdx]->forceAux[coordIdx] += coeff * dNeighborArea[nIdx][coordIdx];
        }
    }
}
