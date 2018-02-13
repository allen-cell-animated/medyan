#include "VolumeConservationMembraneHarmonic.h"

#include "Vertex.h"

double VolumeConservationMembraneHarmonic::energy(double volume, double kBulk, double eqVolume) {
    double dist = volume - eqVolume;
    return 0.5 * kBulk * dist * dist / eqVolume;
}

double VolumeConservationMembraneHarmonic::energy(double stretchedVolume, double kBulk, double eqVolume, double d) {
    double dist = stretchedVolume - eqVolume;
    return 0.5 * kBulk * dist * dist / eqVolume;
}

void VolumeConservationMembraneHarmonic::forces(
    const std::vector<Vertex*>& vertices,
    double volume,
    const std::vector<std::array<double, 3>>& dVolume,
    double kBulk, double eqVolume)
{
    // F_i = -grad_i U = -k / V_0 * (V - V_0) * grad_i V
    // V(olume) and grad_i V(olume) are obtained as function parameters

    double coeff = -kBulk / eqVolume * (volume - eqVolume);
    size_t n = vertices.size();
    for(size_t vtxIdx = 0; vtxIdx < n; ++vtxIdx) {
        for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
            vertices[vtxIdx]->force[coordIdx] += coeff * dVolume[vtxIdx][coordIdx];
        }
    }
}

void VolumeConservationMembraneHarmonic::forcesAux(
    const std::vector<Vertex*>& vertices,
    double volume,
    const std::vector<std::array<double, 3>>& dVolume,
    double kBulk, double eqVolume)
{
    // Same as force calculation

    double coeff = -kBulk / eqVolume * (volume - eqVolume);
    size_t n = vertices.size();
    for(size_t vtxIdx = 0; vtxIdx < n; ++vtxIdx) {
        for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
            vertices[vtxIdx]->forceAux[coordIdx] += coeff * dVolume[vtxIdx][coordIdx];
        }
    }
}
