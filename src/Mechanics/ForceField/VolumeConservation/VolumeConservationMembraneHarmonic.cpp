#include "VolumeConservationMembraneHarmonic.h"

#include "Structure/SurfaceMesh/Vertex.h"

double VolumeConservationMembraneHarmonic::energy(double volume, double kBulk, double eqVolume) {
    double dist = volume - eqVolume;
    return 0.5 * kBulk * dist * dist / eqVolume;
}

void VolumeConservationMembraneHarmonic::forces(
    double* force,
    double volume, const mathfunc::Vec3& dVolume,
    double kBulk, double eqVolume
) {
    // F_i = -grad_i U = -k / V_0 * (V - V_0) * grad_i V
    // V(olume) and grad_i V(olume) are obtained as function parameters

    double coeff = -kBulk / eqVolume * (volume - eqVolume);
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        force[coordIdx] += coeff * dVolume[coordIdx];
    }
}
