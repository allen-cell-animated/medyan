#ifndef MEDYAN_Mechanics_ForceField_VolumeConservation_MembraneHarmonic_Hpp
#define MEDYAN_Mechanics_ForceField_VolumeConservation_MembraneHarmonic_Hpp

#include "MathFunctions.h"

struct VolumeConservationMembraneHarmonic {

    double energy(double volume, double kBulk, double eqVolume) const {
        double dist = volume - eqVolume;
        return 0.5 * kBulk * dist * dist / eqVolume;
    }
    
    void forces(
        floatingpoint* force,
        double volume, const mathfunc::Vec3& dVolume,
        double kBulk, double eqVolume
    ) const {
        // F_i = -grad_i U = -k / V_0 * (V - V_0) * grad_i V
        // V(olume) and grad_i V(olume) are obtained as function parameters

        double coeff = -kBulk / eqVolume * (volume - eqVolume);
        for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
            force[coordIdx] += coeff * dVolume[coordIdx];
        }
    }
};


#endif
