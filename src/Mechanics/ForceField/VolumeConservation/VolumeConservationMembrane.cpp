#include "VolumeConservationMembrane.h"

#include "SysParams.h"

#include "Membrane.h"
#include "MMembrane.h"
#include "GMembrane.h"
#include "Vertex.h"

#include "VolumeConservationMembraneHarmonic.h"

template<>
double VolumeConservationMembrane<VolumeConservationMembraneHarmonic>::computeEnergy(double d) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {
        U_i = 0;

        if(d == 0.0) {
            double kComp = SysParams::Mechanics().CompConst;
            // TODO:
            // Currently kComp is "Compressibility", while we actually need some "incompressibility".

            double eqVolume = m->getMMembrane()->getEqVolume();

            double volume = m->getGMembrane()->getVolume();

            U_i += _FFTpye.energy(volume, kComp, eqVolume);

        } else {
            double kComp = SysParams::Mechanics().CompConst;
            // TODO:
            // Currently kComp is "Compressibility", while we actually need some "incompressibility".

            double eqVolume = m->getMMembrane()->getEqVolume();

            double stretchedVolume = m->getGMembrane()->getStretchedVolume();

            U_i += _FFTpye.energy(stretchedVolume, kComp, eqVolume, d);
        }

        if(fabs(U_i) == numeric_limits<double>::infinity()
            || U_i != U_i || U_i < -1.0) {
            _membraneCulprit = m;
            return -1;
        } else
            U += U_i;
        
    }

    return U;
}

template<>
void VolumeConservationMembrane<VolumeConservationMembraneHarmonic>::computeForces() {
    
    for (auto m: Membrane::getMembranes()) {
    
        double kComp = SysParams::Mechanics().CompConst;
        // TODO:
        // Currently kComp is "Compressibility", while we actually need some "incompressibility".

        double eqVolume = m->getMMembrane()->getEqVolume();

        double volume = m->getGMembrane()->getVolume();
        std::vector<std::array<double,3>>& dVolume = m->getGMembrane()->getDVolume();

        _FFType.forces(m->getVertexVector(), volume, dVolume, kComp, eqVolume);
    }
}

template<>
void VolumeConservationMembrane<VolumeConservationMembraneHarmonic>::computeForcesAux() {
    
    for (auto m: Membrane::getMembranes()) {
    
        double kComp = SysParams::Mechanics().CompConst;
        // TODO:
        // Currently kComp is "Compressibility", while we actually need some "incompressibility".

        double eqVolume = m->getMMembrane()->getEqVolume();

        double volume = m->getGMembrane()->getVolume();
        std::vector<std::array<double,3>>& dVolume = m->getGMembrane()->getDVolume();

        _FFType.forcesAux(m->getVertexVector(), volume, dVolume, kComp, eqVolume);
    }
}
