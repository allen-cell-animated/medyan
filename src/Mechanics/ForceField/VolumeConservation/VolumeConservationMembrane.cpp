#include "Mechanics/ForceField/VolumeConservation/VolumeConservationMembrane.h"

#include "SysParams.h"

#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/MMembrane.h"
#include "Structure/SurfaceMesh/Vertex.h"

#include "Mechanics/ForceField/VolumeConservation/VolumeConservationMembraneHarmonic.h"

template<>
floatingpoint VolumeConservationMembrane<VolumeConservationMembraneHarmonic>::computeEnergy(const floatingpoint* coord, bool stretched) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {
        U_i = 0;

        const auto& mesh = m->getMesh();

        double kBulk = SysParams::Mechanics().BulkModulus;

        double eqVolume = m->getMMembrane()->getEqVolume();

        double volume = 0.0;
        for(const auto& t : mesh.getTriangles())
            volume += stretched ? t.attr.gTriangleS.coneVolume : t.attr.gTriangle.coneVolume;

        U_i += _FFType.energy(volume, kBulk, eqVolume);

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
void VolumeConservationMembrane<VolumeConservationMembraneHarmonic>::computeForces(const floatingpoint* coord, floatingpoint* force) {
    
    for (auto m: Membrane::getMembranes()) {

        const auto& mesh = m->getMesh();

        double kBulk = SysParams::Mechanics().BulkModulus;

        double eqVolume = m->getMMembrane()->getEqVolume();

        double volume = 0.0;
        for(const auto& t : mesh.getTriangles()) volume += t.attr.gTriangle.coneVolume;

        const size_t numVertices = mesh.getVertices().size();
        for(size_t vi = 0; vi < numVertices; ++vi) {
            Vertex* const v = mesh.getVertexAttribute(vi).vertex;
            const auto& dVolume = mesh.getVertexAttribute(vi).gVertex.dVolume;

            _FFType.forces(force + 3 * v->Bead::getIndex(), volume, dVolume, kBulk, eqVolume);
        }
    }
}
