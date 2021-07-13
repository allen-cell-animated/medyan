#ifndef MEDYAN_Mechanics_ForceField_VolumeConservation_VolConsrvMembrane_hpp
#define MEDYAN_Mechanics_ForceField_VolumeConservation_VolConsrvMembrane_hpp

#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/VolumeConservation/VolConsrvMembraneHarmonic.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"


// The force field of elasticity of volume enclosed by the membrane
struct VolumeConservationMembrane : public ForceField {

    VolumeConservationMembraneHarmonic impl;

    virtual void vectorize(const FFCoordinateStartingIndex& si) override {}

    virtual string getName() override { return "Volume conservation"; }

    virtual floatingpoint computeEnergy(floatingpoint* coord, bool stretched) override {
        using namespace std;

        double en = 0;

        for(auto m: Membrane::getMembranes()) if(m->isClosed()) {

            const auto& mesh = m->getMesh();

            const double kBulk = m->mMembrane.kVolume;
            const double eqVolume = m->mMembrane.eqVolume;

            double volume = 0.0;
            for(const auto& t : mesh.getTriangles())
                volume += stretched ? t.attr.gTriangleS.coneVolume : t.attr.gTriangle.coneVolume;

            const double enMem = impl.energy(volume, kBulk, eqVolume);

            if(!isfinite(enMem)) {
                LOG(ERROR) << "In " << getName() << " energy calculation, "
                    << "a membrane has energy " << enMem;

                return numeric_limits<double>::infinity();
            }
            else {
                en += enMem;
            }
            
        }

        return en;
    }

    virtual void computeForces(floatingpoint* coord, floatingpoint* force) override {
        using namespace std;
        using MT = Membrane::MeshType;

        for (auto m: Membrane::getMembranes()) if(m->isClosed()) {

            const auto& mesh = m->getMesh();
            medyan::assertValidIndexCacheForFF(mesh);

            const double kBulk = m->mMembrane.kVolume;
            const double eqVolume = m->mMembrane.eqVolume;

            double volume = 0.0;
            for(const auto& t : mesh.getTriangles()) volume += t.attr.gTriangle.coneVolume;

            const size_t numVertices = mesh.getVertices().size();
            for(MT::VertexIndex vi {0}; vi < numVertices; ++vi) {
                const auto& dVolume = mesh.attribute(vi).gVertex.dVolume;

                impl.forces(force + mesh.attribute(vi).cachedCoordIndex, volume, dVolume, kBulk, eqVolume);
            }
        }

    }

    // Useless overrides
    virtual void cleanup() override {}
    virtual void whoIsCulprit() override {}
    virtual void computeLoadForces() override {}
    virtual vector<NeighborList*> getNeighborLists() override { return {}; }
    virtual std::vector<std::string> getinteractionnames() override { return { getName() }; }
};

#endif
