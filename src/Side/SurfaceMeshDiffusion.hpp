#ifndef MEDYAN_Side_SurfaceMeshDiffusion_hpp
#define MEDYAN_Side_SurfaceMeshDiffusion_hpp

#include "Chemistry/ChemNRMImpl.h"
#include "GController.h"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/AdaptiveMesh.hpp"
#include "Structure/SurfaceMesh/SurfaceMeshGeneratorPreset.hpp"
#include "SysParams.h"

namespace medyan::side {

inline void surfaceMeshDiffusion() {
    LOG(NOTE) << "Running side procedure: surface mesh diffusion test.";


    // Background and other stuff
    //---------------------------------
    LOG(STEP) << "Initializing.";
    SubSystem subSystem;

    GeoParams geoParams;
    geoParams.compartmentSizeX = 500;
    geoParams.compartmentSizeY = 500;
    geoParams.compartmentSizeZ = 500;
    geoParams.NX = 1;
    geoParams.NY = 1;
    geoParams.NZ = 1;

    GController gController(&subSystem);
    gController.initializeGrid(geoParams);

    // Initialize membrane with surface proteins
    //---------------------------------
    MembraneSetup memSetup;
    memSetup.meshParam.push_back({
        "PLANE",
        "0", "0", "50",         // a point on the plane
        "0", "0", "1",          // facing z direction
        "10", "10", "0",        // box origin
        "480", "480", "100"     // box size
    });

    const auto meshData = mesh_gen::generateMeshViaParams< double >(memSetup.meshParam[0]);
    Membrane flatMembrane(
        &subSystem,
        memSetup,
        meshData.vertexCoordinateList,
        meshData.triangleList
    );

    adaptive_mesh::MembraneMeshAdapter meshAdapter(typename adaptive_mesh::MembraneMeshAdapter::Parameter {});
    meshAdapter.adapt(flatMembrane.getMesh());

    flatMembrane.updateGeometryValueForSystem();
    // The geometry of the membrane should be fixed beyond this point.

    LOG(INFO) << "Adding surface chemistry...";
    {
        MembraneMeshChemistryInfo memChemInfo {
            // names
            {"mem-diffu-test-a", "mem-diffu-test-b", "diffu-potential-test"},
            // diffusion (the 2nd species uses 0th category (custom))
            { { 0, 1.0 }, { 1, 0.5 }, { 2, 1.0, 0 } },
            // internal reactions
            {
                { {}, {0}, 0.055 },
                { {}, {1}, 0.062 },
                { {0, 1, 1}, {1, 1, 1}, 1.0 }
            },
            // energy categories
            { MembraneMeshChemistryInfo::SpeciesEnergyCat::custom },
        };
        flatMembrane.setChemistry(memChemInfo);
    }


    // Setup vertex energies
    for (auto& v : flatMembrane.getMesh().getVertices()) {
        auto& vertex = *v.attr.vertex;
        vertex.cVertex.energies.resize(flatMembrane.getMesh().metaAttribute().chemInfo.speciesEnergyCat.size());

        vertex.cVertex.energies[0] = -5 * std::sin(2 * M_PI / 1000 * vertex.coord[0]);
    }

    // Update and fix reaction rates
    setReactionRates(flatMembrane.getMesh());

    // Set initial species. 10000 copies of species 2 on vertex 0.
    flatMembrane.getMesh().getVertices()[0].attr.vertex->cVertex.species.findSpeciesByIndex(2)->setN(10000);


    // Initialize chem sim
    //---------------------------------
    ChemNRMImpl nrm;

    // Activate all reactions
    medyan::forEachReactionInMesh(
        flatMembrane.getMesh(),
        [&](ReactionDy& r) {
            nrm.addReaction(&r);
            r.activateReaction();
        }
    );



    // Start simulation
    //---------------------------------
    LOG(STEP) << "Starting simulation.";

    for (int i = 0; i < 10; ++i) {
        nrm.run(100);
        LOG(INFO) << "Current time: " << nrm.getTime();
    }

}

} // namespace medyan::side

#endif
