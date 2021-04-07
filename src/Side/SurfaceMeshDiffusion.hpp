#ifndef MEDYAN_Side_SurfaceMeshDiffusion_hpp
#define MEDYAN_Side_SurfaceMeshDiffusion_hpp

#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/SurfaceMeshGeneratorPreset.hpp"
#include "SysParams.h"

namespace medyan::side {

inline void surfaceMeshDiffusion() {
    LOG(NOTE) << "Running side procedure: surface mesh diffusion test.";

    LOG(STEP) << "Initializing.";

    // Background and other stuff
    SubSystem subSystem;

    // Initialize membrane with surface proteins
    MembraneSetup memSetup;
    memSetup.meshParam.push_back({ });

    const auto meshData = mesh_gen::generateMeshViaParams< double >(memSetup.meshParam[0]);
    Membrane flatMembrane(
        &subSystem,
        memSetup,
        meshData.vertexCoordinateList,
        meshData.triangleList
    );

    // int numMembranes = 0;
    // const auto addMembrane = [this, &numMembranes](const MembraneSetup& memSetup, const MembraneParser::MembraneInfo& memData) {
    //     auto newMembrane = _subSystem.addTrackable<Membrane>(
    //         &_subSystem,
    //         memSetup,
    //         memData.vertexCoordinateList,
    //         memData.triangleVertexIndexList
    //     );

    //     // Optimize the mesh for membrane
    //     remeshMembrane(*_meshAdapter, *newMembrane);

    //     // Set up mechanics
    //     newMembrane->initMechanicParams(memSetup);

    //     ++numMembranes;
    // };

    // for(auto& memSetup : membraneSettings.setupVec) {

    //     for(auto& initParams : memSetup.meshParam) {
    //         if(initParams.size() == 2 && initParams[0] == "file") {
    //             // The input looks like this: init file path/to/file
    //             // Read membrane mesh information from an external file.
    //             auto memPath = _inputDirectory / std::filesystem::path(initParams[1]);
    //             std::ifstream ifs(memPath);
    //             if (!ifs.is_open()) {
    //                 LOG(ERROR) << "Cannot open membrane file " << memPath;
    //                 throw std::runtime_error("Cannot open membrane file.");
    //             }

    //             const auto memDataVec = MembraneParser::readMembranes(ifs);

    //             for(auto& memData : memDataVec) {
    //                 addMembrane(memSetup, memData);
    //             }
    //         }
    //         else {
    //             // Forward the input to the membrane mesh initializer
    //             const auto newMesh = mesh_gen::generateMeshViaParams< floatingpoint >(initParams);

    //             addMembrane(memSetup, {newMesh.vertexCoordinateList, newMesh.triangleList});
    //         }
    //     }
    // }

    // LOG(INFO) << "Done. " << numMembranes << " membranes created." << endl;

    // // Create a region inside the membrane
    // LOG(INFO) << "Creating membrane regions...";
    // _regionInMembrane = (
    //     numMembranes == 0 ?
    //     make_unique<MembraneRegion<Membrane>>(_subSystem.getBoundary()) :
    //     MembraneRegion<Membrane>::makeByChildren(MembraneHierarchy< Membrane >::root())
    // );
    // _subSystem.setRegionInMembrane(_regionInMembrane.get());


    // LOG(INFO) << "Adding surface chemistry...";
    // {
    //     MembraneMeshChemistryInfo memChemInfo {
    //         // names
    //         {"mem-diffu-test-a", "mem-diffu-test-b"},
    //         // diffusion
    //         { { 0, 1.0 }, { 1, 0.5 } },
    //         // internal reactions
    //         {
    //             { {}, {0}, 0.055 },
    //             { {}, {1}, 0.062 },
    //             { {0, 1, 1}, {1, 1, 1}, 1.0 }
    //         }
    //     };
    //     for(auto m : Membrane::getMembranes()) {
    //         m->setChemistry(memChemInfo);
    //     }
    // }

    // LOG(INFO) << "Adjusting compartments by membranes...";

    // // Deactivate all the compartments outside membrane, and mark boundaries as interesting
    // for(auto c : _subSystem.getCompartmentGrid()->getCompartments()) {
    //     if(!c->getTriangles().empty()) {
    //         // Contains triangles, so this compartment is at the boundary.
    //         c->boundaryInteresting = true;

    //         // Update partial activate status
    //         c->computeSlicedVolumeArea(Compartment::SliceMethod::membrane);
    //         _cController.updateActivation(c, Compartment::ActivateReason::Membrane);

    //     } else if( ! _regionInMembrane->contains(vector2Vec<3, floatingpoint>(c->coordinates()))) {
    //         // Compartment is outside the membrane
    //         _cController.deactivate(c, true);
    //     }
    // }

}

} // namespace medyan::side

#endif
