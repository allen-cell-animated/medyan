#include "VisualHelper.hpp"

#include <array>
#include <functional> // cref, reference_wrapper
#include <mutex>
#include <numeric> // iota
#include <thread>
#include <utility> // move
#include <vector>

#include "MathFunctions.h"
#include "Structure/Bead.h"
#include "Structure/Cylinder.h"
#include "Structure/Filament.h"
#include "Structure/Linker.h"
#include "Structure/MotorGhost.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Visual/Render/PathExtrude.hpp"
#include "Visual/Render/Sphere.hpp"
#include "Visual/SharedData.hpp"
#include "Visual/VisualElement.hpp"

namespace visual {

namespace {

struct SystemDataForVisual {
    struct MembraneIndex {
        std::vector< size_t > vertexIndices; // Vertex indexing in bead data
        std::vector< std::array< size_t, 3 > > triangleVertexIndices; // Index for this membrane only (start from 0)
    };

    std::mutex me;

    sys_data_update::FlagType updated;

    BeadData copiedBeadData;

    std::vector< MembraneIndex > membraneIndices;
    std::vector< std::vector< size_t > > filamentIndices; // [Filament Idx][Bead position in filament]

    std::vector< std::array< mathfunc::Vec3, 2 > > linkerCoords;
    std::vector< std::array< mathfunc::Vec3, 2 > > motorCoords;
};

// Shared data
SystemDataForVisual sdfv;

// This function tranforms the extracted system data to the actual gl
// compatible data structures according to the settings.
// Note:
//   - This function should not touch values used by the simulation.
//   - sdfv must not be changed when this function is active.
//   - When called from another thread, the shared_ptr must be copied to avoid
//     the underlying element being deleted.
void prepareVisualElement(const std::shared_ptr< VisualElement >& ve) {
    using namespace mathfunc;

    std::lock_guard< std::mutex > guard(ve->me);

    if(!ve->profile.enabled) return;

    // Temporary values
    std::size_t curVertexStart = 0; // current filled vertex index in the final vertex attribute array

    if(ve->profile.flag & Profile::targetMembrane) {
        // Membrane
        if(ve->profile.flag & Profile::displayForce) {
            //-----------------------------------------------------------------
            // Membrane Force
            //-----------------------------------------------------------------
            if(sdfv.updated & sys_data_update::BeadPosition) {
                ve->state.vertexAttribs.clear();
                ve->state.attribChanged = true;
                // if(sdfv.updated & sys_data_update::BeadConnection) {
                //     ve->state.vertexIndices.clear();
                //     ve->state.indexChanged = true;
                // }
                const auto& forces = (ve->profile.flag & Profile::forceUseSearchDir) ? sdfv.copiedBeadData.forces : sdfv.copiedBeadData.forcesAux;

                for(const auto& mi : sdfv.membraneIndices) {
                    ve->state.vertexAttribs.reserve(ve->state.vertexAttribs.size() + 2 * GlState::vaStride * mi.vertexIndices.size());
                    for(size_t i : mi.vertexIndices) {
                        const auto coord = sdfv.copiedBeadData.coords[i];
                        ve->state.vertexAttribs.push_back(coord[0]);
                        ve->state.vertexAttribs.push_back(coord[1]);
                        ve->state.vertexAttribs.push_back(coord[2]);
                        ve->state.vertexAttribs.resize(ve->state.vertexAttribs.size() + GlState::vaNormalSize); // dummy normal
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.x);
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.y);
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.z);

                        const auto force = forces[i];
                        const auto forceTip = force * ve->profile.forceScale + coord;
                        ve->state.vertexAttribs.push_back(forceTip[0]);
                        ve->state.vertexAttribs.push_back(forceTip[1]);
                        ve->state.vertexAttribs.push_back(forceTip[2]);
                        ve->state.vertexAttribs.resize(ve->state.vertexAttribs.size() + GlState::vaNormalSize); // dummy normal
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.x);
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.y);
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.z);
                    }
                }

                // const auto numBeads = ve->state.vertexAttribs.size() / GlState::vaStride;
                // if(sdfv.updated & sys_data_update::BeadConnection) {
                //     ve->state.vertexIndices.resize(numBeads);
                //     std::iota(ve->state.vertexIndices.begin(), ve->state.vertexIndices.end(), 0u);
                // }

            }
            ve->state.eleMode = GL_LINES;

        } else {
            //-----------------------------------------------------------------
            // Membrane Shape
            //-----------------------------------------------------------------
            if(sdfv.updated & sys_data_update::BeadPosition) {
                ve->state.vertexAttribs.clear();
                ve->state.attribChanged = true;
                // if(sdfv.updated & sys_data_update::BeadConnection) {
                //     ve->state.vertexIndices.clear();
                //     ve->state.indexChanged = true;
                // }

                for(const auto& mi : sdfv.membraneIndices) {
                    // Update coords
                    ve->state.vertexAttribs.reserve(ve->state.vertexAttribs.size() + 3 * GlState::vaStride * mi.triangleVertexIndices.size());
                    for(const auto& t : mi.triangleVertexIndices) {
                        const decltype(sdfv.copiedBeadData.coords[0]) coord[] {
                            sdfv.copiedBeadData.coords[mi.vertexIndices[t[0]]],
                            sdfv.copiedBeadData.coords[mi.vertexIndices[t[1]]],
                            sdfv.copiedBeadData.coords[mi.vertexIndices[t[2]]]
                        };
                        const auto un = normalizedVector(cross(coord[1] - coord[0], coord[2] - coord[0]));

                        for(size_t i = 0; i < 3; ++i) {
                            ve->state.vertexAttribs.push_back(coord[i][0]);
                            ve->state.vertexAttribs.push_back(coord[i][1]);
                            ve->state.vertexAttribs.push_back(coord[i][2]);
                            ve->state.vertexAttribs.push_back(un[0]);
                            ve->state.vertexAttribs.push_back(un[1]);
                            ve->state.vertexAttribs.push_back(un[2]);
                            ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.x);
                            ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.y);
                            ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.z);
                        }
                    }

                    // if(sdfv.updated & sys_data_update::BeadConnection) {
                    //     // update indices
                    //     ve->state.vertexIndices.reserve(ve->state.vertexIndices.size() + 3 * mi.triangleVertexIndices.size());
                    //     for(const auto& t : mi.triangleVertexIndices) {
                    //         ve->state.vertexIndices.push_back(t[0] + curVertexStart);
                    //         ve->state.vertexIndices.push_back(t[1] + curVertexStart);
                    //         ve->state.vertexIndices.push_back(t[2] + curVertexStart);
                    //     }
                    // }

                    // curVertexStart += mi.vertexIndices.size();
                }
            }
            ve->state.eleMode = GL_TRIANGLES;
        }
    }
    else if(ve->profile.flag & Profile::targetFilament) {
        // Filament
        if(ve->profile.flag & Profile::displayForce) {
            //-----------------------------------------------------------------
            // Filament Force
            //-----------------------------------------------------------------
            if(sdfv.updated & sys_data_update::BeadPosition) {
                ve->state.vertexAttribs.clear();
                ve->state.attribChanged = true;
                // if(sdfv.updated & sys_data_update::BeadConnection) {
                //     ve->state.vertexIndices.clear();
                //     ve->state.indexChanged = true;
                // }
                const auto& forces = (ve->profile.flag & Profile::forceUseSearchDir) ? sdfv.copiedBeadData.forces : sdfv.copiedBeadData.forcesAux;

                for(const auto& fi : sdfv.filamentIndices) {
                    ve->state.vertexAttribs.reserve(ve->state.vertexAttribs.size() + 2 * GlState::vaStride * fi.size());
                    for(size_t i : fi) {
                        const auto coord = sdfv.copiedBeadData.coords[i];
                        ve->state.vertexAttribs.push_back(coord[0]);
                        ve->state.vertexAttribs.push_back(coord[1]);
                        ve->state.vertexAttribs.push_back(coord[2]);
                        ve->state.vertexAttribs.resize(ve->state.vertexAttribs.size() + GlState::vaNormalSize); // dummy normal
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.x);
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.y);
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.z);

                        const auto force = forces[i];
                        const auto forceTip = force * ve->profile.forceScale + coord;
                        ve->state.vertexAttribs.push_back(forceTip[0]);
                        ve->state.vertexAttribs.push_back(forceTip[1]);
                        ve->state.vertexAttribs.push_back(forceTip[2]);
                        ve->state.vertexAttribs.resize(ve->state.vertexAttribs.size() + GlState::vaNormalSize); // dummy normal
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.x);
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.y);
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.z);
                    }
                }

                // const auto numBeads = ve->state.vertexAttribs.size() / 3; // 3 means coord(xyz);
                // if(sdfv.updated & sys_data_update::BeadConnection) {
                //     ve->state.vertexIndices.resize(numBeads);
                //     std::iota(ve->state.vertexIndices.begin(), ve->state.vertexIndices.end(), 0u);
                // }

            }
            ve->state.eleMode = GL_LINES;

        } else {
            //-----------------------------------------------------------------
            // Filament Shape
            //-----------------------------------------------------------------
            if(sdfv.updated & sys_data_update::BeadPosition) {
                ve->state.vertexAttribs.clear();
                ve->state.attribChanged = true;
                // if(sdfv.updated & sys_data_update::BeadConnection) {
                //     ve->state.vertexIndices.clear();
                //     ve->state.indexChanged = true;
                // }

                for(const auto& fi : sdfv.filamentIndices) {
                    mathfunc::VecArray< 3, float > genVertices;
                    std::vector< unsigned > genTriInd;

                    std::tie(genVertices, genTriInd) = visual::PathExtrude<float>{
                        ve->profile.pathExtrudeRadius,
                        ve->profile.pathExtrudeSides
                    }.generate(sdfv.copiedBeadData.coords, fi);

                    // Update coords
                    ve->state.vertexAttribs.reserve(ve->state.vertexAttribs.size() + GlState::vaStride * genTriInd.size());
                    const auto numTriangles = genTriInd.size() / 3;
                    for(size_t t = 0; t < numTriangles; ++t) {
                        const decltype(genVertices[0]) coord[] {
                            genVertices[genTriInd[3 * t + 0]],
                            genVertices[genTriInd[3 * t + 1]],
                            genVertices[genTriInd[3 * t + 2]]
                        };
                        const auto un = normalizedVector(cross(coord[1] - coord[0], coord[2] - coord[0]));

                        for(size_t i = 0; i < 3; ++i) {
                            ve->state.vertexAttribs.push_back(coord[i][0]);
                            ve->state.vertexAttribs.push_back(coord[i][1]);
                            ve->state.vertexAttribs.push_back(coord[i][2]);
                            ve->state.vertexAttribs.push_back(un[0]);
                            ve->state.vertexAttribs.push_back(un[1]);
                            ve->state.vertexAttribs.push_back(un[2]);
                            ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.x);
                            ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.y);
                            ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.z);
                        }
                    }

                    // if(sdfv.updated & sys_data_update::BeadConnection) {
                    //     // Update indices
                    //     ve->state.vertexIndices.reserve(ve->state.vertexIndices.size() + genTriInd.size());
                    //     for(auto i : genTriInd) {
                    //         ve->state.vertexIndices.push_back(i + curVertexStart);
                    //     }
                    // }

                    // curVertexStart += genVertices.size();
                }
            }
            ve->state.eleMode = GL_TRIANGLES;
        }
    }
    else if(ve->profile.flag & (Profile::targetLinker | Profile::targetMotor)) {
        //-----------------------------------------------------------------
        // Linker Shape
        //-----------------------------------------------------------------
        if(sdfv.updated & sys_data_update::BeadPosition) {
            ve->state.vertexAttribs.clear();
            ve->state.attribChanged = true;
            // if(sdfv.updated & sys_data_update::BeadConnection) {
            //     ve->state.vertexIndices.clear();
            //     ve->state.indexChanged = true;
            // }

            const auto& coords = (ve->profile.flag & Profile::targetLinker) ? sdfv.linkerCoords : sdfv.motorCoords;
            for(const auto& c : coords) {
                mathfunc::VecArray< 3, float > genVertices;
                std::vector< unsigned > genTriInd;

                std::tie(genVertices, genTriInd) = visual::PathExtrude<float>{
                    ve->profile.pathExtrudeRadius,
                    ve->profile.pathExtrudeSides
                }.generate(c, std::array<size_t, 2>{0, 1});

                // Update coords
                ve->state.vertexAttribs.reserve(ve->state.vertexAttribs.size() + GlState::vaStride * genTriInd.size());
                const auto numTriangles = genTriInd.size() / 3;
                for(size_t t = 0; t < numTriangles; ++t) {
                    const decltype(genVertices[0]) coord[] {
                        genVertices[genTriInd[3 * t + 0]],
                        genVertices[genTriInd[3 * t + 1]],
                        genVertices[genTriInd[3 * t + 2]]
                    };
                    const auto un = normalizedVector(cross(coord[1] - coord[0], coord[2] - coord[0]));

                    for(size_t i = 0; i < 3; ++i) {
                        ve->state.vertexAttribs.push_back(coord[i][0]);
                        ve->state.vertexAttribs.push_back(coord[i][1]);
                        ve->state.vertexAttribs.push_back(coord[i][2]);
                        ve->state.vertexAttribs.push_back(un[0]);
                        ve->state.vertexAttribs.push_back(un[1]);
                        ve->state.vertexAttribs.push_back(un[2]);
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.x);
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.y);
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.z);
                    }
                }

                // if(sdfv.updated & sys_data_update::BeadConnection) {
                //     // Update indices
                //     ve->state.vertexIndices.reserve(ve->state.vertexIndices.size() + genTriInd.size());
                //     for(auto i : genTriInd) {
                //         ve->state.vertexIndices.push_back(i + curVertexStart);
                //     }
                // }

                // curVertexStart += genVertices.size();
            }
        } // End if updated bead position
        ve->state.eleMode = GL_TRIANGLES;

    } // End if profile target

} // void prepareVisualElement(...)

void helper() {
    std::lock_guard< std::mutex > sdfvGuard(sdfv.me);
    std::lock_guard< std::mutex > veGuard(shared::veMutex);

    for(const auto& spv : shared::visualElements) {
        // Current executed serially, but could be parallelized
        // If parallelized, the shared_ptr must be copied to the working thread
        prepareVisualElement(spv);
    }
} // void helper(...)

} // namespace (anonymous)

void copySystemDataAndRunHelper(sys_data_update::FlagType update) {
    {
        std::lock_guard< std::mutex > guard(sdfv.me);

        if(update & (sys_data_update::BeadPosition | sys_data_update::BeadConnection)) {
            // Copy bead data
            sdfv.copiedBeadData = Bead::getDbDataConst();
        }

        if(update & (sys_data_update::BeadConnection)) {

            // Extract membrane indexing
            sdfv.membraneIndices.clear();
            for(const Membrane* m : Membrane::getMembranes()) {
                const auto& mesh = m->getMesh();

                sdfv.membraneIndices.emplace_back();
                auto& mi = sdfv.membraneIndices.back();

                mi.vertexIndices.reserve(mesh.numVertices());
                for(const auto& v : mesh.getVertices()) {
                    mi.vertexIndices.push_back(v.attr.vertex->Bead::getIndex());
                }

                mi.triangleVertexIndices.reserve(mesh.numTriangles());
                for(const auto& t : mesh.getTriangles()) {
                    size_t vIndex = 0; // 0, 1, 2
                    mi.triangleVertexIndices.emplace_back();
                    mesh.forEachHalfEdgeInPolygon(t, [&](size_t hei) {
                        mi.triangleVertexIndices.back()[vIndex] = mesh.target(hei);
                        ++vIndex;
                    });
                }
            }

            // Extract filament indexing
            sdfv.filamentIndices.clear();
            for(Filament* f : Filament::getFilaments()) {
                sdfv.filamentIndices.emplace_back();
                auto& fi = sdfv.filamentIndices.back();

                const auto& cylinders = f->getCylinderVector(); // TODO make f const Filament*
                fi.reserve(cylinders.size() + 1);
                for(Cylinder* c : cylinders)
                    fi.push_back(c->getFirstBead()->getIndex()); // TODO make c const Cylinder*
                fi.push_back(cylinders.back()->getSecondBead()->getIndex());
            }

            // Extract motors and linkers
            sdfv.linkerCoords.clear();
            for(Linker* l : Linker::getLinkers()) {
                sdfv.linkerCoords.emplace_back();
                const auto pos0 = l->getFirstPosition();
                sdfv.linkerCoords.back()[0]
                    = l->getFirstCylinder()->getFirstBead()->coordinate() * (1 - pos0)
                    + l->getFirstCylinder()->getSecondBead()->coordinate() * pos0;
                const auto pos1 = l->getSecondPosition();
                sdfv.linkerCoords.back()[1]
                    = l->getSecondCylinder()->getFirstBead()->coordinate() * (1 - pos1)
                    + l->getSecondCylinder()->getSecondBead()->coordinate() * pos1;
            }
            sdfv.motorCoords.clear();
            for(MotorGhost* l : MotorGhost::getMotorGhosts()) {
                sdfv.motorCoords.emplace_back();
                const auto pos0 = l->getFirstPosition();
                sdfv.motorCoords.back()[0]
                    = l->getFirstCylinder()->getFirstBead()->coordinate() * (1 - pos0)
                    + l->getFirstCylinder()->getSecondBead()->coordinate() * pos0;
                const auto pos1 = l->getSecondPosition();
                sdfv.motorCoords.back()[1]
                    = l->getSecondCylinder()->getFirstBead()->coordinate() * (1 - pos1)
                    + l->getSecondCylinder()->getSecondBead()->coordinate() * pos1;
            }

            // Save updated
            sdfv.updated = update;
        }
    }

    // Launch helper thread (may use thread pool)
    std::thread(helper).detach();

}

} // namespace visual
