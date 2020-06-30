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
#include "Structure/BranchingPoint.h"
#include "Structure/Cylinder.h"
#include "Structure/Filament.h"
#include "Structure/Linker.h"
#include "Structure/MotorGhost.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "SysParams.h"
#include "Visual/Render/PathExtrude.hpp"
#include "Visual/Render/Sphere.hpp"
#include "Visual/VisualElement.hpp"
#include "Visual/Window.hpp"

namespace visual {

std::weak_ptr< VisualDisplay > vdWeak;

namespace {

// Shared data
visual::SystemRawData sdfv;

// This function tranforms the extracted system data to the actual gl
// compatible data structures according to the settings.
// Note:
//   - This function should not touch values used by the simulation.
//   - sdfv must not be changed when this function is active.
//   - When called from another thread, the shared_ptr must be copied to avoid
//     the underlying element being deleted.
void prepareVisualElement(const std::shared_ptr< VisualElement >& ve) {
    using namespace std;
    using namespace mathfunc;

    std::lock_guard< std::mutex > guard(ve->me);

    if(!ve->profile.enabled) return;

    // Temporary values
    std::size_t curVertexStart = 0; // current filled vertex index in the final vertex attribute array

    if(ve->profile.flag & Profile::targetMembrane) {
        // Membrane
        if(ve->profile.flag & Profile::displayForce) {
            LOG(WARNING) << "Force display is currently unavaiable.";
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

                for(const auto& mi : sdfv.membraneData) {
                    // Update coords
                    ve->state.vertexAttribs.reserve(ve->state.vertexAttribs.size() + 3 * ve->state.size.vaStride * mi.triangleVertexIndices.size());
                    for(const auto& t : mi.triangleVertexIndices) {
                        const Vec3 coord[] {
                            mi.vertexCoords[t[0]],
                            mi.vertexCoords[t[1]],
                            mi.vertexCoords[t[2]]
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
                }
            }
            ve->state.eleMode = GL_TRIANGLES;
        }
    }
    else if(ve->profile.flag & Profile::targetFilament) {
        // Filament
        if(ve->profile.flag & Profile::displayForce) {
            LOG(WARNING) << "Force display is currently unavailable.";
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

                switch(ve->profile.pathMode) {

                case Profile::PathMode::Line:
                    // TODO implement it
                    break;

                case Profile::PathMode::Extrude:
                    for(const auto& fi : sdfv.filamentData) {
                        std::vector< mathfunc::Vec3f > genVertices;
                        std::vector< mathfunc::Vec3f > genVertexNormals;
                        std::vector< std::array< size_t, 3 > > genTriInd;

                        const auto nb = fi.beadCoords.size();
                        vector< size_t > trivial_indices(nb);
                        iota(begin(trivial_indices), end(trivial_indices), 0);
                        std::tie(genVertices, genVertexNormals, genTriInd) = visual::PathExtrude<float>{
                            ve->profile.pathExtrudeRadius,
                            ve->profile.pathExtrudeSides
                        }.generate(fi.beadCoords, trivial_indices);

                        // Update coords
                        ve->state.vertexAttribs.reserve(ve->state.vertexAttribs.size() + ve->state.size.vaStride * 3 * genTriInd.size());
                        const auto numTriangles = genTriInd.size();
                        for(size_t t = 0; t < numTriangles; ++t) {
                            const auto& triInds = genTriInd[t];
                            const mathfunc::Vec3f coord[] {
                                genVertices[triInds[0]],
                                genVertices[triInds[1]],
                                genVertices[triInds[2]]
                            };
                            const mathfunc::Vec3f un[] {
                                genVertexNormals[triInds[0]],
                                genVertexNormals[triInds[1]],
                                genVertexNormals[triInds[2]]
                            };

                            for(size_t i = 0; i < 3; ++i) {
                                ve->state.vertexAttribs.push_back(coord[i][0]);
                                ve->state.vertexAttribs.push_back(coord[i][1]);
                                ve->state.vertexAttribs.push_back(coord[i][2]);
                                ve->state.vertexAttribs.push_back(un[i][0]);
                                ve->state.vertexAttribs.push_back(un[i][1]);
                                ve->state.vertexAttribs.push_back(un[i][2]);
                                ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.x);
                                ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.y);
                                ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.z);
                            }
                        }

                    } // End loop filaments
                    break;

                case Profile::PathMode::Bead:
                    {
                        const auto sphereGen = visual::SphereUv<float> {
                            ve->profile.beadRadius,
                            ve->profile.beadLongitudeSegs,
                            ve->profile.beadLatitudeSegs
                        };
                        const auto sphereCache = sphereGen.makeCache();

                        for(const auto& fi : sdfv.filamentData) {
                            for(const auto& bc : fi.beadCoords) {
                                std::vector< Vec< 3, float > > genVertices;

                                std::tie(genVertices, std::ignore) = sphereGen.generate(
                                    {
                                        static_cast<float>(bc[0]),
                                        static_cast<float>(bc[1]),
                                        static_cast<float>(bc[2])
                                    },
                                    sphereCache
                                );

                                // Update coords
                                ve->state.vertexAttribs.reserve(ve->state.vertexAttribs.size() + ve->state.size.vaStride * 3 * sphereCache.triInd.size());
                                const auto numTriangles = sphereCache.triInd.size();
                                for(size_t t = 0; t < numTriangles; ++t) {
                                    const typename decltype(genVertices)::value_type coord[] {
                                        genVertices[sphereCache.triInd[t][0]],
                                        genVertices[sphereCache.triInd[t][1]],
                                        genVertices[sphereCache.triInd[t][2]]
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
                            } // End loop beads in a filament

                        } // End loop filaments
                    }
                    break;

                } // End switch path mode
            }
            ve->state.eleMode = GL_TRIANGLES;
        }
    }
    else if(ve->profile.flag & (Profile::targetLinker | Profile::targetMotor | Profile::targetBrancher)) {
        //-----------------------------------------------------------------
        // Linker, motor or brancher shape
        //-----------------------------------------------------------------
        if(sdfv.updated & sys_data_update::BeadPosition) {
            ve->state.vertexAttribs.clear();
            ve->state.attribChanged = true;
            // if(sdfv.updated & sys_data_update::BeadConnection) {
            //     ve->state.vertexIndices.clear();
            //     ve->state.indexChanged = true;
            // }

            const auto& coords = (ve->profile.flag & Profile::targetLinker) ? sdfv.linkerCoords :
                                 (ve->profile.flag & Profile::targetMotor ) ? sdfv.motorCoords  :
                                                                              sdfv.brancherCoords;
            for(const auto& c : coords) {
                std::vector< mathfunc::Vec3f > genVertices;
                std::vector< mathfunc::Vec3f > genVertexNormals;
                std::vector< std::array< size_t, 3 > > genTriInd;

                std::tie(genVertices, genVertexNormals, genTriInd) = visual::PathExtrude<float>{
                    ve->profile.pathExtrudeRadius,
                    ve->profile.pathExtrudeSides
                }.generate(c, std::array<size_t, 2>{0, 1});

                // Update coords
                ve->state.vertexAttribs.reserve(ve->state.vertexAttribs.size() + ve->state.size.vaStride * 3 * genTriInd.size());
                const auto numTriangles = genTriInd.size();
                for(size_t t = 0; t < numTriangles; ++t) {
                    const auto& triInds = genTriInd[t];
                    const mathfunc::Vec3f coord[] {
                        genVertices[triInds[0]],
                        genVertices[triInds[1]],
                        genVertices[triInds[2]]
                    };
                    const mathfunc::Vec3f un[] {
                        genVertexNormals[triInds[0]],
                        genVertexNormals[triInds[1]],
                        genVertexNormals[triInds[2]]
                    };

                    for(size_t i = 0; i < 3; ++i) {
                        ve->state.vertexAttribs.push_back(coord[i][0]);
                        ve->state.vertexAttribs.push_back(coord[i][1]);
                        ve->state.vertexAttribs.push_back(coord[i][2]);
                        ve->state.vertexAttribs.push_back(un[i][0]);
                        ve->state.vertexAttribs.push_back(un[i][1]);
                        ve->state.vertexAttribs.push_back(un[i][2]);
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.x);
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.y);
                        ve->state.vertexAttribs.push_back(ve->profile.colorAmbient.z);
                    }
                }

            }
        } // End if updated bead position
        ve->state.eleMode = GL_TRIANGLES;
    }
    else if(ve->profile.flag & Profile::targetCompartment) {
        //-----------------------------------------------------------------
        // Compartments
        //-----------------------------------------------------------------
        if(sdfv.updated & sys_data_update::Compartment) {
            ve->state.vertexAttribs.clear();
            ve->state.attribChanged = true;

            const auto genLines = [&](size_t fixedAxis, size_t ax1, size_t dax1, size_t ax2, size_t dax2) {
                for(size_t x1 = 0; x1 <= sdfv.compartmentNum[ax1]; x1 += dax1) {
                    for(size_t x2 = 0; x2 <= sdfv.compartmentNum[ax2]; x2 += dax2) {
                        const auto v1 = sdfv.compartmentSize[ax1] * x1;
                        const auto v2 = sdfv.compartmentSize[ax2] * x2;
                        Vec3f coord0, coord1;
                        coord0[fixedAxis] = 0;
                        coord0[ax1] = v1;
                        coord0[ax2] = v2;
                        coord1[fixedAxis] = sdfv.compartmentSize[fixedAxis] * sdfv.compartmentNum[fixedAxis];
                        coord1[ax1] = v1;
                        coord1[ax2] = v2;

                        ve->state.vertexAttribs.push_back(coord0[0]);
                        ve->state.vertexAttribs.push_back(coord0[1]);
                        ve->state.vertexAttribs.push_back(coord0[2]);
                        ve->state.vertexAttribs.push_back(ve->profile.colorDiffuse.x);
                        ve->state.vertexAttribs.push_back(ve->profile.colorDiffuse.y);
                        ve->state.vertexAttribs.push_back(ve->profile.colorDiffuse.z);
                        ve->state.vertexAttribs.push_back(coord1[0]);
                        ve->state.vertexAttribs.push_back(coord1[1]);
                        ve->state.vertexAttribs.push_back(coord1[2]);
                        ve->state.vertexAttribs.push_back(ve->profile.colorDiffuse.x);
                        ve->state.vertexAttribs.push_back(ve->profile.colorDiffuse.y);
                        ve->state.vertexAttribs.push_back(ve->profile.colorDiffuse.z);
                    }
                }
            };

            switch(ve->profile.gridMode) {
            case Profile::GridMode::Boundary:
                genLines(0, 1, sdfv.compartmentNum[1], 2, sdfv.compartmentNum[2]);
                genLines(1, 2, sdfv.compartmentNum[2], 0, sdfv.compartmentNum[0]);
                genLines(2, 0, sdfv.compartmentNum[0], 1, sdfv.compartmentNum[1]);
                break;
            case Profile::GridMode::Mesh:
                genLines(0, 1, 1, 2, 1);
                genLines(1, 2, 1, 0, 1);
                genLines(2, 0, 1, 1, 1);
                break;
            }

        } // End if updated compartment
        ve->state.eleMode = GL_LINES;

    } // End if profile target

} // void prepareVisualElement(...)

void helper() {
    std::lock_guard< std::mutex > sdfvGuard(sdfv.me);

    if(auto vd = vdWeak.lock()) {
        for(auto& vp : vd->vps) {
            std::lock_guard< std::mutex > veGuard(vp.veMutex);

            for(const auto& spv : vp.visualElements) {
                // Current executed serially, but could be parallelized
                // If parallelized, the shared_ptr must be copied to the working thread
                prepareVisualElement(spv);
            }
        }
    }
} // void helper(...)

} // namespace (anonymous)

void copySystemDataAndRunHelper(sys_data_update::FlagType update) {
    visual::copySystemData(sdfv, update);

    // Launch helper thread (may use thread pool)
    std::thread(helper).detach();

}

} // namespace visual
