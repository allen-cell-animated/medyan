#include "VisualHelper.hpp"

#include <functional> // cref, reference_wrapper
#include <mutex>
#include <thread>
#include <utility> // move
#include <vector>

#include "Structure/Bead.h"
#include "Structure/Filament.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
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

    std::vector< MembraneIndex > membraneIndices; // [Membrane Idx][Triangle Idx][Bead 0, 1 or 2]
};

SystemDataForVisual sdfv;

// This function tranforms the extracted system data to the actual gl
// compatible data structures according to the settings.
// Note:
//   - This function should not touch values used by the simulation.
//   - sdfv must not be changed when this function is active.
//   - When called from another thread, the shared_ptr must be copied to avoid
//     the underlying element being deleted.
void prepareVisualElement(const std::shared_ptr< VisualElement >& ve) {
    std::lock_guard< std::mutex > guard(ve->me);

    // Temporary values
    std::size_t curVertexStart = 0;

    if(ve->profile.flag & Profile::targetMembrane) {
        // Membrane
        if(sdfv.updated & sys_data_update::BeadPosition) {
            ve->state.vertexAttribs.clear();
            ve->state.attribChanged = true;
            if(sdfv.updated & sys_data_update::BeadConnection) {
                ve->state.vertexIndices.clear();
                ve->state.indexChanged = true;
            }

            for(const auto& mi : sdfv.membraneIndices) {
                // Update coords
                ve->state.vertexAttribs.reserve(ve->state.vertexAttribs.size() + 3 * mi.vertexIndices.size()); // 3 means coord(xyz)
                for(size_t i : mi.vertexIndices) {
                    const auto coord = sdfv.copiedBeadData.coords[i];
                    ve->state.vertexAttribs.push_back(coord[0]);
                    ve->state.vertexAttribs.push_back(coord[1]);
                    ve->state.vertexAttribs.push_back(coord[2]);
                }

                if(sdfv.updated & sys_data_update::BeadConnection) {
                    // update indices
                    ve->state.vertexIndices.reserve(ve->state.vertexIndices.size() + 3 * mi.triangleVertexIndices.size());
                    for(const auto& t : mi.triangleVertexIndices) {
                        ve->state.vertexIndices.push_back(t[0] + curVertexStart);
                        ve->state.vertexIndices.push_back(t[1] + curVertexStart);
                        ve->state.vertexIndices.push_back(t[2] + curVertexStart);
                    }
                }

                curVertexStart += mi.vertexIndices.size();
            }
        }
        ve->state.eleMode = GL_TRIANGLES;
    } // end target membrane

} // void prepareVisualElement(...)

void helper() {
    std::lock_guard< std::mutex > sdfvGuard(sdfv.me);
    std::lock_guard< std::mutex > veGuard(shared::veMutex);

    for(const auto& spv : shared::visualElements) {
        // Current executed serially, but could be parallelized
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

            // Copy filament data
            // TODO
        }
    }

    // Launch helper thread (may use thread pool)
    std::thread(helper).detach();

}

} // namespace visual
