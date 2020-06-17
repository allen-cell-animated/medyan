#ifndef MEDYAN_Mechancis_ForceField_Membrane_Bending_hpp
#define MEDYAN_Mechancis_ForceField_Membrane_Bending_hpp

#include <algorithm> // transform
#include <functional> // plus
#include <type_traits> // is_same
#include <vector>

#include "common.h" // floatingpoint
#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Membrane/MembraneBendingHelfrich.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Vertex.hpp"

/// Represents a Filament bending interaction
template< typename InteractionType >
struct MembraneBending : public ForceField {

    InteractionType impl;

    virtual floatingpoint computeEnergy(floatingpoint *coord, bool stretched) override {
        using namespace std;

        double en = 0;

        for (auto m : Membrane::getMembranes()) {
            const auto &mesh = m->getMesh();

            for (const auto &v : mesh.getVertices()) {
                double enVertex = 0;

                if (v.numTargetingBorderHalfEdges == 0) {
                    const auto kBending = v.attr.vertex->mVertex.kBending;

                    const auto area = (stretched ? v.attr.gVertexS.astar : v.attr.gVertex.astar) / 3;

                    if constexpr(std::is_same_v< InteractionType, MembraneBendingHelfrich >) {
                        const auto eqCurv = v.attr.vertex->mVertex.eqCurv;
                        const auto curv = stretched ? v.attr.gVertexS.curv : v.attr.gVertex.curv;

                        enVertex += impl.energy(area, curv, kBending, eqCurv);
                    }
                    else {
                        const auto curv2 = stretched ? v.attr.gVertexS.curv2 : v.attr.gVertex.curv2;

                        enVertex += impl.energy(area, curv2, kBending);
                    }
                }

                if (!isfinite(enVertex)) {
                    LOG(ERROR) << "In " << getName() << " energy calculation, "
                        << "a vertex has energy " << enVertex;

                    return numeric_limits<double>::infinity();
                }
                else {
                    en += enVertex;
                }
            }

        } // End for membrane

        return en;
    }

    virtual void computeForces(floatingpoint *coord, floatingpoint *force) override {
        using namespace std;
        using MT = Membrane::MeshType;

        for (auto m : Membrane::getMembranes()) {

            medyan::assertValidIndexCacheForFF(m->getMesh());

            const auto &mesh = m->getMesh();
            const auto &cvt = mesh.metaAttribute().cachedVertexTopo;

            const size_t numVertices = mesh.numVertices();
            for (size_t i = 0; i < numVertices; ++i) {
                MT::VertexIndex vi {i};
                if (!mesh.isVertexOnBorder(vi)) {
                    const auto &va = mesh.attribute(vi);

                    const auto kBending = va.vertex->mVertex.kBending;
                    const auto eqCurv = va.vertex->mVertex.eqCurv;

                    const auto area = va.gVertex.astar / 3;
                    const auto curv = va.gVertex.curv;
                    const auto curv2 = va.gVertex.curv2;

                    if constexpr(std::is_same_v< InteractionType, MembraneBendingHelfrich >) {
                        impl.forces(
                            force + 3 * va.cachedCoordIndex,
                            area, va.gVertex.dAstar / 3,
                            curv, va.gVertex.dCurv,
                            kBending, eqCurv
                        );
                    }
                    else {
                        impl.forces(
                            force + 3 * va.cachedCoordIndex,
                            area, va.gVertex.dAstar / 3,
                            curv2, va.gVertex.dCurv2,
                            kBending
                        );
                    }

                    for (size_t n = 0; n < va.cachedDegree; ++n) {
                        const MT::HalfEdgeIndex hei_o { cvt[mesh.metaAttribute().cachedVertexOffsetLeavingHE(i) + n] };
                        const size_t vn_i = cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(i) + n];
                        const auto &dArea = mesh.attribute(hei_o).gHalfEdge.dNeighborAstar / 3;

                        if constexpr(std::is_same_v< InteractionType, MembraneBendingHelfrich >) {
                            const auto& dCurv = mesh.attribute(hei_o).gHalfEdge.dNeighborCurv;
                            impl.forces(
                                force + 3 * vn_i,
                                area, dArea, curv, dCurv, kBending, eqCurv
                            );
                        }
                        else {
                            const auto& dCurv2 = mesh.attribute(hei_o).gHalfEdge.dNeighborCurv2;
                            impl.forces(
                                force + 3 * vn_i,
                                area, dArea, curv2, dCurv2, kBending
                            );
                        }
                    }
                }
            }

        } // End for membrane

    }

    virtual string getName() override { return "Membrane Bending"; }


    // Useless overrides
    virtual void vectorize(const FFCoordinateStartingIndex& si) override {}
    virtual void cleanup() override {}
    virtual void computeLoadForces() override {}
    virtual void whoIsCulprit() override {}
    virtual std::vector<NeighborList*> getNeighborLists() override { return std::vector<NeighborList*>(); }
    virtual vector<std::string> getinteractionnames() override { return { getName() }; }
};

#endif
