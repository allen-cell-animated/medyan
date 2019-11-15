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
#include "Structure/SurfaceMesh/MVoronoiCell.h"
#include "Structure/SurfaceMesh/Vertex.hpp"

/// Represents a Filament bending interaction
template< typename InteractionType >
class MembraneBending : public ForceField {

private:
    InteractionType bending_;

    // Culprit membrane
    const Membrane* membraneCulprit_ = nullptr;

public:
    virtual floatingpoint computeEnergy(floatingpoint *coord, bool stretched) override {
        double U = 0;
        double U_i;

        for (auto m : Membrane::getMembranes()) {
            const auto &mesh = m->getMesh();

            U_i = 0;

            for (const auto &v : mesh.getVertices())
                if (v.numTargetingBorderHalfEdges == 0) {
                    const auto kBending = v.attr.vertex->getMVoronoiCell()->getBendingModulus();

                    const auto area = (stretched ? v.attr.gVertexS.astar : v.attr.gVertex.astar) / 3;

                    if constexpr(std::is_same_v< InteractionType, MembraneBendingHelfrich >) {
                        const auto eqCurv = v.attr.vertex->getMVoronoiCell()->getEqCurv();
                        const auto curv = stretched ? v.attr.gVertexS.curv : v.attr.gVertex.curv;

                        U_i += bending_.energy(area, curv, kBending, eqCurv);
                    }
                    else {
                        const auto curv2 = stretched ? v.attr.gVertexS.curv2 : v.attr.gVertex.curv2;

                        U_i += bending_.energy(area, curv2, kBending);
                    }
                }

            if (fabs(U_i) == numeric_limits<double>::infinity() || U_i != U_i || U_i < -1.0) {
                membraneCulprit_ = m;
                return -1;
            }
            else
                U += U_i;

        } // End for membrane

        return U;
    }

    virtual void computeForces(floatingpoint *coord, floatingpoint *force) override {

        // Configure force buffer
        constexpr bool useForceBuffer = true;
        const std::size_t dof = Bead::getDbDataConst().coords.size_raw();

        if (useForceBuffer) {
            forceBuffer_.assign(dof, 0.0);
        }
        floatingpoint* const f = useForceBuffer ? forceBuffer_.data() : force;

        for (auto m : Membrane::getMembranes()) {

            Membrane::MembraneMeshAttributeType::cacheIndices(m->getMesh());

            const auto &mesh = m->getMesh();
            const auto &cvt = mesh.getMetaAttribute().cachedVertexTopo;

            const size_t numVertices = mesh.getVertices().size();
            for (size_t vi = 0; vi < numVertices; ++vi)
                if (!mesh.isVertexOnBorder(vi)) {
                    const auto &va = mesh.getVertexAttribute(vi);

                    const auto kBending = va.vertex->getMVoronoiCell()->getBendingModulus();
                    const auto eqCurv = va.vertex->getMVoronoiCell()->getEqCurv();

                    const auto area = va.gVertex.astar / 3;
                    const auto curv = va.gVertex.curv;
                    const auto curv2 = va.gVertex.curv2;

                    if constexpr(std::is_same_v< InteractionType, MembraneBendingHelfrich >) {
                        bending_.forces(
                            f + 3 * va.cachedCoordIndex,
                            area, va.gVertex.dAstar / 3,
                            curv, va.gVertex.dCurv,
                            kBending, eqCurv
                        );
                    }
                    else {
                        bending_.forces(
                            f + 3 * va.cachedCoordIndex,
                            area, va.gVertex.dAstar / 3,
                            curv2, va.gVertex.dCurv2,
                            kBending
                        );
                    }

                    for (size_t i = 0; i < va.cachedDegree; ++i) {
                        const size_t hei_o = cvt[mesh.getMetaAttribute().cachedVertexOffsetLeavingHE(vi) + i];
                        const size_t vn_i = cvt[mesh.getMetaAttribute().cachedVertexOffsetNeighborCoord(vi) + i];
                        const auto &dArea = mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dNeighborAstar / 3;

                        if constexpr(std::is_same_v< InteractionType, MembraneBendingHelfrich >) {
                            const auto& dCurv = mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dNeighborCurv;
                            bending_.forces(
                                f + 3 * vn_i,
                                area, dArea, curv, dCurv, kBending, eqCurv
                            );
                        }
                        else {
                            const auto& dCurv2 = mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dNeighborCurv2;
                            bending_.forces(
                                f + 3 * vn_i,
                                area, dArea, curv2, dCurv2, kBending
                            );
                        }
                    }
                }

        } // End for membrane

        if(useForceBuffer) {
            std::transform(
                force, force + dof, forceBuffer_.begin(),
                force,
                std::plus<>{}
            );
        }
    }

    virtual string getName() override { return "Membrane Bending"; }

    virtual void whoIsCulprit() override {
        if(membraneCulprit_) {
            LOG(INFO) << "Printing culprit membrane...";
            membraneCulprit_->printSelf();
        } else {
            LOG(ERROR) << "Membrane culprit is not set.";
        }
    }

    // Useless overrides
    virtual void vectorize() override {}
    virtual void cleanup() override {}
    virtual void computeLoadForces() override {}
    virtual std::vector<NeighborList*> getNeighborLists() override { return std::vector<NeighborList*>(); }
};

#endif
