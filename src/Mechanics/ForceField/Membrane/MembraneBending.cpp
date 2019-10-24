#include "Mechanics/ForceField/Membrane/MembraneBending.hpp"

#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Vertex.hpp"
#include "Structure/SurfaceMesh/MVoronoiCell.h"

#include "Mechanics/ForceField/Membrane/MembraneBendingHelfrich.hpp"

// Using the Helfrich Hamiltonian of mean curvature in Voronoi cells
template< typename InteractionType >
floatingpoint MembraneBending< InteractionType >::computeEnergy(const floatingpoint* coord, bool stretched) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {
        const auto& mesh = m->getMesh();

        U_i = 0;

        for(const auto& v : mesh.getVertices()) if(v.numTargetingBorderHalfEdges == 0) {
            const auto kBending = v.attr.vertex->getMVoronoiCell()->getBendingModulus();
            const auto eqCurv = v.attr.vertex->getMVoronoiCell()->getEqCurv();

            const auto area = (stretched ? v.attr.gVertexS.astar : v.attr.gVertex.astar) / 3;
            const auto curv = stretched ? v.attr.gVertexS.curv : v.attr.gVertex.curv;

            U_i += _FFType.energy(area, curv, kBending, eqCurv);
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

template< typename InteractionType >
void MembraneBending< InteractionType >::computeForces(const floatingpoint* coord, floatingpoint* force) {
    
    for (auto m: Membrane::getMembranes()) {
    
        Membrane::MembraneMeshAttributeType::cacheIndices(m->getMesh());

        const auto& mesh = m->getMesh();
        const auto& cvt = mesh.getMetaAttribute().cachedVertexTopo;

        const size_t numVertices = mesh.getVertices().size();
        for(size_t vi = 0; vi < numVertices; ++vi) if(!mesh.isVertexOnBorder(vi)) {
            const auto& va = mesh.getVertexAttribute(vi);

            const auto kBending = va.vertex->getMVoronoiCell()->getBendingModulus();
            const auto eqCurv = va.vertex->getMVoronoiCell()->getEqCurv();

            const auto area = va.gVertex.astar / 3;
            const auto curv = va.gVertex.curv;

            _FFType.forces(
                force + 3 * va.cachedCoordIndex,
                area, va.gVertex.dAstar / 3,
                curv, va.gVertex.dCurv,
                kBending, eqCurv
            );

            for(size_t i = 0; i < va.cachedDegree; ++i) {
                const size_t hei_o = cvt[mesh.getMetaAttribute().cachedVertexOffsetLeavingHE(vi) + i];
                const size_t vn_i = cvt[mesh.getMetaAttribute().cachedVertexOffsetNeighborCoord(vi) + i];
                const auto& dArea = mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dNeighborAstar / 3;
                const auto& dCurv = mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dNeighborCurv;
                _FFType.forces(
                    force + 3 * vn_i,
                    area, dArea, curv, dCurv, kBending, eqCurv
                );
            }
        }

    }
}

// Explicit instantiations
template class MembraneBending< MembraneBendingHelfrich >;
