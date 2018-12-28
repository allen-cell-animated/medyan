#include "MembraneBending.h"

#include "Membrane.h"
#include "Vertex.h"
#include "GVoronoiCell.h"
#include "MVoronoiCell.h"

#include "MembraneBendingVoronoiHelfrich.h"

// Using the Helfrich Hamiltonian of mean curvature in Voronoi cells
template<>
double MembraneBending<MembraneBendingVoronoiHelfrich>::computeEnergy(double d) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {
        const auto& mesh = m->getMesh();

        U_i = 0;

        if(d == 0.0) {
            for(const auto& v : mesh.getVertices()) {
                const auto kBending = v.vertex->getMVoronoiCell()->getBendingModulus();
                const auto eqCurv = v.vertex->getMVoronoiCell()->getEqCurv();

                const auto area = v.attr.gVertex.area;
                const auto curv = v.attr.gVertex.curv;

                U_i += _FFType.energy(area, curv, kBending, eqCurv);
            }

        } else {
            for(const auto& v : mesh.getVertices()) {
                const auto kBending = v.vertex->getMVoronoiCell()->getBendingModulus();
                const auto eqCurv = v.vertex->getMVoronoiCell()->getEqCurv();

                // The calculation requires that the current stretched area and mean curvature have already been calculated
                // As a result, d is just a dummy variable due to its universality
                const auto sArea = v.attr.gVertex.sArea;
                const auto sCurv = v.attr.gVertex.sCurv;

                U_i += _FFType.energy(sArea, sCurv, kBending, eqCurv, d);
            }
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

template<>
void MembraneBending<MembraneBendingVoronoiHelfrich>::computeForces() {
    
    for (auto m: Membrane::getMembranes()) {
    
        const auto& mesh = m->getMesh();

        const size_t numVertices = mesh.getVertices().size();
        for(size_t vi = 0; vi < numVertices; ++vi) {
            const auto& v = mesh.getVertices()[vi];

            const auto kBending = v.vertex->getMVoronoiCell()->getBendingModulus();
            const auto eqCurv = v.vertex->getMVoronoiCell()->getEqCurv();

            const auto area = v.attr.gVertex.area;
            const auto curv = v.attr.gVertex.curv;
           
            _FFType.forces(
                v.attr.vertex,
                area, v.attr.gVertex.dArea,
                curv, v.attr.gVertex.dCurv,
                kBending, eqCurv
            );

            mesh.forEachHalfEdgeTargetingVertex(vi, [this, &mesh, area, curv, kBending, eqCurv](size_t hei) {
                const size_t hei_o = mesh.opposite(hei);
                auto vt = mesh.getVertexAttribute(mesh.target(hei_o)).vertex;
                const auto& dArea = mesh.getEdgeAttribute(hei_o).gHalfEdge.dNeighborArea;
                const auto& dCurv = mesh.getEdgeAttribute(hei_o).gHalfEdge.dNeighborCurv;
                _FFType.forces(vt, area, dArea, curv, dCurv, kBending, eqCurv);
            });
        }

    }
}

template<>
void MembraneBending<MembraneBendingVoronoiHelfrich>::computeForcesAux() {
    
    for (auto m: Membrane::getMembranes()) {
    
        const auto& mesh = m->getMesh();

        const size_t numVertices = mesh.getVertices().size();
        for(size_t vi = 0; vi < numVertices; ++vi) {
            const auto& v = mesh.getVertices()[vi];

            const auto kBending = v.vertex->getMVoronoiCell()->getBendingModulus();
            const auto eqCurv = v.vertex->getMVoronoiCell()->getEqCurv();

            const auto area = v.attr.gVertex.area;
            const auto curv = v.attr.gVertex.curv;
           
            _FFType.forcesAux(
                v.attr.vertex,
                area, v.attr.gVertex.dArea,
                curv, v.attr.gVertex.dCurv,
                kBending, eqCurv
            );

            mesh.forEachHalfEdgeTargetingVertex(vi, [this, &mesh, area, curv, kBending, eqCurv](size_t hei) {
                const size_t hei_o = mesh.opposite(hei);
                auto vt = mesh.getVertexAttribute(mesh.target(hei_o)).vertex;
                const auto& dArea = mesh.getEdgeAttribute(hei_o).gHalfEdge.dNeighborArea;
                const auto& dCurv = mesh.getEdgeAttribute(hei_o).gHalfEdge.dNeighborCurv;
                _FFType.forcesAux(vt, area, dArea, curv, dCurv, kBending, eqCurv);
            });
        }
    }
}
