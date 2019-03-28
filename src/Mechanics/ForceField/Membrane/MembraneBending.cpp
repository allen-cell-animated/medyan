#include "MembraneBending.h"

#include "Membrane.hpp"
#include "Vertex.h"
#include "MVoronoiCell.h"

#include "MembraneBendingVoronoiHelfrich.h"

// Using the Helfrich Hamiltonian of mean curvature in Voronoi cells
template<>
double MembraneBending<MembraneBendingVoronoiHelfrich>::computeEnergy(bool stretched) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {
        const auto& mesh = m->getMesh();

        U_i = 0;

        for(const auto& v : mesh.getVertices()) {
            const auto kBending = v.attr.vertex->getMVoronoiCell()->getBendingModulus();
            const auto eqCurv = v.attr.vertex->getMVoronoiCell()->getEqCurv();

            const auto area = stretched ? v.attr.gVertexS.area : v.attr.gVertex.area;
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

template<>
void MembraneBending<MembraneBendingVoronoiHelfrich>::computeForces() {
    
    for (auto m: Membrane::getMembranes()) {
    
        const auto& mesh = m->getMesh();

        const size_t numVertices = mesh.getVertices().size();
        for(size_t vi = 0; vi < numVertices; ++vi) {
            const auto& v = mesh.getVertices()[vi];

            const auto kBending = v.attr.vertex->getMVoronoiCell()->getBendingModulus();
            const auto eqCurv = v.attr.vertex->getMVoronoiCell()->getEqCurv();

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
                const auto& dArea = mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dNeighborArea;
                const auto& dCurv = mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dNeighborCurv;
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

            const auto kBending = v.attr.vertex->getMVoronoiCell()->getBendingModulus();
            const auto eqCurv = v.attr.vertex->getMVoronoiCell()->getEqCurv();

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
                const auto& dArea = mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dNeighborArea;
                const auto& dCurv = mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dNeighborCurv;
                _FFType.forcesAux(vt, area, dArea, curv, dCurv, kBending, eqCurv);
            });
        }
    }
}
