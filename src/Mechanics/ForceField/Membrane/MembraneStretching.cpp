#include "Mechanics/ForceField/Membrane/MembraneStretching.hpp"

#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Triangle.h"
#include "Structure/SurfaceMesh/Vertex.h"

// Using the area of the Voronoi cells
template<>
double MembraneStretching<MembraneStretchingAccumulationType::ByVertex>::computeEnergy(const double* coord, bool stretched) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {

        const auto kElastic = m->getMMembrane()->getKElastic();
        const auto eqArea = m->getMMembrane()->getEqArea();

        double area = 0.0;
        for(const auto& v : m->getMesh().getVertices()) if(v.numTargetingBorderHalfEdges == 0)
            area += (stretched ? v.attr.gVertexS.astar : v.attr.gVertex.astar) / 3;

        U_i = _msh.energy(area, kElastic, eqArea); 

        if(fabs(U_i) == numeric_limits<double>::infinity()
            || U_i != U_i || U_i < -1.0) {
            _membraneCulprit = m;
            return -1;
        } else
            U += U_i;
        
    }

    return U;
}

// Force on each vertex is calculated one-time using derivative of vcell area and
// the derivative of neighbor vcell areas on the center vertex.
template<>
void MembraneStretching<MembraneStretchingAccumulationType::ByVertex>::computeForces(const double* coord, double* force) {
    
    for (auto m: Membrane::getMembranes()) {

        Membrane::MembraneMeshAttributeType::cacheIndices(m->getMesh());

        const auto& mesh = m->getMesh();
        const auto& cvt = mesh.getMetaAttribute().cachedVertexTopo;

        const auto kElastic = m->getMMembrane()->getKElastic();
        const auto eqArea = m->getMMembrane()->getEqArea();

        double area = 0.0;
        for(const auto& v : mesh.getVertices()) area += v.attr.gVertex.astar / 3;

        const size_t numVertices = mesh.getVertices().size();
        for(size_t vi = 0; vi < numVertices; ++vi) if(!mesh.isVertexOnBorder(vi)) {
            const auto& va = mesh.getVertexAttribute(vi);
           
            _msh.forces(
                force + 3 * va.cachedCoordIndex,
                area, va.gVertex.dAstar / 3, kElastic, eqArea
            );

            // Position of this vertex also affects neighbor vertex areas
            for(size_t i = 0; i < va.cachedDegree; ++i) {
                const size_t hei = cvt[mesh.getMetaAttribute().cachedVertexOffsetTargetingHE(vi) + i];
                const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dNeighborAstar / 3;
                _msh.forces(
                    force + 3 * va.cachedCoordIndex,
                    area, dArea, kElastic, eqArea
                );
            }
        }
    }
}


// Using the areas of the triangles
template<>
double MembraneStretching<MembraneStretchingAccumulationType::ByTriangle>::computeEnergy(const double* coord, bool stretched) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {

        const auto kElastic = m->getMMembrane()->getKElastic();
        const auto eqArea = m->getMMembrane()->getEqArea();

        double area = 0.0;
        for(const auto& t : m->getMesh().getTriangles())
            area += stretched ? t.attr.gTriangleS.area : t.attr.gTriangle.area;

        U_i = _msh.energy(area, kElastic, eqArea); 

        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) { // (U_i != U_i) => (U_i == NaN)
            
            //set culprit and return
            _membraneCulprit = m;
            
            return -1;
        }
        else
            U += U_i;

    }

    return U;
}

// Currently force calculation using triangles are different with the one using vcells.
// Using triangles, looping through triangles and forces are accumulated on the vertices.
template<>
void MembraneStretching<MembraneStretchingAccumulationType::ByTriangle>::computeForces(const double* coord, double* force) {
    
    for (auto m: Membrane::getMembranes()) {

        Membrane::MembraneMeshAttributeType::cacheIndices(m->getMesh());

        const auto& mesh = m->getMesh();

        const auto kElastic = m->getMMembrane()->getKElastic();
        const auto eqArea = m->getMMembrane()->getEqArea();

        double area = 0.0;
        for(const auto& t : mesh.getTriangles()) area += t.attr.gTriangle.area;

        const size_t numTriangles = mesh.getTriangles().size();
        for(size_t ti = 0; ti < numTriangles; ++ti) {
            const auto& ta = mesh.getTriangleAttribute(ti);

            for(size_t i = 0; i < 3; ++i) {
                const size_t hei = ta.cachedHalfEdgeIndex[i];
                const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dTriangleArea;
                _msh.forces(
                    force + 3 * ta.cachedCoordIndex[i],
                    area, dArea, kElastic, eqArea
                );
            }
        }
    }
}
