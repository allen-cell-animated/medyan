#include "Mechanics/ForceField/Membrane/MembraneStretching.hpp"

#include <type_traits>

#include "Mechanics/ForceField/Membrane/MembraneStretchingImpl.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Triangle.h"
#include "Structure/SurfaceMesh/Vertex.h"

template< typename Impl, MembraneStretchingAccumulationType accuType >
double MembraneStretching< Impl, accuType >::computeEnergy(const double* coord, bool stretched) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {

        double area = 0.0;

        if constexpr (accuType == MembraneStretchingAccumulationType::ByVertex) {
            for(const auto& v : m->getMesh().getVertices()) if(v.numTargetingBorderHalfEdges == 0)
                area += (stretched ? v.attr.gVertexS.astar : v.attr.gVertex.astar) / 3;
        } else {
            for(const auto& t : m->getMesh().getTriangles())
                area += stretched ? t.attr.gTriangleS.area : t.attr.gTriangle.area;
        }

        if constexpr (std::is_same< Impl, MembraneStretchingHarmonic >::value) {
            const auto kElastic = m->getMMembrane()->getKElastic();
            const auto eqArea = m->getMMembrane()->getEqArea();

            U_i = _impl.energy(area, kElastic, eqArea);

        } else {
            const auto tension = m->getMMembrane()->getTension();

            U_i = _impl.energy(area, tension);
        }

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

template< typename Impl, MembraneStretchingAccumulationType accuType >
void MembraneStretching< Impl, accuType >::computeForces(const double* coord, double* force) {
    
    for (auto m: Membrane::getMembranes()) {

        Membrane::MembraneMeshAttributeType::cacheIndices(m->getMesh());

        const auto& mesh = m->getMesh();

        if constexpr (std::is_same< Impl, MembraneStretchingHarmonic >::value) {

            const auto kElastic = m->getMMembrane()->getKElastic();
            const auto eqArea = m->getMMembrane()->getEqArea();

            double area = 0.0;
            if constexpr (accuType == MembraneStretchingAccumulationType::ByVertex) {
                for(const auto& v : mesh.getVertices()) area += v.attr.gVertex.astar / 3;

                const auto& cvt = mesh.getMetaAttribute().cachedVertexTopo;
                const size_t numVertices = mesh.getVertices().size();
                for(size_t vi = 0; vi < numVertices; ++vi) if(!mesh.isVertexOnBorder(vi)) {
                    const auto& va = mesh.getVertexAttribute(vi);
                
                    _impl.forces(
                        force + 3 * va.cachedCoordIndex,
                        area, va.gVertex.dAstar / 3, kElastic, eqArea
                    );

                    // Position of this vertex also affects neighbor vertex areas
                    for(size_t i = 0; i < va.cachedDegree; ++i) {
                        const size_t hei = cvt[mesh.getMetaAttribute().cachedVertexOffsetTargetingHE(vi) + i];
                        const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dNeighborAstar / 3;
                        _impl.forces(
                            force + 3 * va.cachedCoordIndex,
                            area, dArea, kElastic, eqArea
                        );
                    }
                }
            } else { // By triangle
                for(const auto& t : mesh.getTriangles()) area += t.attr.gTriangle.area;

                const size_t numTriangles = mesh.getTriangles().size();
                for(size_t ti = 0; ti < numTriangles; ++ti) {
                    const auto& ta = mesh.getTriangleAttribute(ti);

                    for(size_t i = 0; i < 3; ++i) {
                        const size_t hei = ta.cachedHalfEdgeIndex[i];
                        const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dTriangleArea;
                        _impl.forces(
                            force + 3 * ta.cachedCoordIndex[i],
                            area, dArea, kElastic, eqArea
                        );
                    }
                }
            } // end if (accuType ...)

        } else { // Use linear energy
            const auto tension = m->getMMembrane()->getTension();

            if constexpr (accuType == MembraneStretchingAccumulationType::ByVertex) {

                const auto& cvt = mesh.getMetaAttribute().cachedVertexTopo;
                const size_t numVertices = mesh.getVertices().size();
                for(size_t vi = 0; vi < numVertices; ++vi) if(!mesh.isVertexOnBorder(vi)) {
                    const auto& va = mesh.getVertexAttribute(vi);
                
                    _impl.forces(
                        force + 3 * va.cachedCoordIndex,
                        va.gVertex.dAstar / 3, tension
                    );

                    // Position of this vertex also affects neighbor vertex areas
                    for(size_t i = 0; i < va.cachedDegree; ++i) {
                        const size_t hei = cvt[mesh.getMetaAttribute().cachedVertexOffsetTargetingHE(vi) + i];
                        const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dNeighborAstar / 3;
                        _impl.forces(
                            force + 3 * va.cachedCoordIndex,
                            dArea, tension
                        );
                    }
                }
            } else { // By triangle

                const size_t numTriangles = mesh.getTriangles().size();
                for(size_t ti = 0; ti < numTriangles; ++ti) {
                    const auto& ta = mesh.getTriangleAttribute(ti);

                    for(size_t i = 0; i < 3; ++i) {
                        const size_t hei = ta.cachedHalfEdgeIndex[i];
                        const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dTriangleArea;
                        _impl.forces(
                            force + 3 * ta.cachedCoordIndex[i],
                            dArea, tension
                        );
                    }
                }
            } // end if (accuType ...)

        } // end if (Impl ...)
    }
}


// Explicit instantiations
template class MembraneStretching< MembraneStretchingHarmonic, MembraneStretchingAccumulationType::ByTriangle >;
template class MembraneStretching< MembraneStretchingHarmonic, MembraneStretchingAccumulationType::ByVertex   >;
template class MembraneStretching< MembraneStretchingLinear,   MembraneStretchingAccumulationType::ByTriangle >;
template class MembraneStretching< MembraneStretchingLinear,   MembraneStretchingAccumulationType::ByVertex   >;
