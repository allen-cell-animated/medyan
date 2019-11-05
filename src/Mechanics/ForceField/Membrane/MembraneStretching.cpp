#include "Mechanics/ForceField/Membrane/MembraneStretching.hpp"

#include <type_traits>

#include "Mechanics/ForceField/Membrane/MembraneStretchingImpl.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Triangle.hpp"
#include "Structure/SurfaceMesh/Vertex.hpp"

#ifdef __cpp_if_constexpr
    #define IF_CONSTEXPR if constexpr
#else
    #define IF_CONSTEXPR if

namespace {

template< typename Impl, std::enable_if_t< std::is_same< Impl, MembraneStretchingHarmonic >::value >* = nullptr >
double implEnergy(const Impl& impl, double area, double param0, double param1) {
    return impl.energy(area, param0, param1);
}
template< typename Impl, std::enable_if_t< std::is_same< Impl, MembraneStretchingLinear >::value >* = nullptr >
double implEnergy(const Impl& impl, double area, double param0, double param1) {
    return impl.energy(area, param0);
}

template< typename Impl, std::enable_if_t< std::is_same< Impl, MembraneStretchingHarmonic >::value >* = nullptr >
void implForces(const Impl& impl, floatingpoint* force, double area, const mathfunc::Vec3& dArea, double param0, double param1) {
    impl.forces(force, area, dArea, param0, param1);
}
template< typename Impl, std::enable_if_t< std::is_same< Impl, MembraneStretchingLinear >::value >* = nullptr >
void implForces(const Impl& impl, floatingpoint* force, double area, const mathfunc::Vec3& dArea, double param0, double param1) {
    impl.forces(force, dArea, param0);
}

} // namespace

#endif // ifdef __cpp_if_constexpr

template< typename Impl, MembraneStretchingAccumulationType accuType >
floatingpoint MembraneStretching< Impl, accuType >::computeEnergy(const floatingpoint* coord, bool stretched) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {

        double area = 0.0;

        IF_CONSTEXPR (accuType == MembraneStretchingAccumulationType::ByVertex) {
            for(const auto& v : m->getMesh().getVertices()) if(v.numTargetingBorderHalfEdges == 0)
                area += (stretched ? v.attr.gVertexS.astar : v.attr.gVertex.astar) / 3;
        } else {
            for(const auto& t : m->getMesh().getTriangles())
                area += stretched ? t.attr.gTriangleS.area : t.attr.gTriangle.area;
        }

        IF_CONSTEXPR (std::is_same< Impl, MembraneStretchingHarmonic >::value) {
            const auto kElastic = m->getMMembrane()->getKElastic();
            const auto eqArea = m->getMMembrane()->getEqArea();

#ifdef __cpp_if_constexpr
            U_i = _impl.energy(area, kElastic, eqArea);
#else
            U_i = implEnergy(_impl, area, kElastic, eqArea);
#endif

        } else {
            const auto tension = m->getMMembrane()->getTension();

#ifdef __cpp_if_constexpr
            U_i = _impl.energy(area, tension);
#else
            U_i = implEnergy(_impl, area, tension, 0.0);
#endif
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
void MembraneStretching< Impl, accuType >::computeForces(const floatingpoint* coord, floatingpoint* force) {

    // Configure force buffer
    constexpr bool useForceBuffer = true;
    const std::size_t dof = Bead::getDbDataConst().coords.size_raw();

    if (useForceBuffer) {
        forceBuffer_.assign(dof, 0.0);
    }
    floatingpoint* const f = useForceBuffer ? forceBuffer_.data() : force;

    for (auto m: Membrane::getMembranes()) {

        Membrane::MembraneMeshAttributeType::cacheIndices(m->getMesh());

        const auto& mesh = m->getMesh();

        IF_CONSTEXPR (std::is_same< Impl, MembraneStretchingHarmonic >::value) {

            const auto kElastic = m->getMMembrane()->getKElastic();
            const auto eqArea = m->getMMembrane()->getEqArea();

            double area = 0.0;
            IF_CONSTEXPR (accuType == MembraneStretchingAccumulationType::ByVertex) {
                for(const auto& v : mesh.getVertices()) area += v.attr.gVertex.astar / 3;

                const auto& cvt = mesh.getMetaAttribute().cachedVertexTopo;
                const size_t numVertices = mesh.getVertices().size();
                for(size_t vi = 0; vi < numVertices; ++vi) if(!mesh.isVertexOnBorder(vi)) {
                    const auto& va = mesh.getVertexAttribute(vi);

#ifdef __cpp_if_constexpr
                    _impl.forces(
                        f + 3 * va.cachedCoordIndex,
                        area, va.gVertex.dAstar / 3, kElastic, eqArea
                    );
#else
                    implForces(_impl,
                        f + 3 * va.cachedCoordIndex,
                        area, va.gVertex.dAstar / 3, kElastic, eqArea
                    );
#endif

                    // Position of this vertex also affects neighbor vertex areas
                    for(size_t i = 0; i < va.cachedDegree; ++i) {
                        const size_t hei = cvt[mesh.getMetaAttribute().cachedVertexOffsetTargetingHE(vi) + i];
                        const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dNeighborAstar / 3;
#ifdef __cpp_if_constexpr
                        _impl.forces(
                            f + 3 * va.cachedCoordIndex,
                            area, dArea, kElastic, eqArea
                        );
#else
                        implForces(_impl,
                            f + 3 * va.cachedCoordIndex,
                            area, dArea, kElastic, eqArea
                        );
#endif
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
#ifdef __cpp_if_constexpr
                        _impl.forces(
                            f + 3 * ta.cachedCoordIndex[i],
                            area, dArea, kElastic, eqArea
                        );
#else
                        implForces(_impl,
                            f + 3 * ta.cachedCoordIndex[i],
                            area, dArea, kElastic, eqArea
                        );
#endif
                    }
                }
            } // end if (accuType ...)

        } else { // Use linear energy
            const auto tension = m->getMMembrane()->getTension();

            IF_CONSTEXPR (accuType == MembraneStretchingAccumulationType::ByVertex) {

                const auto& cvt = mesh.getMetaAttribute().cachedVertexTopo;
                const size_t numVertices = mesh.getVertices().size();
                for(size_t vi = 0; vi < numVertices; ++vi) if(!mesh.isVertexOnBorder(vi)) {
                    const auto& va = mesh.getVertexAttribute(vi);
                
#ifdef __cpp_if_constexpr
                    _impl.forces(
                        f + 3 * va.cachedCoordIndex,
                        va.gVertex.dAstar / 3, tension
                    );
#else
                    implForces(_impl,
                        f + 3 * va.cachedCoordIndex,
                        0.0, va.gVertex.dAstar / 3, tension, 0.0
                    );
#endif

                    // Position of this vertex also affects neighbor vertex areas
                    for(size_t i = 0; i < va.cachedDegree; ++i) {
                        const size_t hei = cvt[mesh.getMetaAttribute().cachedVertexOffsetTargetingHE(vi) + i];
                        const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dNeighborAstar / 3;
#ifdef __cpp_if_constexpr
                        _impl.forces(
                            f + 3 * va.cachedCoordIndex,
                            dArea, tension
                        );
#else
                        implForces(_impl,
                            f + 3 * va.cachedCoordIndex,
                            0.0, dArea, tension, 0.0
                        );
#endif
                    }
                }
            } else { // By triangle

                const size_t numTriangles = mesh.getTriangles().size();
                for(size_t ti = 0; ti < numTriangles; ++ti) {
                    const auto& ta = mesh.getTriangleAttribute(ti);

                    for(size_t i = 0; i < 3; ++i) {
                        const size_t hei = ta.cachedHalfEdgeIndex[i];
                        const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dTriangleArea;
#ifdef __cpp_if_constexpr
                        _impl.forces(
                            f + 3 * ta.cachedCoordIndex[i],
                            dArea, tension
                        );
#else
                        implForces(_impl,
                            f + 3 * ta.cachedCoordIndex[i],
                            0.0, dArea, tension, 0.0
                        );
#endif
                    }
                }
            } // end if (accuType ...)

        } // end if (Impl ...)
    }

    if(useForceBuffer) {
        std::transform(
            force, force + dof, forceBuffer_.begin(),
            force,
            std::plus<>{}
        );
    }

}


// Explicit instantiations
template class MembraneStretching< MembraneStretchingHarmonic, MembraneStretchingAccumulationType::ByTriangle >;
template class MembraneStretching< MembraneStretchingHarmonic, MembraneStretchingAccumulationType::ByVertex   >;
template class MembraneStretching< MembraneStretchingLinear,   MembraneStretchingAccumulationType::ByTriangle >;
template class MembraneStretching< MembraneStretchingLinear,   MembraneStretchingAccumulationType::ByVertex   >;
