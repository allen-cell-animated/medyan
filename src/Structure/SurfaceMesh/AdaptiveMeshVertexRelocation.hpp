#ifndef MEDYAN_Structure_SurfaceMesh_AdaptiveMeshVertexRelocation_hpp
#define MEDYAN_Structure_SurfaceMesh_AdaptiveMeshVertexRelocation_hpp

#include "MathFunctions.h"
#include "Structure/SurfaceMesh/AdaptiveMeshGeometryManager.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"

namespace adaptive_mesh {

// Different methods can be used in vertex relocation
//   - (Removed) Vertex relaxations are vertex movement via moving along assigned force fields
//   - Direct vertex relocations are to move the vertex to the optimal location


// Direct vertex relocation as method of vertex relocation
enum class OptimalVertexLocationMethod {
    barycenter // of points that forms equilateral triangles
};
template< OptimalVertexLocationMethod > struct OptimalVertexLocation;
template<> struct OptimalVertexLocation< OptimalVertexLocationMethod::barycenter > {
    // Given vertex v1, v2, and unit normal un,
    // Find the vertex v0 such that triangle (v0, v1, v2) is equilateral
    // and un is the normal direction of the triangle.
    template< typename Float >
    static mathfunc::Vec< 3, Float > findEquilateralTriangle(
        const mathfunc::Vec< 3, Float >& v1,
        const mathfunc::Vec< 3, Float >& v2,
        const mathfunc::Vec< 3, Float >& un
    ) {
        return (Float)0.5 * (v1 + v2) + mathfunc::cross(un, v2 - v1) * (Float)(std::sqrt(3.0) * 0.5);
    }

    template< typename Mesh >
    auto operator()(const Mesh& mesh, typename Mesh::VertexIndex vi) const {
        using namespace mathfunc;
        using CoordinateType = typename Mesh::AttributeType::CoordinateType;

        CoordinateType target {};
        mesh.forEachHalfEdgeTargetingVertex(vi, [&](auto hei) {
            const CoordinateType cn (mesh.attribute(mesh.target(mesh.next(hei))).getCoordinate());
            const CoordinateType cp (mesh.attribute(mesh.target(mesh.prev(hei))).getCoordinate());
            const CoordinateType un (mesh.attribute(mesh.triangle(hei)).gTriangle.unitNormal);
            target += findEquilateralTriangle< typename CoordinateType::float_type >(cp, cn, un);
        });
        target *= (1.0 / mesh.degree(vi));

        // project onto tangent plane
        const CoordinateType ci (mesh.attribute(vi).getCoordinate());
        const auto& un = mesh.attribute(vi).aVertex.unitNormal;
        target -= un * dot(un, target - ci);

        return target;
    }
};

template<
    OptimalVertexLocationMethod opt
> class DirectVertexRelocationManager {
public:
    using MeshType = Membrane::MeshType;
    using OptimalVertexLocationType = OptimalVertexLocation< opt >;
    using GeometryManagerType = GeometryManager< MeshType >;

private:
    size_t maxIterRelocation_;
    size_t maxIter_; // each iteration: (relocation + edge flipping)

    // Returns number of edges flipped
    template< typename EdgeFlipManagerType >
    size_t edgeFlipping_(MeshType& mesh, const EdgeFlipManagerType& efm) const {
        // Edge flipping does not change edge id or total number of edges
        // Also the preferred length does not need to be changed
        size_t res = 0;
        const size_t numEdges = mesh.numEdges();
        for(MeshType::EdgeIndex i {0}; i < numEdges; ++i) {
            if(efm.tryFlip(mesh, i) == EdgeFlipManagerType::State::Success) ++res;
        }
        return res;
    }

    template< typename VT >
    void resetVertexCoordinate_(
        MeshType&             mesh,
        MeshType::VertexIndex vi,
        const VT&             newCoord
    ) const {
        // Record old properties
        double oldTotalEqArea = 0;
        if(mesh.metaAttribute().isMechParamsSet) {
            mesh.forEachHalfEdgeTargetingVertex(vi, [&](MeshType::HalfEdgeIndex hei) {
                if(mesh.isInTriangle(hei)) {
                    const auto ti = mesh.triangle(hei);
                    oldTotalEqArea += mesh.attribute(ti).triangle->mTriangle.eqArea;
                }
            });
        }

        // Set vertex new coordinate
        mesh.attribute(vi).getCoordinate() = newCoord;

        // Set new properties
        if(mesh.metaAttribute().isMechParamsSet) {
            double newTotalArea = 0.0;
            mesh.forEachHalfEdgeTargetingVertex(vi, [&](MeshType::HalfEdgeIndex hei) {
                if(mesh.isInTriangle(hei)) {
                    newTotalArea += medyan::area(mesh, mesh.triangle(hei));
                }
            });
            mesh.forEachHalfEdgeTargetingVertex(vi, [&](MeshType::HalfEdgeIndex hei) {
                if(mesh.isInTriangle(hei)) {
                    mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea
                        = oldTotalEqArea * medyan::area(mesh, mesh.triangle(hei)) / newTotalArea;
                }
            });
        }
    }

public:
    // Constructor
    DirectVertexRelocationManager(size_t maxIterRelocation, size_t maxIter)
        : maxIterRelocation_(maxIterRelocation), maxIter_(maxIter) {}

    template< typename EdgeFlipManagerType >
    void operator()(MeshType& mesh, const EdgeFlipManagerType& efm) const {
        using namespace mathfunc;

        const size_t numVertices = mesh.getVertices().size();

        // Main loop
        size_t flippingCount = 0;
        size_t iter = 0;
        do {
            ++iter;

            // Move vertices
            for(size_t iterRelo = 0; iterRelo < maxIterRelocation_; ++iterRelo) {
                for(MeshType::VertexIndex i {0}; i < numVertices; ++i) if(!mesh.isVertexOnBorder(i)) {
                    // const auto coordOriginal = mesh.attribute(i).getCoordinate();
                    const auto target = OptimalVertexLocationType{}(mesh, i);
                    /*const auto diff = target - coordOriginal;
                    const auto magDiff = magnitude(diff);*/
                    resetVertexCoordinate_(mesh, i, target);
                }
            }

            // Smooth out the folded geometry, without updating normals
            /*while(!membrane_mesh_check::MembraneMeshDihedralCheck{0.0, 0.0}(mesh, true)) {
                // Apply Laplacian smoothing
                _laplacianSmoothing(targets, mesh); 
            }*/

            // Update normals
            GeometryManagerType::computeAllTriangleNormals(mesh);
            GeometryManagerType::computeAllAngles(mesh);
            GeometryManagerType::computeAllVertexNormals(mesh);

            // Try flipping
            flippingCount = edgeFlipping_(mesh, efm);
        } while(flippingCount && iter < maxIter_);
    }
};

//-----------------------------------------------------------------------------
// Implements the vertex smoothing algorithm by flow toward smaller surface
// area.
//
// Reference:
//     Mark Meyer et al. (2003)
//     Discrete Differential-Geometry Operators for Triangulated 2-Manifolds
//     Page 18: "An Anisotropic Smoothing Technique"
//
// Since the surface mesh is mainly for membranes, which should not preserve
// any sharp features, we use iostropic mesh smoothing.
// Therefore, the velocity of a vertex i is
//     ∂x_i/∂t = − Δx_i
// where Δ is the (integrated) discretized Laplace operator (the cotangent
// formula).
//
// Note:
//   - This function helps denoising the initial mesh, but should not be used
//     during simulation, because this WILL change energetics.
//   - As a result, mechanical/chemical attributes will not be updated during
//     the process.

// Note:
//   - Only non-border vertices take part in curvature flow.
//   - dt is dimensionless
inline void meshSmoothingCurvatureFlowStep(
    Membrane::MeshType& mesh,
    double              dt
) {
    using namespace std;
    using namespace mathfunc;
    using MT = Membrane::MeshType;
    using GM = GeometryManager< MT >;
    using VV = vector< Vec< 3, floatingpoint > >;

    // Step 1.1. Calculate all angles
    GM::computeAllAngles(mesh);

    // Step 1.2. Calculate curvature and normal direction
    VV curv(mesh.numVertices());

    for(MT::VertexIndex vi {0}; vi < mesh.numVertices(); ++vi) {
        if(!mesh.isVertexOnBorder(vi)) {

            const auto& ci = mesh.attribute(vi).vertex->coord;

            Vec< 3, floatingpoint > d2x {};
            mesh.forEachHalfEdgeTargetingVertex(vi, [&](MT::HalfEdgeIndex hei) {
                const auto ti0    = mesh.triangle(hei);
                const auto hei_o  = mesh.opposite(hei);
                const auto hei_n  = mesh.next(hei);
                const auto hei_on = mesh.next(hei_o);
                const auto& cn      = mesh.attribute(mesh.target(hei_o )).vertex->coord;
                const auto& c_right = mesh.attribute(mesh.target(hei_on)).vertex->coord;

                const auto sumCotTheta =
                    mesh.attribute(hei_n).gHalfEdge.cotTheta
                    + mesh.attribute(hei_on).gHalfEdge.cotTheta;

                const auto diff = ci - cn;

                d2x += (0.5 * sumCotTheta) * (ci - cn);
            });

            curv[vi.index] = d2x;
        }
    }

    // Step 2. Move vertices along curvature vector (times dt)
    for(MT::VertexIndex vi {0}; vi < mesh.numVertices(); ++vi) {
        mesh.attribute(vi).vertex->coord -= dt * curv[vi.index];
    }
}

inline void meshSmoothing(
    Membrane::MeshType& mesh,
    double              dt,
    unsigned            numIter
) {
    for(unsigned iter = 0; iter < numIter; ++iter) {
        meshSmoothingCurvatureFlowStep(mesh, dt);
    }
}


} // namespace adaptive_mesh

#endif
