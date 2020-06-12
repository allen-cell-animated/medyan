#ifndef MEDYAN_Structure_SurfaceMesh_AdaptiveMeshVertexRelocation_hpp
#define MEDYAN_Structure_SurfaceMesh_AdaptiveMeshVertexRelocation_hpp

#include "MathFunctions.h"
#include "Structure/SurfaceMesh/AdaptiveMeshGeometryManager.hpp"
#include "Structure/SurfaceMesh/MembraneMeshCheck.hpp"
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
        using CoordinateType = typename Mesh::AttributeType::coordinate_type;

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
    typename Mesh,
    OptimalVertexLocationMethod opt
> class DirectVertexRelocationManager {
public:
    using OptimalVertexLocationType = OptimalVertexLocation< opt >;
    using GeometryManagerType = GeometryManager< Mesh >;

private:
    size_t maxIterRelocation_;
    size_t maxIter_; // each iteration: (relocation + edge flipping)

    // Returns number of edges flipped
    template< typename EdgeFlipManagerType >
    size_t edgeFlipping_(Mesh& mesh, const EdgeFlipManagerType& efm) const {
        // Edge flipping does not change edge id or total number of edges
        // Also the preferred length does not need to be changed
        size_t res = 0;
        const size_t numEdges = mesh.numEdges();
        for(typename Mesh::EdgeIndex i {0}; i < numEdges; ++i) {
            if(efm.tryFlip(mesh, i) == EdgeFlipManagerType::State::Success) ++res;
        }
        return res;
    }

public:
    // Constructor
    DirectVertexRelocationManager(size_t maxIterRelocation, size_t maxIter)
        : maxIterRelocation_(maxIterRelocation), maxIter_(maxIter) {}

    template< typename EdgeFlipManagerType >
    void operator()(Mesh& mesh, const EdgeFlipManagerType& efm) const {
        using namespace mathfunc;

        const size_t numVertices = mesh.getVertices().size();

        // Main loop
        size_t flippingCount = 0;
        size_t iter = 0;
        do {
            ++iter;

            // Move vertices
            for(size_t iterRelo = 0; iterRelo < maxIterRelocation_; ++iterRelo) {
                for(typename Mesh::VertexIndex i {}; i < numVertices; ++i) if(!mesh.isVertexOnBorder(i)) {
                    // const auto coordOriginal = mesh.attribute(i).getCoordinate();
                    const auto target = OptimalVertexLocationType{}(mesh, i);
                    /*const auto diff = target - coordOriginal;
                    const auto magDiff = magnitude(diff);*/
                    mesh.attribute(i).getCoordinate() = target;
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

} // namespace adaptive_mesh

#endif
