#ifndef MEDYAN_Structure_SurfaceMesh_AdaptiveMeshVertexRelocation_hpp
#define MEDYAN_Structure_SurfaceMesh_AdaptiveMeshVertexRelocation_hpp

#include "MathFunctions.h"
#include "Structure/SurfaceMesh/AdaptiveMeshGeometryManager.hpp"
#include "Structure/SurfaceMesh/MembraneMeshCheck.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"

namespace adaptive_mesh {

// Different methods can be used in vertex relocation
//   - (Deprecated) Vertex relaxations are vertex movement via moving along assigned force fields
//   - Direct vertex relocations are to move the vertex to the optimal location

// Vertex relaxation as method of vertex relocation
enum class VertexRelaxationType {
    GlobalElastic // E = (l - l_0)^2 / (2 l_0), with const elastic modulus 1.0
};
template< VertexRelaxationType > struct RelaxationForceField;
template<> struct RelaxationForceField< VertexRelaxationType::GlobalElastic > {

    // The function requires the vertex unit normal information
    template< typename Mesh, typename VecType >
    void computeForces(std::vector<VecType>& forces, const Mesh& mesh, const std::vector<VecType>& coords) {
        // The size of forces should be the same as the size of coords.
        // Resizing of forces should be done by the caller.
        const size_t numVertices = coords.size();
        for(size_t i = 0; i < numVertices; ++i) {
            VecType f {};
            mesh.forEachHalfEdgeTargetingVertex(i, [&](size_t hei) {
                const auto l0 = mesh.getEdgeAttribute(mesh.edge(hei)).aEdge.eqLength;
                const auto r = coords[mesh.target(mesh.opposite(hei))] - coords[i];
                const auto mag = mathfunc::magnitude(r);
                f += r * (1.0 / l0 - 1.0 / mag);
            });

            const auto& un = mesh.getVertexAttribute(i).aVertex.unitNormal;
            f -= un * mathfunc::dot(un, f); // Remove normal component

            forces[i] = f;
        }
    }
};

template<
    typename Mesh,
    VertexRelaxationType r
> class [[deprecated]] GlobalRelaxationManager {
public:
    using RelaxationForceFieldType = RelaxationForceField< r >;
    using GeometryManagerType = GeometryManager< Mesh >;

private:
    double _epsilon2; // Square of force tolerance
    double _dt; // Step size for Runge Kutta method
    size_t _maxIterRelocation;
    size_t _maxIterRelaxation; // each iteration: (relocation + edge flipping)

    // Utility for max force magnitude squared
    template< typename VecType > static auto _maxMag2(std::vector<VecType>& v) {
        double res = 0.0;
        for(const auto& i : v) {
            double mag2 = mathfunc::magnitude2(i);
            if(mag2 > res) res = mag2;
        }
        return res;
    }

    // Relocates vertices using 2nd order Runge Kutta method of overdamped dynamics
    // Returns whether the final forces are below threshold.
    template< typename VecType > bool _vertexRelocation(
        std::vector<VecType>& coords,
        std::vector<VecType>& forces,
        std::vector<VecType>& coordsHalfway,
        std::vector<VecType>& forcesHalfway,
        const Mesh& mesh
    ) const {
        const size_t numVertices = coords.size();

        RelaxationForceFieldType().computeForces(forces, mesh, coords);
        auto maxMag2F = _maxMag2(forces);

        size_t iter = 0;
        while(maxMag2F >= _epsilon2 && iter < _maxIterRelocation) {
            ++iter;

            // Test move halfway
            for(size_t i = 0; i < numVertices; ++i)
                coordsHalfway[i] = coords[i] + (0.5 * _dt) * forces[i];

            // Force at halfway
            RelaxationForceFieldType().computeForces(forcesHalfway, mesh, coordsHalfway);

            // Real move
            for(size_t i = 0; i < numVertices; ++i)
                coords[i] += forcesHalfway[i] * _dt;

            // Compute new forces
            RelaxationForceFieldType().computeForces(forces, mesh, coords);
            maxMag2F = _maxMag2(forces);
        }

        if(maxMag2F >= _epsilon2) return false;
        else return true;
    }

    template< typename VecType > void _laplacianSmoothing(
        std::vector< VecType >& targets,
        Mesh& mesh
    ) const {
        using namespace mathfunc;

        const size_t numVertices = mesh.getVertices().size();

        for(size_t i = 0; i < numVertices; ++i) {
            targets[i] = Vec3 {};
            double weightTot = 0.0;
            const Vec3 ci = mesh.getVertexAttribute(i).getCoordinate();
            mesh.forEachHalfEdgeTargetingVertex(i, [&](size_t hei) {
                const size_t vn = mesh.target(mesh.next(hei));
                const size_t vp = mesh.target(mesh.prev(hei));
                const Vec3 cn = mesh.getVertexAttribute(vn).getCoordinate();
                const Vec3 cp = mesh.getVertexAttribute(vp).getCoordinate();
                const auto centroid = (ci + cn + cp) / 3.0;
                const auto len = distance(ci, centroid);
                targets[i] += len * centroid;
                weightTot += len;
            });
            targets[i] /= weightTot;

            // Project to tangent plane
            const auto& un = mesh.getVertexAttribute(i).aVertex.unitNormal;
            targets[i] -= un * dot(targets[i] - ci, un);
        }

        for(size_t i = 0; i < numVertices; ++i) {
            mesh.getVertexAttribute(i).getCoordinate() = targets[i];
        }
    } // End function _laplacianSmoothing(Mesh&)

    // Returns number of edges flipped
    template< typename EdgeFlipManagerType >
    size_t _edgeFlipping(Mesh& mesh, const EdgeFlipManagerType& efm) const {
        // Edge flipping does not change edge id or total number of edges
        // Also the preferred length does not need to be changed
        size_t res = 0;
        const size_t numEdges = mesh.getEdges().size();
        for(size_t i = 0; i < numEdges; ++i) {
            if(efm.tryFlip(mesh, i) == EdgeFlipManagerType::State::Success) ++res;
        }
        return res;
    }

};

// Direct vertex relocation as method of vertex relocation
enum class OptimalVertexLocationMethod {
    Barycenter // of points that forms equilateral triangles
};
template< OptimalVertexLocationMethod > struct OptimalVertexLocation;
template<> struct OptimalVertexLocation< OptimalVertexLocationMethod::Barycenter > {
    // Given vertex v1, v2, and unit normal un,
    // Find the vertex v0 such that triangle (v0, v1, v2) is equilateral
    // and un is the normal direction of the triangle.
    static mathfunc::Vec3 findEquilateralTriangle(
        const mathfunc::Vec3& v1,
        const mathfunc::Vec3& v2,
        const mathfunc::Vec3& un
    ) {
        return 0.5 * (v1 + v2) + mathfunc::cross(un, v2 - v1) * (std::sqrt(3.0) * 0.5);
    }

    template< typename Mesh >
    mathfunc::Vec3 operator()(const Mesh& mesh, size_t vi) const {
        using namespace mathfunc;

        Vec3 target {};
        mesh.forEachHalfEdgeTargetingVertex(vi, [&](size_t hei) {
            const Vec3 cn = mesh.getVertexAttribute(mesh.target(mesh.next(hei))).getCoordinate();
            const Vec3 cp = mesh.getVertexAttribute(mesh.target(mesh.prev(hei))).getCoordinate();
            const auto& un = mesh.getTriangleAttribute(mesh.triangle(hei)).gTriangle.unitNormal;
            target += findEquilateralTriangle(cp, cn, un);
        });
        target *= (1.0 / mesh.degree(vi));

        // project onto tangent plane
        const Vec3 ci = mesh.getVertexAttribute(vi).getCoordinate();
        const auto& un = mesh.getVertexAttribute(vi).aVertex.unitNormal;
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
    size_t _maxIterRelocation;
    size_t _maxIter; // each iteration: (relocation + edge flipping)

    // Returns number of edges flipped
    template< typename EdgeFlipManagerType >
    size_t _edgeFlipping(Mesh& mesh, const EdgeFlipManagerType& efm) const {
        // Edge flipping does not change edge id or total number of edges
        // Also the preferred length does not need to be changed
        size_t res = 0;
        const size_t numEdges = mesh.getEdges().size();
        for(size_t i = 0; i < numEdges; ++i) {
            if(efm.tryFlip(mesh, i) == EdgeFlipManagerType::State::Success) ++res;
        }
        return res;
    }

public:
    // Constructor
    DirectVertexRelocationManager(size_t maxIterRelocation, size_t maxIter)
        : _maxIterRelocation(maxIterRelocation), _maxIter(maxIter) {}

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
            for(size_t iterRelo = 0; iterRelo < _maxIterRelocation; ++iterRelo) {
                for(size_t i = 0; i < numVertices; ++i) if(!mesh.isVertexOnBorder(i)) {
                    // const auto coordOriginal = mesh.getVertexAttribute(i).getCoordinate();
                    const auto target = OptimalVertexLocationType{}(mesh, i);
                    /*const auto diff = target - coordOriginal;
                    const auto magDiff = magnitude(diff);*/
                    mesh.getVertexAttribute(i).getCoordinate() = target;
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
            flippingCount = _edgeFlipping(mesh, efm);
        } while(flippingCount && iter < _maxIter);
    }
};

} // namespace adaptive_mesh

#endif
