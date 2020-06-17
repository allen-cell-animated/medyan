#ifndef MEDYAN_Structure_SurfaceMesh_MembraneMeshGeometry_hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneMeshGeometry_hpp

#include <tuple>

#include "MathFunctions.h"
#include "Mechanics/ForceField/Types.hpp"
#include "Structure/SurfaceMesh/MembraneMeshAttribute.hpp"
#include "Structure/SurfaceMesh/SurfaceMesh.hpp"

namespace medyan {

enum class SurfaceCurvaturePolicy {
    // Computes the signed curvature.
    //
    // This is generally needed when curvature is used in linear cases,
    // such as the bending energy when spontaneous curvature is non-zero.
    //
    // The default implementation is
    //         ∇ Area ⋅ ∇ Vol
    //   H = ------------------
    //        2 ∇ Vol ⋅ ∇ Vol
    //
    // where ∇ acts on the vertex of interest.
    //
    // The signed curvature will be stored in the curv variable, and their
    // derivatives will be stored in curv related places.
    withSign,

    // Computes the squared curvature.
    //
    // This is generally used when curvature is used in quadratic cases,
    // such as the bending energy when spontaneous curvature is zero.
    //
    // The default implementation is
    //           (∇ Area)^2
    //   H^2 = -------------
    //          4 (∇ Vol)^2
    //
    // The square of the curvature will be stored in the curv2 variable, and
    // their derivatives will be stored related to curv2.
    squared
};


//-------------------------------------------------------------------------
// Index caching
//-------------------------------------------------------------------------

// Mesh index caching
//
// The purpose of this function is to reduce pointer chasing while
// traversing the mesh structure.
// The cache is only valid thru force field calculation without mesh
// connection changes.
inline void cacheIndicesForFF(
    MembraneMeshAttribute::MeshType& mesh,
    const FFCoordinateStartingIndex& si
) {
    using MT = MembraneMeshAttribute::MeshType;

    const auto& vertices = mesh.getVertices();
    const auto& halfEdges = mesh.getHalfEdges();
    const auto& edges = mesh.getEdges();
    const auto& triangles = mesh.getTriangles();

    const size_t numVertices = vertices.size();
    const size_t numHalfEdges = halfEdges.size();
    const size_t numEdges = edges.size();
    const size_t numTriangles = triangles.size();

    // Note that the coord index means the index of x-coordinate of the vertex
    // in the vectorized coordinate array.
    const auto vci = [&](MT::VertexIndex vi) {
        return mesh.attribute(vi).vertex->getIndex() * 3 + si.vertex;
    };

    for(MT::HalfEdgeIndex hei {0}; hei < numHalfEdges; ++hei) {
        // The angle is (v0, v1, v2)
        const auto vi0 = mesh.target(mesh.prev(hei));
        const auto vi1 = mesh.target(hei);
        const auto vi2 = mesh.target(mesh.next(hei));

        auto& hea = mesh.attribute(hei);
        hea.cachedCoordIndex[0] = vci(vi0);
        hea.cachedCoordIndex[1] = vci(vi1);
        hea.cachedCoordIndex[2] = vci(vi2);
    }

    for(MT::TriangleIndex ti {0}; ti < numTriangles; ++ti) {
        const auto hei = mesh.halfEdge(ti);
        const auto hei_n = mesh.next(hei);
        const auto hei_p = mesh.prev(hei);
        const auto vi0 = mesh.target(hei);
        const auto vi1 = mesh.target(hei_n);
        const auto vi2 = mesh.target(hei_p);

        auto& ta = mesh.attribute(ti);
        ta.cachedCoordIndex[0] = vci(vi0);
        ta.cachedCoordIndex[1] = vci(vi1);
        ta.cachedCoordIndex[2] = vci(vi2);
        ta.cachedHalfEdgeIndex[0] = hei;
        ta.cachedHalfEdgeIndex[1] = hei_n;
        ta.cachedHalfEdgeIndex[2] = hei_p;
        ta.cachedEdgeIndex[0] = mesh.edge(hei);
        ta.cachedEdgeIndex[1] = mesh.edge(hei_n);
        ta.cachedEdgeIndex[2] = mesh.edge(hei_p);
    }

    for(MT::EdgeIndex ei {0}; ei < numEdges; ++ei) {
        const auto hei = mesh.halfEdge(ei);
        const auto hei_o = mesh.opposite(hei);
        const auto vi0 = mesh.target(hei);
        const auto vi1 = mesh.target(mesh.prev(hei));

        auto& ea = mesh.attribute(ei);
        ea.cachedCoordIndex[0]   = vci(vi0);
        ea.cachedCoordIndex[1]   = vci(vi1);
        ea.cachedPolygonType[0]  = mesh.polygonType(hei);
        ea.cachedPolygonType[1]  = mesh.polygonType(hei_o);
        ea.cachedPolygonIndex[0] = mesh.polygon(hei);
        ea.cachedPolygonIndex[1] = mesh.polygon(hei_o);
    }

    // Determine vertex max number of neighbors
    mesh.metaAttribute().vertexMaxDegree = 0;
    for(MT::VertexIndex vi {0}; vi < numVertices; ++vi) {
        mesh.metaAttribute().vertexMaxDegree = std::max(
            mesh.metaAttribute().vertexMaxDegree,
            mesh.degree(vi)
        );
    }
    mesh.metaAttribute().cachedVertexTopo.resize(
        mesh.metaAttribute().cachedVertexTopoSize() * numVertices);
    // Cache indices around vertices
    for(MT::VertexIndex vi {0}; vi < numVertices; ++vi) {
        auto& va = mesh.attribute(vi);

        va.cachedDegree = mesh.degree(vi);
        va.cachedCoordIndex = vci(vi);

        size_t i = 0; // Neighbor index

        mesh.forEachHalfEdgeTargetingVertex(vi, [&](MT::HalfEdgeIndex hei) {
            const auto hei_o = mesh.opposite(hei);
            const auto vn = mesh.target(hei_o);

            mesh.metaAttribute().cachedVertexTopo[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi) + i]
                = vci(vn);
            mesh.metaAttribute().cachedVertexTopo[mesh.metaAttribute().cachedVertexOffsetTargetingHE  (vi) + i]
                = hei.index;
            mesh.metaAttribute().cachedVertexTopo[mesh.metaAttribute().cachedVertexOffsetLeavingHE    (vi) + i]
                = hei_o.index;
            mesh.metaAttribute().cachedVertexTopo[mesh.metaAttribute().cachedVertexOffsetOuterHE      (vi) + i]
                = mesh.prev(hei).index;
            mesh.metaAttribute().cachedVertexTopo[mesh.metaAttribute().cachedVertexOffsetPolygon      (vi) + i]
                = mesh.polygon(hei);

            ++i;
        });
    }

    mesh.metaAttribute().indexCacheForFFValid = true;
} // void cacheIndicesForFF(...)
inline void invalidateIndexCacheForFF(MembraneMeshAttribute::MeshType& mesh) {
    mesh.metaAttribute().indexCacheForFFValid = false;
}
inline void assertValidIndexCacheForFF(const MembraneMeshAttribute::MeshType& mesh) {
    if(!mesh.metaAttribute().indexCacheForFFValid) {
        LOG(ERROR) << "Mesh index cache valid assertion failed.";
        throw std::runtime_error("Mesh index cache for force field is invalid");
    }
}


//-------------------------------------------------------------------------
// Geometries
//-------------------------------------------------------------------------

// Get vertex indices
inline std::array< MembraneMeshAttribute::MeshType::VertexIndex, 3 >
vertexIndices(
    const MembraneMeshAttribute::MeshType&           mesh,
    const MembraneMeshAttribute::MeshType::Triangle& t
) {
    const auto hei0 = t.halfEdgeIndex;
    const auto hei1 = mesh.next(hei0);
    const auto hei2 = mesh.next(hei1);
    return { mesh.target(hei0), mesh.target(hei1), mesh.target(hei2) };
}

inline std::array< MembraneMeshAttribute::MeshType::VertexIndex, 3 >
vertexIndices(
    const MembraneMeshAttribute::MeshType&         mesh,
    MembraneMeshAttribute::MeshType::TriangleIndex ti
) {
    return vertexIndices(mesh, mesh.element(ti));
}

// Calculate area of triangle
template< typename VT >
inline double area(const VT& c0, const VT& c1, const VT& c2) {
    const auto cp = mathfunc::cross(c1 - c0, c2 - c0);
    return mathfunc::magnitude(cp) * 0.5;
}
inline double area(
    const MembraneMeshAttribute::MeshType&         mesh,
    MembraneMeshAttribute::MeshType::TriangleIndex ti
) {
    using namespace std;
    using namespace mathfunc;
    const auto vis = vertexIndices(mesh, ti);
    const auto& c0 = mesh.attribute(vis[0]).getCoordinate();
    const auto& c1 = mesh.attribute(vis[1]).getCoordinate();
    const auto& c2 = mesh.attribute(vis[2]).getCoordinate();

    return area(c0, c1, c2);
}
// Returns the area of triangle and its derivatives on all vertices
template< typename VT >
inline auto areaAndDerivative(const VT& c0, const VT& c1, const VT& c2) {
    const auto area = area(c0, c1, c2);
    const auto inv4A = 0.25 / area;
    const auto r01 = c1 - c0;
    const auto r02 = c2 - c0;
    const auto dot_01_02 = mathfunc::dot(r01, r02);
    const auto dot_01_01 = mathfunc::dot(r01, r01);
    const auto dot_02_02 = mathfunc::dot(r02, r02);
    return std::tuple {
        area,
        inv4A * ((dot_01_02 - dot_01_01) * r02 + (dot_01_02 - dot_02_02) * r01),
        inv4A * (dot_02_02 * r01 - dot_01_02 * r02),
        inv4A * (dot_01_01 * r02 - dot_01_02 * r01)
    };
}

// This function updates geometries necessary for computing membrane energy
// and signed distance.
// This function uses cached indexing to enhance performance, so a valid
// cache is needed in this function.
template< bool stretched >
inline void updateGeometryValue(
    MembraneMeshAttribute::MeshType& mesh,
    const floatingpoint*             coord,
    SurfaceCurvaturePolicy           curvPol = SurfaceCurvaturePolicy::withSign
) {
    using namespace mathfunc;
    using MT = MembraneMeshAttribute::MeshType;

    assertValidIndexCacheForFF(mesh);

    const auto& vertices = mesh.getVertices();
    const auto& halfEdges = mesh.getHalfEdges();
    const auto& triangles = mesh.getTriangles();

    const size_t numVertices = vertices.size();
    const size_t numHalfEdges = halfEdges.size();
    const size_t numTriangles = triangles.size();

    // Calculate angles stored in half edges
    for(MT::HalfEdgeIndex hei {}; hei < numHalfEdges; ++hei) {
        // The angle is (c0, c1, c2)
        auto& hea = mesh.attribute(hei);
        auto& heag = hea.template getGHalfEdge<stretched>();
        const auto c0 = makeRefVec<3>(coord + hea.cachedCoordIndex[0]);
        const auto c1 = makeRefVec<3>(coord + hea.cachedCoordIndex[1]);
        const auto c2 = makeRefVec<3>(coord + hea.cachedCoordIndex[2]);

        const auto cp = cross(c0 - c1, c2 - c1);
        const auto dp =   dot(c0 - c1, c2 - c1);
        heag.cotTheta = dp / magnitude(cp);
    }

    // Calculate triangle area and cone volume
    for(MT::TriangleIndex ti {}; ti < numTriangles; ++ti) {
        auto& ta = mesh.attribute(ti);
        auto& tag = ta.template getGTriangle<stretched>();
        const auto c0 = makeRefVec<3>(coord + ta.cachedCoordIndex[0]);
        const auto c1 = makeRefVec<3>(coord + ta.cachedCoordIndex[1]);
        const auto c2 = makeRefVec<3>(coord + ta.cachedCoordIndex[2]);

        const auto cp = cross(c1 - c0, c2 - c0);

        // area
        tag.area = magnitude(cp) * 0.5;

        // cone volume
        tag.coneVolume = dot(c0, cp) / 6;
    }

    const auto& cvt = mesh.metaAttribute().cachedVertexTopo;

    // Calculate vertex 1-ring area and local curvature
    for(MT::VertexIndex vi {}; vi < numVertices; ++vi) if(!mesh.isVertexOnBorder(vi)) {
        auto& va = mesh.attribute(vi);
        auto& vag = va.template getGVertex<stretched>();
        const auto ci = makeRefVec<3>(coord + va.cachedCoordIndex);

        // clearing
        vag.astar = 0.0;
        vag.dAstar = {0.0, 0.0, 0.0};
        vag.dVolume = {0.0, 0.0, 0.0};

        // K = 2 * H * n is the result of LB operator.
        // K * A = d AStar (on center vertex), where the choice of A could be different.
        // Here the we use Vector Area for the A above, i.e. AVec = |d Vol (on center vertex)|
        //
        // We choose the following implementation (see star_perp_sq_mean_curvature of Surface Evolver):
        //
        //        (1/2) d AStar  dot  d Vol
        //   H = ----------------------------
        //               d Vol   dot  d Vol

        for(size_t i = 0; i < va.cachedDegree; ++i) {
            const MT::TriangleIndex ti0    { cvt[mesh.metaAttribute().cachedVertexOffsetPolygon(vi.index) + i] };
            const MT::HalfEdgeIndex hei_n  { cvt[mesh.metaAttribute().cachedVertexOffsetLeavingHE(vi.index) + (i + va.cachedDegree - 1) % va.cachedDegree] };
            const MT::HalfEdgeIndex hei_on { cvt[mesh.metaAttribute().cachedVertexOffsetOuterHE(vi.index) + (i + 1) % va.cachedDegree] };
            const auto cn      = makeRefVec<3>(coord + cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi.index) + i]);
            const auto c_right = makeRefVec<3>(coord + cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi.index) + (i + 1) % va.cachedDegree]);

            const auto sumCotTheta =
                mesh.attribute(hei_n).template getGHalfEdge<stretched>().cotTheta
                + mesh.attribute(hei_on).template getGHalfEdge<stretched>().cotTheta;

            const auto diff = ci - cn;

            vag.astar += mesh.attribute(ti0).template getGTriangle<stretched>().area;

            vag.dAstar += 0.5 * sumCotTheta * diff;

            // Added to derivative of sum of cone volume
            vag.dVolume += cross(cn, c_right) * (1.0 / 6);
        }

        const auto dVolume2 = magnitude2(vag.dVolume);

        if(curvPol == SurfaceCurvaturePolicy::withSign) {
            vag.curv = 0.5 * dot(vag.dAstar, vag.dVolume) / dVolume2;
        }
        else {
            vag.curv2 = 0.25 * magnitude2(vag.dAstar) / dVolume2;
        }
    }
} // void updateGeometryValue(...)

// This function updates the geometry value with derivatives necessary in
// the membrane force calculation.
// This function uses cached indexing to enhance performance, so a valid
// cache is needed in this function.
inline void updateGeometryValueWithDerivative(
    MembraneMeshAttribute::MeshType& mesh,
    const floatingpoint*             coord,
    SurfaceCurvaturePolicy           curvPol = SurfaceCurvaturePolicy::withSign
) {
    using namespace mathfunc;
    using MT = MembraneMeshAttribute::MeshType;
    using CT = MembraneMeshAttribute::CoordinateType;

    assertValidIndexCacheForFF(mesh);

    const auto& vertices = mesh.getVertices();
    const auto& triangles = mesh.getTriangles();

    const size_t numVertices = vertices.size();
    const size_t numTriangles = triangles.size();

    // Calculate angles and triangle areas with derivative
    // Calculate triangle cone volumes
    for(MT::TriangleIndex ti {}; ti < numTriangles; ++ti) {
        auto& ta = mesh.attribute(ti);

        const auto& hei = ta.cachedHalfEdgeIndex;
        auto& tag = mesh.attribute(ti).gTriangle;

        const CT c[] {
            static_cast<CT>(makeRefVec<3>(coord + ta.cachedCoordIndex[0])),
            static_cast<CT>(makeRefVec<3>(coord + ta.cachedCoordIndex[1])),
            static_cast<CT>(makeRefVec<3>(coord + ta.cachedCoordIndex[2]))
        };

        const double l2[] {
            distance2(c[2], c[0]),
            distance2(c[0], c[1]),
            distance2(c[1], c[2])
        };

        const double dots[] {
            dot(c[1] - c[0], c[2] - c[0]),
            dot(c[2] - c[1], c[0] - c[1]),
            dot(c[0] - c[2], c[1] - c[2])
        };

        const auto cp = cross(c[1] - c[0], c[2] - c[0]); // Pointing outward

        const auto r0 = cross(c[1], c[2]); // Used in cone volume

        // Calculate area
        const auto area = tag.area = magnitude(cp) * 0.5;
        const auto invA = 1.0 / area;

        // Calculate area gradients
        {
            const auto r01 = c[1] - c[0];
            const auto r02 = c[2] - c[0];
            mesh.attribute(hei[0]).gHalfEdge.dTriangleArea = (-l2[1]* r02 - l2[0]* r01 + dots[0]*(r01 + r02)) * (invA * 0.25);
            mesh.attribute(hei[1]).gHalfEdge.dTriangleArea = (l2[0]* r01 - dots[0]* r02) * (invA * 0.25);
            mesh.attribute(hei[2]).gHalfEdge.dTriangleArea = (l2[1]* r02 - dots[0]* r01) * (invA * 0.25);
        }

        // Calculate cot thetas and gradients
        for(size_t ai = 0; ai < 3; ++ai) {
            auto& heag = mesh.attribute(hei[ai]).gHalfEdge;

            heag.cotTheta = dots[ai] * invA * 0.5;

            const size_t ai_n = (ai + 1) % 3;
            const size_t ai_p = (ai + 2) % 3;
            auto& heag_n = mesh.attribute(hei[ai_n]).gHalfEdge;
            auto& heag_p = mesh.attribute(hei[ai_p]).gHalfEdge;

            const auto r01 = c[ai_n] - c[ai];
            const auto r02 = c[ai_p] - c[ai];

            heag.dCotTheta[1] =
                -(r01 + r02) * (invA * 0.5)
                -(dots[ai] * invA * invA * 0.5) * heag.dTriangleArea;
            heag.dCotTheta[2] = r02 * (invA * 0.5) - (dots[ai] * invA * invA * 0.5) * heag_n.dTriangleArea;
            heag.dCotTheta[0] = r01 * (invA * 0.5) - (dots[ai] * invA * invA * 0.5) * heag_p.dTriangleArea;
        }

        // Calculate cone volume and derivative
        tag.coneVolume = dot(c[0], r0) / 6.0;
        // The derivative of cone volume will be accumulated to each vertex

    }

    const auto& cvt = mesh.metaAttribute().cachedVertexTopo;

    // Calculate vertex 1-ring area and local curvature with derivative
    // Calculate derivative of volume on vertices
    for(MT::VertexIndex vi {0}; vi < numVertices; ++vi) if(!mesh.isVertexOnBorder(vi)) {
        auto& va = mesh.attribute(vi);
        auto& vag = va.gVertex;
        const auto ci = makeRefVec<3>(coord + va.cachedCoordIndex);

        // clearing
        vag.astar = 0.0;
        vag.dAstar = {0.0, 0.0, 0.0};
        vag.dVolume = {0.0, 0.0, 0.0};

        // K = 2 * H * n is the result of LB operator.
        // K * A = d AStar (on center vertex), where the choice of A could be different.
        // Here the we use Vector Area for the A above, i.e. AVec = |d Vol (on center vertex)|
        //
        //-----------------------------------------------------------------
        // Signed curvature implementation
        //
        // We choose the following implementation (see star_perp_sq_mean_curvature of Surface Evolver):
        //
        //        (1/2) d AStar  dot  d Vol
        //   H = ----------------------------
        //               d Vol   dot  d Vol
        //
        // Therefore (here D can operate on neighbor vertices, while d only on center vertex),
        //
        //         (1/2) D(d AStar) d Vol + (1/2) D(d Vol) d AStar - H * 2 D(d Vol) d Vol
        //   DH = -------------------------------------------------------------------------
        //                                  d Vol  dot  d Vol
        //
        // For intermediate variables, let
        //   t1 = D(d AStar) d Vol    (D on both central and neighbor vertices)
        //   t2 = D(d Vol) d AStar    (D on only neighbor vertices because D(d Vol) on center vertex is 0)
        //   t3 = D(d Vol) d Vol      (D on only neighbor vertices because D(d Vol) on center vertex is 0)
        // then
        //
        //          (1/2) t1 + (1/2) t2 - 2 H t3
        //   DH = --------------------------------
        //                d Vol  dot  d Vol
        //
        //-----------------------------------------------------------------
        // Squared curvature implementation
        //
        //             (1/2) D(∇ A_star) ∇ A_star - 2 H^2 D(∇ Vol) ∇ Vol
        //   D(H^2) = ---------------------------------------------------
        //                                (∇ Vol)^2
        //
        // And D can operate on center or neighbor vertices.
        //
        // An additional intermediate variable would be
        //   t4 = D(∇ A_star) ∇ A_star  (D on both central and neighbor vertices)
        // then
        //
        //             (1/2) t4 - 2 H^2 t3
        //   D(H^2) = ---------------------
        //                  (∇ Vol)^2

        // derivative of curvature will be calculated in the next loop
        for(size_t i = 0; i < va.cachedDegree; ++i) {
            const MT::HalfEdgeIndex hei_o  { cvt[mesh.metaAttribute().cachedVertexOffsetLeavingHE(vi.index) + i] };
            const MT::TriangleIndex ti0    { cvt[mesh.metaAttribute().cachedVertexOffsetPolygon(vi.index) + i] };
            const MT::HalfEdgeIndex hei_n  { cvt[mesh.metaAttribute().cachedVertexOffsetLeavingHE(vi.index) + (i + va.cachedDegree - 1) % va.cachedDegree] };
            const MT::HalfEdgeIndex hei_p  { cvt[mesh.metaAttribute().cachedVertexOffsetOuterHE(vi.index) + i] };
            const MT::HalfEdgeIndex hei_on { cvt[mesh.metaAttribute().cachedVertexOffsetOuterHE(vi.index) + (i + 1) % va.cachedDegree] };
            const auto cn      = makeRefVec<3>(coord + cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi.index) + i]);
            const auto c_right = makeRefVec<3>(coord + cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi.index) + (i + 1) % va.cachedDegree]);

            const auto sumCotTheta = mesh.attribute(hei_n).gHalfEdge.cotTheta + mesh.attribute(hei_on).gHalfEdge.cotTheta;

            const auto diff = ci - cn;

            vag.astar += mesh.attribute(ti0).gTriangle.area;

            // Accumulate dAstar
            vag.dAstar += 0.5 * sumCotTheta * diff;

            // Calculate dAstar_n (used in bending force calculation)
            mesh.attribute(hei_o).gHalfEdge.dNeighborAstar
                = mesh.attribute(hei_o).gHalfEdge.dTriangleArea
                + mesh.attribute(hei_p).gHalfEdge.dTriangleArea;

            // Added to derivative of sum of cone volume
            const auto cp = cross(cn, c_right);
            vag.dVolume += cp * (1.0 / 6);
        }

        const auto dVolume2 = magnitude2(vag.dVolume);

        if(curvPol == SurfaceCurvaturePolicy::withSign) {
            vag.curv = 0.5 * dot(vag.dAstar, vag.dVolume) / dVolume2;
        }
        else {
            vag.curv2 = 0.25 * magnitude2(vag.dAstar) / dVolume2;
        }
        // Derivative will be processed later.

        // Calculate derivative of curvature
        // Using another loop because d Vol, d AStar and H are needed for curvature derivative
        std::array<Vec3, 3> dDAstar {}; // On center vertex, indexed by [k1x, k1y, k1z]
        for(size_t i = 0; i < va.cachedDegree; ++i) {
            const MT::HalfEdgeIndex hei_o  { cvt[mesh.metaAttribute().cachedVertexOffsetLeavingHE(vi.index) + i] };
            const MT::HalfEdgeIndex hei_n  { cvt[mesh.metaAttribute().cachedVertexOffsetLeavingHE(vi.index) + (i + va.cachedDegree - 1) % va.cachedDegree] };
            const MT::HalfEdgeIndex hei_p  { cvt[mesh.metaAttribute().cachedVertexOffsetOuterHE(vi.index) + i] };
            const MT::HalfEdgeIndex hei_on { cvt[mesh.metaAttribute().cachedVertexOffsetOuterHE(vi.index) + (i + 1) % va.cachedDegree] };
            const auto cn      = makeRefVec<3>(coord + cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi.index) + i]);
            const auto c_left  = makeRefVec<3>(coord + cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi.index) + (i + va.cachedDegree - 1) % va.cachedDegree]);
            const auto c_right = makeRefVec<3>(coord + cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi.index) + (i + 1) % va.cachedDegree]);

            const auto sumCotTheta = mesh.attribute(hei_n).gHalfEdge.cotTheta + mesh.attribute(hei_on).gHalfEdge.cotTheta;
            const auto& dCotThetaLeft = mesh.attribute(hei_n).gHalfEdge.dCotTheta;
            const auto& dCotThetaRight = mesh.attribute(hei_on).gHalfEdge.dCotTheta;
            const auto sumDCotThetaCenter = dCotThetaLeft[0] + dCotThetaRight[2];
            const auto sumDCotThetaNeighbor = dCotThetaLeft[2] + dCotThetaRight[0];

            const auto diff = ci - cn;
            // Accumulate dDAstar on the center vertex vi
            dDAstar[0] += (0.5 * sumDCotThetaCenter[0]) * diff;
            dDAstar[1] += (0.5 * sumDCotThetaCenter[1]) * diff;
            dDAstar[2] += (0.5 * sumDCotThetaCenter[2]) * diff;
            dDAstar[0][0] += 0.5 * sumCotTheta;
            dDAstar[1][1] += 0.5 * sumCotTheta;
            dDAstar[2][2] += 0.5 * sumCotTheta; // dDAstar += 0.5 * I * sumCotTheta, where I is gradient of diff (identity)

            // Calculate dDAstar and derivative of curvature on neighbor vertex vn
            std::array<Vec3, 3> dDAstar_n {};
            // As direct target
            dDAstar_n[0] = (0.5 * sumDCotThetaNeighbor[0]) * diff;
            dDAstar_n[1] = (0.5 * sumDCotThetaNeighbor[1]) * diff;
            dDAstar_n[2] = (0.5 * sumDCotThetaNeighbor[2]) * diff;
            dDAstar_n[0][0] -= 0.5 * sumCotTheta;
            dDAstar_n[1][1] -= 0.5 * sumCotTheta;
            dDAstar_n[2][2] -= 0.5 * sumCotTheta; // dK1 += -0.5 * I * sumCotTheta

            // As target for left and right
            const auto diff_left = ci - c_left;
            const auto diff_right = ci - c_right;
            const auto& dCotThetaOfLeft = mesh.attribute(hei_p).gHalfEdge.dCotTheta[1];
            const auto& dCotThetaOfRight = mesh.attribute(hei_o).gHalfEdge.dCotTheta[1];
            dDAstar_n[0] += (0.5 * dCotThetaOfLeft[0]) * diff_left;
            dDAstar_n[1] += (0.5 * dCotThetaOfLeft[1]) * diff_left;
            dDAstar_n[2] += (0.5 * dCotThetaOfLeft[2]) * diff_left;
            dDAstar_n[0] += (0.5 * dCotThetaOfRight[0]) * diff_right;
            dDAstar_n[1] += (0.5 * dCotThetaOfRight[1]) * diff_right;
            dDAstar_n[2] += (0.5 * dCotThetaOfRight[2]) * diff_right;

            // D_n (d Vol) = (1/6) D_n (c_left x cn + cn x c_right)
            //             = (1/6) D_n (cn x (c_right - c_left))
            // Then for any vector v,
            // [D_n (d Vol)] v = (1/6) (c_right - c_left) x v
            const auto vec_lr = c_right - c_left;

            // Compute t1_n, t2_n and t3_n
            const Vec3 t3_n = (1.0 / 6) * cross(vec_lr, vag.dVolume);

            if(curvPol == SurfaceCurvaturePolicy::withSign) {
                const Vec3 t1_n {
                    dot(dDAstar_n[0], vag.dVolume),
                    dot(dDAstar_n[1], vag.dVolume),
                    dot(dDAstar_n[2], vag.dVolume)
                };
                const Vec3 t2_n = (1.0 / 6) * cross(vec_lr, vag.dAstar);

                // Derivative of curvature
                mesh.attribute(hei_o).gHalfEdge.dNeighborCurv = (0.5 * (t1_n + t2_n) - (2 * vag.curv) * t3_n) / dVolume2;
            }
            else {
                const Vec3 t4_n {
                    dot(dDAstar_n[0], vag.dAstar),
                    dot(dDAstar_n[1], vag.dAstar),
                    dot(dDAstar_n[2], vag.dAstar)
                };

                // Derivative of curvature squared
                mesh.attribute(hei_o).gHalfEdge.dNeighborCurv2 = (0.5 * t4_n - (2 * vag.curv2) * t3_n) / dVolume2;
            }
        } // End loop neighbor vertices

        if(curvPol == SurfaceCurvaturePolicy::withSign) {
            const Vec3 t1 {
                dot(dDAstar[0], vag.dVolume),
                dot(dDAstar[1], vag.dVolume),
                dot(dDAstar[2], vag.dVolume)
            };

            // Also the derivative of curvature on central vertex
            vag.dCurv = t1 * (0.5 / dVolume2);
        }
        else {
            const Vec3 t4 {
                dot(dDAstar[0], vag.dAstar),
                dot(dDAstar[1], vag.dAstar),
                dot(dDAstar[2], vag.dAstar)
            };

            // Derivative of curvature squared on central vertex
            vag.dCurv2 = t4 * (0.5 / dVolume2);
        }

    } // End loop vertices (V cells)

} // updateGeometryValueWithDerivative(...)

// This function updates geometries necessary for MEDYAN system. Currently
// the system needs
//   - (pseudo) unit normals
//   - triangle areas
// This function uses cached indexing to enhance performance, so a valid
// cache is needed in this function.
inline void updateGeometryValueForSystem(MembraneMeshAttribute::MeshType& mesh) {
    using namespace mathfunc;
    using MT = MembraneMeshAttribute::MeshType;

    const auto& vertices = mesh.getVertices();
    const auto& edges = mesh.getEdges();
    const auto& triangles = mesh.getTriangles();

    const size_t numVertices = vertices.size();
    const size_t numEdges = edges.size();
    const size_t numTriangles = triangles.size();

    // Calculate triangle unit normal and area
    for(MT::TriangleIndex ti {}; ti < numTriangles; ++ti) {
        auto& ta = mesh.attribute(ti);
        auto& tag = ta.gTriangle;
        const auto vis = vertexIndices(mesh, ti);
        const auto& c0 = mesh.attribute(vis[0]).vertex->coord;
        const auto& c1 = mesh.attribute(vis[1]).vertex->coord;
        const auto& c2 = mesh.attribute(vis[2]).vertex->coord;

        const auto cp = cross(c1 - c0, c2 - c0);

        // unit normal
        tag.unitNormal = normalizedVector(cp);

        // area
        tag.area = magnitude(cp) * 0.5;
    }

    // Calculate edge pesudo unit normal
    for(MT::EdgeIndex ei {}; ei < numEdges; ++ei) {
        auto& ea = mesh.attribute(ei);
        const auto hei = mesh.halfEdge(ea);
        const auto hei_o = mesh.opposite(hei);

        // pseudo unit normal
        if(
            mesh.isInTriangle(hei  ) &&
            mesh.isInTriangle(hei_o)
        ) {
            ea.gEdge.pseudoUnitNormal = normalizedVector(
                mesh.attribute(mesh.triangle(hei  )).gTriangle.unitNormal +
                mesh.attribute(mesh.triangle(hei_o)).gTriangle.unitNormal
            );
        }
        else if(mesh.isInTriangle(hei)) {
            ea.gEdge.pseudoUnitNormal = mesh.attribute(mesh.triangle(hei)).gTriangle.unitNormal;
        }
        else if(mesh.isInTriangle(hei_o)) {
            ea.gEdge.pseudoUnitNormal = mesh.attribute(mesh.triangle(hei_o)).gTriangle.unitNormal;
        }
    }

    const auto& cvt = mesh.metaAttribute().cachedVertexTopo;

    // Calculate vertex pseudo unit normal
    for(MT::VertexIndex vi {}; vi < numVertices; ++vi) {
        auto& va = mesh.attribute(vi);
        auto& vag = va.gVertex;

        // clearing
        vag.pseudoUnitNormal = {0.0, 0.0, 0.0};

        mesh.forEachHalfEdgeTargetingVertex(vi, [&](MT::HalfEdgeIndex hei) {
            if(mesh.isInTriangle(hei)) {
                const auto ti0 = mesh.triangle(hei);

                const auto theta = mesh.attribute(hei).gHalfEdge.theta;

                vag.pseudoUnitNormal += theta * mesh.attribute(ti0).gTriangle.unitNormal;
            }
        });

        normalize(vag.pseudoUnitNormal);
    }
} // void updateGeometryValue(...)

// Signed distance using geometric attributes (the inefficient way)
/**************************************************************************
The function works in the following procedure:

- Iterate through the triangles, and for each triangle
    - Find the projection of the point on the triangle plane, and determine
        which element is responsible for being the closest to the point
        (which vertex/ which edge/ this triangle).
    - Find the unsigned distance with the closest element, and then find
        the signed distance using the normal or pseudo normal. Record the
        value with the smallest unsigned distance.

Before this function is used, the following must be calculated:
    - The positions of all the elements are updated
    - The normal and pseudo normal at the triangles, edges and vertices

Note: this method only works if the mesh is closed. This must be ensured by
        the caller of the function.

In fact, the signed distance field serves as a good candidate for membrane
boundary potential. However, this field is not C1-continuous everywhere,
which is detrimental to conjugate gradient methods.
**************************************************************************/
template< typename VecType, std::enable_if_t< VecType::vec_size == 3 >* = nullptr >
inline double signedDistance(const MembraneMeshAttribute::MeshType& mesh, const VecType& p) {
    using namespace mathfunc;
    using MT = MembraneMeshAttribute::MeshType;
    using CT = MembraneMeshAttribute::CoordinateType;

    const size_t numTriangles = mesh.getTriangles().size();

    double minAbsDistance = numeric_limits<double>::infinity();
    for(MT::TriangleIndex ti {}; ti < numTriangles; ++ti) {
        /**********************************************************************
        Calculate the barycentric coordinate of the projection point p'

        See Heidrich 2005, Computing the Barycentric Coordinates of a Projected
        Point.
        **********************************************************************/
        const auto hei0 = mesh.halfEdge(ti);
        const auto hei1 = mesh.next(hei0);
        const auto hei2 = mesh.next(hei1);
        const MT::VertexIndex vi[] {
            mesh.target(hei0), mesh.target(hei1), mesh.target(hei2)
        };
        const CT c[] {
            mesh.attribute(vi[0]).vertex->coord,
            mesh.attribute(vi[1]).vertex->coord,
            mesh.attribute(vi[2]).vertex->coord
        };

        const auto r01 = c[1] - c[0];
        const auto r02 = c[2] - c[0];
        const auto r0p = p - c[0];
        const auto cp = cross(r01, r02);
        const auto oneOver4AreaSquared = 1.0 / magnitude2(cp);

        const auto b1 = dot(cross(r0p, r02), cp) * oneOver4AreaSquared;
        const auto b2 = dot(cross(r01, r0p), cp) * oneOver4AreaSquared;
        const auto b0 = 1.0 - b1 - b2;

        // Now p' = b0*v0 + b1*v1 + b2*v2
        // which is the projection of p in the plane of the triangle

        double d = numeric_limits<double>::infinity();
        if(b0 >= 0 && b1 >= 0 && b2 >= 0) {
            // p' is inside the triangle
            d = dot(mesh.attribute(ti).gTriangle.unitNormal, r0p);
        } else {
            // p' is outside the triangle
            const Vec< 3, typename CT::float_type > r2 {
                distance2(c[1], c[2]),
                distance2(c[2], c[0]),
                distance2(c[0], c[1])
            };
            const auto r1p = p - c[1];
            const auto r2p = p - c[2];
            const auto r12 = c[2] - c[1];
            const auto dot_1p_12 = dot(r1p, r12);
            const auto dot_2p_20 = -dot(r2p, r02);
            const auto dot_0p_01 = dot(r0p, r01);

            if(b0 < 0 && dot_1p_12 >= 0 && dot_1p_12 <= r2[0]) {
                // On edge 12
                d = magnitude(cross(r1p, r12)) / std::sqrt(r2[0]);
                if(dot(mesh.attribute(mesh.edge(hei2)).gEdge.pseudoUnitNormal, r1p) < 0) d = -d;
            } else if(b1 < 0 && dot_2p_20 >= 0 && dot_2p_20 <= r2[1]) {
                // On edge 20
                d = magnitude(cross(r2p, r02)) / std::sqrt(r2[1]);
                if(dot(mesh.attribute(mesh.edge(hei0)).gEdge.pseudoUnitNormal, r2p) < 0) d = -d;
            } else if(b2 < 0 && dot_0p_01 >= 0 && dot_0p_01 <= r2[2]) {
                // On edge 01
                d = magnitude(cross(r0p, r01)) / std::sqrt(r2[2]);
                if(dot(mesh.attribute(mesh.edge(hei1)).gEdge.pseudoUnitNormal, r0p) < 0) d = -d;
            } else if(dot_0p_01 < 0 && dot_2p_20 > r2[1]) {
                // On vertex 0
                d = distance(c[0], p);
                if(dot(mesh.attribute(vi[0]).gVertex.pseudoUnitNormal, r0p) < 0) d = -d;
            } else if(dot_1p_12 < 0 && dot_0p_01 > r2[2]) {
                // On vertex 1
                d = distance(c[1], p);
                if(dot(mesh.attribute(vi[1]).gVertex.pseudoUnitNormal, r1p) < 0) d = -d;
            } else if(dot_2p_20 < 0 && dot_1p_12 > r2[0]) {
                // On vertex 2
                d = distance(c[2], p);
                if(dot(mesh.attribute(vi[2]).gVertex.pseudoUnitNormal, r2p) < 0) d = -d;
            } else {
                // The program should never come here
                throw logic_error("Unknown case of point projection on the plane of triangle.");
            }
        }

        // Update with distance with less absolute value
        if(abs(d) < abs(minAbsDistance)) minAbsDistance = d;
    }
    
    return minAbsDistance;
}
template< typename VecType, std::enable_if_t< VecType::vec_size == 3 >* = nullptr >
inline bool contains(const MembraneMeshAttribute::MeshType& mesh, const VecType& p) {
    return signedDistance(mesh, p) < 0.0;
}

// Attribute computation in adaptive remeshing algorithms
//-------------------------------------------------------------------------

// Triangle normal (Geometric attribute) used in adaptive remeshing
inline void adaptiveComputeTriangleNormal(
    MembraneMeshAttribute::MeshType& mesh,
    MembraneMeshAttribute::MeshType::TriangleIndex ti
) {
    const auto hei = mesh.halfEdge(ti);
    const auto vi0 = mesh.target(hei);
    const auto vi1 = mesh.target(mesh.next(hei));
    const auto vi2 = mesh.target(mesh.prev(hei));
    const auto& c0 = mesh.attribute(vi0).vertex->coord;
    const auto& c1 = mesh.attribute(vi1).vertex->coord;
    const auto& c2 = mesh.attribute(vi2).vertex->coord;
    auto& tag = mesh.attribute(ti).gTriangle;

    const auto cp = mathfunc::cross(c1 - c0, c2 - c0);

    // unit normal
    tag.unitNormal = mathfunc::normalizedVector(cp);
}

// Triangle angles (Geometric attribute, halfedge) used in adaptive remeshing
inline void adaptiveComputeAngle(
    MembraneMeshAttribute::MeshType& mesh,
    MembraneMeshAttribute::MeshType::HalfEdgeIndex hei
) {
    // The angle is (v0, v1, v2)
    const auto vi0 = mesh.target(mesh.prev(hei));
    const auto vi1 = mesh.target(hei);
    const auto vi2 = mesh.target(mesh.next(hei));
    const auto& c0 = mesh.attribute(vi0).vertex->coord;
    const auto& c1 = mesh.attribute(vi1).vertex->coord;
    const auto& c2 = mesh.attribute(vi2).vertex->coord;
    auto& heag = mesh.attribute(hei).gHalfEdge;

    const auto cp = mathfunc::cross(c0 - c1, c2 - c1);
    const auto dp = mathfunc::  dot(c0 - c1, c2 - c1);
    const auto ct = heag.cotTheta = dp / mathfunc::magnitude(cp);
    heag.theta = M_PI_2 - std::atan(ct);
}

// Vertex unit normals (Adaptive attribute) used in adaptive remeshing
// Requires
//   - Unit normals in triangles (geometric)
//   - Angles in halfedges (geometric)
inline void adaptiveComputeVertexNormal(
    MembraneMeshAttribute::MeshType& mesh,
    MembraneMeshAttribute::MeshType::VertexIndex vi
) {
    using MT = MembraneMeshAttribute::MeshType;

    // Using pseudo normal (weighted by angles)
    auto& vaa = mesh.attribute(vi).aVertex;

    // clearing
    vaa.unitNormal = {0.0, 0.0, 0.0};

    mesh.forEachHalfEdgeTargetingVertex(vi, [&](MT::HalfEdgeIndex hei) {
        if(mesh.polygonType(hei) == MT::PolygonType::triangle) {
            const auto ti0 = mesh.triangle(hei);
            const auto theta = mesh.attribute(hei).gHalfEdge.theta;
            vaa.unitNormal += theta * mesh.attribute(ti0).gTriangle.unitNormal;
        }
    });

    mathfunc::normalize(vaa.unitNormal);
}

//-------------------------------------------------------------------------
// Mesh modifiers
//-------------------------------------------------------------------------

// The HalfEdgeMesh already provides basic mesh operations which handles
// the mesh element connections.
// In these operations, some attributes of the mesh need to be handled
// specifically. For example, in material coordinate system, the
// equilibrium area of neighboring triangles after a edge flip need to be
// re-distributed.
// However, these operations cannot be simply achieved by injecting
// behavior (such as passing functions) to the mesh operators of the mesh
// connection class (HalfEdgeMesh), because only with full attribute info
// can one decide what to do/access/modify before/after the mesh
// reconnection.
inline auto insertVertexOnEdge(
    MembraneMeshAttribute::MeshType&             mesh,
    MembraneMeshAttribute::MeshType::EdgeIndex   ei,
    const MembraneMeshAttribute::CoordinateType& newPos
) {
    using MT = MembraneMeshAttribute::MeshType;

    const auto ohei       = mesh.halfEdge(ei);
    const auto ohei_o     = mesh.opposite(ohei);

    const auto opt0       = mesh.polygonType(ohei);
    const bool ist0       = opt0 == MT::PolygonType::triangle;
    const auto opt2       = mesh.polygonType(ohei_o);
    const bool ist2       = opt2 == MT::PolygonType::triangle;

    // Do the vertex insertion
    const auto change = MT::VertexInsertionOnEdge{}(mesh, ei);

    const auto redisEqArea = [](
        double& eqArea1, double a1,
        double& eqArea2, double a2,
        double totalEqArea
    ) {
        // Postcondition: sum of eq areas is equal to the total eq area.
        eqArea1 = totalEqArea * a1 / (a1 + a2);
        eqArea2 = totalEqArea * a2 / (a1 + a2);
    };

    // Set attributes
    mesh.attribute(change.viNew).vertex->coord = newPos;

    // Redistribute properties
    if(ist0) {
        const double a1 = area(mesh, change.tiNew[0]);
        const double a2 = area(mesh, change.tiNew[1]);

        auto& eqArea1 = mesh.attribute(change.tiNew[0]).triangle->mTriangle.eqArea;
        auto& eqArea2 = mesh.attribute(change.tiNew[1]).triangle->mTriangle.eqArea;
        redisEqArea(eqArea1, a1, eqArea2, a2, eqArea1);
    }
    if(ist2) {
        const double a1 = area(mesh, change.tiNew[2]);
        const double a2 = area(mesh, change.tiNew[3]);

        auto& eqArea1 = mesh.attribute(change.tiNew[2]).triangle->mTriangle.eqArea;
        auto& eqArea2 = mesh.attribute(change.tiNew[3]).triangle->mTriangle.eqArea;
        redisEqArea(eqArea1, a1, eqArea2, a2, eqArea1);
    }

    return change;
}

inline auto collapseEdge(
    MembraneMeshAttribute::MeshType&               mesh,
    MembraneMeshAttribute::MeshType::EdgeIndex     ei,
    const std::optional< MembraneMeshAttribute::CoordinateType >& newPos
) {
    using MT = MembraneMeshAttribute::MeshType;

    const auto hei = mesh.halfEdge(ei);
    const auto hei_o = mesh.opposite(hei);
    const auto vi0 = mesh.target(hei);
    const auto vi1 = mesh.target(hei_o);

    // Record properties
    double oldTotalEqArea = 0.0;
    // Accumulate around two vertices, and subtract triangles on the middle edge
    const auto addOldTotalEqArea = [&](MT::HalfEdgeIndex hei) {
        if(mesh.isInTriangle(hei)) {
            oldTotalEqArea += mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea;
        }
    };
    const auto subtractOldTotalEqArea = [&](MT::HalfEdgeIndex hei) {
        if(mesh.isInTriangle(hei)) {
            oldTotalEqArea -= mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea;
        }
    };
    mesh.forEachHalfEdgeTargetingVertex(vi0, addOldTotalEqArea);
    mesh.forEachHalfEdgeTargetingVertex(vi1, addOldTotalEqArea);
    mesh.forEachHalfEdgeInEdge(ei, subtractOldTotalEqArea);

    // Do the edge collapse
    const auto change = MT::EdgeCollapse{}(mesh, ei);

    // Set attributes
    if(newPos.has_value()) mesh.attribute(change.viTo).vertex->coord = *newPos;

    // Redistribute properties
    double newTotalArea = 0.0;
    mesh.forEachHalfEdgeTargetingVertex(change.viTo, [&](MT::HalfEdgeIndex hei) {
        if(mesh.isInTriangle(hei)) {
            newTotalArea += area(mesh, mesh.triangle(hei));
        }
    });
    mesh.forEachHalfEdgeTargetingVertex(change.viTo, [&](MT::HalfEdgeIndex hei) {
        if(mesh.isInTriangle(hei)) {
            mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea
                = oldTotalEqArea * area(mesh, mesh.triangle(hei)) / newTotalArea;
        }
    });

    return change;
}

inline void flipEdge(
    MembraneMeshAttribute::MeshType&               mesh,
    MembraneMeshAttribute::MeshType::EdgeIndex     ei,
) {
    using MT = MembraneMeshAttribute::MeshType;

    // Precondition:
    //   - ei is not on the border

    // Record old attributes
    double oldTotalEqArea = 0.0;
    mesh.forEachHalfEdgeInEdge(ei, [&](MT::HalfEdgeIndex hei) {
        oldTotalEqArea += mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea;
    });

    // Do edge flip
    MT::EdgeFlip{}(mesh, ei);

    // Redistribute attributes
    double newTotalArea = 0.0;
    mesh.forEachHalfEdgeInEdge(ei, [&](MT::HalfEdgeIndex hei) {
        newTotalArea += area(mesh, mesh.triangle(hei));
    });
    mesh.forEachHalfEdgeInEdge(ei, [&](MT::HalfEdgeIndex hei) {
        mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea
            = oldTotalEqArea * area(mesh, mesh.triangle(hei)) / newTotalArea;
    });

}


} // namespace medyan

#endif
