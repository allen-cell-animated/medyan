#ifndef MEDYAN_Structure_SurfaceMesh_MembraneMeshAttribute_hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneMeshAttribute_hpp

#include <algorithm> // max
#include <array>
#include <limits> // numeric_limits
#include <memory> // unique_ptr
#include <stdexcept> // logic_error
#include <vector>

#include "MathFunctions.h"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/AdaptiveMeshAttribute.hpp"
#include "Structure/SurfaceMesh/Edge.hpp"
#include "Structure/SurfaceMesh/GeometricMeshAttribute.hpp"
#include "Structure/SurfaceMesh/SurfaceMesh.hpp"
#include "Structure/SurfaceMesh/Triangle.hpp"
#include "Structure/SurfaceMesh/Vertex.hpp"
#include "Util/Io/Log.hpp"

/******************************************************************************
Implements the attributes of the meshwork used by the membrane, mainly
geometric attribute and adaptive mesh attributes.

Geometric attributes are mainly used in force field computations.
Adaptive attributes provides additional attributes mainly used in adaptive mesh
algorithm.

This struct can connect the topological meshwork with the actual system objects.
******************************************************************************/
struct MembraneMeshAttribute {

    using MeshType = HalfEdgeMesh< MembraneMeshAttribute >;

    struct VertexAttribute {
        using coordinate_type      = Vertex::coordinate_type;
        using coordinate_ref_type  = Vertex::coordinate_ref_type;
        using coordinate_cref_type = Vertex::coordinate_cref_type;

        std::unique_ptr< Vertex > vertex;

        GVertex gVertex;
        GVertex gVertexS; // stretched version (temporary)
        AdaptiveMeshAttribute::VertexAttribute aVertex;

        size_t cachedDegree;
        size_t cachedCoordIndex;

        coordinate_ref_type getCoordinate() const { return vertex->coordinate(); }

        template< bool stretched > const GVertex& getGVertex() const { return stretched ? gVertexS : gVertex; }
        template< bool stretched >       GVertex& getGVertex()       { return stretched ? gVertexS : gVertex; }

        void setIndex(size_t index) {
            vertex->setTopoIndex(index);
        }

    };
    struct EdgeAttribute {
        std::unique_ptr<Edge> edge;

        GEdge gEdge;
        GEdge gEdgeS; // stretched version (temporary)
        AdaptiveMeshAttribute::EdgeAttribute aEdge;

        // The index in Bead::getDbData()
        // [target of half edge, opposite target]
        size_t                                   cachedCoordIndex[2];
        typename MeshType::HalfEdge::PolygonType cachedPolygonType[2];
        // [polygon of half edge, opposite polygon]
        size_t                                   cachedPolygonIndex[2];

        template< bool stretched > const GEdge& getGEdge() const { return stretched ? gEdgeS : gEdge; }
        template< bool stretched >       GEdge& getGEdge()       { return stretched ? gEdgeS : gEdge; }

        void setIndex(size_t index) {
            edge->setTopoIndex(index);
        }
    };
    struct HalfEdgeAttribute {
        GHalfEdge gHalfEdge;
        GHalfEdge gHalfEdgeS; // stretched version (temporary)
        AdaptiveMeshAttribute::HalfEdgeAttribute aHalfEdge;

        // The index in Bead::getDbData()
        // [source, target, next target]
        size_t cachedCoordIndex[3];

        template< bool stretched > const GHalfEdge& getGHalfEdge() const { return stretched ? gHalfEdgeS : gHalfEdge; }
        template< bool stretched >       GHalfEdge& getGHalfEdge()       { return stretched ? gHalfEdgeS : gHalfEdge; }

        void setIndex(size_t index) {}
    };
    struct TriangleAttribute {
        std::unique_ptr<Triangle> triangle;

        GTriangle gTriangle;
        GTriangle gTriangleS; // stretched version (temporary);
        AdaptiveMeshAttribute::TriangleAttribute aTriangle;

        // The index in Bead::getDbData()
        // [target of half edge, next target, prev target]
        size_t cachedCoordIndex[3];
        // [half edge, next, prev]
        size_t cachedHalfEdgeIndex[3];
        // [edge of half edge, next edge, prev edge]
        size_t cachedEdgeIndex[3];

        template< bool stretched > const GTriangle& getGTriangle() const { return stretched ? gTriangleS : gTriangle; }
        template< bool stretched >       GTriangle& getGTriangle()       { return stretched ? gTriangleS : gTriangle; }

        void setIndex(size_t index) {
            triangle->setTopoIndex(index);
        }
    };
    struct BorderAttribute {
        void setIndex(size_t index) {}
    };
    struct MetaAttribute {
        SubSystem *s;
        Membrane *m;

        bool cacheValid = false;

        size_t vertexMaxDegree;

        // A vertex has undetermined number of neighbors, so the cache structure needs to be determined at run-time.
        std::vector< size_t > cachedVertexTopo;
        size_t cachedVertexTopoSize() const { return vertexMaxDegree * 5; }
        size_t cachedVertexOffsetNeighborCoord(size_t idx) const { return cachedVertexTopoSize() * idx; }
        size_t cachedVertexOffsetTargetingHE  (size_t idx) const { return cachedVertexTopoSize() * idx + vertexMaxDegree; }
        size_t cachedVertexOffsetLeavingHE    (size_t idx) const { return cachedVertexTopoSize() * idx + vertexMaxDegree * 2; }
        size_t cachedVertexOffsetOuterHE      (size_t idx) const { return cachedVertexTopoSize() * idx + vertexMaxDegree * 3; }
        size_t cachedVertexOffsetPolygon      (size_t idx) const { return cachedVertexTopoSize() * idx + vertexMaxDegree * 4; }
    };

    using coordinate_type = typename VertexAttribute::coordinate_type;

    struct AttributeInitializerInfo {
        std::vector< coordinate_type > vertexCoordinateList;
    };

    enum class GeometryCurvaturePolicy {
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

    // Mesh element modification (not used in initialization/finalization)
    template< typename VertexInsertionMethod >
    static void newVertex(MeshType& mesh, MeshType::VertexIndex v, const VertexInsertionMethod& op) {
        mesh.attribute(v).vertex = std::make_unique< Vertex >(
            op.coordinate(mesh, v), v);
    }
    template< typename Operation >
    static void newEdge(MeshType& mesh, MeshType::EdgeIndex e, const Operation& op) {
        mesh.attribute(e).edge.reset(
            mesh.metaAttribute().s->template addTrackable<Edge>(mesh.metaAttribute().m, e));
    }
    template< typename Operation >
    static void newHalfEdge(MeshType& mesh, MeshType::HalfEdgeIndex he, const Operation& op) {
        // Do nothing
    }
    template< typename Operation >
    static void newTriangle(MeshType& mesh, MeshType::TriangleIndex t, const Operation& op) {
        mesh.attribute(t).triangle.reset(mesh.metaAttribute().s->template addTrackable<Triangle>(mesh.metaAttribute().m, t));
    }
    template< typename Operation >
    static void newBorder(MeshType& mesh, MeshType::BorderIndex, const Operation& op) {
        // Do nothing
    }

    static void removeElement(MeshType& mesh, MeshType::VertexIndex i) {
        // Do nothing
    }
    static void removeElement(MeshType& mesh, MeshType::EdgeIndex i) {
        mesh.metaAttribute().s->template removeTrackable<Edge>(mesh.attribute(i).edge.get());
    }
    static void removeElement(MeshType& mesh, MeshType::HalfEdgeIndex i) {
        // Do nothing
    }
    static void removeElement(MeshType& mesh, MeshType::TriangleIndex i) {
        mesh.metaAttribute().s->template removeTrackable<Triangle>(mesh.attribute(i).triangle.get());
    }
    static void removeElement(MeshType& mesh, MeshType::BorderIndex i) {
        // Do nothing
    }

    // Mesh index caching
    //
    // The purpose of this function is to reduce pointer chasing while
    // traversing the mesh structure.
    template< bool forceUpdate = false >
    static void cacheIndices(MeshType& mesh) {

        if(forceUpdate || !mesh.metaAttribute().cacheValid) {
            const auto& vertices = mesh.getVertices();
            const auto& halfEdges = mesh.getHalfEdges();
            const auto& edges = mesh.getEdges();
            const auto& triangles = mesh.getTriangles();

            const size_t numVertices = vertices.size();
            const size_t numHalfEdges = halfEdges.size();
            const size_t numEdges = edges.size();
            const size_t numTriangles = triangles.size();

            for(size_t hei = 0; hei < numHalfEdges; ++hei) {
                // The angle is (v0, v1, v2)
                const size_t vi0 = mesh.target(mesh.prev(hei));
                const size_t vi1 = mesh.target(hei);
                const size_t vi2 = mesh.target(mesh.next(hei));

                auto& hea = mesh.getHalfEdgeAttribute(hei);
                hea.cachedCoordIndex[0] = vertices[vi0].attr.vertex->Bead::getStableIndex();
                hea.cachedCoordIndex[1] = vertices[vi1].attr.vertex->Bead::getStableIndex();
                hea.cachedCoordIndex[2] = vertices[vi2].attr.vertex->Bead::getStableIndex();
            }

            for(size_t ti = 0; ti < numTriangles; ++ti) {
                const size_t hei = triangles[ti].halfEdgeIndex;
                const size_t hei_n = mesh.next(hei);
                const size_t hei_p = mesh.prev(hei);
                const size_t vi0 = mesh.target(hei);
                const size_t vi1 = mesh.target(hei_n);
                const size_t vi2 = mesh.target(hei_p);

                auto& ta = mesh.getTriangleAttribute(ti);
                ta.cachedCoordIndex[0] = vertices[vi0].attr.vertex->Bead::getStableIndex();
                ta.cachedCoordIndex[1] = vertices[vi1].attr.vertex->Bead::getStableIndex();
                ta.cachedCoordIndex[2] = vertices[vi2].attr.vertex->Bead::getStableIndex();
                ta.cachedHalfEdgeIndex[0] = hei;
                ta.cachedHalfEdgeIndex[1] = hei_n;
                ta.cachedHalfEdgeIndex[2] = hei_p;
                ta.cachedEdgeIndex[0] = mesh.edge(hei);
                ta.cachedEdgeIndex[1] = mesh.edge(hei_n);
                ta.cachedEdgeIndex[2] = mesh.edge(hei_p);
            }

            for(size_t ei = 0; ei < numEdges; ++ei) {
                const size_t hei = edges[ei].halfEdgeIndex;
                const size_t hei_o = mesh.opposite(hei);
                const size_t vi0 = mesh.target(hei);
                const size_t vi1 = mesh.target(mesh.prev(hei));

                auto& ea = mesh.getEdgeAttribute(ei);
                ea.cachedCoordIndex[0]   = vertices[vi0].attr.vertex->Bead::getStableIndex();
                ea.cachedCoordIndex[1]   = vertices[vi1].attr.vertex->Bead::getStableIndex();
                ea.cachedPolygonType[0]  = halfEdges[hei].polygonType;
                ea.cachedPolygonType[1]  = halfEdges[hei_o].polygonType;
                ea.cachedPolygonIndex[0] = mesh.polygon(hei);
                ea.cachedPolygonIndex[1] = mesh.polygon(hei_o);
            }

            // Determine vertex max number of neighbors
            mesh.metaAttribute().vertexMaxDegree = 0;
            for(size_t vi = 0; vi < numVertices; ++vi) {
                mesh.metaAttribute().vertexMaxDegree = std::max(
                    mesh.metaAttribute().vertexMaxDegree,
                    mesh.degree(vi)
                );
            }
            mesh.metaAttribute().cachedVertexTopo.resize(mesh.metaAttribute().cachedVertexTopoSize() * numVertices);
            // Cache indices around vertices
            for(size_t vi = 0; vi < numVertices; ++vi) {
                auto& va = mesh.getVertexAttribute(vi);

                va.cachedDegree = mesh.degree(vi);
                va.cachedCoordIndex = vertices[vi].attr.vertex->Bead::getStableIndex();

                size_t i = 0; // Neighbor index

                mesh.forEachHalfEdgeTargetingVertex(vi, [&](size_t hei) {
                    const size_t hei_o = mesh.opposite(hei);
                    const size_t vn = mesh.target(hei_o);

                    mesh.metaAttribute().cachedVertexTopo[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi) + i] = vertices[vn].attr.vertex->Bead::getStableIndex();
                    mesh.metaAttribute().cachedVertexTopo[mesh.metaAttribute().cachedVertexOffsetTargetingHE  (vi) + i] = hei;
                    mesh.metaAttribute().cachedVertexTopo[mesh.metaAttribute().cachedVertexOffsetLeavingHE    (vi) + i] = hei_o;
                    mesh.metaAttribute().cachedVertexTopo[mesh.metaAttribute().cachedVertexOffsetOuterHE      (vi) + i] = mesh.prev(hei);
                    mesh.metaAttribute().cachedVertexTopo[mesh.metaAttribute().cachedVertexOffsetPolygon      (vi) + i] = mesh.polygon(hei);

                    ++i;
                });
            }

            mesh.metaAttribute().cacheValid = true;
        }

    } // void cacheIndices(...)

    // Mesh attribute initializing and extracting
    // These operations do not follow the RAII idiom.
    // Initialization should happen only once, as it allocates resources
    static void init(MeshType& mesh, const AttributeInitializerInfo& info) {
        const MetaAttribute& meta = mesh.metaAttribute();
        for(size_t i = 0; i < mesh.getVertices().size(); ++i) {
            mesh.attribute(MeshType::VertexIndex{i}).vertex = std::make_unique(
                info.vertexCoordinateList[i], i);
        }
        for(size_t i = 0; i < mesh.getEdges().size(); ++i) {
            mesh.attribute(MeshType::EdgeIndex{i}).edge.reset(meta.s->template addTrackable<Edge>(meta.m, i));
        }
        for(size_t i = 0; i < mesh.getTriangles().size(); ++i) {
            mesh.attribute(MeshType::TriangleIndex{i}).triangle.reset(meta.s->template addTrackable<Triangle>(meta.m, i));
        }
    }
    // Extraction can be done multiple times without allocating/deallocating
    static auto extract(const MeshType& mesh) {
        AttributeInitializerInfo info;

        const size_t numVertices = mesh.getVertices().size();

        info.vertexCoordinateList.reserve(numVertices);
        for(size_t i = 0; i < numVertices; ++i) {
            info.vertexCoordinateList.push_back(
                static_cast<coordinate_type>(mesh.attribute(MeshType:VertexIndex{i}).vertex->coord));
        }

        return info;
    }

    //-------------------------------------------------------------------------
    // Geometries
    //-------------------------------------------------------------------------

    // This function updates geometries necessary for computing membrane energy
    // and signed distance.
    // This function uses cached indexing to enhance performance, so a valid
    // cache is needed in this function.
    template<
        bool stretched,
        GeometryCurvaturePolicy curvPol = GeometryCurvaturePolicy::withSign
    > static void updateGeometryValue(MeshType& mesh) {
        using namespace mathfunc;

        cacheIndices(mesh);

        const auto& vertices = mesh.getVertices();
        const auto& halfEdges = mesh.getHalfEdges();
        const auto& triangles = mesh.getTriangles();

        const size_t numVertices = vertices.size();
        const size_t numHalfEdges = halfEdges.size();
        const size_t numTriangles = triangles.size();

        const auto& coords = stretched ? Bead::getDbDataConst().coordsStr : Bead::getDbDataConst().coords;

        // Calculate angles stored in half edges
        for(MeshType::HalfEdgeIndex hei {}; hei < numHalfEdges; ++hei) {
            // The angle is (c0, c1, c2)
            auto& hea = mesh.attribute(hei);
            auto& heag = hea.template getGHalfEdge<stretched>();
            const auto& c0 = coords[hea.cachedCoordIndex[0]];
            const auto& c1 = coords[hea.cachedCoordIndex[1]];
            const auto& c2 = coords[hea.cachedCoordIndex[2]];

            const auto cp = cross(c0 - c1, c2 - c1);
            const auto dp =   dot(c0 - c1, c2 - c1);
            heag.cotTheta = dp / magnitude(cp);
        }

        // Calculate triangle area and cone volume
        for(MeshType::TriangleIndex ti {}; ti < numTriangles; ++ti) {
            auto& ta = mesh.attribute(ti);
            auto& tag = ta.template getGTriangle<stretched>();
            const auto& c0 = coords[ta.cachedCoordIndex[0]];
            const auto& c1 = coords[ta.cachedCoordIndex[1]];
            const auto& c2 = coords[ta.cachedCoordIndex[2]];

            const auto cp = cross(c1 - c0, c2 - c0);

            // area
            tag.area = magnitude(cp) * 0.5;

            // cone volume
            tag.coneVolume = dot(c0, cp) / 6;
        }

        const auto& cvt = mesh.metaAttribute().cachedVertexTopo;

        // Calculate vertex 1-ring area and local curvature
        for(MeshType::VertexIndex vi {}; vi < numVertices; ++vi) if(!mesh.isVertexOnBorder(vi)) {
            auto& va = mesh.attribute(vi);
            auto& vag = va.template getGVertex<stretched>();
            const coordinate_type ci (coords[va.cachedCoordIndex]);

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
                const size_t ti0 = cvt[mesh.metaAttribute().cachedVertexOffsetPolygon(vi) + i];
                const size_t hei_n = cvt[mesh.metaAttribute().cachedVertexOffsetLeavingHE(vi) + (i + va.cachedDegree - 1) % va.cachedDegree];
                const size_t hei_on = cvt[mesh.metaAttribute().cachedVertexOffsetOuterHE(vi) + (i + 1) % va.cachedDegree];
                const coordinate_type cn      (coords[cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi) + i]]);
                const coordinate_type c_right (coords[cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi) + (i + 1) % va.cachedDegree]]);

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

            if(curvPol == GeometryCurvaturePolicy::withSign) {
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
    template<
        GeometryCurvaturePolicy curvPol = GeometryCurvaturePolicy::withSign
    > static void updateGeometryValueWithDerivative(MeshType& mesh) {
        using namespace mathfunc;

        cacheIndices(mesh);

        const auto& vertices = mesh.getVertices();
        const auto& triangles = mesh.getTriangles();

        const size_t numVertices = vertices.size();
        const size_t numTriangles = triangles.size();

        const auto& coords = Bead::getDbData().coords;

        // Calculate angles and triangle areas with derivative
        // Calculate triangle cone volumes
        for(MeshType::TriangleIndex ti {}; ti < numTriangles; ++ti) {
            auto& ta = mesh.attribute(ti);

            const auto& hei = ta.cachedHalfEdgeIndex;
            auto& tag = mesh.attribute(ti).gTriangle;

            const coordinate_type c[] {
                static_cast<coordinate_type>(coords[ta.cachedCoordIndex[0]]),
                static_cast<coordinate_type>(coords[ta.cachedCoordIndex[1]]),
                static_cast<coordinate_type>(coords[ta.cachedCoordIndex[2]])
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
        for(MeshType::VertexIndex vi {}; vi < numVertices; ++vi) if(!mesh.isVertexOnBorder(vi)) {
            auto& va = mesh.attribute(vi);
            auto& vag = va.gVertex;
            const coordinate_type ci (coords[va.cachedCoordIndex]);

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
                const size_t hei_o = cvt[mesh.metaAttribute().cachedVertexOffsetLeavingHE(vi) + i];
                const size_t ti0 = cvt[mesh.metaAttribute().cachedVertexOffsetPolygon(vi) + i];
                const size_t hei_n = cvt[mesh.metaAttribute().cachedVertexOffsetLeavingHE(vi) + (i + va.cachedDegree - 1) % va.cachedDegree];
                const size_t hei_p = cvt[mesh.metaAttribute().cachedVertexOffsetOuterHE(vi) + i];
                const size_t hei_on = cvt[mesh.metaAttribute().cachedVertexOffsetOuterHE(vi) + (i + 1) % va.cachedDegree];
                const coordinate_type cn      (coords[cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi) + i]]);
                const coordinate_type c_right (coords[cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi) + (i + 1) % va.cachedDegree]]);

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

            if(curvPol == GeometryCurvaturePolicy::withSign) {
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
                const size_t hei_o = cvt[mesh.metaAttribute().cachedVertexOffsetLeavingHE(vi) + i];
                const size_t hei_n = cvt[mesh.metaAttribute().cachedVertexOffsetLeavingHE(vi) + (i + va.cachedDegree - 1) % va.cachedDegree];
                const size_t hei_p = cvt[mesh.metaAttribute().cachedVertexOffsetOuterHE(vi) + i];
                const size_t hei_on = cvt[mesh.metaAttribute().cachedVertexOffsetOuterHE(vi) + (i + 1) % va.cachedDegree];
                const coordinate_type cn      (coords[cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi) + i]]);
                const coordinate_type c_left  (coords[cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi) + (i + va.cachedDegree - 1) % va.cachedDegree]]);
                const coordinate_type c_right (coords[cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(vi) + (i + 1) % va.cachedDegree]]);

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

                if(curvPol == GeometryCurvaturePolicy::withSign) {
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

            if(curvPol == GeometryCurvaturePolicy::withSign) {
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
    static void updateGeometryValueForSystem(MeshType& mesh) {
        using namespace mathfunc;

        cacheIndices(mesh);

        const auto& vertices = mesh.getVertices();
        const auto& edges = mesh.getEdges();
        const auto& triangles = mesh.getTriangles();

        const size_t numVertices = vertices.size();
        const size_t numEdges = edges.size();
        const size_t numTriangles = triangles.size();

        const auto& coords = Bead::getDbDataConst().coords;

        // Calculate triangle unit normal and area
        for(MeshType::TriangleIndex ti {}; ti < numTriangles; ++ti) {
            auto& ta = mesh.attribute(ti);
            auto& tag = ta.gTriangle;
            const auto& c0 = coords[ta.cachedCoordIndex[0]];
            const auto& c1 = coords[ta.cachedCoordIndex[1]];
            const auto& c2 = coords[ta.cachedCoordIndex[2]];

            const auto cp = cross(c1 - c0, c2 - c0);

            // unit normal
            tag.unitNormal = normalizedVector(cp);

            // area
            tag.area = magnitude(cp) * 0.5;
        }

        // Calculate edge pesudo unit normal
        for(MeshType::EdgeIndex ei {}; ei < numEdges; ++ei) {
            auto& ea = mesh.attribute(ei);

            // pseudo unit normal
            using PolygonType = typename MeshType::HalfEdge::PolygonType;
            if(
                ea.cachedPolygonType[0] == PolygonType::triangle &&
                ea.cachedPolygonType[1] == PolygonType::triangle
            ) {
                ea.gEdge.pseudoUnitNormal = normalizedVector(
                    mesh.attribute(ea.cachedPolygonIndex[0]).gTriangle.unitNormal +
                    mesh.attribute(ea.cachedPolygonIndex[1]).gTriangle.unitNormal
                );
            }
            else if(ea.cachedPolygonType[0] == PolygonType::triangle) {
                ea.gEdge.pseudoUnitNormal = mesh.attribute(ea.cachedPolygonIndex[0]).gTriangle.unitNormal;
            }
            else if(ea.cachedPolygonType[1] == PolygonType::triangle) {
                ea.gEdge.pseudoUnitNormal = mesh.attribute(ea.cachedPolygonIndex[1]).gTriangle.unitNormal;
            }
        }

        const auto& cvt = mesh.metaAttribute().cachedVertexTopo;

        // Calculate vertex pseudo unit normal
        for(MeshType::VertexIndex vi {}; vi < numVertices; ++vi) {
            auto& va = mesh.attribute(vi);
            auto& vag = va.gVertex;

            // clearing
            vag.pseudoUnitNormal = {0.0, 0.0, 0.0};

            for(size_t i = 0; i < va.cachedDegree; ++i) {
                const auto hei = cvt[mesh.metaAttribute().cachedVertexOffsetTargetingHE(vi) + i];
                if(mesh.getHalfEdges()[hei].polygonType == MeshType::HalfEdge::PolygonType::triangle) {
                    const auto ti0 = cvt[mesh.metaAttribute().cachedVertexOffsetPolygon(vi) + i];

                    const auto theta = mesh.attribute(hei).gHalfEdge.theta;

                    vag.pseudoUnitNormal += theta * mesh.attribute(ti0).gTriangle.unitNormal;
                }
            }

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
    static double signedDistance(const MeshType& mesh, const VecType& p) {
        using namespace mathfunc;

        const size_t numTriangles = mesh.getTriangles().size();

        double minAbsDistance = numeric_limits<double>::infinity();
        for(MeshType::TriangleIndex ti {}; ti < numTriangles; ++ti) {
            /**********************************************************************
            Calculate the barycentric coordinate of the projection point p'

            See Heidrich 2005, Computing the Barycentric Coordinates of a Projected
            Point.
            **********************************************************************/
            const auto hei0 = mesh.getTriangles()[ti].halfEdgeIndex;
            const auto hei1 = mesh.next(hei0);
            const auto hei2 = mesh.next(hei1);
            const MeshType::VertexIndex vi[] {
                mesh.target(hei0), mesh.target(hei1), mesh.target(hei2)
            };
            const coordinate_type c[] {
                static_cast<coordinate_type>(mesh.attribute(vi[0]).vertex->coord),
                static_cast<coordinate_type>(mesh.attribute(vi[1]).vertex->coord),
                static_cast<coordinate_type>(mesh.attribute(vi[2]).vertex->coord)
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
                const Vec< 3, typename coordinate_type::float_type > r2 {
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
    static bool contains(const MeshType& mesh, const VecType& p) {
        return signedDistance(mesh, p) < 0.0;
    }

    // Attribute computation in adaptive remeshing algorithms
    //-------------------------------------------------------------------------

    // Triangle normal (Geometric attribute) used in adaptive remeshing
    static void adaptiveComputeTriangleNormal(MeshType& mesh, MeshType::TriangleIndex ti) {
        const auto hei = mesh.getTriangles()[ti].halfEdgeIndex;
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
    static void adaptiveComputeAngle(MeshType& mesh, MeshType::HalfEdgeIndex hei) {
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
    static void adaptiveComputeVertexNormal(MeshType& mesh, MeshType::VertexIndex vi) {
        // Using pseudo normal (weighted by angles)
        auto& vaa = mesh.attribute(vi).aVertex;

        // clearing
        vaa.unitNormal = {0.0, 0.0, 0.0};

        mesh.forEachHalfEdgeTargetingVertex(vi, [&](MeshType::HalfEdgeIndex hei) {
            if(mesh.polygonType(hei) == MeshType::HalfEdge::PolygonType::triangle) {
                const size_t ti0 = mesh.triangle(hei);
                const auto theta = mesh.attribute(hei).gHalfEdge.theta;
                vaa.unitNormal += theta * mesh.attribute(ti0).gTriangle.unitNormal;
            }
        });

        mathfunc::normalize(vaa.unitNormal);
    }

};

#endif
