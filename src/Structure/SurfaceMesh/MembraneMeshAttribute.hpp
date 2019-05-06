#ifndef MEDYAN_Structure_SurfaceMesh_MembraneMeshAttribute_hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneMeshAttribute_hpp

#include <array>
#include <limits> // numeric_limits
#include <memory> // unique_ptr
#include <stdexcept> // logic_error
#include <vector>

#include "MathFunctions.h"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/AdaptiveMeshAttribute.hpp"
#include "Structure/SurfaceMesh/Edge.h"
#include "Structure/SurfaceMesh/GeometricMeshAttribute.hpp"
#include "Structure/SurfaceMesh/Triangle.h"
#include "Structure/SurfaceMesh/Vertex.h"
#include "util/io/log.h"

/******************************************************************************
Implements the attributes of the meshwork used by the membrane, mainly
geometric attribute and adaptive mesh attributes.

Geometric attributes are mainly used in force field computations.
Adaptive attributes provides additional attributes mainly used in adaptive mesh
algorithm.

This struct can connect the topological meshwork with the actual system objects.
******************************************************************************/
template< template< typename > class MeshTopology > // C++17: class -> typename
struct MembraneMeshAttribute {

    using MeshType = MeshTopology< MembraneMeshAttribute >;

    struct VertexAttribute {
        using coordinate_type      = Vertex::coordinate_type;
        using coordinate_ref_type  = Vertex::coordinate_ref_type;
        using coordinate_cref_type = Vertex::coordinate_cref_type;

        Vertex* vertex;

        GVertex gVertex;
        GVertex gVertexS; // stretched version (temporary)
        AdaptiveMeshAttribute::VertexAttribute aVertex;

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
    };

    using coordinate_type = typename VertexAttribute::coordinate_type;

    struct AttributeInitializerInfo {
        std::vector< coordinate_type > vertexCoordinateList;
    };

    // Mesh element modification (not used in initialization/finalization)
    template< typename VertexInsertionMethod >
    static void newVertex(MeshType& mesh, size_t v, const VertexInsertionMethod& op) {
        const auto& meta = mesh.getMetaAttribute();
        mesh.getVertexAttribute(v).vertex = meta.s->template addTrackable<Vertex>(op.coordinate(mesh, v), meta.m, v);
    }
    template< typename Operation >
    static void newEdge(MeshType& mesh, size_t e, const Operation& op) {
        mesh.getEdgeAttribute(e).edge.reset(mesh.getMetaAttribute().s->template addTrackable<Edge>(mesh.getMetaAttribute().m, e));
    }
    template< typename Operation >
    static void newHalfEdge(MeshType& mesh, size_t he, const Operation& op) {
        // Do nothing
    }
    template< typename Operation >
    static void newTriangle(MeshType& mesh, size_t t, const Operation& op) {
        mesh.getTriangleAttribute(t).triangle.reset(mesh.getMetaAttribute().s->template addTrackable<Triangle>(mesh.getMetaAttribute().m, t));
    }
    template< typename Operation >
    static void newBorder(MeshType& mesh, size_t, const Operation& op) {
        // Do nothing
    }

    template< typename Element, std::enable_if_t<std::is_same<Element, typename MeshType::Vertex>::value, void>* = nullptr >
    static void removeElement(MeshType& mesh, size_t i) {
        mesh.getMetaAttribute().s->template removeTrackable<Vertex>(mesh.getVertexAttribute(i).vertex);
        mesh.getMetaAttribute().m->removeChild(mesh.getVertexAttribute(i).vertex);
    }
    template< typename Element, std::enable_if_t<std::is_same<Element, typename MeshType::Edge>::value, void>* = nullptr >
    static void removeElement(MeshType& mesh, size_t i) {
        mesh.getMetaAttribute().s->template removeTrackable<Edge>(mesh.getEdgeAttribute(i).edge.get());
    }
    template< typename Element, std::enable_if_t<std::is_same<Element, typename MeshType::HalfEdge>::value, void>* = nullptr >
    static void removeElement(MeshType& mesh, size_t i) {
        // Do nothing
    }
    template< typename Element, std::enable_if_t<std::is_same<Element, typename MeshType::Triangle>::value, void>* = nullptr >
    static void removeElement(MeshType& mesh, size_t i) {
        mesh.getMetaAttribute().s->template removeTrackable<Triangle>(mesh.getTriangleAttribute(i).triangle.get());
    }
    template< typename Element, std::enable_if_t<std::is_same<Element, typename MeshType::Border>::value, void>* = nullptr >
    static void removeElement(MeshType& mesh, size_t i) {
        // Do nothing
    }

    // Mesh index caching
    template< bool forceUpdate = false >
    static void cacheIndices(MeshType& mesh) {

        if(forceUpdate || !mesh.getMetaAttribute().cacheValid) {
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
                hea.cachedCoordIndex[0] = vertices[vi0].attr.vertex->Bead::getIndex();
                hea.cachedCoordIndex[1] = vertices[vi1].attr.vertex->Bead::getIndex();
                hea.cachedCoordIndex[2] = vertices[vi2].attr.vertex->Bead::getIndex();
            }

            mesh.getMetaAttribute().cacheValid = true;
        }

    } // void cacheIndices(...)

    // Mesh attribute initializing and extracting
    // These operations do not follow the RAII idiom.
    // Initialization should happen only once, as it allocates resources
    static void init(MeshType& mesh, const AttributeInitializerInfo& info) {
        const MetaAttribute& meta = mesh.getMetaAttribute();
        for(size_t i = 0; i < mesh.getVertices().size(); ++i) {
            mesh.getVertexAttribute(i).vertex = meta.s->template addTrackable<Vertex>(info.vertexCoordinateList[i], meta.m, i);
        }
        for(size_t i = 0; i < mesh.getEdges().size(); ++i) {
            mesh.getEdgeAttribute(i).edge.reset(meta.s->template addTrackable<Edge>(meta.m, i));
        }
        for(size_t i = 0; i < mesh.getTriangles().size(); ++i) {
            mesh.getTriangleAttribute(i).triangle.reset(meta.s->template addTrackable<Triangle>(meta.m, i));
        }
    }
    // Extraction can be done multiple times without allocating/deallocating
    static auto extract(const MeshType& mesh) {
        AttributeInitializerInfo info;

        const size_t numVertices = mesh.getVertices().size();

        info.vertexCoordinateList.reserve(numVertices);
        for(size_t i = 0; i < numVertices; ++i) {
            info.vertexCoordinateList.push_back(mesh.getVertexAttribute(i).vertex->coordinate());
        }

        return info;
    }

    // Geometries
    template< bool stretched > static void updateGeometryValue(MeshType& mesh) {
        using namespace mathfunc;

        if(!mesh.getMetaAttribute().cacheValid) {
            LOG(ERROR) << "Updating membrane mesh geometry values without a valid index cache."
            throw std::runtime_error("Membrane mesh index cache invalid");
        }

        const auto& vertices = mesh.getVertices();
        const auto& halfEdges = mesh.getHalfEdges();
        const auto& edges = mesh.getEdges();
        const auto& triangles = mesh.getTriangles();

        const size_t numVertices = vertices.size();
        const size_t numHalfEdges = halfEdges.size();
        const size_t numEdges = edges.size();
        const size_t numTriangles = triangles.size();

        // Calculate angles stored in half edges
        for(size_t hei = 0; hei < numHalfEdges; ++hei) {
            // The angle is (c0, c1, c2)
            auto& hea = mesh.getHalfEdgeAttribute(hei);
            auto& heag = mesh.getHalfEdgeAttribute(hei).template getGHalfEdge<stretched>();
            const auto& c0 = Bead::getDbData().coords[hea.cachedCoordIndex[0]];
            const auto& c1 = Bead::getDbData().coords[hea.cachedCoordIndex[1]];
            const auto& c2 = Bead::getDbData().coords[hea.cachedCoordIndex[2]];

            const auto cp = cross(c0 - c1, c2 - c1);
            const auto dp =   dot(c0 - c1, c2 - c1);
            const auto ct = heag.cotTheta = dp / magnitude(cp);
            heag.theta = M_PI_2 - atan(ct);
        }

        // Calculate triangle area, unit normal and cone volume
        for(size_t ti = 0; ti < numTriangles; ++ti) {
            const size_t hei = triangles[ti].halfEdgeIndex;
            const size_t vi0 = mesh.target(hei);
            const size_t vi1 = mesh.target(mesh.next(hei));
            const size_t vi2 = mesh.target(mesh.prev(hei));
            const auto& c0 = vertices[vi0].attr.vertex->template getCoordinate<stretched>();
            const auto& c1 = vertices[vi1].attr.vertex->template getCoordinate<stretched>();
            const auto& c2 = vertices[vi2].attr.vertex->template getCoordinate<stretched>();
            auto& tag = mesh.getTriangleAttribute(ti).template getGTriangle<stretched>();

            const auto cp = cross(c1 - c0, c2 - c0);

            // area
            tag.area = magnitude(cp) * 0.5;

            // unit normal
            tag.unitNormal = normalizedVector(cp);

            // cone volume
            tag.coneVolume = dot(c0, cp) / 6;
        }

        // Calculate edge length and pesudo unit normal
        for(size_t ei = 0; ei < numEdges; ++ei) {
            const size_t hei = edges[ei].halfEdgeIndex;
            const size_t hei_o = mesh.opposite(hei);
            const size_t vi0 = mesh.target(hei);
            const size_t vi1 = mesh.target(mesh.prev(hei));

            // length
            mesh.getEdgeAttribute(ei).template getGEdge<stretched>().length = distance(
                vertices[vi0].attr.vertex->template getCoordinate<stretched>(),
                vertices[vi1].attr.vertex->template getCoordinate<stretched>()
            );

            // pseudo unit normal
            using PolygonType = typename MeshType::HalfEdge::PolygonType;
            if(
                halfEdges[hei].polygonType   == PolygonType::Triangle &&
                halfEdges[hei_o].polygonType == PolygonType::Triangle
            ) {
                const size_t ti0 = mesh.triangle(hei);
                const size_t ti1 = mesh.triangle(mesh.opposite(hei));
                mesh.getEdgeAttribute(ei).template getGEdge<stretched>().pseudoUnitNormal = normalizedVector(
                    mesh.getTriangleAttribute(ti0).template getGTriangle<stretched>().unitNormal
                    + mesh.getTriangleAttribute(ti1).template getGTriangle<stretched>().unitNormal
                );
            }
        }

        // Calculate vcell area, curvature and vertex pseudo unit normal
        for(size_t vi = 0; vi < numVertices; ++vi) if(!mesh.isVertexOnBorder(vi)) {
            auto& vag = mesh.getVertexAttribute(vi).template getGVertex<stretched>();

            // clearing
            vag.area = 0.0;
            vag.pseudoUnitNormal = {0.0, 0.0, 0.0};

            // k1 = 2A * k, where k is the result of LB operator
            Vec3 k1 {};

            mesh.forEachHalfEdgeTargetingVertex(vi, [&mesh, &vertices, &vi, &vag, &k1](size_t hei) {
                const size_t hei_o = mesh.opposite(hei);
                const size_t ti0 = mesh.triangle(hei);
                const size_t vn = mesh.target(hei_o);
                const size_t hei_n = mesh.next(hei);
                const size_t hei_on = mesh.next(hei_o);
                const Vec3 ci = vertices[vi].attr.vertex->template getCoordinate<stretched>();
                const Vec3 cn = vertices[vn].attr.vertex->template getCoordinate<stretched>();

                const auto sumCotTheta =
                    mesh.getHalfEdgeAttribute(hei_n).template getGHalfEdge<stretched>().cotTheta
                    + mesh.getHalfEdgeAttribute(hei_on).template getGHalfEdge<stretched>().cotTheta;

                const auto theta = mesh.getHalfEdgeAttribute(hei).template getGHalfEdge<stretched>().theta;

                const auto diff = ci - cn;
                const auto dist2 = magnitude2(diff);

                vag.area += sumCotTheta * dist2 * 0.125;

                k1 += sumCotTheta * diff;
                vag.pseudoUnitNormal += theta * mesh.getTriangleAttribute(ti0).template getGTriangle<stretched>().unitNormal;
            });

            const double invA = 1 / vag.area;
            const double magK1 = magnitude(k1);

            normalize(vag.pseudoUnitNormal);

            const int flippingCurv = (dot(k1, vag.pseudoUnitNormal) > 0 ? 1 : -1);

            vag.curv = flippingCurv * magK1 * 0.25 * invA;
        }
    }
    static void updateGeometryValueWithDerivative(MeshType& mesh) {
        using namespace mathfunc;

        const auto& vertices = mesh.getVertices();
        const auto& halfEdges = mesh.getHalfEdges();
        const auto& edges = mesh.getEdges();
        const auto& triangles = mesh.getTriangles();

        const size_t numVertices = vertices.size();
        const size_t numHalfEdges = halfEdges.size();
        const size_t numEdges = edges.size();
        const size_t numTriangles = triangles.size();

        // Calculate edge length with deriviative
        for(size_t ei = 0; ei < numEdges; ++ei) {
            // The edge must have an opposite
            const size_t hei = edges[ei].halfEdgeIndex;
            const size_t hei_o = mesh.opposite(hei);
            const size_t vi0 = mesh.target(hei);
            const size_t vi1 = mesh.target(hei_o);
            const Vec3 c0 = vertices[vi0].attr.vertex->coordinate();
            const Vec3 c1 = vertices[vi1].attr.vertex->coordinate();

            const auto length = mesh.getEdgeAttribute(ei).gEdge.length = distance(c0, c1);
            const auto invL = 1.0 / length;
            mesh.getHalfEdgeAttribute(hei).gHalfEdge.dEdgeLength = (c0 - c1) * invL;
            mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dEdgeLength = (c1 - c0) * invL;
        }

        // Calculate angles and triangle areas with derivative
        // Calculate triangle normals and cone volumes
        for(size_t ti = 0; ti < numTriangles; ++ti) {
            size_t hei[3];
            hei[0] = triangles[ti].halfEdgeIndex;
            hei[1] = mesh.next(hei[0]);
            hei[2] = mesh.next(hei[1]);
            auto& tag = mesh.getTriangleAttribute(ti).gTriangle;

            const size_t vi[] {mesh.target(hei[0]), mesh.target(hei[1]), mesh.target(hei[2])};
            const Vec3 c[] {
                vertices[vi[0]].attr.vertex->coordinate(),
                vertices[vi[1]].attr.vertex->coordinate(),
                vertices[vi[2]].attr.vertex->coordinate()
            };

            const double l[] {
                edges[mesh.edge(hei[0])].attr.gEdge.length,
                edges[mesh.edge(hei[1])].attr.gEdge.length,
                edges[mesh.edge(hei[2])].attr.gEdge.length
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
                mesh.getHalfEdgeAttribute(hei[0]).gHalfEdge.dTriangleArea = (-l[1]*l[1]* r02 - l[0]*l[0]* r01 + dots[0]*(r01 + r02)) * (invA * 0.25);
                mesh.getHalfEdgeAttribute(hei[1]).gHalfEdge.dTriangleArea = (l[0]*l[0]* r01 - dots[0]* r02) * (invA * 0.25);
                mesh.getHalfEdgeAttribute(hei[2]).gHalfEdge.dTriangleArea = (l[1]*l[1]* r02 - dots[0]* r01) * (invA * 0.25);
            }

            // Calculate thetas and gradients
            for(size_t ai = 0; ai < 3; ++ai) {
                auto& heag = mesh.getHalfEdgeAttribute(hei[ai]).gHalfEdge;

                const auto ct = heag.cotTheta = dots[ai] * invA * 0.5;
                heag.theta = M_PI_2 - atan(ct);

                const size_t ai_n = (ai + 1) % 3;
                const size_t ai_p = (ai + 2) % 3;
                auto& heag_n = mesh.getHalfEdgeAttribute(hei[ai_n]).gHalfEdge;
                auto& heag_p = mesh.getHalfEdgeAttribute(hei[ai_p]).gHalfEdge;

                const auto r01 = c[ai_n] - c[ai];
                const auto r02 = c[ai_p] - c[ai];

                heag.dCotTheta[1] =
                    -(r01 + r02) * (invA * 0.5)
                    -(dots[ai] * invA * invA * 0.5) * heag.dTriangleArea;
                heag.dCotTheta[2] = r02 * (invA * 0.5) - (dots[ai] * invA * invA * 0.5) * heag_n.dTriangleArea;
                heag.dCotTheta[0] = r01 * (invA * 0.5) - (dots[ai] * invA * invA * 0.5) * heag_p.dTriangleArea;
            }

            // Calculate unit normal
            tag.unitNormal = normalizedVector(cp);

            // Calculate cone volume and derivative
            tag.coneVolume = dot(c[0], r0) / 6.0;
            // The derivative of cone volume will be accumulated to each vertex

        }

        // Clear derivative of vcell area on neighbors
        for(size_t hei = 0; hei < numHalfEdges; ++hei) {
            mesh.getHalfEdgeAttribute(hei).gHalfEdge.dNeighborArea = {0.0, 0.0, 0.0};
        }

        // Calculate vcell area, curvature with derivative
        // Calculate vertex pseudo unit normal
        // Calculate derivative of volume on vertices
        for(size_t vi = 0; vi < numVertices; ++vi) if(!mesh.isVertexOnBorder(vi)) {
            auto& vag = mesh.getVertexAttribute(vi).gVertex;

            // clearing
            vag.area = 0.0;
            vag.dArea = {0.0, 0.0, 0.0};
            vag.pseudoUnitNormal = {0.0, 0.0, 0.0};
            vag.dVolume = {0.0, 0.0, 0.0};

            // K = 2*H*n is the result of LB operator
            // And let k1 = 2*A*k (as an intermediate variable)
            Vec3 k1 {};

            // derivative of k1 and curvature will be calculated in the next loop
            mesh.forEachHalfEdgeTargetingVertex(vi, [&mesh, &vertices, vi, &vag, &k1](size_t hei) {
                const size_t hei_o = mesh.opposite(hei);
                const size_t ti0 = mesh.triangle(hei);
                const size_t vn = mesh.target(hei_o);
                const size_t hei_n = mesh.next(hei);
                const size_t hei_on = mesh.next(hei_o);
                const size_t hei_right = mesh.opposite(mesh.next(hei_on)); // hei_left is hei_n
                const size_t vi_right = mesh.target(hei_right);
                const Vec3 ci = vertices[vi].attr.vertex->coordinate();
                const Vec3 cn = vertices[vn].attr.vertex->coordinate();
                const Vec3 c_right = vertices[vi_right].attr.vertex->coordinate();

                const auto sumCotTheta = mesh.getHalfEdgeAttribute(hei_n).gHalfEdge.cotTheta + mesh.getHalfEdgeAttribute(hei_on).gHalfEdge.cotTheta;
                const auto& dCotThetaLeft = mesh.getHalfEdgeAttribute(hei_n).gHalfEdge.dCotTheta;
                const auto& dCotThetaRight = mesh.getHalfEdgeAttribute(hei_on).gHalfEdge.dCotTheta;

                const auto theta = mesh.getHalfEdgeAttribute(hei).gHalfEdge.theta;

                const auto diff = ci - cn;
                const auto dist2 = magnitude2(diff);

                vag.area += sumCotTheta * dist2 * 0.125;

                // Area derivative
                vag.dArea +=
                    (dCotThetaLeft[0] + dCotThetaRight[2]) * (dist2 * 0.125)
                    + (sumCotTheta * 0.25) * diff; // d(dist2) / dx = 2 * diff
                mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dNeighborArea +=
                    (dCotThetaLeft[2] + dCotThetaRight[0]) * (dist2 * 0.125)
                    - (sumCotTheta * 0.25) * diff; // d(dist2) / dx = -2 * diff
                mesh.getHalfEdgeAttribute(hei_n).gHalfEdge.dNeighborArea +=
                    dCotThetaLeft[1] * (dist2 * 0.125);
                mesh.getHalfEdgeAttribute(hei_right).gHalfEdge.dNeighborArea +=
                    dCotThetaRight[1] * (dist2 * 0.125);

                // Accumulate k1
                k1 += sumCotTheta * diff;

                // Accumulate pseudo unit normal
                vag.pseudoUnitNormal += theta * mesh.getTriangleAttribute(ti0).gTriangle.unitNormal;

                // Added to derivative of sum of cone volume
                const auto cp = cross(cn, c_right);
                vag.dVolume += cp * (1.0 / 6);
            });

            const auto invA = 1.0 / vag.area;
            const auto magK1 = magnitude(k1);

            // Calculate pseudo unit normal
            normalize(vag.pseudoUnitNormal);

            // Calculate mean curvature H = |k1| / 4A
            // dH = (dK1)K1 / 4A|K1| - |K1|dA / 4A^2
            const int flippingCurv = (dot(k1, vag.pseudoUnitNormal) > 0 ? 1 : -1);
            const auto dCurvFac1 = 0.25 * invA * flippingCurv / magK1;
            const auto dCurvFac2 = -0.25 * invA * invA * magK1 * flippingCurv;

            vag.curv = flippingCurv * magK1 * 0.25 * invA;
            // Derivative will be processed later.

            // Calculate derivative of k1 and curvature
            // Using another loop because k1 is needed for curvature derivative
            std::array<Vec3, 3> dK1 {}; // On center vertex, indexed by [k1x, k1y, k1z]
            mesh.forEachHalfEdgeTargetingVertex(vi, [&mesh, &vertices, vi, &k1, &dK1, dCurvFac1, dCurvFac2](size_t hei) {
                const size_t hei_o = mesh.opposite(hei);
                const size_t vn = mesh.target(hei_o);
                const size_t hei_n = mesh.next(hei);
                const size_t hei_on = mesh.next(hei_o);
                const size_t hei_right = mesh.opposite(mesh.next(hei_on)); // hei_left is hei_n
                const size_t vi_left = mesh.target(hei_n);
                const size_t vi_right = mesh.target(hei_right);
                const Vec3 ci = vertices[vi].attr.vertex->coordinate();
                const Vec3 cn = vertices[vn].attr.vertex->coordinate();
                const Vec3 c_left = vertices[vi_left].attr.vertex->coordinate();
                const Vec3 c_right = vertices[vi_right].attr.vertex->coordinate();

                const auto sumCotTheta = mesh.getHalfEdgeAttribute(hei_n).gHalfEdge.cotTheta + mesh.getHalfEdgeAttribute(hei_on).gHalfEdge.cotTheta;
                const auto& dCotThetaLeft = mesh.getHalfEdgeAttribute(hei_n).gHalfEdge.dCotTheta;
                const auto& dCotThetaRight = mesh.getHalfEdgeAttribute(hei_on).gHalfEdge.dCotTheta;
                const auto sumDCotThetaCenter = dCotThetaLeft[0] + dCotThetaRight[2];
                const auto sumDCotThetaNeighbor = dCotThetaLeft[2] + dCotThetaRight[0];

                const auto diff = ci - cn;
                // Accumulate dK1 on the center vertex vi
                dK1[0] += sumDCotThetaCenter[0] * diff;
                dK1[1] += sumDCotThetaCenter[1] * diff;
                dK1[2] += sumDCotThetaCenter[2] * diff;
                dK1[0][0] += sumCotTheta;
                dK1[1][1] += sumCotTheta;
                dK1[2][2] += sumCotTheta; // dK1 += I * sumCotTheta, where I is gradient of diff (identity)

                // Calculate dK1 and derivative of curvature on neighbor vertex vn
                std::array<Vec3, 3> dK1_n {};
                // As direct target
                dK1_n[0] = sumDCotThetaNeighbor[0] * diff;
                dK1_n[1] = sumDCotThetaNeighbor[1] * diff;
                dK1_n[2] = sumDCotThetaNeighbor[2] * diff;
                dK1_n[0][0] -= sumCotTheta;
                dK1_n[1][1] -= sumCotTheta;
                dK1_n[2][2] -= sumCotTheta; // dK1 += (-I) * sumCotTheta

                // As target for left and right
                const auto diff_left = ci - c_left;
                const auto diff_right = ci - c_right;
                const auto& dCotThetaOfLeft = mesh.getHalfEdgeAttribute(mesh.next(hei_n)).gHalfEdge.dCotTheta[1];
                const auto& dCotThetaOfRight = mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dCotTheta[1];
                dK1_n[0] += dCotThetaOfLeft[0] * diff_left;
                dK1_n[1] += dCotThetaOfLeft[1] * diff_left;
                dK1_n[2] += dCotThetaOfLeft[2] * diff_left;
                dK1_n[0] += dCotThetaOfRight[0] * diff_right;
                dK1_n[1] += dCotThetaOfRight[1] * diff_right;
                dK1_n[2] += dCotThetaOfRight[2] * diff_right;

                // Derivative of curvature
                const Vec3 mp {{{
                    dot(dK1_n[0], k1),
                    dot(dK1_n[1], k1),
                    dot(dK1_n[2], k1)
                }}}; // A matrix product dK1_n * k1
                mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dNeighborCurv =
                    dCurvFac1 * mp + dCurvFac2 * mesh.getHalfEdgeAttribute(hei_o).gHalfEdge.dNeighborArea;
            });

            // Also the derivative of curvature on central vertex
            vag.dCurv =
                dCurvFac1 * Vec3{ dot(dK1[0], k1), dot(dK1[1], k1), dot(dK1[2], k1) }
                + dCurvFac2 * vag.dArea;

        } // End loop vertices (V cells)
    }

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
        - The length of edges

    Note: this method only works if the mesh is closed. This must be ensured by
          the caller of the function.
    
    In fact, the signed distance field serves as a good candidate for membrane
    boundary potential. However, this field is not C1-continuous everywhere,
    which is detrimental to conjugate gradient methods.
    **************************************************************************/
    static double signedDistance(const MeshType& mesh, const mathfunc::Vec3& p) {
        using namespace mathfunc;

        const size_t numTriangles = mesh.getTriangles().size();

        double minAbsDistance = numeric_limits<double>::infinity();
        for(size_t ti = 0; ti < numTriangles; ++ti) {
            /**********************************************************************
            Calculate the barycentric coordinate of the projection point p'

            See Heidrich 2005, Computing the Barycentric Coordinates of a Projected
            Point.
            **********************************************************************/
            const size_t hei0 = mesh.getTriangles()[ti].halfEdgeIndex;
            const size_t hei1 = mesh.next(hei0);
            const size_t hei2 = mesh.next(hei1);
            const size_t vi[] {
                mesh.target(hei0), mesh.target(hei1), mesh.target(hei2)
            };
            const Vec3 c[] {
                mesh.getVertexAttribute(vi[0]).vertex->coordinate(),
                mesh.getVertexAttribute(vi[1]).vertex->coordinate(),
                mesh.getVertexAttribute(vi[2]).vertex->coordinate()
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
                d = dot(mesh.getTriangleAttribute(ti).gTriangle.unitNormal, r0p);
            } else {
                // p' is outside the triangle
                const Vec3 r {
                    mesh.getEdgeAttribute(mesh.edge(hei2)).gEdge.length, // 1->2
                    mesh.getEdgeAttribute(mesh.edge(hei0)).gEdge.length, // 2->0
                    mesh.getEdgeAttribute(mesh.edge(hei1)).gEdge.length  // 0->1
                };
                const auto r1p = p - c[1];
                const auto r2p = p - c[2];
                const auto r12 = c[2] - c[1];
                const auto dot_1p_12 = dot(r1p, r12);
                const auto dot_2p_20 = -dot(r2p, r02);
                const auto dot_0p_01 = dot(r0p, r01);

                if(b0 < 0 && dot_1p_12 >= 0 && dot_1p_12 <= r[0]*r[0]) {
                    // On edge 12
                    d = magnitude(cross(r1p, r12)) / r[0];
                    if(dot(mesh.getEdgeAttribute(mesh.edge(hei2)).gEdge.pseudoUnitNormal, r1p) < 0) d = -d;
                } else if(b1 < 0 && dot_2p_20 >= 0 && dot_2p_20 <= r[1]*r[1]) {
                    // On edge 20
                    d = magnitude(cross(r2p, r02)) / r[1];
                    if(dot(mesh.getEdgeAttribute(mesh.edge(hei0)).gEdge.pseudoUnitNormal, r2p) < 0) d = -d;
                } else if(b2 < 0 && dot_0p_01 >= 0 && dot_0p_01 <= r[2]*r[2]) {
                    // On edge 01
                    d = magnitude(cross(r0p, r01)) / r[2];
                    if(dot(mesh.getEdgeAttribute(mesh.edge(hei1)).gEdge.pseudoUnitNormal, r0p) < 0) d = -d;
                } else if(dot_0p_01 < 0 && dot_2p_20 > r[1]*r[1]) {
                    // On vertex 0
                    d = distance(c[0], p);
                    if(dot(mesh.getVertexAttribute(vi[0]).gVertex.pseudoUnitNormal, r0p) < 0) d = -d;
                } else if(dot_1p_12 < 0 && dot_0p_01 > r[2]*r[2]) {
                    // On vertex 1
                    d = distance(c[1], p);
                    if(dot(mesh.getVertexAttribute(vi[1]).gVertex.pseudoUnitNormal, r1p) < 0) d = -d;
                } else if(dot_2p_20 < 0 && dot_1p_12 > r[0]*r[0]) {
                    // On vertex 2
                    d = distance(c[2], p);
                    if(dot(mesh.getVertexAttribute(vi[2]).gVertex.pseudoUnitNormal, r2p) < 0) d = -d;
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
    static bool contains(const MeshType& mesh, const mathfunc::Vec3& p) {
        return signedDistance(mesh, p) < 0.0;
    }

    // Attribute computation in adaptive remeshing algorithms

    // Triangle normal (Geometric attribute) used in adaptive remeshing
    static void adaptiveComputeTriangleNormal(MeshType& mesh, size_t ti) {
        const size_t hei = mesh.getTriangles()[ti].halfEdgeIndex;
        const size_t vi0 = mesh.target(hei);
        const size_t vi1 = mesh.target(mesh.next(hei));
        const size_t vi2 = mesh.target(mesh.prev(hei));
        const auto c0 = mesh.getVertexAttribute(vi0).vertex->coordinate();
        const auto c1 = mesh.getVertexAttribute(vi1).vertex->coordinate();
        const auto c2 = mesh.getVertexAttribute(vi2).vertex->coordinate();
        auto& tag = mesh.getTriangleAttribute(ti).gTriangle;

        const auto cp = mathfunc::cross(c1 - c0, c2 - c0);

        // unit normal
        tag.unitNormal = mathfunc::normalizedVector(cp);
    }

    // Triangle angles (Geometric attribute, halfedge) used in adaptive remeshing
    static void adaptiveComputeAngle(MeshType& mesh, size_t hei) {
        // The angle is (v0, v1, v2)
        const size_t vi0 = mesh.target(mesh.prev(hei));
        const size_t vi1 = mesh.target(hei);
        const size_t vi2 = mesh.target(mesh.next(hei));
        const auto c0 = mesh.getVertexAttribute(vi0).vertex->coordinate();
        const auto c1 = mesh.getVertexAttribute(vi1).vertex->coordinate();
        const auto c2 = mesh.getVertexAttribute(vi2).vertex->coordinate();
        auto& heag = mesh.getHalfEdgeAttribute(hei).gHalfEdge;

        const auto cp = mathfunc::cross(c0 - c1, c2 - c1);
        const auto dp = mathfunc::  dot(c0 - c1, c2 - c1);
        const auto ct = heag.cotTheta = dp / mathfunc::magnitude(cp);
        heag.theta = M_PI_2 - std::atan(ct);
    }

    // Vertex unit normals (Adaptive attribute) used in adaptive remeshing
    // Requires
    //   - The vertex is not on the border
    //   - Unit normals in triangles (geometric)
    //   - Angles in halfedges (geometric)
    static void adaptiveComputeVertexNormal(MeshType& mesh, size_t vi) {
        // Using pseudo normal (weighted by angles)
        auto& vaa = mesh.getVertexAttribute(vi).aVertex;

        // clearing
        vaa.unitNormal = {0.0, 0.0, 0.0};

        mesh.forEachHalfEdgeTargetingVertex(vi, [&](size_t hei) {
            const size_t ti0 = mesh.triangle(hei);
            const auto theta = mesh.getHalfEdgeAttribute(hei).gHalfEdge.theta;
            vaa.unitNormal += theta * mesh.getTriangleAttribute(ti0).gTriangle.unitNormal;
        });

        mathfunc::normalize(vaa.unitNormal);
    }

};

#endif
