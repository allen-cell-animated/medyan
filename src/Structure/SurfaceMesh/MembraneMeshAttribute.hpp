#ifndef MEDYAN_Structure_SurfaceMesh_MembraneMeshAttribute_hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneMeshAttribute_hpp

#include <algorithm> // max
#include <array>
#include <limits> // numeric_limits
#include <memory> // unique_ptr
#include <stdexcept> // logic_error
#include <tuple>
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

    // If vertices represent fixed coordinates in a material coordinate system,
    // a vertex represents the location of a specific lipid molecule.
    // The local area elasticity is usually a necessity in this case.
    // An exception is the vertices on the border connected to a lipid
    // reservoir, where the vertices stand for the boundary locations which are
    // usually fixed in the ambient space.
    //
    // If vertices represent fixed coordinates in a normal coordinate system,
    // a vertex is only a representative point on the surface, where the motion
    // of the vertex must be in the local normal direction.
    // Local area elasticity cannot be directly defined on mesh elements.
    enum class VertexSystem { material, normal };

    struct VertexAttribute {
        using CoordinateType      = Vertex::CoordinateType;
        using CoordinateRefType   = Vertex::CoordinateType&;
        using CoordinateCrefType  = const Vertex::CoordinateType&;

        std::unique_ptr< Vertex > vertex;

        GVertex gVertex;
        GVertex gVertexS; // stretched version (temporary)
        AdaptiveMeshAttribute::VertexAttribute aVertex;

        size_t cachedDegree;
        size_t cachedCoordIndex;

        CoordinateCrefType getCoordinate() const { return vertex->coord; }

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

        // The index of x coordinate of vertices in the vectorized dof array
        // [target of half edge, opposite target]
        size_t                          cachedCoordIndex[2];
        MeshType::HalfEdge::PolygonType cachedPolygonType[2];
        // [polygon of half edge, opposite polygon]
        size_t                          cachedPolygonIndex[2];

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

        // The index of x coordinate of vertices in the vectorized dof array
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

        // The index of x coordinate of vertices in the vectorized dof array
        // [target of half edge, next target, prev target]
        size_t                  cachedCoordIndex[3];
        // [half edge, next, prev]
        MeshType::HalfEdgeIndex cachedHalfEdgeIndex[3];
        // [edge of half edge, next edge, prev edge]
        MeshType::EdgeIndex     cachedEdgeIndex[3];

        template< bool stretched > const GTriangle& getGTriangle() const { return stretched ? gTriangleS : gTriangle; }
        template< bool stretched >       GTriangle& getGTriangle()       { return stretched ? gTriangleS : gTriangle; }

        void setIndex(size_t index) {
            triangle->setTopoIndex(index);
        }
    };
    struct BorderAttribute {
        // Whether this border touches a reservoir.
        //
        // If the border is connected to a lipid reservoir, a surface tension
        // would be applied at this border, and the value is supplied as a
        // membrane mechanical parameter.
        // The surface tension must be the same for all borders of a membrane,
        // to make energy minimization possible.
        //
        // If the border is not connected to a lipid reservoir, then it simply
        // represents the free boundary of the membrane.
        bool reservoir = false;

        void setIndex(size_t index) {}
    };
    struct MetaAttribute {
        SubSystem *s;
        Membrane *m;

        VertexSystem vertexSystem = VertexSystem::material;

        // Cache related
        //-----------------------------
        bool indexCacheForFFValid = false;

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

    using CoordinateType = typename VertexAttribute::CoordinateType;

    struct AttributeInitializerInfo {
        std::vector< CoordinateType > vertexCoordinateList;
    };

    //-------------------------------------------------------------------------
    // Basic mesh operations
    //-------------------------------------------------------------------------

    // Mesh element modification (not used in initialization/finalization)
    static void newVertex(MeshType& mesh, MeshType::VertexIndex v) {
        mesh.attribute(v).vertex = std::make_unique< Vertex >(CoordinateType{}, v.index);
    }
    static void newEdge(MeshType& mesh, MeshType::EdgeIndex e) {
        mesh.attribute(e).edge.reset(
            mesh.metaAttribute().s->template addTrackable<Edge>(mesh.metaAttribute().m, e));
    }
    static void newHalfEdge(MeshType& mesh, MeshType::HalfEdgeIndex he) {
        // Do nothing
    }
    static void newTriangle(MeshType& mesh, MeshType::TriangleIndex tp) {
        mesh.attribute(t).triangle.reset(mesh.metaAttribute().s->template addTrackable<Triangle>(mesh.metaAttribute().m, t));
    }
    static void newBorder(MeshType& mesh, MeshType::BorderIndex) {
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
                static_cast<CoordinateType>(mesh.attribute(MeshType:VertexIndex{i}).vertex->coord));
        }

        return info;
    }

    // Mesh modifiers
    //
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
    static auto insertVertexOnEdge(
        MeshType&             mesh,
        MeshType::EdgeIndex   ei,
        const CoordinateType& newPos
    ) {
        const auto ohei       = mesh.halfEdge(edgeIndex);
        const auto ohei_o     = mesh.opposite(ohei);

        const auto opt0       = mesh.polygonType(ohei);
        const bool ist0       = opt0 == HalfEdge::PolygonType::triangle;
        const auto opt2       = mesh.polygonType(ohei_o);
        const bool ist2       = opt2 == HalfEdge::PolygonType::triangle;

        // Do the vertex insertion
        const auto change = MeshType::VertexInsertionOnEdge{}(mesh, ei);

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
        if(opt0) {
            const double a1 = area(mesh, change.tiNew[0]);
            const double a2 = area(mesh, change.tiNew[1]);

            auto& eqArea1 = mesh.attribute(change.tiNew[0]).triangle->mTriangle.eqArea;
            auto& eqArea2 = mesh.attribute(change.tiNew[1]).triangle->mTriangle.eqArea;
            redisEqArea(eqArea1, a1, eqArea2, a2, eqArea1);
        }
        if(opt2) {
            const double a1 = area(mesh, change.tiNew[2]);
            const double a2 = area(mesh, change.tiNew[3]);

            auto& eqArea1 = mesh.attribute(change.tiNew[2]).triangle->mTriangle.eqArea;
            auto& eqArea2 = mesh.attribute(change.tiNew[3]).triangle->mTriangle.eqArea;
            redisEqArea(eqArea1, a1, eqArea2, a2, eqArea1);
        }

        return change;
    }

    static auto collapseEdge(
        MeshType&               mesh,
        MeshType::EdgeIndex     ei,
        const std::optional< CoordinateType >& newPos
    ) {
        const auto hei = mesh.halfEdge(ei);
        const auto hei_o = mesh.opposite(hei);
        const auto vi0 = mesh.target(hei);
        const auto vi1 = mesh.target(hei_o);

        // Record properties
        double oldTotalEqArea = 0.0;
        // Accumulate around two vertices, and subtract triangles on the middle edge
        const auto addOldTotalEqArea = [&](MeshType::HalfEdgeIndex hei) {
            if(mesh.isInTriangle(hei)) {
                oldTotalEqArea += mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea;
            }
        };
        const auto subtractOldTotalEqArea = [&](MeshType::HalfEdgeIndex hei) {
            if(mesh.isInTriangle(hei)) {
                oldTotalEqArea -= mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea;
            }
        };
        mesh.forEachHalfEdgeTargetingVertex(vi0, addOldTotalEqArea);
        mesh.forEachHalfEdgeTargetingVertex(vi1, addOldTotalEqArea);
        mesh.forEachHalfEdgeInEdge(ei, subtractOldTotalEqArea);

        // Do the edge collapse
        const auto change = MeshType::EdgeCollapse{}(mesh, hei);

        // Set attributes
        if(newPos.has_value()) mesh.attribute(change.viTo).vertex->coord = *newPos;

        // Redistribute properties
        double newTotalArea = 0.0;
        mesh.forEachHalfEdgeTargetingVertex(change.viTo, [&](MeshType::HalfEdgeIndex hei) {
            if(mesh.isInTriangle(hei)) {
                newTotalArea += area(mesh, mesh.triangle(hei));
            }
        });
        mesh.forEachHalfEdgeTargetingVertex(change.viTo, [&](MeshType::HalfEdgeIndex hei) {
            if(mesh.isInTriangle(hei)) {
                mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea
                    = oldTotalEqArea * area(mesh, mesh.triangle(hei)) / newTotalArea;
            }
        });

        return change;
    }

    static void flipEdge(
        MeshType&               mesh,
        MeshType::EdgeIndex     ei,
    ) {
        // Precondition:
        //   - ei is not on the border

        // Record old attributes
        double oldTotalEqArea = 0.0;
        mesh.forEachHalfEdgeInEdge(ei, [&](MeshType::HalfEdgeIndex hei) {
            oldTotalEqArea += mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea;
        });

        // Do edge flip
        MeshType::EdgeFlip{}(mesh, ei);

        // Redistribute attributes
        double newTotalArea = 0.0;
        mesh.forEachHalfEdgeInEdge(ei, [&](MeshType::HalfEdgeIndex hei) {
            newTotalArea += area(mesh, mesh.triangle(hei));
        });
        mesh.forEachHalfEdgeInEdge(ei, [&](MeshType::HalfEdgeIndex hei) {
            mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea
                = oldTotalEqArea * area(mesh, mesh.triangle(hei)) / newTotalArea;
        });

    }

};

#endif
