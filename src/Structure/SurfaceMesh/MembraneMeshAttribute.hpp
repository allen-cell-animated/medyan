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
#include "Structure/SurfaceMesh/HalfEdge.hpp"
#include "Structure/SurfaceMesh/SurfaceMesh.hpp"
#include "Structure/SurfaceMesh/Triangle.hpp"
#include "Structure/SurfaceMesh/Vertex.hpp"
#include "Util/Io/Log.hpp"

// Forward declarations
class Membrane;

// Membrane mesh chemistry data
struct MembraneMeshChemistryInfo {
    // Note:
    //   - indices in this struct must correspond to the actual index as is in
    //     the diffusing species names vector

    struct DiffusionInfo {
        unsigned speciesIndex = 0;

        // Diffusion coeffecient, with dimension L^2 T^-1
        double   diffusionCoeff = 0;
    };
    struct InternalReactionInfo {
        std::vector< unsigned > reactantSpeciesIndices;
        std::vector< unsigned > productSpeciesIndices;

        // Rate constant, with dimension L^(2*(numReactants - 1)) T^-1
        double rateConstant = 0;
    };

    // diffusing species
    std::vector< std::string > diffusingSpeciesNames;

    // diffusion reactions
    std::vector< DiffusionInfo > diffusion;

    // Internal reactions involving species
    std::vector< InternalReactionInfo > internalReactions;
};

// Returns whether the chemisty info is valid
inline bool validate(const MembraneMeshChemistryInfo& ci) {
    bool valid = true;

    const auto ns = ci.diffusingSpeciesNames.size();

    // Check diffusion
    for(const auto& di : ci.diffusion) {
        if(di.speciesIndex >= ns) {
            LOG(ERROR) << "Diffusion species index " << di.speciesIndex
                << " is out of range.";
            valid = false;
        }
        if(di.diffusionCoeff < 0) {
            LOG(ERROR) << "Diffusion coefficient " << di.diffusionCoeff
                << " is unacceptable.";
            valid = false;
        }
    }

    // Check internal reaction
    for(const auto& iri : ci.internalReactions) {
        for(auto i : iri.reactantSpeciesIndices) {
            if(i >= ns) {
                LOG(ERROR) << "Internal reaction reactant species index " << i
                    << " is out of range.";
                valid = false;
            }
        }
        for(auto i : iri.productSpeciesIndices) {
            if(i >= ns) {
                LOG(ERROR) << "Internal reaction product species index " << i
                    << " is out of range.";
                valid = false;
            }
        }
        if(iri.rateConstant < 0) {
            LOG(ERROR) << "Internal reaction rate constant " << iri.rateConstant
                << " is unacceptable.";
            valid = false;
        }
    }

    return valid;
}


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
    //
    // If vertices represent fixed coordinates in a normal coordinate system,
    // then the system is similar to the "normal" coordinate case, but without
    // the requirement for vertices to move only in the normal directions.
    enum class VertexSystem { material, normal, general };

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

        CoordinateRefType  getCoordinate()       { return vertex->coord; }
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
        size_t                              cachedCoordIndex[2];
        HalfEdgeMeshConnection::PolygonType cachedPolygonType[2];
        // [polygon of half edge, opposite polygon]
        size_t                              cachedPolygonIndex[2];

        template< bool stretched > const GEdge& getGEdge() const { return stretched ? gEdgeS : gEdge; }
        template< bool stretched >       GEdge& getGEdge()       { return stretched ? gEdgeS : gEdge; }

        void setIndex(size_t index) {
            edge->setTopoIndex(index);
        }
    };
    struct HalfEdgeAttribute {
        std::unique_ptr< medyan::HalfEdge > halfEdge;

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
        size_t                                cachedCoordIndex[3];
        // [half edge, next, prev]
        HalfEdgeMeshConnection::HalfEdgeIndex cachedHalfEdgeIndex[3];
        // [edge of half edge, next edge, prev edge]
        HalfEdgeMeshConnection::EdgeIndex     cachedEdgeIndex[3];

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

        // About mech and chem
        //-----------------------------
        bool isMechParamsSet = false;
        bool isChemParamsSet = false;

        // Chemistry information
        //-----------------------------
        MembraneMeshChemistryInfo chemInfo;

        // Cache related
        //-----------------------------
        bool indexCacheForFFValid = false;

        size_t vertexMaxDegree;

        // A vertex has undetermined number of neighbors, so the cache structure needs to be determined at run-time.
        //
        // Notes:
        //   - The outer half edge targets the corresponding "neighbor vertex", and is the "prev" of the targeting half edge.
        //   - The polygon corresponds to the targeting half edge.
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
    static void newVertex(MeshType& mesh, HalfEdgeMeshConnection::VertexIndex v) {
        mesh.attribute(v).vertex = std::make_unique< Vertex >(CoordinateType{}, v.index);
    }
    static void newEdge(MeshType& mesh, HalfEdgeMeshConnection::EdgeIndex e) {
        mesh.attribute(e).edge.reset(
            mesh.metaAttribute().s->addTrackable<Edge>(mesh.metaAttribute().m, e.index));
    }
    static void newHalfEdge(MeshType& mesh, HalfEdgeMeshConnection::HalfEdgeIndex he) {
        mesh.attribute(he).halfEdge = std::make_unique< medyan::HalfEdge >();
    }
    static void newTriangle(MeshType& mesh, HalfEdgeMeshConnection::TriangleIndex t) {
        mesh.attribute(t).triangle.reset(
            mesh.metaAttribute().s->addTrackable<Triangle>(mesh.metaAttribute().m, t.index));
    }
    static void newBorder(MeshType& mesh, HalfEdgeMeshConnection::BorderIndex) {
        // Do nothing
    }

    static void removeElement(MeshType& mesh, HalfEdgeMeshConnection::VertexIndex i) {
        // Do nothing
    }
    static void removeElement(MeshType& mesh, HalfEdgeMeshConnection::EdgeIndex i) {
        mesh.metaAttribute().s->template removeTrackable<Edge>(mesh.attribute(i).edge.get());
    }
    static void removeElement(MeshType& mesh, HalfEdgeMeshConnection::HalfEdgeIndex i) {
        // Do nothing
    }
    static void removeElement(MeshType& mesh, HalfEdgeMeshConnection::TriangleIndex i) {
        mesh.metaAttribute().s->template removeTrackable<Triangle>(mesh.attribute(i).triangle.get());
    }
    static void removeElement(MeshType& mesh, HalfEdgeMeshConnection::BorderIndex i) {
        // Do nothing
    }

    // Mesh attribute initializing and extracting
    // These operations do not follow the RAII idiom.
    // Initialization should happen only once, as it allocates resources
    static void init(MeshType& mesh, const AttributeInitializerInfo& info) {
        const MetaAttribute& meta = mesh.metaAttribute();
        for(size_t i = 0; i < mesh.getVertices().size(); ++i) {
            mesh.attribute(HalfEdgeMeshConnection::VertexIndex{i}).vertex = std::make_unique<Vertex>(
                info.vertexCoordinateList[i], i);
        }
        for(size_t i = 0; i < mesh.getHalfEdges().size(); ++i) {
            mesh.attribute(HalfEdgeMeshConnection::HalfEdgeIndex{i}).halfEdge
                = std::make_unique< medyan::HalfEdge >();
        }
        for(size_t i = 0; i < mesh.getEdges().size(); ++i) {
            mesh.attribute(HalfEdgeMeshConnection::EdgeIndex{i}).edge.reset(
                meta.s->template addTrackable<Edge>(meta.m, i));
        }
        for(size_t i = 0; i < mesh.getTriangles().size(); ++i) {
            mesh.attribute(HalfEdgeMeshConnection::TriangleIndex{i}).triangle.reset(
                meta.s->template addTrackable<Triangle>(meta.m, i));
        }
    }
    // Extraction can be done multiple times without allocating/deallocating
    static auto extract(const MeshType& mesh) {
        AttributeInitializerInfo info;

        const size_t numVertices = mesh.getVertices().size();

        info.vertexCoordinateList.reserve(numVertices);
        for(size_t i = 0; i < numVertices; ++i) {
            info.vertexCoordinateList.push_back(
                static_cast<CoordinateType>(mesh.attribute(HalfEdgeMeshConnection::VertexIndex{i}).vertex->coord));
        }

        return info;
    }

};

#endif
