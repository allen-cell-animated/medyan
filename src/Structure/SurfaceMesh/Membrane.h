#ifndef MEDYAN_Membrane_h
#define MEDYAN_Membrane_h

#include <array>
#include <tuple>
#include <vector>
#include <memory>

#include "Database.h"
#include "Geometric.h"
#include "Trackable.h"
#include "Composite.h"

#include "GMembrane.h"
#include "MathFunctions.h"
#include "MMembrane.h"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/Edge.h"
#include "Structure/SurfaceMesh/GEdge.h"
#include "Structure/SurfaceMesh/GHalfEdge.h"
#include "Structure/SurfaceMesh/GTriangle.h"
#include "Structure/SurfaceMesh/GVertex.h"
#include "Structure/SurfaceMesh/SurfaceMesh.hpp"
#include "Structure/SurfaceMesh/Triangle.h"
#include "Structure/SurfaceMesh/Vertex.h"

// Mesh type specification
struct MembraneMeshAttribute {
    struct VertexAttribute {
        using coordinate_type = decltype(Vertex::coordinate);
        Vertex* vertex;

        GVertex gVertex;

        coordinate_type& getCoordinate() { return vertex->coordinate; }

        void setIndex(size_t index) {
            vertex->setTopoIndex(index);
        }

        // TODO adaptive
    };
    struct EdgeAttribute {
        Edge* edge;

        GEdge gEdge;

        void setIndex(size_t index) {
            edge->setTopoIndex(index);
        }
        // TODO geometry / adaptive
    };
    struct HalfEdgeAttribute {
        GHalfEdge gHalfEdge;

        void setIndex(size_t index) {}
    };
    struct TriangleAttribute {
        Triangle* triangle;

        GTriangle gTriangle;

        void setIndex(size_t index) {
            triangle->setTopoIndex(index);
        }
        // TODO geometry / adaptive
    };
    struct MetaAttribute {
        SubSystem *s;
        Membrane *m;
    };

    using coordinate_type = VertexAttribute::coordinate_type;

    template< typename Mesh > static void newVertex(const MetaAttribute& meta, Mesh& mesh, size_t v, const typename Mesh::VertexInsertionOnEdge::InsertMid& op) {
        coordinate_type c0 = mesh.getVertexAttribute(op.v0).getCoordinate();
        coordinate_type c1 = mesh.getVertexAttribute(op.v1).getCoordinate();
        coordinate_type c = mathfunc::midPointCoordinate(c0, c1, 0.5);
        mesh.getVertexAttribute(v).vertex = meta.s->addTrackable<Vertex>(c, meta.m, v);
    }
    template< typename Mesh > static void newVertex(const MetaAttribute& meta, Mesh& mesh, size_t v, const coordinate_type& coord, const typename Mesh::GeometricVertexInit& op) {
        mesh.getVertexAttribute(v).vertex = meta.s->addTrackable<Vertex>(coord, meta.m, v);
    }
    template< typename Mesh, typename Operation > static void newEdge(const MetaAttribute& meta, Mesh& mesh, size_t e, const Operation& op) {
        mesh.getEdgeAttribute(e).edge = meta.s->addTrackable<Edge>(meta.m, e);
    }
    template< typename Mesh, typename Operation > static void newHalfEdge(const MetaAttribute& meta, Mesh& mesh, size_t he, const Operation& op) {
    }
    template< typename Mesh, typename Operation > static void newTriangle(const MetaAttribute& meta, Mesh& mesh, size_t t, const Operation& op) {
        mesh.getTriangleAttribute(t).triangle = meta.s->addTrackable<Triangle>(meta.m, t);
    }

    template< typename Mesh > static void removeVertex(const MetaAttribute& meta, Mesh& mesh, size_t v) {
        meta.s->removeTrackable<Vertex>(mesh.getVertexAttribute(v).vertex);
    }
    template< typename Mesh > static void removeEdge(const MetaAttribute& meta, Mesh& mesh, size_t e) {
        meta.s->removeTrackable<Edge>(mesh.getEdgeAttribute(e).edge);
    }
    template< typename Mesh > static void removeHalfEdge(const MetaAttribute& meta, Mesh& mesh, size_t he) {
    }
    template< typename Mesh > static void removeTriangle(const MetaAttribute& meta, Mesh& mesh, size_t t) {
        meta.s->removeTrackable<Triangle>(mesh.getEdgeAttribute(t).triangle);
    }

};

/******************************************************************************
Topologically, the membrane is represented by a 2d surface with 2 sides (which
means no Klein bottles are allowed!). The surface is constructed by
interconnected vertices, edges and triangles.

The Membrane class is a manager for constructing the mesh and computing the
geometry. It contains a meshwork instance that's responsible for adding and
removing vertices, edges (halfedges) and triangles to/from the SubSystem.
However, the ownership of all elements is in this Membrane class through
inheriting Composite.
******************************************************************************/

class Membrane: public Composite, public Trackable, public Geometric {

private:

    SurfaceTriangularMesh< MembraneMeshAttribute > _mesh;

    vector<Triangle*> _triangleVector; // collection of triangles
    vector<Edge*> _edgeVector; // collection of edges
    vector<Vertex*> _vertexVector; // collection of vertices

    unique_ptr<GMembrane> _gMembrane; // pointer to geometric membrane object
    unique_ptr<MMembrane> _mMembrane; // pointer to mechanical membrane object

    short _memType; // Membrane type

    SubSystem* _subSystem; // SubSystem pointer

    static Database<Membrane*> _membranes; // Collection in SubSystem
    int _id; // Unique integer id of this membrane

    bool _isClosed = true; // Whether the membrane is topologically closed, regardless of genus
    int _genus = 0; // Genus of the surface. Normally 0, as for a topologically spherical shape

public:
    // Constructors
    // This constructor creates a membrane according to vertex and neighbor data
    Membrane(
        SubSystem* s,
        short membraneType,
        const std::vector< coordinate_type >& vertexCoordinateList,
        const std::vector< std::array< size_t, 3 > >& triangleVertexIndexList
    );

    /// Get vector of triangles/edges/vertices that this membrane contains.
    const SurfaceTriangularMesh< MembraneMeshAttribute >& getMesh() const { return _mesh; }
    vector<Triangle*>& getTriangleVector() { return _triangleVector; }
    vector<Edge*>& getEdgeVector() { return _edgeVector; }
    vector<Vertex*>& getVertexVector() { return _vertexVector; }

    // Get Id
    int getId()const { return _id; }
    
    
    // SubSystem management, inherited from Trackable
    virtual void addToSubSystem()override { _membranes.addElement(this); }
    virtual void removeFromSubSystem()override { _membranes.removeElement(this); }
    
    /// Get all instances of this class from the SubSystem
    static const vector<Membrane*>& getMembranes() {
        return _membranes.getElements();
    }
    /// Get the number of membranes in this system
    static int numMembranes() {
        return _membranes.countElements();
    }

    //@{
    /// Implements Component
    // Get type
    int getType()override { return _memType; }
    // Print self information
    virtual void printSelf()const override;
    //@}

    /**************************************************************************
    Geometric
    **************************************************************************/
    /**
     * Implements Geometric.
     * If calcDerivative is true, most implementation would assume d is zero
     *   regardless of the actual passed d value. If calcDerivative is false,
     *   most implementation would store the result in member variables with
     *   "stretched" in their name.
     */
    void updateGeometryValue();
    virtual void updateGeometry(bool calcDerivative=false, double d=0.0)override;

    // Get geo membrane
    GMembrane* getGMembrane() { return _gMembrane.get(); }

    /**
     * Use pseudo normal signed distance field method to get the signed distance to a point.
     * If the point is outside, the result is positive and vice versa.
     * Throws an exception if the membrane is not closed.
     * The function will search through the whole meshwork, so it might not be efficient.
     * However, if the "safe" flag is turned off and neighboring compartments happen to contain membrane elements,
     *   search space will be limited to those compartments to save time. In this case meshwork size should be
     *   much smaller than the compartment size.
     */
    double signedDistance(const std::array<double, 3>& p, bool safe=false)const;
    /**
     * Use signed distance or other methods to judge whether a point is inside membrane.
     * Throws an exception if the membrane is not closed.
     */
    bool contains(const std::array<double, 3>& point)const;

    // Function to monitor the quality of the meshwork
    double meshworkQuality()const; // Must be used after updating the geometry
                                   // Returns a value between 0 and 1,
                                   // 1 being best and 0 being worst.

    /**************************************************************************
    Topological
    **************************************************************************/
    bool isClosed()const { return _isClosed; }

    /**************************************************************************
    Mechanics
    **************************************************************************/
    // Get mech membrane
    MMembrane* getMMembrane() { return _mMembrane.get(); }


};





#endif
