#ifndef MEDYAN_SurfaceMesh_hpp
#define MEDYAN_SurfaceMesh_hpp

#include <array>
#include <vector>

/******************************************************************************
The data structure for an orientable, manifold 2d triangular meshwork in 3d
space.

The connectivity is similar to CGAL library, which is halfedge based.

target:     The vertex that the halfedge points to.
opposite:   The halfedge on the same edge but in the opposite direction.
face:       The triangle that the halfedge is circling ccw.
next:       The next halfedge in the current face circling ccw.
prev:       The previous halfedge in the current face circling ccw.
edge:       The undirected edge of this halfedge.

All the other elements must have at least one index pointing to a halfedge.
******************************************************************************/

class SurfaceTriangularMesh {
public:

    // The elements should be trivially copyable.
    struct Vertex {
        bool markedAsDeleted = false;
        size_t halfEdgeIndex; // Only one HalfEdge targeting the vertex is needed.
    };
    struct HalfEdge {
        bool markedAsDeleted = false;
        bool hasOpposite;
        size_t triangleIndex;
        size_t targetVertexIndex;
        size_t oppositeHalfEdgeIndex;
        size_t edgeIndex;
    };
    struct Edge {
        bool markedAsDeleted = false;
        size_t halfEdgeIndex; // Only one HalfEdge is needed.
    };
    struct Triangle {
        bool markedAsDeleted = false;
        size_t halfEdgeIndex; // Only one HalfEdge is needed.
    };

private:

    std::vector<Triangle> _triangles; // collection of triangles
    std::vector<HalfEdge> _halfEdges; // collection of halfedges
    std::vector<Edge> _edges; // collection of edges
    std::vector<Vertex> _vertices; // collection of vertices

    bool _isClosed = true; // Whether the meshwork is topologically closed
    int _genus = 0; // Genus of the surface. Normally 0, as for a topologically spherical shape

public:
    // Constructors
    SurfaceMesh(SubSystem* s, short membraneType,
        const vector<tuple<array<double, 3>, vector<size_t>>>& membraneData);

    // This destructor is called when a membrane is to be removed from the system.
    // Removes all vertices, edges and triangles associated with the membrane.
    ~SurfaceMesh();

    // Initialize the meshwork using triangle vertex index lists. Throws on error.
    void init(const std::vector< std::array< size_t, 3 > >& triangleVertexIndexList);

    /// Get vector of triangles/edges/vertices that this membrane contains.
    vector<Triangle*>& getTriangleVector() {return _triangleVector;}
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
