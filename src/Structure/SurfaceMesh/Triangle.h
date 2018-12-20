#ifndef MEDYAN_Triangle_h
#define MEDYAN_Triangle_h

#include <array>

#include "common.h"

#include "Database.h"
#include "Trackable.h"
#include "Movable.h"
//#include "Reactable.h"
#include "DynamicNeighbor.h"
#include "Component.h"

#include "Vertex.h"
#include "Edge.h"
#include "GTriangle.h"
#include "MTriangle.h"

// Forward declarations
class Edge;
class Vertex;
class Compartment;

/******************************************************************************
Triangles are the only element in the meshwork that has area and act as patches
that construct the surface.

The triangle patches have geometric and mechanical properties.

The Triangle class has pointers to the vertices and edges.
******************************************************************************/
class Triangle:
    public Component,
    public Trackable,
    public Movable,
    // public Reactable,
    public DynamicNeighbor {

private:
    // Pointers to 3 vertices. The rotation of vertices must be pointing outside.
    array<Vertex*, 3> _v;

    unique_ptr<GTriangle> _gTriangle; // pointer to geo triangle
    unique_ptr<MTriangle> _mTriangle; // pointer to mech triangle

    // ptr to 3 edges. Order is (b0, b1), (b1, b2), (b2, b0)
    std::array<Edge*, 3> _edges;
    std::array<size_t, 3> _edgeHead; // The index of edge head vertex.
                                     // e.g. {1, 0, 1} means _edges' heads (Edge->v[0]) are at v1, v1, v0.

    static Database<Triangle*> _triangles; // Collection of triangles in SubSystem
    int _id; // Unique integer id of this triangle

    void updateCoordinate(); // helper function to update coordiante of this triangle

    Compartment* _compartment = nullptr; // The compartment containing this triangle

public:
    std::array<double, 3> coordinate; // Coordinate of the center point, updated with updateCoordiante()

    Triangle(Composite *parent, Vertex *v1, Vertex *v2, Vertex *v3);

    // Get Beads
    array<Vertex*, 3>& getVertices() { return _v; }
    array<Edge*, 3>& getEdges() { return _edges; }
    array<size_t, 3>& getEdgeHead() { return _edgeHead; }
    
    // Get geo triangle
    GTriangle* getGTriangle() { return _gTriangle.get(); }
    // Get mech triangle
    MTriangle* getMTriangle() { return _mTriangle.get(); }

    /// Get all instances of this class from the SubSystem
    static const vector<Triangle*>& getTriangles() {
        return _triangles.getElements();
    }
    /// Get ID
    int getId()const { return _id; }

    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem()override { _triangles.addElement(this); }
    virtual void removeFromSubSystem()override { _triangles.removeElement(this); }
    //@}

    //@{
    /// Implements Component
    virtual int getType() override { return getParent()->getType(); }
    virtual void printSelf()const override;
    //@}

    //@{
    /// Implements Movable
    virtual void updatePosition() override;
    //@}
    Compartment* getCompartment() { return _compartment; }


};

#endif
