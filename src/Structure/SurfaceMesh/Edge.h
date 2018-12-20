#ifndef MEDYAN_Edge_h
#define MEDYAN_Edge_h

/*
 
An unordered edge contains 2 vertices.
 
*/

#include <array>

#include "common.h"

#include "Trackable.h"
#include "DynamicNeighbor.h"
#include "Component.h"
#include "Movable.h"

#include "Vertex.h"
#include "GEdge.h"

// Forward declarations
class Vertex;
class Triangle;
class Compartment;

class Edge:
    public Component,
    public Movable,
    public Trackable {

private:
    // Pointers to the vertices.
    std::array<Vertex*, 2> _v;

    std::array<Triangle*, 2> _triangles; // Neighbor triangles (v0, v1, other1) and (v1, v0, other2)

    unique_ptr<GEdge> _gEdge; // pointer to mech edge

    static Database<Edge*> _edges; // Collection of edges in SubSystem
    int _id; // Unique integer id of this edge

    void updateCoordinate(); // helper function to update coordiante of this edge

    Compartment* _compartment = nullptr; // The compartment containing this edge

public:
    Edge(Composite *parent, Vertex* v1, Vertex* v2);

    std::array<double, 3> coordinate; // Coordinate of the mid point, updated with updateCoordiante()

    // Get vertices
    std::array<Vertex*, 2>& getVertices() { return _v; }
    // Get neighbor triangles
    std::array<Triangle*, 2>& getTriangles() { return _triangles; }

    // Get mech edge
    GEdge* getGEdge() { return _gEdge.get(); }

    /// Get all instances of this class from the SubSystem
    static const vector<Edge*>& getEdges() {
        return _edges.getElements();
    }
    /// Get Id
    int getId()const { return _id; }

    //@{
    /// Implements Component
    virtual int getType() override { return getParent()->getType(); }
    virtual void printSelf()const override;
    //@}

    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem()override { _edges.addElement(this); }
    virtual void removeFromSubSystem()override { _edges.removeElement(this); }
    //@}

    //@{
    /// Implements Movable
    virtual void updatePosition() override;
    //@}
    Compartment* getCompartment() { return _compartment; }

};

#endif
