#ifndef MEDYAN_Edge_h
#define MEDYAN_Edge_h

/*
 
 The edge containes two halfedges in the opposite direction.
 
 By using collection of edges, the 
 
*/

#include <array>

#include "common.h"

#include "Trackable.h"
#include "Movable.h"
#include "DynamicNeighbor.h"
#include "Component.h"

#include "Vertex.h"
#include "MEdge.h"

// Forward declarations
class Vertex;

class Edge:
    public Component,
    public Trackable,
    public Movable,
    public DynamicNeighbor {

private:
    // Pointers to the vertices.
    std::array<Vertex*, 2> _v;
    // For simplicity, currently the neighbor triangles are not recorded.

    unique_ptr<MEdge> _mEdge; // pointer to mech edge

    static Database<Edge*> _edges; // Collection of edges in SubSystem
    int _Id; // Unique integer id of this edge

public:
    Edge(Composite *parent, Vertex* v1, Vertex* v2);

    // Get vertices
    std::array<Vertex*, 2>& getVertices() { return _v; }

    // Get mech edge
    MEdge* getMEdge() { return _mEdge.get(); }

    /// Get all instances of this class from the SubSystem
    static const vector<Edge*>& getEdges() {
        return _edges.getElements();
    }
    /// Get ID
    int getId()const { return _Id; }

    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem()override { _edges.addElement(this); }
    virtual void removeFromSubSystem()override { _edges.removeElement(this); }
    //@}


};


#endif
