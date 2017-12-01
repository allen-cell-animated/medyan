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

class Edge:
    public Component,
    public Trackable,
    public Movable,
    public DynamicNeighbor {

    friend class MEdge;

private:
    // Pointers to the vertices.
    std::array<Vertex*, 2> _v;

    unique_ptr<MEdge> _mEdge; // pointer to mech edge


public:
    Edge(Composite *parent, Vertex* v1, Vertex* v2);

    // Get vertices
    std::array<Vertex*, 2>& getVertices() { return _v; }

    // Get mech edge
    MEdge* getMEdge() {return _mEdge.get();}

};


#endif
