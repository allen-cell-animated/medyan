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

#include "HalfEdge.h"

class Edge:
    public Component,
    public Trackable,
    public Movable,
    public DynamicNeighbor {

private:
    // Pointers to the halfedges.
    std::array<HalfEdge*, 2> _he;

    unique_ptr<MEdge> _mEdge; // pointer to mech edge


public:
    Edge(Composite *parent, HalfEdge* he1, HalfEdge* he2);

    // Get Beads
    std::array<HalfEdge*, 2>& getHalfEdges() { return _he; }

};


#endif
