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

class HalfEdge:
    public Component,
    public Trackable,
    public Movable,
    public DynamicNeighbor {

private:
    // Pointers to the halfedges.
    HalfEdge* _he1;
	HalfEdge* _he2;

public:
    Edge(Composite *parent, HalfEdge* he1, HalfEdge* he2);

    // Get Beads
    Bead* getFirstHalfEdge() { return _he1; }
	Bead* getSecondHalfEdge() { return _he2; }

};


#endif
