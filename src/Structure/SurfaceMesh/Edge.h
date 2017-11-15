#ifndef MEDYAN_Edge_h
#define MEDYAN_Edge_h

/*
 
 This structure defines the shared edge between two
 triangular surface patches.
 
 The structure is similar to cylinder
 
*/

#include <array>

#include "common.h"

#include "Trackable.h"
#include "Movable.h"
#include "DynamicNeighbor.h"
#include "Component.h"

#include "Bead.h"
#include "Triangle.h"
#include "MTriangle.h"

class Edge:
    public Component,
    public Trackable,
    public Movable,
    public DynamicNeighbor {

private:
    // Pointers to the beads. The order is irrelevant.
    Bead* _b1;
	Bead* _b2;

    // ptr to neighbor triangles (_b1, _b2, other1) and (_b2, _b1, other2) <- counter-clockwise dir
    std::array<Triangle*, 2> _neighborTriangle;

public:
    Edge(Composite *parent, Bead *b1, Bead *b2);

    // Get Beads
    Bead* getBead1() { return _b1; }
	Bead* getBead2() { return _b2; }
    
	// TODO: Add getter and setter for neighbor triangles

};


#endif
