#ifndef MEDYAN_Triangle_h
#define MEDYAN_Triangle_h

#include <array>

#include "common.h"

#include "Trackable.h"
#include "Movable.h"
#include "Reactable.h"
#include "DynamicNeighbor.h"
#include "Component.h"

#include "Bead.h"
#include "MTriangle.h"

class Triangle:
    public Component,
    public Trackable,
    public Movable,
    public Reactable, // TODO: is it reactable?
    public DynamicNeighbor {

private:
    // Pointers to 3 beads. The rotation of beads must be pointing outside.
    array<Bead*, 3> _b;

    unique_ptr<MTriangle> _mTriangle; // pointer to mech triangle

    // ptr to neighbor triangles
    // For each vertex, the neighbor is the nearest triangle in the counter-clockwise direction
    std::array<Triangle*, 3> _neighborTriangle;

public:
    Triangle(Composite *parent, Bead *b1, Bead *b2, Bead *b3);

    // Get Beads
    array<Bead*, 3>& getBeads() { return _b; }
    
    // Get mech triangle
    MTriangle* getMTriangle() {return _mTriangle.get();}


};


#endif
