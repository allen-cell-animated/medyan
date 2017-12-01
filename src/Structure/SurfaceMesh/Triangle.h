#ifndef MEDYAN_Triangle_h
#define MEDYAN_Triangle_h

#include <array>

#include "common.h"

#include "Trackable.h"
#include "Movable.h"
#include "Reactable.h"
#include "DynamicNeighbor.h"
#include "Component.h"

#include "Vertex.h"
#include "HalfEdge.h"
#include "MTriangle.h"

class Triangle:
    public Component,
    public Trackable,
    public Movable,
    public Reactable, // TODO: is it reactable?
    public DynamicNeighbor {

private:
    // Pointers to 3 vertices. The rotation of vertices must be pointing outside.
    array<Vertex*, 3> _v;

    unique_ptr<MTriangle> _mTriangle; // pointer to mech triangle

    // ptr to 3 edges. Order is (b0, b1), (b1, b2), (b2, b0)
    std::array<HalfEdge*, 3> _halfEdges;

public:
    Triangle(Composite *parent, Vertex *v1, Vertex *v2, Vertex *v3);

    // Get Beads
    array<Vertex*, 3>& getVertices() { return _v; }
    
    // Get mech triangle
    MTriangle* getMTriangle() {return _mTriangle.get();}


};


#endif
