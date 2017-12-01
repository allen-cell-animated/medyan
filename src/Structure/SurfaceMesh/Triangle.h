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
    std::array<Edge*, 3> _edges;
    std::array<size_t, 3> _edgeHead; // The index of edge head vertex.
                                     // e.g. {1, 0, 1} means _edges' heads (Edge->v[0]) are at v1, v1, v0.

public:
    Triangle(Composite *parent, Vertex *v1, Vertex *v2, Vertex *v3);

    // Get Beads
    array<Vertex*, 3>& getVertices() { return _v; }
    array<Edge*, 3>& getEdges() { return _edges; }
    array<size_t, 3>& getEdgeHead() { return _edgeHead; }
    
    // Get mech triangle
    MTriangle* getMTriangle() {return _mTriangle.get();}


};


#endif
