#include "Edge.h"

Edge::Edge(Composite* parent, HalfEdge* he1, HalfEdge* he2):
    _he{he1, he2} {
    
    parent -> addChild(unique_ptr<Component>(this));

}