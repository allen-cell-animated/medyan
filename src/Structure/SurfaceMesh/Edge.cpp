#include "Edge.h"

Edge::Edge(Composite* parent, HalfEdge* he1, HalfEdge* he2):
    _he1(b1), _he2(he2) {
    
    parent -> addChild(unique_ptr<Component>(this));

}