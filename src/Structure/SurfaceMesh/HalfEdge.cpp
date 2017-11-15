#include "HalfEdge.h"

HalfEdge::HalfEdge(Composite* parent, Bead* b1, Bead* b2):
    _b1(b1), _b2(b2), _neighborTriangle{nullptr, nullptr} {
    
    parent -> addChild(unique_ptr<Component>(this));

}