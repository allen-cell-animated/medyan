#include "Edge.h"

Edge::Edge(Composite* parent, Bead* b1, Bead* b2):
    _b1(b1), _b2(b2), _neighborTriangle{nullptr, nullptr} {
    
    parent -> addChild(unique_ptr<Component>(this));

}