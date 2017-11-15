#include "Triangle.h"

Triangle::Triangle(Composite* parent, Bead* b1, Bead* b2, Bead* b3):
    _b{b1, b2, b3}, _halfEdges{nullptr, nullptr, nullptr} {
    
    parent -> addChild(unique_ptr<Component>(this));

}