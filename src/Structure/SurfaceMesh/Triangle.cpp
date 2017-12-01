#include "Triangle.h"

Triangle::Triangle(Composite* parent, Vertex* v1, Vertex* v2, Vertex* v3):
    _v{v1, v2, v3}, _halfEdges{nullptr, nullptr, nullptr} {
    
    parent -> addChild(unique_ptr<Component>(this));

}