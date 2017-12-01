#include "Edge.h"

Edge::Edge(Composite* parent, Vertex* v1, Vertex* v2):
    _v{v1, v2} {
    
    parent -> addChild(unique_ptr<Component>(this));

}