#include "Vertex.h"

Vertex::Vertex(Composite* parent, Bead* b):
    _b(b) {
    
    parent -> addChild(unique_ptr<Component>(this));

}