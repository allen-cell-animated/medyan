#include "Edge.h"

Database<Edge*> Edge::_edges;

Edge::Edge(Composite* parent, Vertex* v1, Vertex* v2):
    _v{v1, v2}, _Id(_edges.getID()) {
    
    parent -> addChild(unique_ptr<Component>(this));

#ifdef MECHANICS
    _mEdge = unique_ptr<MEdge>(new MEdge);
    _mEdge->setEdge(this);
#endif

}