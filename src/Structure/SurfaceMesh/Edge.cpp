#include "Edge.h"

Database<Edge*> Edge::_edges;

Edge::Edge(Composite* parent, Vertex* v1, Vertex* v2):
    Trackable(),
    _v{v1, v2}, _id(_edges.getID()) {
    
    parent -> addChild(unique_ptr<Component>(this));

#ifdef MECHANICS
    _mEdge = unique_ptr<MEdge>(new MEdge);
    _mEdge->setEdge(this);
#endif

}

void Edge::printSelf() {
    cout << endl;
    
    cout << "Edge: ptr = " << this << endl;
    cout << "Edge ID = " << _Id << endl;
    cout << "Parent ptr = " << getParent() << endl;
        
    cout << endl;

    cout << "Vertex information..." << endl;
    
    _v[0]->printSelf();
    _v[1]->printSelf();
    
    cout << endl;
}

void Edge::updateGeometry(bool calcDerivative, double d) {

#ifdef MECHANICS
    _mEdge->updateGeometry(calcDerivative, d);
#endif

}
