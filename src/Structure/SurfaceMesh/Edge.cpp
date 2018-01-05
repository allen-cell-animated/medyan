#include "Edge.h"

Database<Edge*> Edge::_edges;

Edge::Edge(Composite* parent, Vertex* v1, Vertex* v2):
    Trackable(),
    _v{v1, v2}, _id(_edges.getID()) {
    
    parent -> addChild(unique_ptr<Component>(this));

    _gEdge = unique_ptr<GEdge>(new GEdge);
    _gEdge->setEdge(this);

}

void Edge::updateCoordinate() {
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        coordinate[coordIdx] = (_v[0]->coordinate[coordIdx] + _v[1]->coordinate[coordIdx]) / 2;
    }
}

void Edge::updatePosition() {
    updateCoordinate();
    
    // If compartments are implemented, they are also updated here.
}

void Edge::printSelf() {
    cout << endl;
    
    cout << "Edge: ptr = " << this << endl;
    cout << "Edge ID = " << _id << endl;
    cout << "Parent ptr = " << getParent() << endl;
        
    cout << endl;

    cout << "Vertex information..." << endl;
    
    _v[0]->printSelf();
    _v[1]->printSelf();
    
    cout << endl;
}

void Edge::updateGeometry(bool calcDerivative, double d) {

#ifdef MECHANICS
    _gEdge->updateGeometry(calcDerivative, d);
#endif

}
