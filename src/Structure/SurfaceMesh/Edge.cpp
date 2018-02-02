#include "Edge.h"

#include "Compartment.h"
#include "GController.h"
#include "MathFunctions.h"

Database<Edge*> Edge::_edges;

Edge::Edge(Composite* parent, Vertex* v1, Vertex* v2):
    Trackable(),
    _v{v1, v2}, _id(_edges.getID()) {
    
    parent -> addChild(unique_ptr<Component>(this));

    _gEdge = unique_ptr<GEdge>(new GEdge);
    _gEdge->setEdge(this);

    // Set coordinate and add to compartment
    updateCoordinate();
    try { _compartment = GController::getCompartment(mathfunc::array2Vector(coordinate)); }
    catch (exception& e) {
        cout << e.what() << endl;
        printSelf();
        exit(EXIT_FAILURE);
    }
   _compartment->addEdge(this);
}

void Edge::updateCoordinate() {
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        coordinate[coordIdx] = (_v[0]->coordinate[coordIdx] + _v[1]->coordinate[coordIdx]) / 2;
    }
}

void Edge::updatePosition() {
    updateCoordinate();
    
    // Get the current compartment
    Compartment *c;
    try { c = GController::getCompartment(mathfunc::array2Vector(coordinate)); }
    catch (exception& e) {
        cout << e.what() << endl;
        printSelf();
        exit(EXIT_FAILURE);
    }

    // Things to do if the comparment changes
    if(c != _compartment) {
     
        //remove from old compartment, add to new
        _compartment->removeEdge(this);
        _compartment = c;
        _compartment->addEdge(this);
    }
}

void Edge::printSelf()const {
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
