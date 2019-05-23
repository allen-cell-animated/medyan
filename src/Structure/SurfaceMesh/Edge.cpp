#include "Structure/SurfaceMesh/Edge.h"

#include "Core/Globals.hpp"
#include "Compartment.h"
#include "GController.h"
#include "MathFunctions.h"
#include "Structure/SurfaceMesh/Membrane.hpp"

Edge::Edge(Membrane* parent, size_t topoIndex):
    Trackable(),
    _parent(parent), _topoIndex{topoIndex} {
    
    // Set coordinate and add to compartment
    updateCoordinate();
    if(medyan::global().mode == medyan::GlobalVar::RunMode::Simulation) {
        try { _compartment = GController::getCompartment(mathfunc::vec2Vector(coordinate)); }
        catch (exception& e) {
            cout << e.what() << endl;
            printSelf();
            exit(EXIT_FAILURE);
        }
        _compartment->addEdge(this);
    }
}

Edge::~Edge() {
    _compartment->removeEdge(this);
}

void Edge::updateCoordinate() {
    const auto& mesh = _parent->getMesh();
    const size_t hei0 = mesh.getEdges()[_topoIndex].halfEdgeIndex;
    const size_t hei1 = mesh.opposite(hei0);
    const size_t v0 = mesh.target(hei0);
    const size_t v1 = mesh.target(hei1);

    coordinate = 0.5 * (mesh.getVertexAttribute(v0).getCoordinate() + mesh.getVertexAttribute(v1).getCoordinate());
}

void Edge::updatePosition() {
    updateCoordinate();
    
    // Get the current compartment
    Compartment *c;
    try { c = GController::getCompartment(mathfunc::vec2Vector(coordinate)); }
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

int Edge::getType()const {
    return getParent()->getType();
}

void Edge::printSelf()const {
    cout << endl;
    
    cout << "Edge: ptr = " << this << endl;
    cout << "Edge ID = " << getId() << endl;
    cout << "Parent ptr = " << getParent() << endl;
        
    cout << endl;
}
