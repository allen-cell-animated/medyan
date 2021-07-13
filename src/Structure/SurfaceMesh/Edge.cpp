#include "Structure/SurfaceMesh/Edge.hpp"

#include "Core/Globals.hpp"
#include "Compartment.h"
#include "GController.h"
#include "MathFunctions.h"
#include "Structure/SurfaceMesh/Membrane.hpp"

Edge::Edge(Membrane* parent, size_t topoIndex):
    Trackable(),
    parent_(parent), topoIndex_{topoIndex} {
    
    // Set coordinate and add to compartment
    updateCoordinate();

    Compartment* compartment;
    try { compartment = GController::getCompartment(mathfunc::vec2Vector(coordinate)); }
    catch (exception& e) {
        cout << e.what() << endl;
        printSelf();
        exit(EXIT_FAILURE);
    }
    _cellElement.manager = compartment->edgeCell.manager;
    _cellElement.manager->addElement(this, _cellElement, compartment->edgeCell);
}

Edge::~Edge() {
    _cellElement.manager->removeElement(_cellElement);
}

void Edge::updateCoordinate() {
    const auto& mesh = parent_->getMesh();
    const auto hei0 = mesh.halfEdge(Membrane::MeshType::EdgeIndex{topoIndex_});
    const auto hei1 = mesh.opposite(hei0);
    const auto v0 = mesh.target(hei0);
    const auto v1 = mesh.target(hei1);

    coordinate = 0.5 * (mesh.attribute(v0).getCoordinate() + mesh.attribute(v1).getCoordinate());
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
    Compartment* curCompartment = getCompartment();
    if(c != curCompartment) {
     
        //remove from old compartment, add to new
        _cellElement.manager->updateElement(_cellElement, c->edgeCell);
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
