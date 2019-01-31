#include "Edge.h"

#include "core/globals.h"
#include "Compartment.h"
#include "core/controller/GController.h"
#include "MathFunctions.h"
#include "Structure/SurfaceMesh/Membrane.hpp"

Database<Edge*> Edge::_edges;

Edge::Edge(Composite* parent, size_t topoIndex):
    Trackable(),
    _topoIndex{topoIndex}, _id(_edges.getID()) {
    
    parent -> addChild(unique_ptr<Component>(this));

    // Set coordinate and add to compartment
    updateCoordinate();
    if(medyan::Global::readGlobal().mode == medyan::GlobalVar::RunMode::Simulation) {
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
    const auto& mesh = static_cast<Membrane*>(getParent())->getMesh();
    const size_t hei0 = mesh.getEdges()[_topoIndex].halfEdgeIndex;
    const size_t hei1 = mesh.opposite(hei0);
    const size_t v0 = mesh.target(hei0);
    const size_t v1 = mesh.target(hei1);

    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        coordinate[coordIdx] = (
            mesh.getVertexAttribute(v0).vertex->coordinate[coordIdx]
            + mesh.getVertexAttribute(v1).vertex->coordinate[coordIdx]
        ) / 2;
    }
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

void Edge::printSelf()const {
    cout << endl;
    
    cout << "Edge: ptr = " << this << endl;
    cout << "Edge ID = " << _id << endl;
    cout << "Parent ptr = " << getParent() << endl;
        
    cout << endl;
}
