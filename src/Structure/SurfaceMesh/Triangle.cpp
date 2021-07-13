#include "Triangle.hpp"

#include "Core/Globals.hpp"
#include "Compartment.h"
#include "MathFunctions.h"
#include "GController.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Util/Io/Log.hpp"

Triangle::Triangle(Membrane* parent, size_t topoIndex):
    Trackable(true, false, true, false),
    _parent(parent), _topoIndex{topoIndex} {

    // Set coordinate and add to compartment
    updateCoordinate();

    Compartment* compartment;
    try { compartment = GController::getCompartment(mathfunc::vec2Vector(coordinate)); }
    catch (exception& e) {
        cout << e.what() << endl;
        printSelf();
        exit(EXIT_FAILURE);
    }
    _cellElement.manager = compartment->triangleCell.manager;
    _cellElement.manager->addElement(this, _cellElement, compartment->triangleCell);
   
}

Triangle::~Triangle() {
    _cellElement.manager->removeElement(_cellElement);
}

void Triangle::updateCoordinate() {
    const auto& mesh = _parent->getMesh();
    const auto hei0 = mesh.halfEdge(Membrane::MeshType::TriangleIndex{_topoIndex});
    const auto hei1 = mesh.next(hei0);
    const auto hei2 = mesh.next(hei1);
    const auto v0 = mesh.target(hei0);
    const auto v1 = mesh.target(hei1);
    const auto v2 = mesh.target(hei2);

    coordinate = (
        mesh.attribute(v0).getCoordinate()
        + mesh.attribute(v1).getCoordinate()
        + mesh.attribute(v2).getCoordinate()
    ) / 3;
}

void Triangle::updatePosition() {
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

/* TODO:
#ifdef CHEMISTRY
        auto oldCompartment = _compartment;
        auto newCompartment = c;
#endif
*/
        
        //remove from old compartment, add to new
        _cellElement.manager->updateElement(_cellElement, c->triangleCell);

/* TODO:
#ifdef CHEMISTRY
        auto oldCCylinder = _cCylinder.get();
        
        //Remove old ccylinder from binding managers
        for(auto &manager : oldCompartment->getFilamentBindingManagers())
            manager->removePossibleBindings(oldCCylinder);
        
        //clone and set new ccylinder
        CCylinder* clone = _cCylinder->clone(c);
        setCCylinder(clone);
        
        auto newCCylinder = _cCylinder.get();
        
        //Add new ccylinder to binding managers
        for(auto &manager : newCompartment->getFilamentBindingManagers())
            manager->addPossibleBindings(newCCylinder);
#endif
*/
    }

}

int Triangle::getType()const {
    return _parent->getType();
}

void Triangle::printSelf()const {
    
    cout << endl;
    
    cout << "Triangle: ptr = " << this << endl;
    cout << "Triangle ID = " << getId() << endl;
    cout << "Parent ptr = " << getParent() << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    cout << endl;
    
    /*
#ifdef CHEMISTRY
    cout << "Chemical composition of triangle:" << endl;
    _cTriangle->printCTriangle();
#endif
    
    cout << endl;
    */
}
