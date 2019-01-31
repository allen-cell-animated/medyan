#include "Triangle.h"

#include "core/globals.h"
#include "Compartment.h"
#include "MathFunctions.h"
#include "core/controller/GController.h"
#include "Structure/SurfaceMesh/Membrane.hpp"

Database<Triangle*> Triangle::_triangles;

Triangle::Triangle(Composite* parent, size_t topoIndex):
    Trackable(true, false, true, false),
    _topoIndex{topoIndex}, _id(_triangles.getID()) {
    
    parent -> addChild(unique_ptr<Component>(this));

#ifdef MECHANICS
    // eqArea cannot be obtained at this moment
    _mTriangle = std::make_unique<MTriangle>(getType());
#endif

    // Set coordinate and add to compartment
    updateCoordinate();
    if(medyan::Global::readGlobal().mode == medyan::GlobalVar::RunMode::Simulation) {
        try { _compartment = GController::getCompartment(mathfunc::vec2Vector(coordinate)); }
        catch (exception& e) {
            cout << e.what() << endl;
            printSelf();
            exit(EXIT_FAILURE);
        }
        _compartment->addTriangle(this);
    }
   
}

Triangle::~Triangle() {
    _compartment->removeTriangle(this);
}

void Triangle::updateCoordinate() {
    const auto& mesh = static_cast<Membrane*>(getParent())->getMesh();
    const size_t hei0 = mesh.getTriangles()[_topoIndex].halfEdgeIndex;
    const size_t hei1 = mesh.next(hei0);
    const size_t hei2 = mesh.next(hei1);
    const size_t v0 = mesh.target(hei0);
    const size_t v1 = mesh.target(hei1);
    const size_t v2 = mesh.target(hei2);

    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        coordinate[coordIdx] = (
            mesh.getVertexAttribute(v0).vertex->coordinate[coordIdx]
            + mesh.getVertexAttribute(v1).vertex->coordinate[coordIdx]
            + mesh.getVertexAttribute(v2).vertex->coordinate[coordIdx]
        ) / 3;
    }
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
    if(c != _compartment) {

/* TODO:
#ifdef CHEMISTRY
        auto oldCompartment = _compartment;
        auto newCompartment = c;
#endif
*/
        
        //remove from old compartment, add to new
        _compartment->removeTriangle(this);
        _compartment = c;
        _compartment->addTriangle(this);

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

void Triangle::printSelf()const {
    
    cout << endl;
    
    cout << "Triangle: ptr = " << this << endl;
    cout << "Triangle ID = " << _id << endl;
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
