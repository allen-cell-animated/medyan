#include "Triangle.h"

#include "Compartment.h"
#include "MathFunctions.h"
#include "GController.h"

Database<Triangle*> Triangle::_triangles;

Triangle::Triangle(Composite* parent, Vertex* v1, Vertex* v2, Vertex* v3):
    Trackable(true, false, true, false),
    _v{v1, v2, v3}, _edges{nullptr, nullptr, nullptr}, _id(_triangles.getID()) {
    
    parent -> addChild(unique_ptr<Component>(this));

    _gTriangle = unique_ptr<GTriangle>(new GTriangle);
    _gTriangle->setTriangle(this);
#ifdef MECHANICS
    // eqArea cannot be obtained at this moment
    _mTriangle = unique_ptr<MTriangle>(new MTriangle(getType()));
    _mTriangle->setTriangle(this);
#endif

}

void Triangle::updateCoordinate() {
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        coordinate[coordIdx] = (_v[0]->coordinate[coordIdx] + _v[1]->coordinate[coordIdx] + _v[2]->coordinate[coordIdx]) / 3;
    }
}

void Triangle::updatePosition() {
    updateCoordinate();
    
    // Get the current compartment
    Compartment *c;
    try { c = GController::getCompartment(mathfunc::array2Vector(coordinate)); }
    catch (exception& e) {
        cout << e.what();
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

void Triangle::printSelf() {
    
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
    
    cout << "Vertex information..." << endl;
    
    _v[0]->printSelf();
    _v[1]->printSelf();
    _v[2]->printSelf();
    
    cout << endl;
}
