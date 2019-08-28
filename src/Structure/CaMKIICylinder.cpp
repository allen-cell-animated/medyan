
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "Cylinder.h"
#include "CaMKIICylinder.h"

#include "SubSystem.h"
#include "CController.h"
#include "ChemManager.h"
#include "ChemRNode.h"

#include "Filament.h"
#include "Bead.h"

#include "GController.h"
#include "MathFunctions.h"

using namespace mathfunc;

CaMKIICylinder::CaMKIICylinder(CaMKIIingPoint *camkiiPoint, Bead* b1, short type, int position):
Cylinder(nullptr, b1, b1, type, position, false, false, false), _camkiiPoint(camkiiPoint){
    _camkiiPoint = camkiiPoint;
    for(auto &manager : _compartment->getFilamentBindingManagers())
      manager->addPossibleBindings(_cCylinder.get());
};

CaMKIICylinder::~CaMKIICylinder() noexcept {
    for(auto &manager : _compartment->getFilamentBindingManagers())
      manager->removePossibleBindings(_cCylinder.get());
    //remove from compartment
    _compartment->removeCylinder(this);

}

void CaMKIICylinder::updateCoordinate(){
  coordinate = _b1->coordinate;
}

//TODO override this for CaMKII
void CaMKIICylinder::updatePosition() {

    //check if were still in same compartment, set new position
    updateCoordinate();

    Compartment* c;
    
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }
    
    if(c != _compartment) {

#ifdef CHEMISTRY
        auto oldCompartment = _compartment;
        auto newCompartment = c;
#endif
        
        //remove from old compartment, add to new
        _compartment->removeCylinder(this);
        _compartment = c;
        _compartment->addCylinder(this);
        
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
    }
#endif
    
#ifdef MECHANICS
    //update length
    _mCylinder->setLength(twoPointDistance(this->_b1->coordinate,
                                           this->_b2->coordinate));
#endif

}

/// @note -  The function uses the bead load force to calculate this changed rate.
/// If there is no force on the beads the reaction rates are set to the bare.

//void Cylinder::updateReactionRates() {
//
//    double force;
//    //if no rate changer was defined, skip
//    if(_polyChanger.empty()) return;
//
//    //load force from front (affects plus end polymerization)
//    if(_plusEnd) {
//
//        //get force of front bead
//        force = _b2->getLoadForcesP();
//
//        //change all plus end polymerization rates
//        for(auto &r : _cCylinder->getInternalReactions()) {
//
//            if(r->getReactionType() == ReactionType::POLYMERIZATIONPLUSEND) {
//
//                float newRate = _polyChanger[_type]->changeRate(r->getBareRate(), force);
//
//                r->setRate(newRate);
//                r->updatePropensity();
//            }
//        }
//    }
//
//    //load force from back (affects minus end polymerization)
//    if(_minusEnd) {
//
//        //get force of front bead
//        force = _b1->getLoadForcesM();
//
//        //change all plus end polymerization rates
//        for(auto &r : _cCylinder->getInternalReactions()) {
//
//            if(r->getReactionType() == ReactionType::POLYMERIZATIONMINUSEND) {
//
//                float newRate =  _polyChanger[_type]->changeRate(r->getBareRate(), force);
//
//                r->setRate(newRate);
//                r->updatePropensity();
//            }
//        }
//    }
//}

//bool Cylinder::isFullLength() {
//
//#ifdef MECHANICS
//    return areEqual(_mCylinder->getEqLength(), SysParams::Geometry().cylinderSize[_type]);
//#else
//    return true;
//#endif
//}

//void Cylinder::printSelf() {
//
//    cout << endl;
//
//    cout << "Cylinder: ptr = " << this << endl;
//    cout << "Cylinder ID = " << _ID << endl;
//    cout << "Parent ptr = " << getParent() << endl;
//    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
//
//    if(_plusEnd) cout << "Is a plus end." << endl;
//    if(_minusEnd) cout << "Is a minus end." << endl;
//
//    if(_branchingCylinder != nullptr) cout << "Has a branching cylinder." << endl;
//    if(_camkiiingCylinder != nullptr) cout << "Has a camkiiing cylinder." << endl;
//
//    cout << "Position = " << _position << endl;
//
//    cout << endl;
//
//#ifdef CHEMISTRY
//    cout << "Chemical composition of cylinder:" << endl;
//    _cCylinder->printCCylinder();
//#endif
//
//    cout << endl;
//
//    cout << "Bead information..." << endl;
//
//    _b1->printSelf();
//    _b2->printSelf();
//
//    cout << endl;
//}

bool CaMKIICylinder::within(Cylinder* other, double dist) {

    //check midpoints
    if(twoPointDistance(coordinate, other->coordinate) <= dist)
        return true;

    //briefly check endpoints of other
    if(twoPointDistance(coordinate, other->_b1->coordinate) <= dist ||
       twoPointDistance(coordinate, other->_b2->coordinate) <= dist)
        return true;

    return false;
}

//vector<FilamentRateChanger*> CaMKIICylinder::_polyChanger;
//ChemManager* CaMKIICylinder::_chemManager = 0;

//Database<CaMKIICylinder*> CaMKIICylinder::_cylinders;
