
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
};



CaMKIICylinder::~CaMKIICylinder() noexcept {

	removeFromFilamentBindingManagers();

    //remove from compartment
    _compartment->removeCylinder(this);

}

void CaMKIICylinder::addToFilamentBindingManagers() {
    for (auto &manager : _compartment->getFilamentBindingManagers())
        manager->addPossibleBindings(_cCylinder.get());
}

void CaMKIICylinder::removeFromFilamentBindingManagers() {
    for(auto &manager : _compartment->getFilamentBindingManagers())
        manager->removePossibleBindings(_cCylinder.get());
}

void CaMKIICylinder::updateCoordinate(){
  coordinate = _b1->coordinate;
}

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
        
        // Creating a copy of CaMKII cylinder in case the compartment is changed.
        CCylinder* clone = _cCylinder->clone(c);
        setCCylinder(clone);
        assert(getCaMKIIPointParent()->getCoordinationNumber() > 0);
        const auto cpoint = getCaMKIIPointParent()->getCCaMKIIingPoint();

        // Copy the correct off reaction depending on the coordination number
        if(getCaMKIIPointParent()->getCoordinationNumber() == 1L) { // unbinding
            cpoint->setOffRxnBinding(cpoint->getOffReaction());
            cpoint->setOffRxnBundling(nullptr);
        } else { // unbundling
            cpoint->setOffRxnBinding(nullptr);
            cpoint->setOffRxnBundling(cpoint->getOffReaction());
        }

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

