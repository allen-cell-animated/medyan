#ifdef CAMKII
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

#include "Camkii.h"

#include "SubSystem.h"
#include "CController.h"
#include "ChemManager.h"
#include "ChemRNode.h"

#include <Cylinder.h>
#include "Filament.h"
#include "Bead.h"

#include "GController.h"
#include "MathFunctions.h"

using namespace mathfunc;

Database<Camkii*> _camkiis; ///< Collection in SubSystem
void Camkii::updateCoordinate() {
    auto itr1 = _cylinders.begin();
    auto itr2 = _cylinders.begin() + 3;
    vector<double> avgMidpoint = {0.0, 0.0, 0.0};
    for(auto i=0; i<_cylinders.size()/2; i++){
        auto mp = midPointCoordinate(itr1[i]->coordinate, itr2[i]->coordinate, 0.5);
        avgMidpoint[0] += mp[0];
        avgMidpoint[1] += mp[1];
        avgMidpoint[2] += mp[2];
    }

    coordinate[0] = avgMidpoint[0]/3.0;
    coordinate[1] = avgMidpoint[1]/3.0;
    coordinate[2] = avgMidpoint[2]/3.0;
}

Camkii::Camkii(Composite* parent, int position, bool initialization)
    : Trackable(true, true, true, false), _cylinders(), _position(position), _ID(_camkiiDB.getID()) {
    // TODO init cylinders

    if (initialization)
        cout<<"dddd";
    // TODO do we need this?
    parent->addChild(unique_ptr<Component>(this));
          
    //Set coordinate
    updateCoordinate();

    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what() << endl;
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }
// TODO verify with James this isn't needed
//   _compartment->addCamkii(this);
//
//#ifdef CHEMISTRY
//    _cCamkii = unique_ptr<CCamkii>(new CCamkii(_compartment, this));
//    _cCamkii->setCamkii(this);
//
//    //init using chem manager
//    _chemManager->initializeCCamkii(_cCamkii.get(), extensionFront,
//                                      extensionBack, initialization);
//#endif
//
//#ifdef MECHANICS
//    //set eqLength according to cylinder size
//    double eqLength  = twoPointDistance(b1->coordinate, b2->coordinate);
//
//    _mCamkii = unique_ptr<MCylinder>(new MCamkii(_type, eqLength));
//    _mCamkii->setCamkii(this);
//#endif
        
}

Camkii::~Camkii() {
    //remove from compartment
    //_compartment->removeCamkii(this);
    // TODO verify with James we need this
    for (auto &c: _cylinders){
        c->getCompartment()->removeCylinder(c);
    }
}

// TODO verify that's all we need with James.
void Camkii::updatePosition() {

    //check if were still in same compartment, set new position
    updateCoordinate();

    // TODO check if I need to update the binding manager
//
//    Compartment* c;
//
//    try {c = GController::getCompartment(coordinate);}
//    catch (exception& e) {
//        cout << e.what();
//
//        printSelf();
//
//        exit(EXIT_FAILURE);
//    }
//
//    if(c != _compartment) {
//
//#ifdef CHEMISTRY
//        auto oldCompartment = _compartment;
//        auto newCompartment = c;
//#endif
//
//        //remove from old compartment, add to new
//        _compartment->removeCamkii(this);
//        _compartment = c;
//        _compartment->addCamkii(this);
//
//#ifdef CHEMISTRY
//        auto oldCCylinder = _cCamkii.get();
//
//        //TODO camkii
//        //Remove old ccylinder from binding managers
//        for(auto &manager : oldCompartment->getFilamentBindingManagers())
//            manager->removePossibleBindings(oldCCylinder);
//
//        //clone and set new ccylinder
//        CCylinder* clone = _cCylinder->clone(c);
//        setCCylinder(clone);
//
//        auto newCCylinder = _cCylinder.get();
//
//        //Add new ccylinder to binding managers
//        for(auto &manager : newCompartment->getFilamentBindingManagers())
//            manager->addPossibleBindings(newCCylinder);
//    }
//#endif
//
//#ifdef MECHANICS
//    //update length
//    // TODO what to do here
//    _mCylinder->setLength(twoPointDistance(_b1->coordinate,
//                                           _b2->coordinate));
//#endif

}

/// @note -  The function uses the bead load force to calculate this changed rate.
/// If there is no force on the beads the reaction rates are set to the bare.

void Camkii::updateReactionRates() {
    // TODO double check if should be left empty
}

void Camkii::printSelf() {
    
    cout << endl;
    
    cout << "Camkii: ptr = " << this << endl;
    cout << "Camkii ID = " << _ID << endl;
    cout << "Parent ptr = " << getParent() << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    cout << "Position = " << _position << endl;
    cout << endl;

    cout << "Print cylinders..." << endl;
    for (const auto& c: _cylinders)
        c->printSelf();

    cout << endl;
}

// TODO
bool Camkii::isConsistent() {
    return true;
}

// TODO ????
// ChemManager* Camkii::_chemManager = 0; TODO what about that


#endif //CAMKII