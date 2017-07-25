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

// TODO how to implement Camkii::updateCoordinate?
void Camkii::updateCoordinate() {
    coordinate = midPointCoordinate(_b1->coordinate, _b2->coordinate, 0.5);
}

Camkii::Camkii(Composite* parent, int position, bool initialization = false)
    : Trackable(true, true, true, false),
      _cylinders(), _position(position), _ID(_camkiiDB.getID()) {
    
    // TODO init cylinders

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
                   
   // TODO implement addCamkii                
   //add to compartment TODO add camkii or all the cylinders indivdually? 
   _compartment->addCamkii(this);
    
#ifdef CHEMISTRY
    _cCamkii = unique_ptr<CCamkii>(new CCamkii(_compartment, this));
    _cCamkii->setCamkii(this);
          
    //init using chem manager
    _chemManager->initializeCCamkii(_cCamkii.get(), extensionFront,
                                      extensionBack, initialization);
#endif

#ifdef MECHANICS
    //set eqLength according to cylinder size
    double eqLength  = twoPointDistance(b1->coordinate, b2->coordinate);
        
    _mCamkii = unique_ptr<MCylinder>(new MCamkii(_type, eqLength));
    _mCamkii->setCamkii(this);
#endif
        
}

Camkii::~Camkii() {
    
    //remove from compartment
    _compartment->removeCamkii(this);
    // TODO free all the cylinders
    
}

void Camkii::updatePosition() {

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
        _compartment->removeCamkii(this);
        _compartment = c;
        _compartment->addCamkii(this);
        
#ifdef CHEMISTRY
        auto oldCCylinder = _cCamkii.get();
        
        //TODO camkii
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
    // TODO what to do here
    _mCylinder->setLength(twoPointDistance(_b1->coordinate,
                                           _b2->coordinate));
#endif

}

/// @note -  The function uses the bead load force to calculate this changed rate.
/// If there is no force on the beads the reaction rates are set to the bare.

void Camkii::updateReactionRates() {
    // TODO what to do here???
    double force;
    
    //if no rate changer was defined, skip
    if(_polyChanger.empty()) return;
    
    //load force from front (affects plus end polymerization)
    if(_plusEnd) {
        
        //get force of front bead
        force = _b2->getLoadForcesP();
        
        //change all plus end polymerization rates
        for(auto &r : _cCylinder->getInternalReactions()) {
            
            if(r->getReactionType() == ReactionType::POLYMERIZATIONPLUSEND) {
            
                float newRate = _polyChanger[_type]->changeRate(r->getBareRate(), force);
                
                r->setRate(newRate);
                r->updatePropensity();
            }
        }
    }
    
    //load force from back (affects minus end polymerization)
    if(_minusEnd) {
        
        //get force of front bead
        force = _b1->getLoadForcesM();
        
        //change all plus end polymerization rates
        for(auto &r : _cCylinder->getInternalReactions()) {
            
            if(r->getReactionType() == ReactionType::POLYMERIZATIONMINUSEND) {
                
                float newRate =  _polyChanger[_type]->changeRate(r->getBareRate(), force);
                
                r->setRate(newRate);
                r->updatePropensity();
            }
        }
    }
}

// TODO what's that
bool Camkii::isFullLength() {
    
#ifdef MECHANICS
    return areEqual(_mCylinder->getEqLength(), SysParams::Geometry().cylinderSize[_type]);
#else
    return true;
#endif
}

void Camkii::printSelf() {
    
    cout << endl;
    
    cout << "Camkii: ptr = " << this << endl;
    cout << "Camkii ID = " << _ID << endl;
    cout << "Parent ptr = " << getParent() << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    cout << "Position = " << _position << endl;
    
    cout << endl;
    
#ifdef CHEMISTRY
    cout << "Chemical composition of cylinder:" << endl;
    _cCylinder->printCCamkii();
#endif
    
    cout << endl;

    cout << "Print cylinders..." << endl;
    for (const auto& c: _cylinders)
        c->printSelf();

    cout << "Bead information..." << endl;
    
    _b1->printSelf();
    _b2->printSelf();
    
    cout << endl;
}

bool Camkii::within(Filament* other, double dist) {
    
    // TODO
    //check midpoints
    if(twoPointDistance(coordinate, other->coordinate) <= dist)
        return true;
    
    //briefly check endpoints of other
    if(twoPointDistance(coordinate, other->_b1->coordinate) <= dist ||
       twoPointDistance(coordinate, other->_b2->coordinate) <= dist)
        return true;
    
    return false;
}

// TODO
bool Camkii::isConsistent() {
    return false;
}

// TODO ????
vector<FilamentRateChanger*> Cylinder::_polyChanger;
ChemManager* Camkii::_chemManager = 0;
Database<Cylinder*> Cylinder::_cylinders;
#endif //CAMKII