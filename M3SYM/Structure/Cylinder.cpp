
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "Cylinder.h"

#include "SubSystem.h"
#include "CController.h"
#include "ChemManager.h"
#include "ChemRNode.h"

#include "Bead.h"

#include "GController.h"
#include "MathFunctions.h"

using namespace mathfunc;

FilamentRateChanger* Cylinder::_polyChanger = 0;
ChemManager* Cylinder::_chemManager = 0;

Database<Cylinder*> Cylinder::_cylinders;

void Cylinder::updateCoordinate() {
    
    coordinate = midPointCoordinate(_b1->coordinate, _b2->coordinate, 0.5);
}


Cylinder::Cylinder(Filament* f, Bead* b1, Bead* b2, int positionFilament,
                   bool extensionFront,
                   bool extensionBack,
                   bool initialization)

    : Trackable(true, true, true, false),
      _b1(b1), _b2(b2), _pFilament(f),
      _positionFilament(positionFilament),
      _ID(_cylinders.getID()) {
    
    //Set coordinate
    updateCoordinate();

    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what() << endl;
        
        printInfo();
        
        exit(EXIT_FAILURE);
    }
                   
   //add to compartment
   _compartment->addCylinder(this);
    
#ifdef CHEMISTRY
    _cCylinder = unique_ptr<CCylinder>(new CCylinder(_compartment, this));
    _cCylinder->setCylinder(this);
          
    //init using chem manager
    _chemManager->initializeCCylinder(_cCylinder.get(),
                                      extensionFront,
                                      extensionBack,
                                      initialization);
#endif

#ifdef MECHANICS
    //set eqLength according to cylinder size
    double eqLength  = twoPointDistance(b1->coordinate, b2->coordinate);
        
    _mCylinder = unique_ptr<MCylinder>(new MCylinder(eqLength));
    _mCylinder->setCylinder(this);
#endif
        
}

Cylinder::~Cylinder() noexcept {
    
    //remove from compartment
    _compartment->removeCylinder(this);
}

void Cylinder::updatePosition() {

    //check if were still in same compartment, set new position
    updateCoordinate();

    Compartment* c;
    
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printInfo();
        
        exit(EXIT_FAILURE);
    }
    
    if(c != _compartment) {

        //remove from old compartment, add to new
        _compartment->removeCylinder(this);
        _compartment = c;
        _compartment->addCylinder(this);
        
#ifdef CHEMISTRY
        CCylinder* clone = _cCylinder->clone(c);
        setCCylinder(clone);
#endif
    }
    
#ifdef MECHANICS
    //update length
    _mCylinder->setLength(twoPointDistance(_b1->coordinate,
                                           _b2->coordinate));
#endif

}

/// @note -  The function uses the bead load force to calculate this changed rate.
/// If there is no force on the beads the reaction rates are set to the bare.

void Cylinder::updateReactionRates() {
    
    double force;
    
    //if no rate changer was defined, skip
    if(_polyChanger == nullptr) return;
    
    //load force from front (affects plus end polymerization)
    if(_plusEnd) {
        
        //get force of front bead
        force = _b2->loadForce;
        
        //change all plus end polymerization rates
        for(auto &r : _cCylinder->getInternalReactions()) {
            
            if(r->getReactionType() == ReactionType::POLYMERIZATIONPLUSEND) {
            
                float newRate = _polyChanger->changeRate(r->getBareRate(), force);
                
                r->setRate(newRate);
                r->updatePropensity();
            }
        }
    }
    
    //load force from back (affects minus end polymerization)
    else if(_minusEnd) {
        
        //get force of front bead
        force = _b1->loadForce;
        
        //change all plus end polymerization rates
        for(auto &r : _cCylinder->getInternalReactions()) {
            
            if(r->getReactionType() == ReactionType::POLYMERIZATIONMINUSEND) {
                
                float newRate =  _polyChanger->changeRate(r->getBareRate(), force);

                r->setRate(newRate);
                r->updatePropensity();
            }
        }
    }
}

bool Cylinder::isFullLength() {
    
#ifdef MECHANICS
    return _mCylinder->getEqLength() ==
            SysParams::Geometry().cylinderSize;
#else
    return true;
#endif
}

void Cylinder::printInfo() {
    
    cout << endl;
    
    cout << "Cylinder: ptr = " << this << endl;
    cout << "Cylinder ID = " << _ID << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    if(_plusEnd) cout << "Is a plus end." << endl;
    if(_minusEnd) cout << "Is a minus end." << endl;
    
    if(_branchingCylinder != nullptr) cout << "Has a branching cylinder." << endl;
    
    cout << "Position on filament = " << _positionFilament << endl;
    
    cout << endl;
    
#ifdef CHEMISTRY
    cout << "Chemical composition of cylinder:" << endl;
    _cCylinder->printCCylinder();
#endif
    
    cout << endl;
    
    cout << "Bead information..." << endl;
    
    _b1->printInfo();
    _b2->printInfo();
    
    cout << endl;
}

