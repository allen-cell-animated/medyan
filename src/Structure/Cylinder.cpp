
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

#include "SubSystem.h"
#include "CController.h"
#include "ChemManager.h"
#include "ChemRNode.h"

#include "Filament.h"
#include "Bead.h"

#include "GController.h"
#include "MathFunctions.h"

using namespace mathfunc;

void Cylinder::updateCoordinate() {
    
    coordinate = midPointCoordinate(_b1->coordinate, _b2->coordinate, 0.5);
}


Cylinder::Cylinder(Composite* parent, Bead* b1, Bead* b2, short type, int position,
                   bool extensionFront, bool extensionBack, bool initialization)

    : Trackable(true, true, true, false),
      _b1(b1), _b2(b2), _type(type), _position(position), _ID(_cylinders.getID()) {
    
    parent->addChild(unique_ptr<Component>(this));
          
    //Set coordinate
    updateCoordinate();

    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what() << endl;

//        cout << "Creating Cylinder: ptr = " << this << endl;
        
        exit(EXIT_FAILURE);
    }
                   
   //add to compartment
   _compartment->addCylinder(this);
          
#ifdef MECHANICS
          //set eqLength according to cylinder size
          
              double eqLength  = twoPointDistance(b1->coordinate, b2->coordinate);
          if(!SysParams::RUNSTATE) //RESTARTPHASE
          {
              int nummonomers = (int) round(eqLength/ SysParams::Geometry().monomerSize[type]);
              double tpd = eqLength;
              
              if(nummonomers ==0){
                  eqLength = SysParams::Geometry().monomerSize[type];
              }
              else{
                  eqLength = (nummonomers) * SysParams::Geometry().monomerSize[type];
                  double mindis = abs(tpd - eqLength);

                  for(auto i=nummonomers-1;i<=min(nummonomers+1, SysParams::Geometry().cylinderNumMon[type]);i++){
                      if(mindis > abs(tpd - i * SysParams::Geometry().monomerSize[type]))
                      {
                          eqLength = i * SysParams::Geometry().monomerSize[type];
                          mindis = abs(tpd - eqLength);
                      }
                  }
              }
              
            
          }
          _mCylinder = unique_ptr<MCylinder>(new MCylinder(_type, eqLength));
          _mCylinder->setCylinder(this);
#endif
    
#ifdef CHEMISTRY
    _cCylinder = unique_ptr<CCylinder>(new CCylinder(_compartment, this));
    _cCylinder->setCylinder(this);
          
    //init using chem manager
    _chemManager->initializeCCylinder(_cCylinder.get(), extensionFront,
                                      extensionBack, initialization);
#endif

        
}

Cylinder::~Cylinder() noexcept {
//    std::cout<<"Removing cylinder from Cmpt "<<this<<endl;
    //remove from compartment
    _compartment->removeCylinder(this);
    
}

/// Get filament type
int Cylinder::getType() {return _type;}

void Cylinder::updatePosition() {

    //check if were still in same compartment, set new position
    updateCoordinate();

    Compartment* c;
//    std::cout<<coordinate[0]<<" "<<coordinate[1]<<" "<<coordinate[2]<<endl;
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
    _mCylinder->setLength(twoPointDistance(_b1->coordinate,
                                           _b2->coordinate));
#endif

}

/// @note -  The function uses the bead load force to calculate this changed rate.
/// If there is no force on the beads the reaction rates are set to the bare.

void Cylinder::updateReactionRates() {
    
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

bool Cylinder::isFullLength() {
    
#ifdef MECHANICS
    return areEqual(_mCylinder->getEqLength(), SysParams::Geometry().cylinderSize[_type]);
#else
    return true;
#endif
}

void Cylinder::printSelf() {
    
    cout << endl;
    
    cout << "Cylinder: ptr = " << this << endl;
    cout << "Cylinder ID = " << _ID << endl;
    cout << "Parent ptr = " << getParent() << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    if(_plusEnd) cout << "Is a plus end." << endl;
    if(_minusEnd) cout << "Is a minus end." << endl;
    
    if(_branchingCylinder != nullptr) cout << "Has a branching cylinder." << endl;
    
    cout << "Position = " << _position << endl;
    
    cout << endl;
    
#ifdef CHEMISTRY
    cout << "Chemical composition of cylinder:" << endl;
    _cCylinder->printCCylinder();
#endif
    
    cout << endl;
    
    cout << "Bead information..." << endl;
    
    _b1->printSelf();
    _b2->printSelf();
    
    cout << endl;
}

bool Cylinder::within(Cylinder* other, double dist) {
    
    //check midpoints
    if(twoPointDistance(coordinate, other->coordinate) <= dist)
        return true;
    
    //briefly check endpoints of other
    if(twoPointDistance(coordinate, other->_b1->coordinate) <= dist ||
       twoPointDistance(coordinate, other->_b2->coordinate) <= dist)
        return true;
    
    return false;
}

vector<FilamentRateChanger*> Cylinder::_polyChanger;
ChemManager* Cylinder::_chemManager = 0;

Database<Cylinder*> Cylinder::_cylinders;
