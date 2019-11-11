
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

#ifdef CAMKII_ENABLED
#include "CaMKIIingPoint.h"
#endif

using namespace mathfunc;

void Cylinder::updateCoordinate() {
    coordinate = midPointCoordinate(_b1->coordinate, _b2->coordinate, 0.5);
}


Cylinder::Cylinder(Composite* parent, Bead* b1, Bead* b2, short type, int position,
                   bool extensionFront, bool extensionBack, bool initialization)

    : Trackable(true, true, true, false),
      _b1(b1), _b2(b2), _type(type), _position(position), _ID(_cylinders.getID()) {
    
	if (parent) {
		parent->addChild(unique_ptr<Component>(this));
		updateCoordinate();
	} else {
		coordinate = _b1->coordinate;
	}
    //Set coordinate

    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what() << endl;
        
        printSelf();
        
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
//              std::cout<<eqLength<<" ";
              
              if(nummonomers ==0){
                  eqLength = SysParams::Geometry().monomerSize[type];
              }
              else{
                  eqLength = (nummonomers) * SysParams::Geometry().monomerSize[type];
                  double mindis = abs(tpd - eqLength);
//                  std::cout<<eqLength<<" ";
                  for(auto i=nummonomers-1;i<=min(nummonomers+1, SysParams::Geometry().cylinderNumMon[type]);i++){
                      if(mindis > abs(tpd - i * SysParams::Geometry().monomerSize[type]))
                      {
                          eqLength = i * SysParams::Geometry().monomerSize[type];
                          mindis = abs(tpd - eqLength);
                      }
                  }
              }
              
              
//              for(auto i=nummonomers ;i<=min(nummonomers+1, SysParams::Geometry().cylinderNumMon[type]);i++){
//                  if(mindis > abs(tpd - i * SysParams::Geometry().monomerSize[type]))
//                  {
//                      eqLength = i * SysParams::Geometry().monomerSize[type];
//                      mindis = abs(tpd - eqLength);
//                  }
//              }
              
//              std::cout<<eqLength<<endl;
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
    
	//loop through camkii points
//	for camkii in
//	get bond
//	get cylinder
//	if cylinder==this
//	/error message
//	for b in _bonds

#if 0
    cerr << "=======================BEFORE REMOVE CYLINDER==============" << endl;
	for (auto camkii : CaMKIIingPoint::getCaMKIIingPoints()) {
        cerr << "=== CaMKII " << camkii->getID() << ": (Coord number = " << camkii->getCoordinationNumber() << ") " << endl;
		for(auto v : camkii->getBonds()) {
//			vector<tuple<Cylinder*, double>>
			Cylinder *c = get<0>(v);
			if(c == this) {
			  //TODO: Solve this error caused by cylinder missing after filament depolymerization
				cerr << "We find camkii cylinder in removal.\n";
				for(int i=0;i<10;i++) {
				    for(int j=0;j<6;j++) {
				        auto x = c->getCCylinder()->getCMonomer(i)->speciesBound(j);
                        cerr << "For monomer index " << i << " speciesBound index " << j << " : "
                             << x->getFullName() << "  N: " << x->getN() << endl;
                    }
					}
				exit(1);
			}
		}
	}
	cerr << "==========================================================" << endl;
#endif

    //remove from compartment
    _compartment->removeCylinder(this);

#if 0
	cerr << "=======================AFTER REMOVE CYLINDER==============" << endl;
	for (auto camkii : CaMKIIingPoint::getCaMKIIingPoints()) {
		cerr << "=== CaMKII " << camkii->getID() << ": (Coord number = " << camkii->getCoordinationNumber() << ") " << endl;
		for(auto v : camkii->getBonds()) {
//			vector<tuple<Cylinder*, double>>
			Cylinder *c = get<0>(v);
			if(c == this) {
				//TODO: Solve this error caused by cylinder missing after filament depolymerization
				cerr << "We find camkii cylinder in removal.\n";
				for(int i=0;i<10;i++) {
					for(int j=0;j<6;j++) {
						auto x = c->getCCylinder()->getCMonomer(i)->speciesBound(j);
						cerr << "For monomer index " << i << " speciesBound index " << j << " : "
							 << x->getFullName() << "  N: " << x->getN() << endl;
					}
				}
				exit(1);
			}
		}
	}
	cerr << "==========================================================" << endl;
#endif
}

/// Get filament type
int Cylinder::getType() {return _type;}

void Cylinder::updatePosition() {

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
        for(auto &manager : newCompartment->getFilamentBindingManagers()) {
            if(dynamic_cast<CaMKIIBundlingManager*>(manager.get()) != nullptr) {
                int count = Cylinder::getCylinders().size();
                int bef = manager->numBindingSites();
                manager->addPossibleBindings(newCCylinder);
                int aft = manager->numBindingSites();

                cerr << "========MILLAD: Cylinder::UpdatePosition  " << newCCylinder->getCylinder() << "  " << newCCylinder << "  " <<
                     (dynamic_cast<CaMKIICylinder *>(newCCylinder->getCylinder()) == NULL ? "normal" : "camkii_cylinder") << " count: "
                     << count << "  before: " << bef << "  after: " << aft << endl;
            } else {
                manager->addPossibleBindings(newCCylinder);
            }
        }
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
#ifdef CAMKII_ENABLED
    if(_camkiiingCylinder != nullptr) cout << "Has a camkiiing cylinder." << endl;
#endif

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
