
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "Cylinder.h"

#include "Bead.h"
#include "NeighborListDB.h"
#include "ChemManager.h"
#include "ChemRNode.h"

#include "GController.h"
#include "MathFunctions.h"

using namespace mathfunc;

Cylinder::Cylinder(Filament* f, Bead* b1, Bead* b2, int positionFilament,
                   bool extensionFront, bool extensionBack, bool creation)
    : _pFilament(f), _b1(b1), _b2(b2), _positionFilament(positionFilament) {
    
    //Add to cylinder DB
    CylinderDB::instance()->addCylinder(this);
    _ID = CylinderDB::instance()->getCylinderID();
                       
    //Set coordinate
    coordinate = midPointCoordinate(_b1->coordinate, _b2->coordinate, 0.5);

    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
                   
   //add to compartment
   _compartment->addCylinder(this);
    
#ifdef CHEMISTRY
    _cCylinder = unique_ptr<CCylinder>(new CCylinder(_compartment));
    _cCylinder->setCylinder(this);
    ChemManager::initializeCCylinder(
        _cCylinder.get(), f, extensionFront, extensionBack, creation);
    
    if(creation || extensionFront || extensionBack)
        //Update filament reactions, only if not initialization
        ChemManager::updateCCylinder(_cCylinder.get());
    
#endif

#ifdef MECHANICS
    double eqLength;
    
    //set eqLength according to cylinder size
    if(extensionFront || extensionBack)
        eqLength = SystemParameters::Geometry().monomerSize;
    else if(creation)
        eqLength = SystemParameters::Geometry().minCylinderSize;
    else
        eqLength = SystemParameters::Geometry().cylinderSize;
    
    _mCylinder = unique_ptr<MCylinder>(new MCylinder(eqLength));
    _mCylinder->setCylinder(this);
#endif
                       
   //add to neighbor list db
   NeighborListDB::instance()->addDynamicNeighbor(this);
}

Cylinder::~Cylinder() {
    
    //Remove from cylinder DB
    CylinderDB::instance()->removeCylinder(this);
    
    //remove from compartment
    _compartment->removeCylinder(this);
    
    //remove from neighbor lists
    NeighborListDB::instance()->removeDynamicNeighbor(this);
}

void Cylinder::updatePosition() {

    //check if were still in same compartment, set new position
    coordinate = midPointCoordinate(_b1->coordinate, _b2->coordinate, 0.5);
    
    //update length
    _mCylinder->setLength(twoPointDistance(_b1->coordinate, _b2->coordinate));
    
    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
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

}

/// @note - This function updates polymerization rates based on the
/// Elastic Brownian Ratchet Model:
///
///                 k = k_0 * exp(-f * a / kT)
///
/// The function uses the bead load force to calculate this changed rate.
/// If there is no force on the beads the reaction rates are set to the bare.

void Cylinder::updateReactionRates() {
    
    double force;
    
    //characteristic length
    float a = SystemParameters::DynamicRates().FDPLength;
    
    //load force from front (affects plus end polymerization)
    if(_plusEnd) {
        
        //get force of front bead
        force = _b2->loadForce;
        
        //change all plus end polymerization rates
        for(auto &r : _cCylinder->getInternalReactions()) {
            
            if(r->getReactionType() == ReactionType::POLYMERIZATIONPLUSEND) {
            
                float newRate = r->getBareRate() * exp( - force * a / kT);
                r->setRate(newRate);
                r->getRNode()->activateReaction();
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
                
                float newRate = r->getBareRate() * exp( - force * a / kT);
                r->setRate(newRate);
                r->getRNode()->activateReaction();
            }
        }
    }
}

