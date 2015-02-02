
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

#include "MotorGhost.h"

#include "Bead.h"
#include "Cylinder.h"
#include "ChemRNode.h"

#include "GController.h"
#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

MotorGhost::MotorGhost(Cylinder* c1, Cylinder* c2, short motorType,
                       double position1, double position2, bool creation)
    : _c1(c1), _c2(c2),
      _position1(position1), _position2(position2), _motorType(motorType) {
    
    //add to motor ghost db
    MotorGhostDB::instance()->addMotorGhost(this);
    _motorID = MotorGhostDB::instance()->getMotorID();
    _birthTime = tau();
    
    //Find compartment
    auto m1 =
        midPointCoordinate(_c1->getFirstBead()->coordinate,
                           _c1->getSecondBead()->coordinate, _position1);
    auto m2 =
        midPointCoordinate(_c2->getFirstBead()->coordinate,
                           _c2->getSecondBead()->coordinate, _position2);
          
    coordinate = midPointCoordinate(m1, m2, 0.5);
    
    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
#ifdef CHEMISTRY
    _cMotorGhost = unique_ptr<CMotorGhost>(new CMotorGhost(_compartment));
    _cMotorGhost->setMotorGhost(this);
    
    //Find species on cylinder that should be marked. If initialization,
    //this should be done. But, if this is because of a reaction callback,
    //it will have already been done.
    
    int pos1 = int(position1 * SystemParameters::Geometry().cylinderIntSize);
    int pos2 = int(position2 * SystemParameters::Geometry().cylinderIntSize);
    
    SpeciesMotor* sm1 = _c1->getCCylinder()->
                        getCMonomer(pos1)->speciesMotor(_motorType);
    SpeciesMotor* sm2 = _c2->getCCylinder()->
                        getCMonomer(pos2)->speciesMotor(_motorType);
    
    if(!creation) {
        SpeciesBound* se1 =
            _c1->getCCylinder()->getCMonomer(pos1)->speciesBound(0);
        sm1->up();
        se1->down();
        
        SpeciesBound* se2 =
            _c2->getCCylinder()->getCMonomer(pos2)->speciesBound(0);
        sm2->up();
        se2->down();
    }
    
    //attach this motor to the species
    _cMotorGhost->setFirstSpecies(sm1);
    _cMotorGhost->setSecondSpecies(sm2);
    
#endif
    
#ifdef MECHANICS
    _mMotorGhost = unique_ptr<MMotorGhost>(
        new MMotorGhost(motorType, position1, position2,
            _c1->getFirstBead()->coordinate, _c1->getSecondBead()->coordinate,
            _c2->getFirstBead()->coordinate, _c2->getSecondBead()->coordinate));
    _mMotorGhost->setMotorGhost(this);
#endif
    
}

MotorGhost::~MotorGhost() noexcept {
    
    //remove from motor ghost db
    MotorGhostDB::instance()->removeMotorGhost(this);
    
#ifdef CHEMISTRY
    //remove the unbinding reaction
    CCylinder* cc1 = _c1->getCCylinder();
    CCylinder* cc2 = _c2->getCCylinder();
    cc1->removeCrossCylinderReaction(cc2, _cMotorGhost->getOffReaction());
#endif
}

void MotorGhost::updatePosition() {
    
    //check if were still in same compartment
    auto m1 =
        midPointCoordinate(_c1->getFirstBead()->coordinate,
                           _c1->getSecondBead()->coordinate, _position1);
    auto m2 =
        midPointCoordinate(_c2->getFirstBead()->coordinate,
                           _c2->getSecondBead()->coordinate, _position2);
    
    coordinate = midPointCoordinate(m1, m2, 0.5);
    
    _mMotorGhost->setLength(twoPointDistance(m1, m2));
    
    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
    if(c != _compartment) {
        
        _compartment = c;
#ifdef CHEMISTRY
        SpeciesBound* firstSpecies = _cMotorGhost->getFirstSpecies();
        SpeciesBound* secondSpecies = _cMotorGhost->getSecondSpecies();
        
        CMotorGhost* clone = _cMotorGhost->clone(c);
        setCMotorGhost(clone);
        
        _cMotorGhost->setFirstSpecies(firstSpecies);
        _cMotorGhost->setSecondSpecies(secondSpecies);
#endif
    }
}

/// @note - This function updates walking rates based on the
/// following exponential form:
///
///                 k = k_0 * exp(-f * a / kT)
///
/// The function uses the motor's stretching force at the current state
/// to change this rate.

/// @note - This function updates unbinding rates based on the
/// following exponential form:
///
///                 k = k_0 * exp(f * a / kT)
///
/// The function uses the motor's stretching force at the current
/// state to change this rate.

void MotorGhost::updateReactionRates() {

    //current force
    double force = _mMotorGhost->stretchForce;
    
    //characteristic lengths
    double aUnbinding = SystemParameters::DynamicRates().MDULength[_motorType];
    double aWalking = SystemParameters::DynamicRates().MDWLength[_motorType];
    
    //get all walking reactions
    Species* s1 = _cMotorGhost->getFirstSpecies();
    Species* s2 = _cMotorGhost->getSecondSpecies();
    
    for(auto r : s1->getRSpecies().reactantReactions()) {
        
        if(r->getReactionType() == ReactionType::MOTORWALKINGFORWARD) {
            
            float newRate = r->getBareRate() * exp( - force * aWalking / kT);
            r->setRate(newRate);
            r->activateReaction();
        }
        if(r->getReactionType() == ReactionType::MOTORUNBINDING) {
            
            float newRate = r->getBareRate() * exp( force * aUnbinding / kT);
            r->setRate(newRate);
            r->activateReaction();
        }
    }
    for(auto r : s2->getRSpecies().reactantReactions()) {
        
        if(r->getReactionType() == ReactionType::MOTORWALKINGFORWARD) {
            
            float newRate = r->getBareRate() * exp( - force * aWalking / kT);
            r->setRate(newRate);
            r->activateReaction();
        }
        
        if(r->getReactionType() == ReactionType::MOTORUNBINDING) {
            
            float newRate = r->getBareRate() * exp( force * aUnbinding / kT);
            r->setRate(newRate);
            r->activateReaction();
        }
    }
}
