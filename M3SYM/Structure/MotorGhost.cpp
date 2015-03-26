
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

#include <cmath>

#include "MotorGhost.h"

#include "Bead.h"
#include "Cylinder.h"
#include "ChemRNode.h"

#include "GController.h"
#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

vector<MotorRateChanger*> MotorGhost::_unbindingChangers;
vector<MotorRateChanger*> MotorGhost::_walkingChangers;

MotorGhost::MotorGhost(Cylinder* c1, Cylinder* c2, short motorType,
                       double position1, double position2, bool creation)
    : _c1(c1), _c2(c2),
      _position1(position1), _position2(position2), _motorType(motorType) {
    
    //add to motor ghost db
    MotorGhostDB::instance()->addMotorGhost(this);
          
    _motorID = MotorGhostDB::instance()->getMotorID();
    _birthTime = tau();
          
    auto c1b1 = _c1->getFirstBead()->coordinate;
    auto c1b2 = _c1->getSecondBead()->coordinate;
    auto c2b1 = _c2->getFirstBead()->coordinate;
    auto c2b2 = _c2->getSecondBead()->coordinate;
          
    //Find compartment
    auto m1 = midPointCoordinate(c1b1, c1b2, _position1);
    auto m2 = midPointCoordinate(c2b1, c2b2, _position2);
          
    coordinate = midPointCoordinate(m1, m2, 0.5);
    
    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
          
    int pos1 = int(position1 * SysParams::Geometry().cylinderIntSize);
    int pos2 = int(position2 * SysParams::Geometry().cylinderIntSize);
          
    //set number of heads by picking random int between maxheads and minheads
    _numHeads = (int) randomDouble(
        SysParams::Chemistry().motorNumHeadsMin[_motorType],
        SysParams::Chemistry().motorNumHeadsMax[_motorType]);
    
#ifdef CHEMISTRY
    _cMotorGhost = unique_ptr<CMotorGhost>(
        new CMotorGhost(motorType, _compartment,
                        _c1->getCCylinder(), _c2->getCCylinder(), pos1, pos2));
    _cMotorGhost->setMotorGhost(this);
#endif
    
#ifdef MECHANICS
    _mMotorGhost = unique_ptr<MMotorGhost>(
        new MMotorGhost(motorType, _numHeads,
                        position1, position2,
                        c1b1, c1b2, c2b1, c2b2));
    _mMotorGhost->setMotorGhost(this);
#endif
    
}

MotorGhost::~MotorGhost() noexcept {
    
    //remove from motor ghost db
    MotorGhostDB::instance()->removeMotorGhost(this);
}

void MotorGhost::updatePosition() {
    
#ifdef CHEMISTRY
    //update ccylinders
    _cMotorGhost->setFirstCCylinder(_c1->getCCylinder());
    _cMotorGhost->setSecondCCylinder(_c2->getCCylinder());
    
#endif
    
    //check if were still in same compartment
    auto c1b1 = _c1->getFirstBead()->coordinate;
    auto c1b2 = _c1->getSecondBead()->coordinate;
    auto c2b1 = _c2->getFirstBead()->coordinate;
    auto c2b2 = _c2->getSecondBead()->coordinate;
    
    auto m1 = midPointCoordinate(c1b1, c1b2, _position1);
    auto m2 = midPointCoordinate(c2b1, c2b2, _position2);
    
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

/// @note - This function updates forward walking rates using the
/// stetching force in the opposite direction of the forward walk.
/// Does not consider negative forces in this direction.

/// Updates unbinding rates based on the stretch force. Does not
/// consider compression forces, only stretching.
void MotorGhost::updateReactionRates() {

    //current force
    double force = max(0.0, _mMotorGhost->stretchForce);
    
    //walking rate changer
    if(!_walkingChangers.empty()) {
        
        auto c1b1 = _c1->getFirstBead()->coordinate;
        auto c1b2 = _c1->getSecondBead()->coordinate;
        auto c2b1 = _c2->getFirstBead()->coordinate;
        auto c2b2 = _c2->getSecondBead()->coordinate;
        
        //get component of force in direction of forward walk for C1, C2
        vector<double> motorC1Direction =
        twoPointDirection(midPointCoordinate(c1b1, c1b2, _position1),
                          midPointCoordinate(c2b1, c2b2, _position2));
        
        vector<double> motorC2Direction =
        twoPointDirection(midPointCoordinate(c2b1, c2b2, _position2),
                          midPointCoordinate(c1b1, c1b2, _position1));
        
        vector<double> c1Direction = twoPointDirection(c1b2,c1b1);
        vector<double> c2Direction = twoPointDirection(c2b2,c2b1);
        
        double forceDotDirectionC1 =
        max(0.0, force * dotProduct(motorC1Direction, c1Direction));
        double forceDotDirectionC2 =
        max(0.0, force * dotProduct(motorC2Direction, c2Direction));
        
        //WALKING REACTIONS
        Species* s1 = _cMotorGhost->getFirstSpecies();
        Species* s2 = _cMotorGhost->getSecondSpecies();
        
        for(auto r : s1->getRSpecies().reactantReactions()) {
            
            if(r->getReactionType() == ReactionType::MOTORWALKINGFORWARD) {
                
                float newRate =
                    _walkingChangers[_motorType]->
                    changeRate(r->getBareRate(), _numHeads, forceDotDirectionC1);
                
                r->setRate(newRate);
                r->activateReaction();
            }
        }
        for(auto r : s2->getRSpecies().reactantReactions()) {
            
            if(r->getReactionType() == ReactionType::MOTORWALKINGFORWARD) {
                
                float newRate =
                    _walkingChangers[_motorType]->
                    changeRate(r->getBareRate(), _numHeads, forceDotDirectionC2);
                
                r->setRate(newRate);
                r->activateReaction();
            }
        }
    }
    
    //unbinding rate changer
    if(!_unbindingChangers.empty()) {
        
        //get the unbinding reaction
        ReactionBase* offRxn = _cMotorGhost->getOffReaction();
        
        //change the rate
        float newRate =
            _unbindingChangers[_motorType]->
            changeRate(offRxn->getBareRate(), _numHeads, force);
        
        offRxn->setRate(newRate);
        offRxn->activateReaction();
    }
    
    
}
