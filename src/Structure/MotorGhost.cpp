
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include <cmath>

#include "MotorGhost.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "ChemRNode.h"

#include "GController.h"
#include "SysParams.h"
#include "MathFunctions.h"
#include "Rand.h"

using namespace mathfunc;

void MotorGhost::updateCoordinate() {
    
    auto x1 = _c1->getFirstBead()->coordinate;
    auto x2 = _c1->getSecondBead()->coordinate;
    auto x3 = _c2->getFirstBead()->coordinate;
    auto x4 = _c2->getSecondBead()->coordinate;
    
    auto m1 = midPointCoordinate(x1, x2, _position1);
    auto m2 = midPointCoordinate(x3, x4, _position2);
    
    coordinate = midPointCoordinate(m1, m2, 0.5);
}


MotorGhost::MotorGhost(Cylinder* c1, Cylinder* c2, short motorType,
                       double position1, double position2)

    : Trackable(true, true),
      _c1(c1), _c2(c2),
      _position1(position1), _position2(position2),
      _motorType(motorType), _motorID(_motorGhosts.getID()),
      _birthTime(tau()) {
          
    //find compartment
    updateCoordinate();
    
    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }
    short filType = c1->getType();
          
    int pos1 = int(position1 * SysParams::Geometry().cylinderNumMon[filType]);
    int pos2 = int(position2 * SysParams::Geometry().cylinderNumMon[filType]);
          
    //set number of heads by picking random int between maxheads and minheads
    _numHeads = Rand::randInteger(SysParams::Chemistry().motorNumHeadsMin[_motorType],
                                  SysParams::Chemistry().motorNumHeadsMax[_motorType]);
    
#ifdef CHEMISTRY
    _cMotorGhost = unique_ptr<CMotorGhost>(
    new CMotorGhost(motorType, _compartment, _c1->getCCylinder(), _c2->getCCylinder(), pos1, pos2));
    _cMotorGhost->setMotorGhost(this);
#endif
    
#ifdef MECHANICS
    auto x1 = _c1->getFirstBead()->coordinate;
    auto x2 = _c1->getSecondBead()->coordinate;
    auto x3 = _c2->getFirstBead()->coordinate;
    auto x4 = _c2->getSecondBead()->coordinate;
          
    _mMotorGhost = unique_ptr<MMotorGhost>(
    new MMotorGhost(motorType, _numHeads, position1, position2, x1, x2, x3, x4));
    _mMotorGhost->setMotorGhost(this);
#endif
    
}

///@note - nothing for now, but could record data here
MotorGhost::~MotorGhost() noexcept {}

void MotorGhost::updatePosition() {
    
#ifdef CHEMISTRY
    //update ccylinders
    _cMotorGhost->setFirstCCylinder(_c1->getCCylinder());
    _cMotorGhost->setSecondCCylinder(_c2->getCCylinder());
    
#endif
    //check if in same compartment
    updateCoordinate();
    
    Compartment* c;
    
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }
    
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
    
#ifdef MECHANICS
    auto x1 = _c1->getFirstBead()->coordinate;
    auto x2 = _c1->getSecondBead()->coordinate;
    auto x3 = _c2->getFirstBead()->coordinate;
    auto x4 = _c2->getSecondBead()->coordinate;
    
    auto m1 = midPointCoordinate(x1, x2, _position1);
    auto m2 = midPointCoordinate(x3, x4, _position2);
    
    _mMotorGhost->setLength(twoPointDistance(m1, m2));
#endif
    
}

/// @note - This function updates forward walking rates using the
/// stetching force in the opposite direction of the motor walk.
/// Does not consider negative forces in this direction.

/// Updates unbinding rates based on the stretch force. Does not
/// consider compression forces, only stretching.
void MotorGhost::updateReactionRates() {

    //current force
    double force = max(0.0, _mMotorGhost->stretchForce);
    
    //walking rate changer
    if(!_walkingChangers.empty()) {
        
        auto x1 = _c1->getFirstBead()->coordinate;
        auto x2 = _c1->getSecondBead()->coordinate;
        auto x3 = _c2->getFirstBead()->coordinate;
        auto x4 = _c2->getSecondBead()->coordinate;
        
        //get component of force in direction of forward walk for C1, C2
        vector<double> motorC1Direction =
        twoPointDirection(midPointCoordinate(x1, x2, _position1),
                          midPointCoordinate(x3, x4, _position2));
        
        vector<double> motorC2Direction =
        twoPointDirection(midPointCoordinate(x3, x4, _position2),
                          midPointCoordinate(x1, x2, _position1));
        
        vector<double> c1Direction = twoPointDirection(x2,x1);
        vector<double> c2Direction = twoPointDirection(x4,x3);
        
        double forceDotDirectionC1 = force * dotProduct(motorC1Direction, c1Direction);
        double forceDotDirectionC2 = force * dotProduct(motorC2Direction, c2Direction);
        
        //WALKING REACTIONS
        Species* s1 = _cMotorGhost->getFirstSpecies();
        Species* s2 = _cMotorGhost->getSecondSpecies();
        
        for(auto r : s1->getRSpecies().reactantReactions()) {
            
            if(r->getReactionType() == ReactionType::MOTORWALKINGFORWARD) {
                
                float newRate =
                _walkingChangers[_motorType]->
                changeRate(_cMotorGhost->getOnRate(),
                           _cMotorGhost->getOffRate(),
                           _numHeads, max(0.0, forceDotDirectionC1));
                if(SysParams::RUNSTATE==false){
                    newRate=0.0;}
                r->setRate(newRate);
                r->updatePropensity();

            }
            else if(r->getReactionType() == ReactionType::MOTORWALKINGBACKWARD) {
                float newRate =
                _walkingChangers[_motorType]->
                changeRate(_cMotorGhost->getOnRate(),
                           _cMotorGhost->getOffRate(),
                           _numHeads, max(0.0, -forceDotDirectionC1));
                
                if(SysParams::RUNSTATE==false){
                    newRate=0.0;}
                r->setRate(newRate);

                r->updatePropensity();
            }
        }
        for(auto r : s2->getRSpecies().reactantReactions()) {
            
            if(r->getReactionType() == ReactionType::MOTORWALKINGFORWARD) {
                
                float newRate =
                _walkingChangers[_motorType]->
                changeRate(_cMotorGhost->getOnRate(),
                           _cMotorGhost->getOffRate(),
                           _numHeads, max(0.0, forceDotDirectionC2));
                if(SysParams::RUNSTATE==false)
                { newRate=0.0;}
                
                r->setRate(newRate);
                r->updatePropensity();
            }
            else if(r->getReactionType() == ReactionType::MOTORWALKINGFORWARD) {
                
                float newRate =
                _walkingChangers[_motorType]->
                changeRate(_cMotorGhost->getOnRate(),
                           _cMotorGhost->getOffRate(),
                           _numHeads, max(0.0, -forceDotDirectionC2));
                if(SysParams::RUNSTATE==false)
                { newRate=0.0;}
                r->setRate(newRate);
                r->updatePropensity();
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
        changeRate(_cMotorGhost->getOnRate(), _cMotorGhost->getOffRate(), _numHeads, force);
        
        offRxn->setRate(newRate);
        offRxn->activateReaction();
    }
}

void MotorGhost::moveMotorHead(Cylinder* c,
                               double oldPosition, double newPosition,
                               short boundType, SubSystem* ps) {
    
    //shift the position of one side of the motor
    double shift =  newPosition - oldPosition;
    
    //shift the head
    if(c == _c1) {
        _position1 += shift;
    }
    else {
        _position2 += shift;
    }
    short filType = c->getType();
    
    //record walk length
    _walkLength += shift * SysParams::Geometry().cylinderSize[filType];
    
#ifdef CHEMISTRY
    short oldpos = int (oldPosition * SysParams::Geometry().cylinderNumMon[filType]);
    short newpos = int (newPosition * SysParams::Geometry().cylinderNumMon[filType]);
    
    _cMotorGhost->moveMotorHead(c->getCCylinder(), oldpos, newpos,
                                _motorType, boundType, ps);
#endif
    
}

void MotorGhost::moveMotorHead(Cylinder* oldC, Cylinder* newC,
                               double oldPosition, double newPosition,
                               short boundType, SubSystem* ps) {
    
    //shift the head
    if(oldC == _c1) {
        _position1 = newPosition;
        _c1 = newC;
    }
    else {
        _position2 = newPosition;
        _c2 = newC;
    }
    short filType = _c1->getType();
    
    //record walk length
    _walkLength += (1-oldPosition + newPosition) * SysParams::Geometry().cylinderSize[filType];
    
#ifdef CHEMISTRY
    short oldpos = int (oldPosition * SysParams::Geometry().cylinderNumMon[filType]);
    short newpos = int (newPosition * SysParams::Geometry().cylinderNumMon[filType]);
    
    _cMotorGhost->moveMotorHead(oldC->getCCylinder(), newC->getCCylinder(),
                                oldpos, newpos, _motorType, boundType, ps);
#endif
}


void MotorGhost::printSelf() {
    
    cout << endl;
    
    cout << "MotorGhost: ptr = " << this << endl;
    cout << "Motor type = " << _motorType << ", Motor ID = " << _motorID << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    cout << "Position on first cylinder (double) = " << _position1 << endl;
    cout << "Position on second cylinder (double) = " << _position2 << endl;
    
    cout << "Number of heads = " << _numHeads << endl;
    cout << "Birth time = " << _birthTime << endl;
    
    cout << endl;
    
#ifdef CHEMISTRY
    cout << "Associated species 1 = " << _cMotorGhost->getFirstSpecies()->getName()
    << " , copy number = " << _cMotorGhost->getFirstSpecies()->getN()
    << " , position on first cylinder (int) = " << _cMotorGhost->getFirstPosition() << endl;
    
    cout << "Associated species 2 = " << _cMotorGhost->getSecondSpecies()->getName()
    << " , copy number = " << _cMotorGhost->getSecondSpecies()->getN()
    << " , position on second cylinder (int) = " << _cMotorGhost->getSecondPosition() << endl;
#endif
    
    cout << endl;
    
    cout << "Associated cylinders (one and two): " << endl;
    _c1->printSelf();
    _c2->printSelf();
    
    cout << endl;
}

species_copy_t MotorGhost::countSpecies(const string& name) {
    
    species_copy_t copyNum = 0;
    
    for(auto m : _motorGhosts.getElements()) {
        
        auto s = m->getCMotorGhost()->getFirstSpecies();
        string sname = SpeciesNamesDB::removeUniqueFilName(s->getName());
        
        if(sname == name)
            copyNum += s->getN();
    }
    return copyNum;
}

vector<MotorRateChanger*> MotorGhost::_unbindingChangers;
vector<MotorRateChanger*> MotorGhost::_walkingChangers;

Database<MotorGhost*> MotorGhost::_motorGhosts;
Histogram* MotorGhost::_lifetimes;
Histogram* MotorGhost::_walkLengths;

