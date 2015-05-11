
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

#include "CMotorGhost.h"

#include "ChemCallbacks.h"
#include "CCylinder.h"

CMotorGhost::CMotorGhost(short motorType, Compartment* c,
                         CCylinder* cc1, CCylinder* cc2, int position1, int position2)

    : CBound(c, cc1, cc2, position1, position2) {
    
    //Find species on cylinder that should be marked
    SpeciesBound* sm1 = _cc1->getCMonomer(_position1)->speciesMotor(motorType);
    SpeciesBound* sm2 = _cc2->getCMonomer(_position2)->speciesMotor(motorType);

    SpeciesBound* se1 = _cc1->getCMonomer(_position1)->speciesBound(BOUND_EMPTY);
    SpeciesBound* se2 = _cc2->getCMonomer(_position2)->speciesBound(BOUND_EMPTY);
        
    //mark species
    assert(sm1->getN() == 0 && sm2->getN() == 0 &&
           se1->getN() == 1 && se2->getN() == 1 &&
           "Major bug: Motor binding to an occupied site.");
        
    sm1->up(); sm2->up();
    se1->down(); se2->down();
        
    //attach this motor to the species
    setFirstSpecies(sm1);
    setSecondSpecies(sm2);
}

CMotorGhost::~CMotorGhost() {
    
    //remove the off reaction
    _cc1->removeCrossCylinderReaction(_cc2, _offRxn);
    
}

void CMotorGhost::createOffReaction(ReactionBase* onRxn, SubSystem* ps) {
    
    RSpecies** rs = onRxn->rspecies();
    vector<Species*> os;
    
    //copy into offspecies vector
    os.push_back(_firstSpecies);
    os.push_back(_secondSpecies);
    
    os.push_back(&rs[1]->getSpecies());
    
    Species* empty1 = _cc1->getCMonomer(_position1)->speciesBound(BOUND_EMPTY);
    Species* empty2 = _cc2->getCMonomer(_position2)->speciesBound(BOUND_EMPTY);
    
    os.push_back(empty1);
    os.push_back(empty2);
    
    ReactionBase* offRxn =
    new Reaction<LMUNBINDINGREACTANTS,LMUNBINDINGPRODUCTS>(os, _offRate);
    offRxn->setReactionType(ReactionType::MOTORUNBINDING);
    
    //Attach the callback to the off reaction, add it
    MotorUnbindingCallback mcallback(_pMotorGhost, ps);
    boost::signals2::shared_connection_block rcb(offRxn->connect(mcallback,false));
    
    _cc1->addCrossCylinderReaction(_cc2, offRxn);
    setOffReaction(offRxn);
}

void CMotorGhost::moveMotorHead(CCylinder* cc,
                                short oldPosition,
                                short newPosition,
                                short motorType,
                                short boundType,
                                SubSystem* ps) {
    
    auto smOld = cc->getCMonomer(oldPosition)->speciesMotor(motorType);
    auto smNew = cc->getCMonomer(newPosition)->speciesMotor(motorType);

    auto seNew = cc->getCMonomer(newPosition)->speciesBound(boundType);
    
    ReactionBase* newOffRxn;
    
    if(getFirstSpecies() == smOld) {
        
        _position1 = newPosition;
        
        setFirstSpecies(smNew);
        
        //change off reaction to include new species
        Species* smOther = _secondSpecies;
        Species* seOther = _cc2->getCMonomer(_position2)->speciesBound(BOUND_EMPTY);
        
        Species* sbd = &(_offRxn->rspecies()[2]->getSpecies());
        
        //create new reaction
        newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                          ({smNew, smOther, sbd, seNew, seOther}, _offRate);
    }
    else {
        _position2 = newPosition;
        
        setSecondSpecies(smNew);
        
        //change off reaction to include new species
        Species* smOther = _firstSpecies;
        Species* seOther = _cc1->getCMonomer(_position1)->speciesBound(BOUND_EMPTY);
        
        Species* sbd = &(_offRxn->rspecies()[2]->getSpecies());
        
        //create new reaction
        newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                          ({smOther, smNew, sbd, seNew, seOther}, _offRate);
    }
    //set new reaction type
    newOffRxn->setReactionType(ReactionType::MOTORUNBINDING);
    
    //attach signal
    MotorUnbindingCallback mcallback(_pMotorGhost, ps);
    boost::signals2::shared_connection_block rcb(newOffRxn->connect(mcallback,false));
    
    //remove old reaction, add new one
    _cc1->removeCrossCylinderReaction(_cc2, _offRxn);
    _cc1->addCrossCylinderReaction(_cc2, newOffRxn);
    
    //set new unbinding reaction
    setOffReaction(newOffRxn);
    
}

void CMotorGhost::moveMotorHead(CCylinder* oldCC,
                                CCylinder* newCC,
                                short oldPosition,
                                short newPosition,
                                short motorType,
                                short boundType,
                                SubSystem* ps) {

    
    auto smOld = oldCC->getCMonomer(oldPosition)->speciesMotor(motorType);
    auto smNew = newCC->getCMonomer(newPosition)->speciesMotor(motorType);
    
    auto seNew = newCC->getCMonomer(newPosition)->speciesBound(boundType);
    
    ReactionBase* newOffRxn;
    
    if(getFirstSpecies() == smOld) {
        
        _position1 = newPosition;
        
        setFirstSpecies(smNew);
        
        //change off reaction to include new species
        Species* smOther = _secondSpecies;
        Species* seOther = _cc2->getCMonomer(_position2)->speciesBound(BOUND_EMPTY);
        
        Species* sbd = &(_offRxn->rspecies()[2]->getSpecies());
        
        //create new reaction
        newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                          ({smNew, smOther, sbd, seNew, seOther}, _offRate);
        
        //remove old off reaction
        _cc1->removeCrossCylinderReaction(_cc2, _offRxn);
        
        //set new ccylinders
        setFirstCCylinder(newCC);
    }
    else {
        _position2 = newPosition;
        
        setSecondSpecies(smNew);
        
        //change off reaction to include new species
        Species* smOther = _firstSpecies;
        Species* seOther = _cc1->getCMonomer(_position1)->speciesBound(BOUND_EMPTY);
        
        Species* sbd = &(_offRxn->rspecies()[2]->getSpecies());
        
        //create new reaction
        newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                         ({smOther, smNew, sbd, seNew, seOther}, _offRate);
        
        //remove old off reaction
        _cc1->removeCrossCylinderReaction(_cc2, _offRxn);
        
        //set new ccylinders
        setSecondCCylinder(newCC);
    }
    
    //set new reaction type
    newOffRxn->setReactionType(ReactionType::MOTORUNBINDING);
    
    //attach signal
    MotorUnbindingCallback mcallback(_pMotorGhost, ps);
    boost::signals2::shared_connection_block rcb(newOffRxn->connect(mcallback,false));

    //add new
    _cc1->addCrossCylinderReaction(_cc2, newOffRxn);
    
    //set new unbinding reaction
    setOffReaction(newOffRxn);
    
}