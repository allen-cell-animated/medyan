
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
                         CCylinder* cc1, CCylinder* cc2, int pos1, int pos2)

    : CBound(c, cc1, cc2) {
    
    //Find species on cylinder that should be marked
    SpeciesBound* sm1 = _cc1->getCMonomer(pos1)->speciesMotor(motorType);
    SpeciesBound* sm2 = _cc2->getCMonomer(pos2)->speciesMotor(motorType);
        
    //attach this motor to the species
    setFirstSpecies(sm1);
    setSecondSpecies(sm2);
}

CMotorGhost::~CMotorGhost() {
    
    //remove the off reaction
    _cc1->removeCrossCylinderReaction(_cc2, _offRxn);
    
}

void CMotorGhost::createOffReaction(ReactionBase* onRxn, float offRate, SubSystem* ps) {
    
    RSpecies** rs = onRxn->rspecies();
    vector<Species*> os;
    
    //copy into offspecies vector in opposite order
    for(int i = LMBINDINGREACTANTS; i < LMBINDINGREACTANTS+LMBINDINGPRODUCTS; i++)
        os.push_back(&rs[i]->getSpecies());
    
    for(int i = 0; i < LMBINDINGREACTANTS; i++)
        os.push_back(&rs[i]->getSpecies());
    
    ReactionBase* offRxn =
    new Reaction<LMUNBINDINGREACTANTS,LMUNBINDINGPRODUCTS>(os, offRate);
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
    
    auto sm1 = cc->getCMonomer(oldPosition)->speciesMotor(motorType);
    auto sm2 = cc->getCMonomer(newPosition)->speciesMotor(motorType);
    auto sb2 = cc->getCMonomer(newPosition)->speciesBound(boundType);
    
    ReactionBase* newOffRxn;
    
    if(getFirstSpecies() == sm1) {
        
        setFirstSpecies(sm2);
        
        //change off reaction to include new species
        
        //take the first, third, and fourth species
        Species* s1 = &(_offRxn->rspecies()[1]->getSpecies());
        Species* s3 = &(_offRxn->rspecies()[3]->getSpecies());
        Species* s4 = &(_offRxn->rspecies()[4]->getSpecies());
        
        //create new reaction
        newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                          ({sm2, s1, sb2, s3, s4}, _offRxn->getBareRate());
    }
    else {
        setSecondSpecies(sm2);
        
        //change off reaction to include new species
        
        //take the zeroth, second, and fourth species
        Species* s0 = &(_offRxn->rspecies()[0]->getSpecies());
        Species* s2 = &(_offRxn->rspecies()[2]->getSpecies());
        Species* s4 = &(_offRxn->rspecies()[4]->getSpecies());
        
        //create new reaction
        newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                          ({s0, sm2, s2, sb2, s4}, _offRxn->getBareRate());
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

    
    auto sm1 =oldCC->getCMonomer(oldPosition)->speciesMotor(motorType);
    auto sm2 =newCC->getCMonomer(newPosition)->speciesMotor(motorType);
    auto sb2 =newCC->getCMonomer(newPosition)->speciesBound(boundType);
    
    ReactionBase* newOffRxn;
    
    if(getFirstSpecies() == sm1) {
        
        setFirstSpecies(sm2);
        
        //change off reaction to include new species
        
        //take the first, third, and fourth species
        Species* s1 = &(_offRxn->rspecies()[1]->getSpecies());
        Species* s3 = &(_offRxn->rspecies()[3]->getSpecies());
        Species* s4 = &(_offRxn->rspecies()[4]->getSpecies());
        
        //create new reaction
        newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                          ({sm2, s1, sb2, s3, s4}, _offRxn->getBareRate());
        
        //remove old off reaction
        _cc1->removeCrossCylinderReaction(_cc2, _offRxn);
        
        //set new ccylinders
        setFirstCCylinder(newCC);
    }
    else {
        setSecondSpecies(sm2);
        
        //change off reaction to include new species
        
        //take the zeroth, second, and fourth species
        Species* s0 = &(_offRxn->rspecies()[0]->getSpecies());
        Species* s2 = &(_offRxn->rspecies()[2]->getSpecies());
        Species* s4 = &(_offRxn->rspecies()[4]->getSpecies());
        
        //create new reaction
        newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                          ({s0, sm2, s2, sb2, s4}, _offRxn->getBareRate());
        
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