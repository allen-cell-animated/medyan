
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

#include "CLinker.h"

#include "ChemCallbacks.h"
#include "CCylinder.h"
#include "Cylinder.h"

CLinker::CLinker(short linkerType, Compartment* c,
                 CCylinder* cc1, CCylinder* cc2, int position1, int position2)

    : CBound(cc1->getType(), c, cc1, cc2, position1, position2) {
        
    //Find species on cylinder that should be marked
    SpeciesBound* sl1 = _cc1->getCMonomer(_position1)->speciesLinker(linkerType);
    SpeciesBound* sl2 = _cc2->getCMonomer(_position2)->speciesLinker(linkerType);
        
    SpeciesBound* se1 = _cc1->getCMonomer(_position1)->speciesBound(
                        SysParams::Chemistry().linkerBoundIndex[_filamentType]);
    SpeciesBound* se2 = _cc2->getCMonomer(_position2)->speciesBound(
                        SysParams::Chemistry().linkerBoundIndex[_filamentType]);
    
    //mark species
    assert(areEqual(sl1->getN(), 0.0) && areEqual(sl2->getN(), 0.0) &&
           areEqual(se1->getN(), 1.0) && areEqual(se2->getN(), 1.0) &&
           "Major bug: Linker binding to an occupied site.");
    
    sl1->up(); sl2->up();
    se1->down(); se2->down();
    
    //attach this linker to the species
    setFirstSpecies(sl1);
    setSecondSpecies(sl2);
}

CLinker::~CLinker() {

    //remove the off reaction
    _cc1->removeCrossCylinderReaction(_cc2, _offRxn);

}

void CLinker::createOffReaction(ReactionBase* onRxn, SubSystem* ps) {
    
    RSpecies** rs = onRxn->rspecies();
    vector<Species*> os;
    
    //copy into offspecies vector
    os.push_back(_firstSpecies);
    os.push_back(_secondSpecies);
    
    os.push_back(&rs[SPECIESL_BINDING_INDEX]->getSpecies());
    
    Species* empty1 = _cc1->getCMonomer(_position1)->speciesBound(
                      SysParams::Chemistry().linkerBoundIndex[_filamentType]);
    Species* empty2 = _cc2->getCMonomer(_position2)->speciesBound(
                      SysParams::Chemistry().linkerBoundIndex[_filamentType]);
    
    os.push_back(empty1);
    os.push_back(empty2);
    
    ReactionBase* offRxn =
    new Reaction<LMUNBINDINGREACTANTS,LMUNBINDINGPRODUCTS>(os, _offRate);
    offRxn->setReactionType(ReactionType::LINKERUNBINDING);
    
    //Attach the callback to the off reaction, add it
    LinkerUnbindingCallback lcallback(_pLinker, ps);
    ConnectionBlock rcb(offRxn->connect(lcallback,false));
    
    _cc1->addCrossCylinderReaction(_cc2, offRxn);
    setOffReaction(offRxn);
}