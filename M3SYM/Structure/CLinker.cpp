
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

#include "CLinker.h"

#include "ChemCallbacks.h"
#include "CCylinder.h"

CLinker::CLinker(short linkerType, Compartment* c,
                 CCylinder* cc1, CCylinder* cc2, int pos1, int pos2)

    : CBound(c), _cc1(cc1), _cc2(cc2) {
    
        //Find species on cylinder that should be marked
        SpeciesBound* sl1 = _cc1->getCMonomer(pos1)->speciesLinker(linkerType);
        SpeciesBound* sl2 = _cc2->getCMonomer(pos2)->speciesLinker(linkerType);
        
        //attach this linker to the species
        setFirstSpecies(sl1);
        setSecondSpecies(sl2);
}

CLinker::~CLinker() {

    //remove the off reaction
    _cc1->removeCrossCylinderReaction(_cc2, _offRxn);

}

void CLinker::createOffReaction(ReactionBase* onRxn, float offRate, SubSystem* ps) {
    
    RSpecies** rs = onRxn->rspecies();
    vector<Species*> os;
    
    //copy into offspecies vector in opposite order
    for(int i = LMBINDINGREACTANTS; i < LMBINDINGREACTANTS+LMBINDINGPRODUCTS; i++)
        os.push_back(&rs[i]->getSpecies());
    
    for(int i = 0; i < LMBINDINGREACTANTS; i++)
        os.push_back(&rs[i]->getSpecies());
    
    ReactionBase* offRxn =
    new Reaction<LMUNBINDINGREACTANTS,LMUNBINDINGPRODUCTS>(os, offRate);
    offRxn->setReactionType(ReactionType::LINKERUNBINDING);
    
    //Attach the callback to the off reaction, add it
    LinkerUnbindingCallback lcallback(_pLinker, ps);
    boost::signals2::shared_connection_block rcb(offRxn->connect(lcallback,false));
    
    _cc1->addCrossCylinderReaction(_cc2, offRxn);
    setOffReaction(offRxn);

}