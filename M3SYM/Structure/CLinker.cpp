
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

void CLinker::createOffReaction(ReactionBase* onRxn, float offRate, SubSystem* ps) {
    
    RSpecies** rSpecies = onRxn->rspecies();
    vector<Species*> offSpecies;
    
    //copy into offspecies vector in opposite order
    for(int i = LMBINDINGREACTANTS; i < LMBINDINGREACTANTS + LMBINDINGPRODUCTS; i++)
        offSpecies.push_back(&rSpecies[i]->getSpecies());
    
    for(int i = 0; i < LMBINDINGREACTANTS; i++)
        offSpecies.push_back(&rSpecies[i]->getSpecies());
    
    ReactionBase* offRxn =
    new Reaction<LMUNBINDINGREACTANTS,LMUNBINDINGPRODUCTS>(offSpecies, offRate);
    offRxn->setReactionType(ReactionType::LINKERUNBINDING);
    
    //Attach the callback to the off reaction, add it
    LinkerUnbindingCallback lcallback(_pLinker, ps);
    boost::signals2::shared_connection_block
    rcb(offRxn->connect(lcallback,false));
    
    _cc1->addCrossCylinderReaction(_cc2, offRxn);
    setOffReaction(offRxn);

}