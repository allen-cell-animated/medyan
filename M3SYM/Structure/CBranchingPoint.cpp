
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

#include "CBranchingPoint.h"

#include "ChemCallbacks.h"

ReactionBase* CBranchingPoint::createOffReaction(vector<Species*> species,
                                                 float rate,
                                                 SubSystem* s) {
    
    ReactionBase* offRxn =
    new Reaction<BUNBINDINGREACTANTS,BUNBINDINGPRODUCTS>(species, rate);
    
    offRxn->setReactionType(ReactionType::BRANCHUNBINDING);
    
    BranchingPointUnbindingCallback bcallback(_pBranchingPoint, s);
    boost::signals2::shared_connection_block
        rcb(offRxn->connect(bcallback,false));
    
    setOffReaction(offRxn);

    return offRxn;
}