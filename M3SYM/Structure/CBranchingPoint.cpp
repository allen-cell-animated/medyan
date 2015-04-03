
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

#include "CBranchingPoint.h"

#include "ChemCallbacks.h"
#include "CompartmentContainer.h"
#include "CCylinder.h"
#include "CMonomer.h"

CBranchingPoint::CBranchingPoint(short branchType, Compartment* c,
                                 CCylinder* cc1, CCylinder* cc2, int pos)
    : CBound(c, cc1, cc2), _pos(pos), _branchType(branchType) {

    //Find species on cylinder that should be marked
    SpeciesBound* sb1 =
    _cc1->getCMonomer(pos)->speciesBrancher(branchType);
    
    //attach this branchpoint to the species
    setFirstSpecies(sb1);
}

CBranchingPoint::~CBranchingPoint() {
    
    //remove the unbinding reaction
    _cc1->removeInternalReaction(_offRxn);
    
}

void CBranchingPoint::createOffReaction(ReactionBase* onRxn, SubSystem* ps){
    
    //first, find the correct diffusing or bulk species
    RSpecies** rs = onRxn->rspecies();
    Species* sfb = &(rs[0]->getSpecies());
    
    //create the reaction species
    CMonomer* m = _cc1->getCMonomer(_pos);
    vector<Species*> os = {m->speciesBrancher(_branchType),
                           m->speciesBound(0), sfb};
    
    //create reaction, add to cylinder
    ReactionBase* offRxn =
    new Reaction<BUNBINDINGREACTANTS,BUNBINDINGPRODUCTS>(os, _offRate);
    
    offRxn->setReactionType(ReactionType::BRANCHUNBINDING);
    
    //add the unbinding reaction and callback
    BranchingPointUnbindingCallback bcallback(_pBranchingPoint, ps);
    boost::signals2::shared_connection_block
    rcb(offRxn->connect(bcallback,false));
    
    setOffReaction(offRxn);
    _cc1->addInternalReaction(offRxn);
    
}