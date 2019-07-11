
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "CCaMKIIingPoint.h"

#include "ChemCallbacks.h"
#include "CompartmentGrid.h"
#include "CCylinder.h"
#include "CMonomer.h"

CCaMKIIingPoint::CCaMKIIingPoint(short camkiiType, Compartment* c,
                                 CCylinder* cc1, CCylinder* cc2, int position)

    : CBound(cc1->getType(), c, cc1, cc2, position, 0), _camkiiType(camkiiType) {

    //Find species on cylinder that should be marked
    SpeciesBound* sb1 = _cc1->getCMonomer(_position1)->speciesCaMKIIer(camkiiType);
    SpeciesBound* se1 = _cc1->getCMonomer(_position1)->speciesBound(
                        SysParams::Chemistry().camkiierBoundIndex[_filamentType]);

    //mark species
    //TODO assert fails, need to check whether neighbor list update is correct (remove bound binding site from the neighbor list)
    assert(areEqual(sb1->getN(), 0) && areEqual(se1->getN(), 1) &&
           "Major bug: CaMKIIer binding to an occupied site.");
        
    sb1->up(); se1->down();

	assert(areEqual(sb1->getN(), 1) && areEqual(se1->getN(), 0) &&
		   "Major bug: CaMKIIer didn't bind to the site.");

    //attach this camkiipoint to the species
    setFirstSpecies(sb1);

	string str = getFirstSpecies()->getName();
	str = str + string("hello");

}

CCaMKIIingPoint::~CCaMKIIingPoint() {

    // TODO: July 1st, 2019 - Check whether removeInternalReaction is appropriate here!
    //remove the unbinding reaction
    if(!_offRxn) {
		cerr << "Off reaction in CCaMKIIingPoint " << __FILE__ << " (" << __LINE__ << ") is NULL." << endl;
    	exit(1);
    }
	_cc1->removeInternalReaction(_offRxn);

}

void CCaMKIIingPoint::createOffReaction(ReactionBase* onRxn, SubSystem* ps){
    
    //first, find the correct diffusing or bulk species
    RSpecies** rs = onRxn->rspecies();
    Species* sfb = &(rs[SPECIESCaMKII_BINDING_INDEX]->getSpecies());
    
    //create the reaction species
    CMonomer* m = _cc1->getCMonomer(_position1);
    vector<Species*> os = {m->speciesCaMKIIer(_camkiiType),
                           m->speciesBound(SysParams::Chemistry().camkiierBoundIndex[_filamentType]), sfb};
    
    //create reaction, add to cylinder
    ReactionBase* offRxn =
    new Reaction<CAMKIIUNBINDINGREACTANTS,CAMKIIUNBINDINGPRODUCTS>(os, _offRate);
    
    offRxn->setReactionType(ReactionType::CAMKIIUNBINDING);
    
    //add the unbinding reaction and callback
    CaMKIIingPointUnbindingCallback bcallback(_pCaMKIIingPoint, ps);
    ConnectionBlock rcb(offRxn->connect(bcallback,false));
    
    setOffReaction(offRxn);
    _cc1->addInternalReaction(offRxn);
    
}
