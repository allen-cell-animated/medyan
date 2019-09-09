
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
#include "Reaction.h"

CCaMKIIingPoint::CCaMKIIingPoint(short camkiiType, Compartment* c,
                                 CCylinder* cc1, CCylinder* cc2, int position)

    : CBound(cc1->getType(), c, cc1, cc2, position, 0), _camkiiType(camkiiType),
    		_offRxnBinding(nullptr),
    		_offRxnBundling(nullptr) {
//    Code block moved to addBond
//    //Find species on cylinder that should be marked
    SpeciesBound* sb1 = _cc1->getCMonomer(_position1)->speciesCaMKIIer(camkiiType);
    SpeciesBound* se1 = _cc1->getCMonomer(_position1)->speciesBound(
                        SysParams::Chemistry().camkiierBindingBoundIndex[_filamentType]);
//
//    //mark species
    assert(areEqual(sb1->getN(), 0) && areEqual(se1->getN(), 1) && "Major bug: CaMKIIer binding to an occupied site.");
//
//    sb1->up(); se1->down();
//
//	assert(areEqual(sb1->getN(), 1) && areEqual(se1->getN(), 0) &&
//		   "Major bug: CaMKIIer didn't bind to the site.");
//
//    //attach this camkiipoint to the species
    setFirstSpecies(sb1);

}

void CCaMKIIingPoint::addBond(CCylinder* cc, short pos){

  //Find species on cylinder that should be marked
  SpeciesBound* sb1 = cc->getCMonomer(pos)->speciesCaMKIIer(_camkiiType);
  SpeciesBound* se1 = cc->getCMonomer(pos)->speciesBound(
          SysParams::Chemistry().camkiierBindingBoundIndex[_filamentType]);

  //mark species
  assert(areEqual(sb1->getN(), 0) && areEqual(se1->getN(), 1) &&
         "Major bug: CaMKIIer binding to an occupied site.");

  sb1->up(); se1->down();

  assert(areEqual(sb1->getN(), 1) && areEqual(se1->getN(), 0) &&
         "Major bug: CaMKIIer didn't bind to the site.");
  
  
};

CCaMKIIingPoint::~CCaMKIIingPoint() {

   // TODO: July 1st, 2019 - Check whether removeInternalReaction is appropriate here!
	_cc1->removeInternalReaction(_offRxn);

}

void CCaMKIIingPoint::createOffReactionBinding(ReactionBase *onRxn, SubSystem *ps) {
	//first, find the correct diffusing or bulk species
	RSpecies** rs = onRxn->rspecies();
	Species* sfb = &(rs[SPECIESCaMKII_BINDING_INDEX]->getSpecies());

	//create the reaction species
	CMonomer* m = _cc1->getCMonomer(_position1);
	vector<Species*> os = {m->speciesCaMKIIer(_camkiiType),
						   m->speciesBound(SysParams::Chemistry().camkiierBindingBoundIndex[_filamentType]), sfb};

	//create reaction, add to cylinder
	ReactionBase* offRxn =
			new Reaction<CAMKIIUNBINDINGREACTANTS,CAMKIIUNBINDINGPRODUCTS>(os, _offRate);

	offRxn->setReactionType(ReactionType::CAMKIIUNBINDING);

	//add the unbinding reaction and callback
	CaMKIIingPointUnbindingCallback bcallback(_pCaMKIIingPoint, ps);
	ConnectionBlock rcb(offRxn->connect(bcallback,false));

	setOffReaction(offRxn);
	_cc1->addInternalReaction(offRxn);
	_offRxnBinding = offRxn;
}

void CCaMKIIingPoint::createOffReactionBundling(ReactionBase *onRxn, SubSystem *ps) {
	//first, find the correct diffusing or bulk species
	RSpecies** rs = onRxn->rspecies();

	//create the reaction species
	CMonomer* m = _cc1->getCMonomer(_position1);
	Species* sfb = m->speciesCaMKIIer(_camkiiType);
	vector<Species*> os = {sfb,
						   m->speciesBound(SysParams::Chemistry().camkiierBindingBoundIndex[_filamentType]), sfb};

	//create reaction, add to cylinder
	ReactionBase* offRxn =
			new Reaction<CAMKIIUNBUNDLINGGREACTANTS,CAMKIIUNBUNDLINGPRODUCTS>(os, _offRate);

	offRxn->setReactionType(ReactionType::CAMKIIUNBUNDLING);

	//add the unbinding reaction and callback
	CaMKIIingPointUnbundlingCallback bcallback(_pCaMKIIingPoint, ps);
	ConnectionBlock rcb(offRxn->connect(bcallback,false));

	_offRxnBinding->passivateReaction();

	auto bonds = _pCaMKIIingPoint->getBonds();
	for(int i=0;i<bonds.size();i++) {
		auto cc = get<0>(bonds[i])->getCCylinder();
		cc->removeInternalReaction(_offRxn);
		cc->addInternalReaction(offRxn);
	}

	setOffReaction(offRxn);

	_offRxnBundling = offRxn;
}

void CCaMKIIingPoint::createOffReaction(ReactionBase* onRxn, SubSystem* ps){
    // TODO: This reaction needs to be implemented for unbinding and unbundling mix.
    // Currently it only implements unbinding.

    assert(_pCaMKIIingPoint->getCoordinationNumber() != 0);

    if(_pCaMKIIingPoint->getCoordinationNumber() == 1) {
		    createOffReactionBinding(onRxn, ps);
    } else if(_pCaMKIIingPoint->getCoordinationNumber() == 2){
		    createOffReactionBundling(onRxn, ps);
    } else {
        auto bonds = _pCaMKIIingPoint->getBonds();
        for(int i=0;i<bonds.size();i++) {
          auto cc = get<0>(bonds[i])->getCCylinder();
          cc->addInternalReaction(_offRxnBundling);
        }
    }

    _pCaMKIIingPoint->updateReactionRates();

    assert(_offRxn->isPassivated() == false);

    cout << "========== CaMKII OffReaction created " <<  "  _offRxnBinding: "<< _offRxnBinding << "  _offRxnBundling: " << _offRxnBundling <<"  _offRxn: " <<_offRxn << "  isPassivated: " << _offRxnBinding->isPassivated() << " " << ((_offRxnBundling != NULL) ? _offRxnBundling->isPassivated() : true) << " " << _offRxn->isPassivated() << endl;
}