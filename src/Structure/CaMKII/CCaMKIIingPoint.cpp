
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

	/*
	 * _pCaMKIIingPoint is initialized on the constructor of CaMKIIingPoint class.
	 *
	 */


	// Find species on cylinder that should be marked
	SpeciesBound* sb1 = _cc1->getCMonomer(_position1)->speciesCaMKIIer(camkiiType);
    SpeciesBound* se1 = _cc1->getCMonomer(_position1)->speciesBound(
                        SysParams::Chemistry().camkiierBindingBoundIndex[_filamentType]);

    // mark species
    assert(areEqual(sb1->getN(), 0) && areEqual(se1->getN(), 1) && "Major bug: CaMKIIer binding to an occupied site.");
}

void CCaMKIIingPoint::addBond(CCylinder* cc, short pos){

	//Find species on cylinder that should be marked
	SpeciesBound *camkiier = cc->getCMonomer(pos)->speciesCaMKIIer(_camkiiType);
	SpeciesBound *camkii_bindingsite = cc->getCMonomer(pos)->speciesBound(
			SysParams::Chemistry().camkiierBindingBoundIndex[_filamentType]);

	SpeciesBound *camkii_species_cylinder = getSpeciesCaMKIICylinder();

	//mark species
	assert(areEqual(camkiier->getN(), 0) && areEqual(camkii_bindingsite->getN(), 1) &&
		   "Major bug: CaMKIIer binding to an occupied site.");

	camkiier->up();
	camkii_bindingsite->down();

	// Increasing N for the cylinder species
	camkii_species_cylinder->up();

	assert(areEqual(camkiier->getN(), 1) && areEqual(camkii_bindingsite->getN(), 0) &&
		   "Major bug: CaMKIIer didn't bind to the site.");

};

void CCaMKIIingPoint::removeBond(CCylinder* cc, short pos) {
	auto m = cc->getCMonomer(pos);
	auto camkiier = m->speciesCaMKIIer(_camkiiType);
	auto camkii_bindingsite =
			m->speciesBound(SysParams::Chemistry().camkiierBundlingBoundIndex[cc->getCylinder()->getType()]);

	assert(areEqual(camkiier->getN(), 1) && areEqual(camkii_bindingsite->getN(), 0)
		   && "Major bug: CaMKIIer unbundling on a free site.");

	camkiier->down();
	camkii_bindingsite->up();

}

void CCaMKIIingPoint::removeBond(tuple<Cylinder*, short> input) {
	const auto cc = ((Cylinder*)get<0>(input))->getCCylinder();
	const auto pos = get<1>(input);
	removeBond(cc, pos);
}

CCaMKIIingPoint::~CCaMKIIingPoint() {

   // The removeInternalReaction has been moved to this callback: CaMKIIingPointUnbindingCallback.

}

void CCaMKIIingPoint::createOffReactionBinding(ReactionBase* onRxn, SubSystem *ps) {
	SpeciesBound *scad = getSpeciesCaMKIICylinder();

	// Recording the corresponding diffusing species name
	// to recreate the unbinding reaction in the future
	Species *sfb;
	if(onRxn != nullptr) {
		RSpecies **rs = onRxn->rspecies();
		sfb = &(rs[SPECIESCaMKII_BINDING_INDEX]->getSpecies());
		_diffusingSpeciesName = sfb->getName();
	} else {
		assert(_diffusingSpeciesName.size() != 0);
		sfb = this->getCompartment()->findSpeciesByName(_diffusingSpeciesName);
	}


	//create the reaction species
	vector<Species*> os = {scad,
							sfb};

	//create reaction, add to cylinder
	ReactionBase* offRxn =
			new Reaction<1,1>(os, _offRate);

	offRxn->setReactionType(ReactionType::CAMKIIUNBINDING);

	//add the unbinding reaction and callback
	CaMKIIingPointUnbindingCallback bcallback(_pCaMKIIingPoint, ps);
	ConnectionBlock rcb(offRxn->connect(bcallback,false));

	// Setting the current off reaction to the binding off reaction pointer
	setOffReaction(offRxn);
	_pCaMKIIingPoint->getCaMKIICylinder()->getCCylinder()->addInternalReaction(offRxn);
	_offRxnBinding = offRxn;
}

void CCaMKIIingPoint::createOffReactionBundling(SubSystem *ps, FilamentBindingManager *fm) {
	SpeciesBound *scad = getSpeciesCaMKIICylinder();

	//create the reaction species
	vector<Species*> os = {scad};

	//create reaction, add to cylinder
	ReactionBase* offRxn =
			new Reaction<1,0>(os, _offRate);

	offRxn->setReactionType(ReactionType::CAMKIIUNBUNDLING);

	_offRxnBinding->passivateReaction();

	auto cc = _pCaMKIIingPoint->getCaMKIICylinder()->getCCylinder();
	cc->removeInternalReaction(_offRxn);
	cc->addInternalReaction(offRxn);

	_offRxnBundling = offRxn;
	setOffReaction(_offRxnBundling);

	_offRxnBinding = nullptr;

	//add the unbinding reaction and callback
	CaMKIIBundlingManager *bundlingManager = dynamic_cast<CaMKIIBundlingManager*>(fm);
	assert(bundlingManager != nullptr);
	CaMKIIingPointUnbundlingCallback bcallback(_pCaMKIIingPoint, ps, bundlingManager);
	ConnectionBlock rcb(offRxn->connect(bcallback,false));

}

SpeciesBound *CCaMKIIingPoint::getSpeciesCaMKIICylinder() const {
	// This variable is the position of species CaMKII cylinder (speciesCaMKIICylinder) on CaMKII cylinder.
	short positionOfCaMKIIC = 0;
	CCylinder* camkii_cc = _pCaMKIIingPoint->getCaMKIICylinder()->getCCylinder();
	SpeciesBound *scad = camkii_cc->getCMonomer(positionOfCaMKIIC)->speciesCaMKIICylinder(_camkiiType);
	return scad;
}

void CCaMKIIingPoint::createOffReactionCaMKII(ReactionBase* onRxn, SubSystem* ps, FilamentBindingManager *fm) {
	assert(_pCaMKIIingPoint->getCoordinationNumber() != 0L);

	if(_pCaMKIIingPoint->getCoordinationNumber() == 1L)
		createOffReactionBinding(onRxn, ps);
	else if(_pCaMKIIingPoint->getCoordinationNumber() == 2L)
		createOffReactionBundling(ps, fm);


	_pCaMKIIingPoint->updateReactionRates();

	assert(_offRxn->isPassivated() == false);
}

void CCaMKIIingPoint::createOffReaction(ReactionBase* onRxn, SubSystem* ps){
	cerr << "Undefined behavior for createOffReaction on CCaMKIIingPoint. createOffReaction was called "
		 "without bindingManager. Please use the createOffReactionCaMKII method." << endl;
	exit(1);
}
