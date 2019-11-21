
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
//    setFirstSpecies(sfs);

}

void CCaMKIIingPoint::addBond(CCylinder* cc, short pos){

	//Find species on cylinder that should be marked
	SpeciesBound *sb1 = cc->getCMonomer(pos)->speciesCaMKIIer(_camkiiType);
	SpeciesBound *se1 = cc->getCMonomer(pos)->speciesBound(
			SysParams::Chemistry().camkiierBindingBoundIndex[_filamentType]);

	SpeciesBound *scad = getSpeciesCaMKIIDummyCylinder();

	//mark species
	assert(areEqual(sb1->getN(), 0) && areEqual(se1->getN(), 1) &&
		   "Major bug: CaMKIIer binding to an occupied site.");

	cout<<"====== JAMES: Before: ";
	cc->getCMonomer(pos)->print();
	cout << endl;
//	cout<<cc->getCylinder()->getID()<<" "<<pos<<" "<<sb1->getN()<<" "<<se1->getN()<<endl;

	sb1->up();
	se1->down();

	// Increasing N for the dummy cylinder species
	scad->up();

	cout<<"====== JAMES: After : ";
	cc->getCMonomer(pos)->print();
	cout << endl;

	assert(areEqual(sb1->getN(), 1) && areEqual(se1->getN(), 0) &&
		   "Major bug: CaMKIIer didn't bind to the site.");



};

void CCaMKIIingPoint::removeBond(CCylinder* cc, short pos) {
	auto m = cc->getCMonomer(pos);
	auto camkiier = m->speciesCaMKIIer(_camkiiType);
	auto camkii_bindingsite =
			m->speciesBound(SysParams::Chemistry().camkiierBundlingBoundIndex[cc->getCylinder()->getType()]);

	assert(areEqual(camkiier->getN(), 1) && areEqual(camkii_bindingsite->getN(), 0)
		   && "Major bug: CaMKIIer unbundling on a free site.");

	cout<<"====== JAMES: Before: ";
	cc->getCMonomer(pos)->print();
	cout << endl;

	camkiier->down();
	camkii_bindingsite->up();

	cout<<"====== JAMES: After : ";
	cc->getCMonomer(pos)->print();
	cout << endl;
}

void CCaMKIIingPoint::removeBond(tuple<Cylinder*, short> input) {
	const auto cc = ((Cylinder*)get<0>(input))->getCCylinder();
	const auto pos = get<1>(input);
	removeBond(cc, pos);
}

CCaMKIIingPoint::~CCaMKIIingPoint() {

   // TODO: July 1st, 2019 - Check whether removeInternalReaction is appropriate here!
//	_pCaMKIIingPoint->getCaMKIICylinder()->getCCylinder()->removeInternalReaction(_offRxn);

}

void CCaMKIIingPoint::createOffReactionBinding(ReactionBase* onRxn, SubSystem *ps) {
	SpeciesBound *scad = getSpeciesCaMKIIDummyCylinder();

	// Recording the corresponding diffusing species name
	// for unbinding reaction re-creation in future
	Species *sfb;
	if(onRxn != nullptr) {
		RSpecies **rs = onRxn->rspecies();
		sfb = &(rs[SPECIESCaMKII_BINDING_INDEX]->getSpecies());
		diffusingSpeciesName = sfb->getName();
	} else {
		assert(diffusingSpeciesName.size() != 0);
		sfb = this->getCompartment()->findSpeciesByName(diffusingSpeciesName);
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

	setOffReaction(offRxn);
	_pCaMKIIingPoint->getCaMKIICylinder()->getCCylinder()->addInternalReaction(offRxn);
	_offRxnBinding = offRxn;
}

void CCaMKIIingPoint::createOffReactionBundling(SubSystem *ps, FilamentBindingManager *fm) {
	SpeciesBound *scad = getSpeciesCaMKIIDummyCylinder();

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

SpeciesBound *CCaMKIIingPoint::getSpeciesCaMKIIDummyCylinder() const {
	// This variable is the position of species CaMKII Dummy cylinder (speciesCaMKIIDummyCylinder) on CaMKII cylinder.
	short positionOfCaMKIIDummyC = 0;
	CCylinder* camkii_cc = _pCaMKIIingPoint->getCaMKIICylinder()->getCCylinder();
	SpeciesBound *scad = camkii_cc->getCMonomer(positionOfCaMKIIDummyC)->speciesCaMKIIDummyCylinder(_camkiiType);
	return scad;
}

void CCaMKIIingPoint::createOffReactionCaMKII(ReactionBase* onRxn, SubSystem* ps, FilamentBindingManager *fm) {
	// TODO: This reaction needs to be implemented for unbinding and unbundling mix.
	// Currently it only implements unbinding.

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
