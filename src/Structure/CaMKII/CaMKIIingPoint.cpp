
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

#include "CaMKIIingPoint.h"

#include "SubSystem.h"
#include "Bead.h"
#include "Cylinder.h"
#include "Filament.h"
#include "ChemRNode.h"
#include "CompartmentGrid.h"

#include "GController.h"
#include "SysParams.h"
#include "MathFunctions.h"
#include "Rand.h"
#include "GController.h"

using namespace mathfunc;

CaMKIIingPoint::CaMKIIingPoint(Cylinder* cylinder, short camkiiType, double position, vector<double> &cp)
    : Trackable(true), _camkiiType(camkiiType), _camkiiID(_camkiiingPoints.getID()), _birthTime(tau()) {

    int pos = int(position * SysParams::Geometry().cylinderNumMon[cylinder->getType()]);

#ifdef CHEMISTRY
    _cCaMKIIingPoint = unique_ptr<CCaMKIIingPoint>(
    new CCaMKIIingPoint(camkiiType, _compartment, cylinder->getCCylinder(), cylinder->getCCylinder(), pos));
    _cCaMKIIingPoint->setCaMKIIingPoint(this);
#endif

#ifdef MECHANICS
    _mCaMKIIingPoint = unique_ptr<MCaMKIIingPoint>(new MCaMKIIingPoint(camkiiType));
    _mCaMKIIingPoint->setCaMKIIingPoint(this);
#endif


    // Creating and initializing CaMKII cylinder
    Bead* b1 = _subSystem->addTrackable<Bead>(cp, nullptr, 0);

    // Always point the coordinate of the first bead
    coordinate = &b1->coordinate;

    // Creating the CaMKII cylinder
    setCaMKIICylinder(_subSystem->addTrackable<CaMKIICylinder>(this, b1, _filType, 0));

    // Add the first bond to the CaMKII
	addBond(cylinder, pos);
    _camkiiCylinder->addToFilamentBindingManagers();

    // Setting the first species to CaMKII cylinder species
	SpeciesBound* sfs = _camkiiCylinder->getCCylinder()->getCMonomer(0)->speciesCaMKIICylinder(camkiiType);
    _cCaMKIIingPoint->setFirstSpecies(sfs);


	//Find compartment
	try {_compartment = GController::getCompartment(*coordinate);}
	catch (exception& e) {
		cout << e.what();
		printSelf();
		exit(EXIT_FAILURE);
	}
	_cCaMKIIingPoint->setCompartment(_compartment);

}

void CaMKIIingPoint::addBond(Cylinder *c, short pos) {
	_cCaMKIIingPoint->addBond(c->getCCylinder(), pos);
	_bonds.push_back(tuple<Cylinder *, short>(c, pos));

}

CaMKIIingPoint::~CaMKIIingPoint() noexcept {
	// removeTrackable of CaMKIICylinder is moved to CaMKIIingPointUnbindingCallback.
}

void CaMKIIingPoint::updatePosition() {

	assert(coordinate != nullptr);

    //Find compartment
    Compartment* c;
    
    try {c = GController::getCompartment(*coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printSelf();

        exit(EXIT_FAILURE);
    }
    
    if(c != _compartment) {
        _compartment = c;
#ifdef CHEMISTRY
        SpeciesBound* firstSpecies = _cCaMKIIingPoint->getFirstSpecies();
        
        CCaMKIIingPoint* clone = _cCaMKIIingPoint->clone(c);
        setCCaMKIIingPoint(clone);

        _cCaMKIIingPoint->setFirstSpecies(firstSpecies);
#endif
    }
}

void CaMKIIingPoint::printSelf() {

    cout << endl;
    
    cout << "CaMKIIingPoint: ptr = " << this << endl;
    cout << "CaMKIIing type = " << _camkiiType << ", CaMKII ID = " << _camkiiID << endl;
    cout << "Coordinates = " << (*coordinate)[0] << ", " << (*coordinate)[1] << ", " << (*coordinate)[2] << endl;
    
    cout << "Position on mother cylinder (short) = " << get<1>(_bonds.at(0)) << endl;
    cout << "Birth time = " << _birthTime << endl;
    
    cout << endl;
    
#ifdef CHEMISTRY
    cout << "Associated species = " << _cCaMKIIingPoint->getFirstSpecies()->getName()
         << " , copy number = " << _cCaMKIIingPoint->getFirstSpecies()->getN()
         << " , position on mother cylinder (int) = " << _cCaMKIIingPoint->getFirstPosition() << endl;
#endif
    
    cout << endl;
    
    cout << "Associated cylinders (mother and camkiiing): " << endl;
    getCylinder(0)->printSelf();
    //getCylinder(1)->printSelf();
    
    cout << endl;
}
            
species_copy_t CaMKIIingPoint::countSpecies(const string& name) {

	/*
	 * This function counts the total number of the CaMKII cylinder species in the system.
	 * In case of multiple type of CaMKII, sname and name as in countSpecies(const string& name) should be implemented.
	 */

	species_copy_t copyNum = 0;

    for(auto b : _camkiiingPoints.getElements()) {

        auto s = b->getCCaMKIIingPoint()->getFirstSpecies();
        string sname = SpeciesNamesDB::removeUniqueFilName(s->getName());

        if(sname == (name + "D"))
			copyNum += 1;
    }
    return copyNum;
}

species_copy_t CaMKIIingPoint::countCaMKIICylinderSpecies(const string &name) {

	/*
	 * This function counts the total number of the CaMKII cylinder species in the system.
	 * In case of multiple type of CaMKII, sname and name as in countSpecies(const string& name) should be implemented.
	 */

	species_copy_t copyNum = 0;

	for(auto b : _camkiiingPoints.getElements()) {

		auto s = b->getCCaMKIIingPoint()->getFirstSpecies();
		string sname = SpeciesNamesDB::removeUniqueFilName(s->getName());

		if(sname == name)
            copyNum += s->getN();
	}
	return copyNum;

}


tuple<Cylinder*, short> CaMKIIingPoint::removeRandomBond(CaMKIIBundlingManager* _bManager) {
	assert(_bonds.size() > 1);

	size_t sz = _bonds.size();
	size_t index = (size_t) Rand::randInteger(0, sz-1);
	tuple<Cylinder*, short> bondToRemove = _bonds[index];

	// Adding back the rest of cylinders of the filament to possible binding
	// before adding back the corresponding cylinder.
	_bManager->addAllCylindersOfFilamentToCaMKIIPossibleBindings(get<0>(bondToRemove), this);

	_cCaMKIIingPoint->removeBond(get<0>(bondToRemove)->getCCylinder(), get<1>(bondToRemove));

	_bonds.erase(_bonds.begin() + index);

	return bondToRemove;
}

void CaMKIIingPoint::updateReactionRates() {

	//get the unbinding reaction
	ReactionBase* offRxn = getCCaMKIIingPoint()->getOffReaction();

	// The propensity of off reaction depends on the coordination number as tracked by species CaMKII cylinder.
	float newRate = offRxn->getBareRate();
	if(SysParams::RUNSTATE==false)
		newRate=0.0;

	offRxn->setRate(newRate);
	offRxn->activateReaction();
	offRxn->updatePropensity();
}


Database<CaMKIIingPoint*> CaMKIIingPoint::_camkiiingPoints;
