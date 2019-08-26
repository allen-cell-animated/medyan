
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
#include "../GController.h"


using namespace mathfunc;

CaMKIIingPoint::CaMKIIingPoint(Cylinder* cylinder, short camkiiType, double position)
    : Trackable(true), _camkiiType(camkiiType), _camkiiID(_camkiiingPoints.getID()), _birthTime(tau()), coordinate(3,0.0) {

	assert(camkiiType == 0);

    int pos = int(position * SysParams::Geometry().cylinderNumMon[cylinder->getType()]);
    //std::cout<<c1->getID()<<" "<<c2->getID()<<" "<<pos<<endl;

#ifdef CHEMISTRY
    _cCaMKIIingPoint = unique_ptr<CCaMKIIingPoint>(
    new CCaMKIIingPoint(camkiiType, _compartment, cylinder->getCCylinder(), cylinder->getCCylinder(), pos));
    _cCaMKIIingPoint->setCaMKIIingPoint(this);
#endif

#ifdef MECHANICS
    _mCaMKIIingPoint = unique_ptr<MCaMKIIingPoint>(new MCaMKIIingPoint(camkiiType));
    _mCaMKIIingPoint->setCaMKIIingPoint(this);
#endif

    addBond(cylinder, pos);

    //Find compartment
    updateCoordinate();
    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) {
      cout << e.what();
      printSelf();
      exit(EXIT_FAILURE);
    }

    Bead* b1 = _subSystem->addTrackable<Bead>(coordinate, nullptr, 0);
  	_camkiiCylinder = unique_ptr<CaMKIICylinder>(new CaMKIICylinder(this, b1, _filType, 0)); // init the dummy cylinder for CaMKII
}

void CaMKIIingPoint::addBond(Cylinder* c, short pos){
  _cCaMKIIingPoint->addBond(c->getCCylinder(), pos);
  _bonds.push_back(tuple<Cylinder*, short>(c, pos));
}

void CaMKIIingPoint::updateCaMKIIingPointCoM(){

	vector<double> temp(3, 0.0);
	for (int i=0; i<_bonds.size(); i++) {
		 Cylinder *bond = get<0>(_bonds[i]);
		 auto pos = get<1>(_bonds[i]);
		 auto mp = midPointCoordinate(bond->_b1->coordinate, bond->_b2->coordinate, pos);
		 temp[0] += mp[0];
		 temp[1] += mp[1];
		 temp[2] += mp[2];
	}
	_coordinate[0] = temp[0]/_bonds.size();
	_coordinate[1] = temp[1]/_bonds.size();
	_coordinate[2] = temp[2]/_bonds.size();
}

CaMKIIingPoint::~CaMKIIingPoint() noexcept {}

void CaMKIIingPoint::updatePosition() {
    
#ifdef CHEMISTRY
    //update ccylinders
    for (int i=0; i<_bonds.size(); i++)
        _cCaMKIIingPoint->setConnectedCCylinder(getCylinder(i)->getCCylinder());
//	_cCaMKIIingPoint->setConnectedCCylinder(get<0>(_bonds[0])->getCCylinder());

//	auto ccyl=get<0>(_bonds[0]);
//	if(ccyl == nullptr){
//	    //A cylinder was removed
//	    cerr << "CC1 is NULL\n";
//		exit(1);
//	}
//	assert(_bonds.size() == 1);

#endif
    //Find compartment
    updateCoordinate();
    
    Compartment* c;
    
    try {c = GController::getCompartment(coordinate);}
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
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    cout << "Position on mother cylinder (double) = " << get<1>(_bonds.at(0)) << endl;
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
    
    species_copy_t copyNum = 0;
    
    for(auto b : _camkiiingPoints.getElements()) {

        auto s = b->getCCaMKIIingPoint()->getFirstSpecies();
        string sname = SpeciesNamesDB::removeUniqueFilName(s->getName());

        if(sname == name)
            copyNum += s->getN();
    }
    return copyNum;
}

void CaMKIIingPoint::removeRandomBond() {
	size_t sz = _bonds.size();
	size_t index = Rand::randInteger(0, sz-1);
	tuple<Cylinder*, double> &bondToRemove = _bonds[index];
	assert(this->getCCaMKIIingPoint()->getOffReaction() == this->getCCaMKIIingPoint()->getOffRxnBundling());
	get<0>(bondToRemove)->getCCylinder()->removeInternalReaction(this->getCCaMKIIingPoint()->getOffReaction());
	_bonds.erase(_bonds.begin() + index);
}

void CaMKIIingPoint::updateReactionRates() {

	//get the unbinding reaction
	ReactionBase* offRxn = getCCaMKIIingPoint()->getOffReaction();

	//change the rate
	float newRate = offRxn->getBareRate() * getCoordinationNumber();
	if(SysParams::RUNSTATE==false)
		newRate=0.0;

	offRxn->setRate(newRate);
	offRxn->updatePropensity();
}


Database<CaMKIIingPoint*> CaMKIIingPoint::_camkiiingPoints;
