
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


using namespace mathfunc;

void CaMKIIingPoint::updateCoordinate() {
    //The coordinate of the Brancher seems to be on the binding site.
	//TODO The coordinate of the CAMKII needs to be on the middle of all the cylinders
	coordinate = midPointCoordinate(get<0>(_bonds.at(0))->getFirstBead()->coordinate,
									get<0>(_bonds.at(0))->getSecondBead()->coordinate,
									get<1>(_bonds.at(0)));
}

CaMKIIingPoint::CaMKIIingPoint(Cylinder* cylinder, short camkiiType, double position)
    : Trackable(true), _camkiiType(camkiiType), _camkiiID(_camkiiingPoints.getID()), _birthTime(tau()), _coordinate(3,0.0) {

    addBond(cylinder, position);
    //Find compartment
    updateCoordinate();

        
    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }
        
    int pos = int(position * SysParams::Geometry().cylinderNumMon[getCylinder(0)->getType()]);
          //std::cout<<c1->getID()<<" "<<c2->getID()<<" "<<pos<<endl;

#ifdef CHEMISTRY
    _cCaMKIIingPoint = unique_ptr<CCaMKIIingPoint>(
    new CCaMKIIingPoint(camkiiType, _compartment, getCylinder(0)->getCCylinder(), getCylinder(0)->getCCylinder(), pos));
    _cCaMKIIingPoint->setCaMKIIingPoint(this);
#endif

#ifdef MECHANICS
    _mCaMKIIingPoint = unique_ptr<MCaMKIIingPoint>(new MCaMKIIingPoint(camkiiType));
    _mCaMKIIingPoint->setCaMKIIingPoint(this);
#endif

    //Composite* parent, Bead* b1, Bead* b2, short type, int position;

   // updatePosition();

    //Dummy Cylinder for CaMKII
    //choose length
    //Composite *Dummy = NULL;

    updateCaMKIIingPointCoM();

    Bead* b1 = _subSystem->addTrackable<Bead>(_coordinate, nullptr, 0);

	double length = 0.0;


	//length = SysParams::Geometry().monomerSize[_filType];
	//length = SysParams::Geometry().monomerSize[_filType];
    //length = SysParams::Geometry().minCylinderSize[_filType];

	//auto pos2 = nextPointProjection(position, length, direction);
	Bead* b2 = _subSystem->addTrackable<Bead>(_coordinate, nullptr, 1);

	//create cylinder
	//Cylinder* c0 = _subSystem->addTrackable<CaMKIICylinder>(nullptr, b1, b2, _filType, 0);

	//c0->setPlusEnd(true);
	//c0->setMinusEnd(true);
	_camkiiCylinder = unique_ptr<CaMKIICylinder>(new CaMKIICylinder(this, b1, _filType, 0)); // init the dummy cylinder for CaMKII
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

        string str = firstSpecies->getName();
        str += "hello";
        
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

Database<CaMKIIingPoint*> CaMKIIingPoint::_camkiiingPoints;
