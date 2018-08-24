
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
    coordinate = midPointCoordinate(_cylinders.at(0)->getFirstBead()->coordinate,
    								_cylinders.at(0)->getSecondBead()->coordinate,
                                    _position);
}

CaMKIIingPoint::CaMKIIingPoint(vector<Cylinder*> cylinders,
                               short camkiiType, double position)

    : Trackable(true), _cylinders(cylinders), _position(position), _coordinationNum(1),
      _camkiiType(camkiiType), _camkiiID(_camkiiingPoints.getID()), _birthTime(tau()) {
    
    //Find compartment
    updateCoordinate();
        
    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }
        
    int pos = int(position * SysParams::Geometry().cylinderNumMon[cylinders.at(0)->getType()]);
          //std::cout<<c1->getID()<<" "<<c2->getID()<<" "<<pos<<endl;
#ifdef CHEMISTRY
    _cCaMKIIingPoint = unique_ptr<CCaMKIIingPoint>(
    new CCaMKIIingPoint(camkiiType, _compartment, _cylinders.at(0)->getCCylinder(),_cylinders.at(1)->getCCylinder(), pos));
    _cCaMKIIingPoint->setCaMKIIingPoint(this);
#endif
    
#ifdef MECHANICS
    _mCaMKIIingPoint = unique_ptr<MCaMKIIingPoint>(new MCaMKIIingPoint(camkiiType));
    _mCaMKIIingPoint->setCaMKIIingPoint(this);
#endif
        
}

CaMKIIingPoint::~CaMKIIingPoint() noexcept {
    
#ifdef MECHANICS
    //offset the camkiiing cylinder's bead by a little for safety
    auto msize = SysParams::Geometry().monomerSize[_cylinders.at(0)->getType()];
    
    vector<double> offsetCoord =
    {(Rand::randInteger(0,1) ? -1 : +1) * Rand::randDouble(msize, 2 * msize),
     (Rand::randInteger(0,1) ? -1 : +1) * Rand::randDouble(msize, 2 * msize),
     (Rand::randInteger(0,1) ? -1 : +1) * Rand::randDouble(msize, 2 * msize)};
    
    auto b = _cylinders.at(1)->getFirstBead();
    
    b->coordinate[0] += offsetCoord[0];
    b->coordinate[1] += offsetCoord[1];
    b->coordinate[2] += offsetCoord[2];
#endif
    
    
#ifdef CHEMISTRY
    //mark the correct species on the minus end of the camkiied
    //filament. If this is a filament species, change it to its
    //corresponding minus end. If a plus end, release a diffusing
    //or bulk species, depending on the initial reaction.
    CMonomer* m = _cylinders.at(1)->getCCylinder()->getCMonomer(0);
    short speciesFilament = m->activeSpeciesFilament();
    
    //there is a filament species, mark its corresponding minus end
    if(speciesFilament != -1) {
        m->speciesMinusEnd(speciesFilament)->up();
        
        //unmark the filament and bound species
        m->speciesFilament(speciesFilament)->down();
        m->speciesBound(SysParams::Chemistry().camkiierBoundIndex[_cylinders.at(0)->getType()])->down();
    }
    //mark the free species instead
    else {
        //find the free species
        Species* speciesFilament = m->speciesFilament(m->activeSpeciesPlusEnd());
        
        string speciesName = SpeciesNamesDB::removeUniqueFilName(speciesFilament->getName());
        string speciesFirstChar = speciesName.substr(0,1);
        
        //find the free monomer, either bulk or diffusing
        Species* freeMonomer = nullptr;
        auto grid = _subSystem->getCompartmentGrid();
        
        Species* dMonomer  = _compartment->findSpeciesByName(speciesName);
        Species* dfMonomer = _compartment->findSpeciesByName(speciesFirstChar);
        
        Species* bMonomer  = grid->findSpeciesBulkByName(speciesName);
        Species* bfMonomer = grid->findSpeciesBulkByName(speciesFirstChar);
        
        //try diffusing
        if(dMonomer != nullptr) freeMonomer = dMonomer;
        // try bulk
        else if(bMonomer  != nullptr) freeMonomer = bMonomer;
        //diffusing, remove all but first char
        else if(dfMonomer != nullptr) freeMonomer = dfMonomer;
        //bulk, remove all but first char
        else if(bfMonomer != nullptr) freeMonomer = bfMonomer;
        //could not find. exit ungracefully
        else {
            cout << "In uncamkiiing reaction, could not find corresponding " <<
                    "diffusing species of filament species " << speciesName <<
                    ". Exiting." << endl;
            exit(EXIT_FAILURE);
        }
            
        //remove the filament from the system
        Filament *bf = (Filament*)(_cylinders.at(1)->getParent());
        _subSystem->removeTrackable<Filament>(bf);
        
        delete bf;
            
        //mark species, update reactions
        freeMonomer->up();
        freeMonomer->updateReactantPropensities();
    }
#endif
    //reset camkiiing cylinder
    _cylinders.at(0)->setCaMKIIingCylinder(nullptr);
}

void CaMKIIingPoint::updatePosition() {
    
#ifdef CHEMISTRY
    //update ccylinders
    _cCaMKIIingPoint->setFirstCCylinder(_cylinders.at(0)->getCCylinder());
    _cCaMKIIingPoint->setSecondCCylinder(_cylinders.at(1)->getCCylinder());
    
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
    
    cout << "Position on mother cylinder (double) = " << _position << endl;
    cout << "Birth time = " << _birthTime << endl;
    
    cout << endl;
    
#ifdef CHEMISTRY
    cout << "Associated species = " << _cCaMKIIingPoint->getFirstSpecies()->getName()
         << " , copy number = " << _cCaMKIIingPoint->getFirstSpecies()->getN()
         << " , position on mother cylinder (int) = " << _cCaMKIIingPoint->getFirstPosition() << endl;
#endif
    
    cout << endl;
    
    cout << "Associated cylinders (mother and camkiiing): " << endl;
    _cylinders.at(0)->printSelf();
    _cylinders.at(1)->printSelf();
    
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
