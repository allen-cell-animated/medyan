
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

#include "BranchingPoint.h"

#include "SubSystem.h"
#include "Bead.h"
#include "Cylinder.h"
#include "Filament.h"
#include "ChemRNode.h"
#include "CompartmentGrid.h"

#include "GController.h"
#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

Database<BranchingPoint*> BranchingPoint::_branchingPoints;

void BranchingPoint::updateCoordinate() {
    
    coordinate = midPointCoordinate(_c1->getFirstBead()->coordinate,
                                    _c1->getSecondBead()->coordinate,
                                    _position);
}

BranchingPoint::BranchingPoint(Cylinder* c1, Cylinder* c2,
                               short branchType, double position)

    : Trackable(true),
      _c1(c1), _c2(c2), _position(position),
      _branchType(branchType), _branchID(_branchingPoints.getID()),
      _birthTime(tau()) {
    
    //Find compartment
    updateCoordinate();
        
    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
        
    int pos = int(position * SysParams::Geometry().cylinderIntSize);
    
#ifdef CHEMISTRY
    _cBranchingPoint = unique_ptr<CBranchingPoint>(
                       new CBranchingPoint(branchType, _compartment,
                       c1->getCCylinder(), c2->getCCylinder(), pos));
    _cBranchingPoint->setBranchingPoint(this);
#endif
    
#ifdef MECHANICS
    _mBranchingPoint = unique_ptr<MBranchingPoint>(
                       new MBranchingPoint(branchType));
    _mBranchingPoint->setBranchingPoint(this);
#endif
        
    //set the branching cylinder
    _c1->setBranchingCylinder(_c2);
}

BranchingPoint::~BranchingPoint() noexcept {
    
#ifdef CHEMISTRY
    //mark the correct species on the minus end of the branched
    //filament. If this is a filament species, change it to its
    //corresponding minus end. If a plus end, release a diffusing
    //or bulk species, depending on the initial reaction.
    CMonomer* m = _c2->getCCylinder()->getCMonomer(0);
    short speciesFilament = m->activeSpeciesFilament();
    
    //there is a filament species, mark its corresponding minus end
    if(speciesFilament != -1) {
        m->speciesMinusEnd(speciesFilament)->up();
    }
    //mark the free species instead
    else {
        //find the free species
        Species* speciesFilament = m->speciesFilament(m->activeSpeciesPlusEnd());
        
        string speciesName = SpeciesNamesDB::
                             removeUniqueFilName(speciesFilament->getName());
        string speciesFirstChar = speciesName.substr(0,1);
        
        //find the free monomer, either bulk or diffusing
        Species* freeMonomer = nullptr;
        
        Species* dMonomer  = _compartment->findSpeciesByName(speciesName);
        Species* dfMonomer = _compartment->findSpeciesByName(speciesFirstChar);
        
        Species* bMonomer  = _subSystem->getCompartmentGrid()->
                             findSpeciesBulkByName(speciesName);
        Species* bfMonomer = _subSystem->getCompartmentGrid()->
                             findSpeciesBulkByName(speciesFirstChar);
        
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
            cout << "In unbranching reaction, could not find corresponding " <<
                    "diffusing species of filament species " << speciesName <<
                    ". Exiting." << endl;
            exit(EXIT_FAILURE);
        }
            
        //remove the filament from the system
        _subSystem->removeTrackable<Filament>(_c2->getFilament());
            
        //update reactions
        freeMonomer->getRSpecies().activateAssocReactantReactions();
    }
#endif
}

void BranchingPoint::updatePosition() {
    
#ifdef CHEMISTRY
    //update ccylinders
    _cBranchingPoint->setFirstCCylinder(_c1->getCCylinder());
    _cBranchingPoint->setSecondCCylinder(_c2->getCCylinder());
    
#endif
    //Find compartment
    updateCoordinate();
    
    Compartment* c;
    
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
    if(c != _compartment) {
        _compartment = c;
#ifdef CHEMISTRY
        SpeciesBound* firstSpecies = _cBranchingPoint->getFirstSpecies();
        
        CBranchingPoint* clone = _cBranchingPoint->clone(c);
        setCBranchingPoint(clone);
        
        _cBranchingPoint->setFirstSpecies(firstSpecies);
#endif
    }
}

