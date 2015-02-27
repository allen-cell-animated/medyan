
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

#include "Bead.h"
#include "Cylinder.h"
#include "Filament.h"
#include "ChemRNode.h"
#include "CompartmentContainer.h"

#include "GController.h"
#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

BranchingPoint::BranchingPoint(Cylinder* c1, Cylinder* c2,
                               short branchType, double position,
                               bool creation)
    : _c1(c1), _c2(c2), _position(position), _branchType(branchType) {
                             
    //Add to branch point db
    BranchingPointDB::instance()->addBranchingPoint(this);
    _branchID = BranchingPointDB::instance()->getBranchID();
    
    _birthTime = tau();
    
    //Find compartment
    coordinate = midPointCoordinate(_c1->getFirstBead()->coordinate,
                                    _c1->getSecondBead()->coordinate, _position);
    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
#ifdef CHEMISTRY
    _cBranchingPoint = unique_ptr<CBranchingPoint>(new CBranchingPoint(_compartment));
    _cBranchingPoint->setBranchingPoint(this);
        
    //Find species on cylinder that should be marked. If initialization,
    //this should be done. But, if this is because of a reaction callback,
    //it will have already been done.
    int pos = int(position * SysParams::Geometry().cylinderIntSize);
    
    SpeciesBrancher* sb1 =
    _c1->getCCylinder()->getCMonomer(pos)->speciesBrancher(branchType);
    SpeciesBrancher* sb2 =
    _c2->getCCylinder()->getCMonomer(pos)->speciesBrancher(branchType);
    
    if(!creation) {
        SpeciesBound* se1 =
        _c1->getCCylinder()->getCMonomer(pos)->speciesBound(0);
        sb1->up();
        se1->down();
        
        SpeciesBound* se2 =
        _c2->getCCylinder()->getCMonomer(pos)->speciesBound(0);
        sb2->up();
        se2->down();
    }
    
    //attach this branchpoint to the species
    _cBranchingPoint->setFirstSpecies(sb1);
    _cBranchingPoint->setSecondSpecies(sb2);
#endif
    
#ifdef MECHANICS
    _mBranchingPoint = unique_ptr<MBranchingPoint>(new MBranchingPoint(branchType));
    _mBranchingPoint->setBranchingPoint(this);
#endif
}

BranchingPoint::~BranchingPoint() noexcept {
    
    //Remove from branch point db
    BranchingPointDB::instance()->removeBranchingPoint(this);
    
#ifdef CHEMISTRY
    //remove the unbinding reaction
    CCylinder* cc = _c1->getCCylinder();
    cc->removeInternalReaction(_cBranchingPoint->getOffReaction());
    
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
        
        string speciesName = SpeciesNamesDB::Instance()->
        removeUniqueName(speciesFilament->getName());
        
        //find the free monomer, either bulk or diffusing
        Species* freeMonomer = _compartment->findSpeciesByName(speciesName);
        if(freeMonomer == nullptr)
            freeMonomer = CompartmentGrid::instance()->
                          findSpeciesBulkByName(speciesName);
        
        //remove the filament from the system
        delete _c2->getFilament();
        
        //update reaction rates
        for(auto &r : freeMonomer->getRSpecies().reactantReactions())
            r->activateReaction();
    }
#endif
    
}

void BranchingPoint::updatePosition() {
    
    //Find compartment
    coordinate = midPointCoordinate(_c1->getFirstBead()->coordinate,
                                    _c1->getSecondBead()->coordinate, _position);

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

