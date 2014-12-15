
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "Linker.h"

#include "Bead.h"
#include "Cylinder.h"

#include "GController.h"
#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

Linker::Linker(Cylinder* c1, Cylinder* c2, short linkerType, double position1, double position2, bool creation)
                : _c1(c1), _c2(c2), _linkerType(linkerType), _position1(position1), _position2(position2) {
                    
    //Add to linker db
    LinkerDB::instance()->addLinker(this);
    _linkerID = LinkerDB::instance()->getLinkerID();
                                            
    _birthTime = tau();
        
    //Find compartment
    auto m1 = midPointCoordinate(_c1->getFirstBead()->coordinate, _c1->getSecondBead()->coordinate, _position1);
    auto m2 = midPointCoordinate(_c2->getFirstBead()->coordinate, _c2->getSecondBead()->coordinate, _position2);
    coordinate = midPointCoordinate(m1, m2, 0.5);

    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
        
#ifdef CHEMISTRY
    _cLinker = unique_ptr<CLinker>(new CLinker(_compartment));
    _cLinker->setLinker(this);
        
    //Find species on cylinder that should be marked. If initialization, this should be done. But,
    //if this is because of a reaction callback, it will have already been done.
    int pos1 = int(position1 * SystemParameters::Geometry().cylinderIntSize);
    int pos2 = int(position2 * SystemParameters::Geometry().cylinderIntSize);
    
    SpeciesLinker* sl1 = _c1->getCCylinder()->getCMonomer(pos1)->speciesLinker(linkerType);
    SpeciesLinker* sl2 = _c2->getCCylinder()->getCMonomer(pos2)->speciesLinker(linkerType);
        
    if(!creation) {
        
        SpeciesBound* se1 = _c1->getCCylinder()->getCMonomer(pos1)->speciesBound(0);
        sl1->getRSpecies().up();
        se1->getRSpecies().down();
        
        SpeciesBound* se2 = _c2->getCCylinder()->getCMonomer(pos2)->speciesBound(0);
        sl2->getRSpecies().up();
        se2->getRSpecies().down();
    }
        
    //attach this linker to the species
    _cLinker->setFirstSpecies(sl1);
    _cLinker->setSecondSpecies(sl2);
#endif
    
#ifdef MECHANICS
    _mLinker = unique_ptr<MLinker>(new MLinker(linkerType, position1, position2,
                                   _c1->getFirstBead()->coordinate, _c1->getSecondBead()->coordinate,
                                   _c2->getFirstBead()->coordinate, _c2->getSecondBead()->coordinate));
    _mLinker->setLinker(this);
#endif
}

Linker::~Linker() {
    
    //Remove from linker db
    LinkerDB::instance()->removeLinker(this);
}


void Linker::updatePosition() {
    
    //check if were still in same compartment
    auto m1 = midPointCoordinate(_c1->getFirstBead()->coordinate, _c1->getSecondBead()->coordinate, _position1);
    auto m2 = midPointCoordinate(_c2->getFirstBead()->coordinate, _c2->getSecondBead()->coordinate, _position2);
    coordinate = midPointCoordinate(m1, m2, 0.5);
    
    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
    if(c != _compartment) {
        
        _compartment = c;
#ifdef CHEMISTRY
        CLinker* clone = _cLinker->clone(c);
        
        clone->setFirstSpecies(_cLinker->getFirstSpecies());
        clone->setSecondSpecies(_cLinker->getSecondSpecies());
        
        setCLinker(clone);
#endif
    }
}

void Linker::updateReactionRates() {}

