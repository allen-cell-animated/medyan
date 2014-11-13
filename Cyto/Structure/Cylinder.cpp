//
//  Cylinder.cpp
//  Cyto
//
//  Created by James Komianos on 7/31/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//
#include <exception>
#include "Cylinder.h"

#include "ChemManager.h"
#include "Bead.h"
#include "Composite.h"
#include "GController.h"
#include "NeighborListDB.h"

#include "MathFunctions.h"


using namespace mathfunc;

Cylinder::Cylinder(Filament* f, Bead* b1, Bead* b2, int ID, bool extensionFront, bool extensionBack, bool creation)
                   : _pFilament(f), _b1(b1), _b2(b2), _ID(ID) {
    
    ////Set coordinate
    coordinate = MidPointCoordinate(_b1->coordinate, _b2->coordinate, 0.5);

    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
                   
   ///add to compartment
   _compartment->addCylinder(this);
   ///add to neighbor lists
   NeighborListDB::instance(NeighborListDBKey())->addNeighbor(this);
    
#ifdef CHEMISTRY
    _cCylinder = unique_ptr<CCylinder>(
        ChemManager::createCCylinder(ChemManagerCylinderKey(), f, _compartment, extensionFront, extensionBack, creation));
    _cCylinder->setCylinder(this);
    
    if(creation || extensionFront || extensionBack) {
        ///Update filament reactions, only if not initialization
        ChemManager::updateCCylinder(ChemManagerCylinderKey(), _cCylinder.get());
    }
    
#endif

#ifdef MECHANICS
    double eqLength;
    
    ///set eqLength according to cylinder size
    if(extensionFront || extensionBack) eqLength = SystemParameters::Geometry().monomerSize;
    else if(creation) eqLength = SystemParameters::Geometry().monomerSize * 2;
    else eqLength = SystemParameters::Geometry().cylinderSize;
    
    _mCylinder = unique_ptr<MCylinder>(new MCylinder(eqLength));
    
    _mCylinder->setCylinder(this);
    
#endif
    

}

Cylinder::~Cylinder() {
    ///remove from compartment
    _compartment->removeCylinder(this);
    ///remove from neighbor lists
    NeighborListDB::instance(NeighborListDBKey())->removeNeighbor(this);
}

void Cylinder::updatePosition() {

    ///check if were still in same compartment, set new position
    coordinate = MidPointCoordinate(_b1->coordinate, _b2->coordinate, 0.5);
    
    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
    if(c != _compartment) {
        
        ///remove from old compartment, add to new
        _compartment->removeCylinder(this);
        _compartment = c;
        _compartment->addCylinder(this);
        
#ifdef CHEMISTRY
        CCylinder* clone = _cCylinder->clone(c);
        setCCylinder(clone);
#endif
    }
    
///Update filament reactions
#ifdef CHEMISTRY
    ChemManager::updateCCylinder(ChemManagerCylinderKey(), _cCylinder.get());
#endif

}


