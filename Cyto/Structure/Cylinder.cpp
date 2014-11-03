//
//  Cylinder.cpp
//  Cyto
//
//  Created by James Komianos on 7/31/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//
#include <exception>
#include "Cylinder.h"

#include "ChemInitializer.h"
#include "Composite.h"
#include "GController.h"
#include "MathFunctions.h"
#include "Bead.h"

using namespace mathfunc;

Cylinder::Cylinder(Filament* pf, Bead* firstBead, Bead* secondBead, bool extensionFront, bool extensionBack, bool creation) {
    
    ///Set beads
    _pFirst = firstBead;
    _pSecond = secondBead;
    setFilament(pf);
    
    ///check if were still in same compartment
    coordinate = MidPointCoordinate(firstBead->coordinate, secondBead->coordinate, 0.5);

    try {_compartment = GController::getCompartment(coordinate);}
    catch (std::exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
    
    
#ifdef CHEMISTRY
    _cCylinder = std::unique_ptr<CCylinder>(
        ChemInitializer::createCCylinder(ChemInitializerCylinderKey(), pf, _compartment, extensionFront, extensionBack, creation));
    _cCylinder->setCylinder(this);
    
    ///Update filament reactions, only if not initialization
//    if(!(!creation && !extensionFront && !extensionBack))
//        ChemInitializer::updateCCylinder(ChemInitializerCylinderKey(), _cCylinder.get());
    
#endif

#ifdef MECHANICS
    double eqLength;
    
    ///set eqLength according to cylinder size
    if(extensionFront || extensionBack) eqLength = SystemParameters::Geometry().monomerSize;
    else if(creation) eqLength = SystemParameters::Geometry().monomerSize * 2;
    else eqLength = SystemParameters::Geometry().cylinderSize;
    
    _mCylinder = std::unique_ptr<MCylinder>(new MCylinder(eqLength));
    
    _mCylinder->setCylinder(this);
    _mCylinder->setCoordinate(coordinate);
#endif
    
    ///add to compartment
    _compartment->addCylinder(this);
}

Cylinder::~Cylinder() {

    ///remove from compartment
    _compartment->removeCylinder(this);
}

bool Cylinder::IfLast(){
    if (_ifLast) return true;
    
    return false;
}

void Cylinder::SetLast(bool b){ _ifLast = b;}

void Cylinder::updatePosition() {

    ///check if were still in same compartment, set new position
    coordinate = MidPointCoordinate(_pFirst->coordinate, _pSecond->coordinate, 0.5);
    
    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (std::exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
    
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
    ChemInitializer::updateCCylinder(ChemInitializerCylinderKey(), _cCylinder.get());
#endif
    
    ///update exvol neighbors list
#ifdef MECHANICS
    _mCylinder->setCoordinate(coordinate);
    
    ///Find mcylinders that are in surrounding compartments
    std::vector<MCylinder*> nearbyMCylinders;
    
    ///Find surrounding compartments
    std::vector<Compartment*> compartments;
    GController::findCompartments(coordinate, _compartment,
        SystemParameters::Mechanics().VolumeCutoff + SystemParameters::Geometry().largestCompartmentSide * 2, compartments);
    
    for(auto &c : compartments) {
        for(auto &cyl : c->getCylinders()) nearbyMCylinders.push_back(cyl->getMCylinder());
    }
    _mCylinder->updateExVolNeighborsList(nearbyMCylinders);
#endif
    
}


