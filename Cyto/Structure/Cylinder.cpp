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

Cylinder::Cylinder(Filament* pf, Bead* firstBead, Bead* secondBead, Compartment* c, bool extensionFront, bool extensionBack) {
    
    setFilament(pf);
    
#ifdef CHEMISTRY
    _cCylinder = std::unique_ptr<CCylinder>(
        ChemInitializer::createCCylinder(ChemInitializerCylinderKey(), pf, c, extensionFront, extensionBack));
    _cCylinder->setCylinder(this);
#endif

#ifdef MECHANICS
    double eqLength;
    if(extensionFront || extensionBack) eqLength = SystemParameters::Geometry().monomerSize;
    else eqLength = SystemParameters::Geometry().cylinderSize;
    
    _mCylinder = std::unique_ptr<MCylinder>(new MCylinder(eqLength));
    _mCylinder->setCylinder(this);
#endif
    
    ///Set beads
    _pFirst = firstBead;
    _pSecond = secondBead;
    
}

Cylinder::~Cylinder() {}

bool Cylinder::IfLast(){
    if (_ifLast) return true;
    
    return false;
}

void Cylinder::SetLast(bool b){ _ifLast = b;}

void Cylinder::updatePosition() {

#ifdef CHEMISTRY
    ///check if were still in same compartment
    auto midpoint = MidPointCoordinate(_pFirst->coordinate, _pSecond->coordinate, 0.5);
    
    Compartment* c;
    try {c = GController::getCompartment(midpoint);}
    catch (std::exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
    
    if(c != _cCylinder->getCompartment()) {
        CCylinder* clone = _cCylinder->clone(c);
        setCCylinder(clone);
    }
    
#endif
}


