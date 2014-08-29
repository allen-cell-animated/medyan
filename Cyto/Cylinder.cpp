//
//  Cylinder.cpp
//  Cyto
//
//  Created by James Komianos on 7/31/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "Cylinder.h"
#include "GController.h"
#include "ChemInitializer.h"

Cylinder::Cylinder(Filament* pf, Bead* firstBead, bool extensionFront, bool extensionBack) {
    
    _mCylinder = std::unique_ptr<MCylinder>(new MCylinder(pf, firstBead));
    
    ///Get compartment that this cylinder should be in
    Compartment* c = GController::getCompartment(0.0,0.0,0.0);
    
    _cCylinder = std::unique_ptr<CCylinder>(
        ChemInitializer::createCCylinder(ChemInitializerCylinderKey() , pf, c, extensionFront, extensionBack));
    
    _mCylinder->setCylinder(this);
    _cCylinder->setCylinder(this);
    
}


bool Cylinder::IfLast(){
    if (_ifLast) return true;
    
    return false;
}

void Cylinder::SetLast(bool b){ _ifLast = b;}
