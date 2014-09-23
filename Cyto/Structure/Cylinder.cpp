//
//  Cylinder.cpp
//  Cyto
//
//  Created by James Komianos on 7/31/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "Cylinder.h"

#include "ChemInitializer.h"
#include "Composite.h"
#include "GController.h"

Cylinder::Cylinder(Filament* pf, Bead* firstBead, Compartment* c, bool extensionFront, bool extensionBack) {
    
    setFilament(pf);
    
    _mCylinder = std::unique_ptr<MCylinder>(new MCylinder(pf, firstBead));
    _mCylinder->setCylinder(this);
    
#ifdef CHEMISTRY
    _cCylinder = std::unique_ptr<CCylinder>(
        ChemInitializer::createCCylinder(ChemInitializerCylinderKey() , pf, c, extensionFront, extensionBack));
    _cCylinder->setCylinder(this);
#endif

}

Cylinder::~Cylinder() {
    
#ifdef CHEMISTRY
    ChemInitializer::removeCCylinder(ChemInitializerCylinderKey(), _pFilament, IfLast() , !IfLast());
#endif
    
}



bool Cylinder::IfLast(){
    if (_ifLast) return true;
    
    return false;
}

void Cylinder::SetLast(bool b){ _ifLast = b;}
