//
//  Cylinder.cpp
//  Cyto
//
//  Created by James Komianos on 7/31/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "Cylinder.h"

Cylinder::Cylinder(Filament* pf, Bead* firstBead, bool extension = false) {
    
    _mCylinder = new MCylinder(pf, firstBead);
    _cCylinder = ChemInitializer::createCCylinder(, pf->getLastCylinder()->getCCylinder(), extension);
    
}


bool Cylinder::IfLast(){
    if (_ifLast) return true;
    
    return false;
}

void Cylinder::SetLast(bool b){ _ifLast = b;}
