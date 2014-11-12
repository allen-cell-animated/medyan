//
//  ChemManager.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "ChemManager.h"

ChemManagerImpl* ChemManager::_pimpl = 0;

void ChemManager::setInstance(ChemManagerInitKey k, ChemManagerImpl *cii)
{
    if(_pimpl != 0) delete _pimpl;
    _pimpl=cii;
}

void ChemManager::initialize(ChemManagerGridKey k, ChemistryData& chem)
{
    _pimpl->initialize(chem);
}

CCylinder* ChemManager::createCCylinder(ChemManagerCylinderKey k, Filament* pf, Compartment* c,
                                        bool extensionFront, bool extensionBack, bool creation)

{
    return _pimpl->createCCylinder(pf, c, extensionFront, extensionBack, creation);
}

void ChemManager::updateCCylinder(ChemManagerCylinderKey k, CCylinder* cc) {
    
    _pimpl->updateCCylinder(cc);
}