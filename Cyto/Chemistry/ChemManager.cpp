//
//  ChemManager.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "ChemManager.h"

ChemManagerImpl* ChemManager::_pimpl = 0;

void ChemManager::setInstance(ChemManagerImpl *cii)
{
    if(_pimpl != 0) delete _pimpl;
    _pimpl=cii;
}

void ChemManager::initialize(ChemistryData& chem)
{
    _pimpl->initialize(chem);
}

void ChemManager::initializeCCylinder(CCylinder* cc, Filament* f, bool extensionFront, bool extensionBack, bool creation)

{
    _pimpl->initializeCCylinder(cc, f, extensionFront, extensionBack, creation);
}

void ChemManager::updateCCylinder(CCylinder* cc) {
    
    _pimpl->updateCCylinder(cc);
}