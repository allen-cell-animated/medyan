//
//  ChemInitializer.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "ChemInitializer.h"

ChemInitializerImpl* ChemInitializer::_pimpl = 0;

void ChemInitializer::setInstance(ChemInitializerInitKey k, ChemInitializerImpl *cii)
{
    if(_pimpl != 0) delete _pimpl;
    _pimpl=cii;
}

void ChemInitializer::initialize(ChemInitializerGridKey k, ChemistrySpeciesAndReactions& chemSR)
{
    _pimpl->initialize(chemSR);
}

CCylinder* ChemInitializer::createCCylinder(ChemInitializerCylinderKey k, Filament* pf, Compartment* c,
                                                                          bool extensionFront, bool extensionBack, bool creation)

{
    return _pimpl->createCCylinder(pf, c, extensionFront, extensionBack, creation);
}

void ChemInitializer::updateCCylinder(ChemInitializerCylinderKey k, CCylinder* cc) {
    
    _pimpl->updateCCylinder(cc);
}