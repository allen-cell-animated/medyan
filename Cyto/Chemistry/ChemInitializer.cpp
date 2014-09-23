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

void ChemInitializer::initialize(ChemInitializerGridKey k, ChemistrySetup& chemSetup)
{
    _pimpl->initialize(chemSetup);
}

CCylinder* ChemInitializer::createCCylinder(ChemInitializerCylinderKey k, Filament* pf, Compartment* c,
                                                                          bool extensionFront, bool extensionBack)

{
    return _pimpl->createCCylinder(pf, c, extensionFront, extensionBack);
}

void ChemInitializer::removeCCylinder(ChemInitializerCylinderKey k, Filament* pf,
                                                                    bool retractionFront, bool retractionBack)
{
    _pimpl->removeCCylinder(pf, retractionFront, retractionBack);
}
