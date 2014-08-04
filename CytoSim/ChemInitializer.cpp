//
//  ChemInitializer.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "ChemInitializer.h"
#include "ChemInitializerImpl.h"

namespace chem {
    
    ChemInitializerImpl* ChemInitializer::_pimpl = 0;
    
    void ChemInitializer::setInstance(ChemInitializerInitKey k, ChemInitializerImpl *cii)
    {
        _pimpl=cii;
    }
    
    void ChemInitializer::initializeGrid(ChemInitializerGridKey k, CompartmentGrid* grid)
    {
        _pimpl->initializeGrid(grid);
    }

    CCylinder* ChemInitializer::createCCylinder(ChemInitializerCylinderKey k, Compartment* c, CCylinder* lastCCylinder, bool extension)

    {
        return _pimpl->createCCylinder(c, lastCCylinder, extension);
    }
    
    void ChemInitializer::removeCCylinder(ChemInitializerCylinderKey k, CCylinder *cylinder)
    {
        _pimpl->removeCCylinder(cylinder);
    }
    
}; ///end namespace chem