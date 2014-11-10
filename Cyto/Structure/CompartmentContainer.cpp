//
//  CompartmentContainer.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 9/5/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include "CompartmentContainer.h"

CompartmentGrid* CompartmentGrid::_instance = 0;

void CompartmentGrid::setInstance(CompartmentGridKey k, int numCompartments)
{
    if(_instance != 0)
        delete _instance;
    _instance = new CompartmentGrid(numCompartments);
}

CompartmentGrid* CompartmentGrid::instance(CompartmentGridKey k) {
    if(_instance==0)
        _instance = new CompartmentGrid(0);
    return _instance;
}

