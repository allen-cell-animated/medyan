
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "CompartmentContainer.h"

CompartmentGrid* CompartmentGrid::_instance = 0;

void CompartmentGrid::setInstance(int numCompartments)
{
    if(_instance != 0)
        delete _instance;
    _instance = new CompartmentGrid(numCompartments);
}

CompartmentGrid* CompartmentGrid::instance() {
    if(_instance==0)
        _instance = new CompartmentGrid(0);
    return _instance;
}

