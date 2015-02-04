
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

#include "BoundaryElementDB.h"

BoundaryElementDB* BoundaryElementDB::_instance = 0;

BoundaryElementDB* BoundaryElementDB::instance() {
    if(_instance==0)
        _instance = new BoundaryElementDB;
    return _instance;
}