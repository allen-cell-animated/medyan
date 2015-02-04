
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

#include "NeighborListDB.h"

NeighborListDB* NeighborListDB::_instance = 0;

NeighborListDB* NeighborListDB::instance() {
    if(_instance==0)
        _instance = new NeighborListDB;
    return _instance;
}