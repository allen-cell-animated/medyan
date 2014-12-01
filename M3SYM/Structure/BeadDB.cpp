
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "BeadDB.h"

BeadDB* BeadDB::_instance = 0;

BeadDB* BeadDB::instance() {
    if(_instance==0)
        _instance = new BeadDB;
    return _instance;
}