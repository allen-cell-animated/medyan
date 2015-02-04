
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

#include "LinkerDB.h"

int LinkerDB::_currentLinkerID = 0;
LinkerDB* LinkerDB::_instance = 0;

LinkerDB* LinkerDB::instance() {
    if(_instance==0)
        _instance = new LinkerDB;
    return _instance;
}