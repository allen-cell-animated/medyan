
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

#include "FilamentDB.h"

int FilamentDB::_currentFilamentID = 0;
FilamentDB* FilamentDB::_instance = 0;

FilamentDB* FilamentDB::instance() {
    if(_instance==0)
        _instance = new FilamentDB;
    return _instance;
}