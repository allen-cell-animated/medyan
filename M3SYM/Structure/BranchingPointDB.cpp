
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

#include "BranchingPointDB.h"

int BranchingPointDB::_currentBranchID = 0;
BranchingPointDB* BranchingPointDB::_instance = 0;

BranchingPointDB* BranchingPointDB::instance() {
    if(_instance==0)
        _instance = new BranchingPointDB;
    return _instance;
}