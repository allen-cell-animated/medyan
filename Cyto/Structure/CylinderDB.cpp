
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

#include "CylinderDB.h"

CylinderDB* CylinderDB::_instance = 0;

CylinderDB* CylinderDB::instance() {
    if(_instance==0)
        _instance = new CylinderDB;
    return _instance;
}