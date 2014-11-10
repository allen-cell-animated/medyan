//
//  CylinderDB.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CylinderDB.h"

CylinderDB* CylinderDB::_instance = 0;

CylinderDB* CylinderDB::instance(CylinderDBKey k) {
    if(_instance==0)
        _instance = new CylinderDB;
    return _instance;
}