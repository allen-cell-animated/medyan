//
//  MMotorGhostDB.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MMotorGhostDB.h"

MotorGhostDB* MotorGhostDB::_instance = 0;

MotorGhostDB* MotorGhostDB::Instance(MotorGhostDBKey k) {
    if(_instance==0)
        _instance = new MotorGhostDB;
    return _instance;
}