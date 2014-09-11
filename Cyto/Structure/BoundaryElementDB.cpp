//
//  BoundaryElementDB.cpp
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BoundaryElementDB.h"

BoundaryElementDB* BoundaryElementDB::_instance = 0;

BoundaryElementDB* BoundaryElementDB::Instance(BoundaryElementDBKey k) {
    if(_instance==0)
        _instance = new BoundaryElementDB;
    return _instance;
}