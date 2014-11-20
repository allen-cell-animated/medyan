//
//  NeighborListDB.cpp
//  Cyto
//
//  Created by James Komianos on 11/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "NeighborListDB.h"

NeighborListDB* NeighborListDB::_instance = 0;

NeighborListDB* NeighborListDB::instance() {
    if(_instance==0)
        _instance = new NeighborListDB;
    return _instance;
}