//
//  MBeadDB.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MBeadDB.h"

BeadDB* BeadDB::_instance = 0;

BeadDB* BeadDB::Instance(BeadDBKey k) {
    if(_instance==0)
        _instance = new BeadDB;
    return _instance;
}