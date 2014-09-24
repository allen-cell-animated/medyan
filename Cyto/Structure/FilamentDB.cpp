//
//  FilamentDB.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "FilamentDB.h"

int FilamentDB::_currentFilamentID = 0;
FilamentDB* FilamentDB::_instance = 0;

FilamentDB* FilamentDB::Instance(FilamentDBKey k) {
    if(_instance==0)
        _instance = new FilamentDB;
    return _instance;
}