
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------
#include "Bin.h"
#include "Cylinder.h"
void Bin::updatecindices(){
    cindicesvector.clear();
    cindicesvector.reserve(_cylinders.size());
    for(auto &c:_cylinders)
        cindicesvector.push_back(c->_dcIndex);
}
