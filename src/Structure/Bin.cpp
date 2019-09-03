
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

void Bin::addCylinder(Cylinder* c) {_cylinders.insert(c);}

//bincylinderdatatype& Bin::getCylinders() {return _cylinders;}

void Bin::removeCylinder(Cylinder* c) {
    auto it = _cylinders.find(c);
    if(it != _cylinders.end()) _cylinders.erase(it);
}
void Bin::updatecindices(){
    cindicesvector.clear();
    cindicesvector.reserve(_cylinders.size());
    for(auto &c:_cylinders)
        cindicesvector.push_back(c->getStableIndex());
}
