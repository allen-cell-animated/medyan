
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
void Bin::vectorize(){
    delete []cindicesvector;
    std::cout<<_cylinders.size()<<endl;
    cindicesvector = new int[_cylinders.size()];
    int idx = 0;
    for (auto &ncylinder : _cylinders) {
        std::cout<<"CID "<<ncylinder->_dcIndex<<endl;
        idx++;
    }
}