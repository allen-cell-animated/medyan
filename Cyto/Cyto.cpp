//
//  Cyto.cpp
//  Cyto
//
//  Created by James Komianos on 7/21/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include <iostream>
#include <vector>
#include "System.h"

int main(int argc, const char * argv[])
{
    
    std::vector<double> v0(3,0.0);
    std::vector<double> v1; v1.push_back(0.0); v1.push_back(0.0); v1.push_back(100.0);
    
    std::vector<std::vector<double>> vv; vv.push_back(v0); vv.push_back(v1);
    std::vector<std::vector<std::vector<double>>> vvv; vvv.push_back(vv);
    
    
    
    
    FullSystem<1> system;
    system.initialize(vvv);
    system.print();
    
}

