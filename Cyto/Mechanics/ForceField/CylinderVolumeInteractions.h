//
//  CylindrExclVolInteractions.h
//  Cyto
//
//  Created by Konstantin Popov on 10/29/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_CylinderVolInteractions_h
#define Cyto_CylinderVolInteractions_h


#include <iostream>
#include "common.h"

class Cylinder;

class CylinderVolumeInteractions {
private:
    std::string _name;
    
public:
    virtual double ComputeEnergy( Cylinder*, Cylinder*,  double d) = 0;
    virtual void ComputeForces(Cylinder*, Cylinder*) = 0;
    virtual void ComputeForcesAux(Cylinder*, Cylinder*) = 0;
    
    // std::string getName() {return _name;}
    const std::string& getName() {return _name;}
    
};


#endif
