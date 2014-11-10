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
    string _name;
    
public:
    virtual double computeEnergy(Cylinder*, Cylinder*,  double d) = 0;
    virtual void computeForces(Cylinder*, Cylinder*) = 0;
    virtual void computeForcesAux(Cylinder*, Cylinder*) = 0;
    
    // string getName() {return _name;}
    const string& getName() {return _name;}
    
};


#endif
