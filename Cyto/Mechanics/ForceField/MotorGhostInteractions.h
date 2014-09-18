//
//  MotorGhostInteractions.h
//  Cyto
//
//  Created by Konstantin Popov on 9/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_MotorGhostInteractions_h
#define Cyto_MotorGhostInteractions_h

#include "common.h"
#include <iostream>

class MotorGhost;

class MotorGhostInteractions {
private:
    std::string _name;
    
public:
    virtual double ComputeEnergy( MotorGhost*,  double d) = 0;
    virtual void ComputeForces(MotorGhost*) = 0;
    virtual void ComputeForcesAux(MotorGhost*) = 0;
    
    // std::string getName() {return _name;}
    const std::string& getName() {return _name;}
};

#endif /* defined(__Cyto__MotorGhostInteractions__) */
