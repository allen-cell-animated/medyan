//
//  MotorGhostInteractions.h
//  Cyto
//
//  Created by Konstantin Popov on 9/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_MotorGhostInteractions_h
#define Cyto_MotorGhostInteractions_h

#include <iostream>

#include "common.h"

//FORWARD DECLARATIONS
class MotorGhost;

class MotorGhostInteractions {
private:
    string _name;
    
public:
    virtual double computeEnergy( MotorGhost*,  double d) = 0;
    virtual void computeForces(MotorGhost*) = 0;
    virtual void computeForcesAux(MotorGhost*) = 0;

    const string& getName() {return _name;}
};

#endif /* defined(__Cyto__MotorGhostInteractions__) */
