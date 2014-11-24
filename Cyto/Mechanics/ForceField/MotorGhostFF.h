//
//  MotorGhostFF.h
//  Cyto
//
//  Created by Konstantin Popov on 9/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MotorGhostFF__
#define __Cyto__MotorGhostFF__

#include <iostream>
#include <vector>

#include "common.h"

#include "ForceField.h"

///FORWARD DECLARATIONS
class MotorGhostInteractions;

class MotorGhostFF : public ForceField {
    
private:
    vector <unique_ptr<MotorGhostInteractions>> _motorGhostInteractionVector;
    
public:
    MotorGhostFF(string& stretching, string& bending, string& twisting);
    
    // Public interfaces to compute forces:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
};




#endif /* defined(__Cyto__MotorGhostFF__) */
