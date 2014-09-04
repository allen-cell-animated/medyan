//
//  MMotorGhostFF.h
//  Cyto
//
//  Created by Konstantin Popov on 9/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MMotorGhostFF__
#define __Cyto__MMotorGhostFF__

#include <iostream>
#include <vector>
#include "ForceField.h"


class MotorGhostInteractions;

class MotorGhostFFFF : private ForceField
{
    
private:
    std::vector <MotorGhostInteractions> _motorGhostInteractionVector;
    
    
public:
    
    
    // Public interfaecs to compute forces:
    
    double ComputeEnergy(double d);
    
    void ComputeForces();
    
    void ComputeForcesAux();
    
};




#endif /* defined(__Cyto__MMotorGhostFF__) */
