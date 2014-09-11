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
#include "ForceField.h"
#include "MotorGhostStretching.h"
#include "MotorGhostStretchingHarmonic.h"
#include "MotorGhostDB.h"

class MotorGhostInteractions;

class MotorGhostFF : public ForceField {
    
private:
    std::vector <std::unique_ptr<MotorGhostInteractions>> _motorGhostInteractionVector;
    
public:
    MotorGhostFF(std::string Stretching, std::string Bending, std::string Twisting );
    
    // Public interfaecs to compute forces:
    double ComputeEnergy(double d);
    void ComputeForces();
    void ComputeForcesAux();
};




#endif /* defined(__Cyto__MotorGhostFF__) */
