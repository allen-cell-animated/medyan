//
//  MMotorGhostInteractions.h
//  Cyto
//
//  Created by Konstantin Popov on 9/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_MMotorGhostInteractions_h
#define Cyto_MMotorGhostInteractions_h


#include <iostream>

class MotorGhost;

class MotorGhostInteractions
{
    
    
private:
    
    std::string _name;
    
    
public:
    virtual double ComputeEnergy( MotorGhost*,  double d) = 0;
    virtual double ComputeForces(MotorGhost*) = 0;
    virtual double ComputeForcesAux(MotorGhost*) = 0;
    
    // std::string getName() {return _name;}
    
    std::string getName() {return _name;}
    
};




#endif
