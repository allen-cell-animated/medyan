//
//  MLinkerFF.h
//  Cyto
//
//  Created by Konstantin Popov on 9/5/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MLinkerFF__
#define __Cyto__MLinkerFF__

#include <iostream>
#include <vector>
#include "ForceField.h"


class LinkerInteractions;

class LinkerFF : private ForceField
{
    
private:
    std::vector <LinkerInteractions> _linkerInteractionVector;
    
    
public:
    
    
    // Public interfaecs to compute forces:
    
    double ComputeEnergy(double d);
    
    void ComputeForces();
    
    void ComputeForcesAux();
    
};






#endif /* defined(__Cyto__MLinkerFF__) */
