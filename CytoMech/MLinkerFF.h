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
#include "MLinkerStretching.h"
#include "MLinkerStretchingHarmonic.h"
#include "MLinkerDB.h"

class LinkerInteractions;

class LinkerFF : public ForceField
{
    
private:
    std::vector<std::unique_ptr<LinkerInteractions>> _linkerInteractionVector;
    
public:
    
    LinkerFF(std::string Stretching, std::string Bending, std::string Twisting );
    
    // Public interfaecs to compute forces:
    
    double ComputeEnergy(double d);
    
    void ComputeForces();
    
    void ComputeForcesAux();
    
};






#endif /* defined(__Cyto__MLinkerFF__) */
