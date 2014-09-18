//
//  FilamentFF.h
//  Cyto
//
//  Created by Konstantin Popov on 8/19/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__FilamentFF__
#define __Cyto__FilamentFF__

#include <iostream>
#include <vector>
#include <stdlib.h>

#include "common.h"
#include "ForceField.h"

class FilamentInteractions;

class FilamentFF : public ForceField {
 
private:
    std::vector<std::unique_ptr<FilamentInteractions>> _filamentInteractionVector;
    
public:
    FilamentFF(std::string& Stretching, std::string& Bending, std::string& Twisting );
    
   // Public interfaecs to compute forces:
    virtual double ComputeEnergy(double d);
    virtual void ComputeForces();
    virtual void ComputeForcesAux();
};


#endif /* defined(__Cyto__FilamentFF__) */
