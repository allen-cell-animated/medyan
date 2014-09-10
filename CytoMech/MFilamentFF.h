//
//  MFilamentFF.h
//  Cyto
//
//  Created by Konstantin Popov on 8/19/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MFilamentFF__
#define __Cyto__MFilamentFF__

#include <iostream>
#include <vector>
#include "ForceField.h"
#include "MFilamentStretching.h"
#include "MFilamentStretchingHarmonic.h"
#include "MFilamentBending.h"
#include "MFilamentBendingHarmonic.h"
#include "MFilamentDB.h"
#include <stdlib.h>

class FilamentInteractions;

class FilamentFF : public ForceField
{
 
private:
    std::vector<std::unique_ptr<FilamentInteractions>> _filamentInteractionVector;
    
public:
    
    FilamentFF(std::string Stretching, std::string Bending, std::string Twisting );
    
   // Public interfaecs to compute forces:
    
    double ComputeEnergy(double d);
    
    void ComputeForces();
    
    void ComputeForcesAux();
    
};


#endif /* defined(__Cyto__MFilamentFF__) */
