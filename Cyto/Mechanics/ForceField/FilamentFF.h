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

///FORWARD DECLARATIONS
class FilamentInteractions;

class FilamentFF : public ForceField {
 
private:
    vector<unique_ptr<FilamentInteractions>> _filamentInteractionVector;
    
public:
    FilamentFF(string& stretching, string& bending, string& twisting);
    
   // Public interfaecs to compute forces:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
};


#endif /* defined(__Cyto__FilamentFF__) */
